"""Harvest ENA metadata for the MMC2 INSDC accessions.

Two stages, both cached and resumable:

  A. portal `filereport` / `read_run` per study accession -> run rows, with the
     curated field set validated against the live `returnFields` listing.
  B. sample XML in batches of 50 -> SAMPLE_ATTRIBUTES (TAG/VALUE), which is the
     only place custom keys such as antibiotic fields live.

Writes MMC2_ena_samples.csv (one row per sample) and MMC2_ena_failures.tsv.

Re-running skips anything already in ena_cache/, so the harvest can be
interrupted and resumed freely.
"""

import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402

ACC_LIST = os.path.join(ROOT, "MMC2_insdc_accessions.txt")
CACHE = os.path.join(ROOT, "ena_cache")
RUN_CACHE = os.path.join(CACHE, "runs")
ATTR_CACHE = os.path.join(CACHE, "sample_attrs.jsonl")
OUT = os.path.join(ROOT, "MMC2_ena_samples.csv")
FAIL = os.path.join(ROOT, "MMC2_ena_failures.tsv")

WORKERS = 10
XML_BATCH = 50
# Samples per study pulled for XML attributes. Presence needs 1; detecting
# per-sample variation needs a spread. 100 evenly-spaced samples is ample and
# avoids the "first N are all controls" ordering bias.
MAX_SAMPLES_PER_STUDY = 100
RETRIES = 3

# Curated read_run fields (ena_portal_fields.md). Intersected with the live
# returnFields listing at runtime -- an unknown name 400s the whole batch.
# `common_name` is deliberately absent: valid for read_sample, not read_run.
WANT_FIELDS = """
run_accession experiment_accession experiment_alias experiment_title
study_accession secondary_study_accession study_alias study_title
sample_accession secondary_sample_accession sample_alias sample_title
sample_description submission_accession run_alias
library_strategy library_source library_selection library_layout library_name
library_construction_protocol instrument_platform instrument_model
sequencing_method target_gene checklist investigation_type
read_count base_count first_public last_updated first_created run_date
collection_date collection_date_start collection_date_end
tax_id tax_lineage host_scientific_name host_tax_id host_sex submitted_host_sex
host_status host_body_site host_genotype host_gravidity host_growth_conditions
host_phenotype age sex dev_stage disease isolate isolation_source host
project_name sample_capture_status sample_collection sample_material
sample_storage sampling_campaign sampling_site sampling_platform status
environment_biome environment_feature environment_material environmental_medium
broad_scale_environmental_context local_environmental_context
country location location_end location_start lat lon altitude depth elevation
ph salinity temperature tissue_lib culture_collection
""".split()

# Attributes present on essentially every ENA sample; they carry no biology and
# would swamp the custom-key survey.
BORING_TAGS = {
    "ena-first-public", "ena-last-update", "ena-checklist", "insdc first public",
    "insdc last update", "insdc center name", "insdc center alias", "insdc status",
    "external id", "submitter id", "biosamplemodel", "ncbi submission package",
    "ncbi submission model", "insdc secondary accession",
}


def fetch(url, timeout=180, data=None):
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(
                url, data=data, headers={"User-Agent": "mmc2-harvest/1.0"}
            )
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:  # noqa: BLE001 - network layer, retry everything
            last = e
            time.sleep(1.5 * (attempt + 1))
    raise last


def valid_fields():
    txt = fetch(
        "https://www.ebi.ac.uk/ena/portal/api/returnFields"
        "?result=read_run&dataPortal=ena&format=tsv"
    )
    have = {l.split("\t")[0].strip() for l in txt.splitlines()[1:] if l.strip()}
    use = [f for f in WANT_FIELDS if f in have]
    dropped = [f for f in WANT_FIELDS if f not in have]
    if dropped:
        print(f"  dropped {len(dropped)} unavailable field(s): {', '.join(dropped)}")
    return use


def fetch_runs(acc, fields):
    """Stage A: all runs for one study accession. Cached as TSV."""
    path = os.path.join(RUN_CACHE, f"{acc}.tsv")
    if os.path.exists(path):
        return acc, open(path, encoding="utf8").read(), None
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(
        {
            "accession": acc,
            "result": "read_run",
            "fields": ",".join(fields),
            "format": "tsv",
            "limit": "0",
        }
    )
    try:
        txt = fetch(url)
    except urllib.error.HTTPError as e:
        return acc, None, f"HTTP {e.code}"
    except Exception as e:  # noqa: BLE001
        return acc, None, type(e).__name__
    if not txt.strip() or not txt.splitlines()[0].startswith("run_accession"):
        return acc, None, "no runs returned"
    with open(path, "w", encoding="utf8") as fh:
        fh.write(txt)
    time.sleep(0.05)
    return acc, txt, None


def parse_runs(txt):
    lines = [l for l in txt.splitlines() if l.strip()]
    if len(lines) < 2:
        return []
    hdr = lines[0].split("\t")
    return [dict(zip(hdr, l.split("\t"))) for l in lines[1:]]


def fetch_attrs(batch):
    """Stage B: SAMPLE_ATTRIBUTES for a batch of sample accessions."""
    url = "https://www.ebi.ac.uk/ena/browser/api/xml/" + ",".join(batch)
    try:
        txt = fetch(url)
        root = ET.fromstring(txt)
    except Exception as e:  # noqa: BLE001
        return [], f"{type(e).__name__} on batch of {len(batch)}"
    out = []
    for sample in root.iter("SAMPLE"):
        ids = sample.find("IDENTIFIERS")
        pid = ids.findtext("PRIMARY_ID") if ids is not None else None
        if not pid:
            continue
        attrs = {}
        node = sample.find("SAMPLE_ATTRIBUTES")
        if node is not None:
            for a in node.findall("SAMPLE_ATTRIBUTE"):
                tag = (a.findtext("TAG") or "").strip()
                val = (a.findtext("VALUE") or "").strip()
                if tag and tag.lower() not in BORING_TAGS:
                    attrs[tag] = val
        out.append({
            "sample_accession": pid,
            "sample_title_xml": (sample.findtext("TITLE") or "").strip(),
            "custom_attributes": attrs,
        })
    time.sleep(0.05)
    return out, None


def evenly_spaced(seq, k):
    """Up to k items spread across seq, preserving order."""
    if len(seq) <= k:
        return list(seq)
    step = len(seq) / k
    return [seq[int(i * step)] for i in range(k)]


def main():
    os.makedirs(RUN_CACHE, exist_ok=True)
    accs = [l.strip() for l in open(ACC_LIST, encoding="utf8") if l.strip()]
    print(f"accessions to harvest: {len(accs):,}")

    print("validating portal fields...")
    fields = valid_fields()
    print(f"  using {len(fields)} read_run fields")

    # ---- Stage A ----------------------------------------------------------
    print("\nStage A: read_run per accession")
    runs_by_acc, failures = {}, []
    done = 0
    with ThreadPoolExecutor(WORKERS) as ex:
        futs = {ex.submit(fetch_runs, a, fields): a for a in accs}
        for f in as_completed(futs):
            acc, txt, err = f.result()
            done += 1
            if err:
                failures.append((acc, "read_run", err))
            else:
                runs_by_acc[acc] = parse_runs(txt)
            if done % 100 == 0:
                print(f"  {done:,}/{len(accs):,}  ok={len(runs_by_acc):,} fail={len(failures):,}",
                      flush=True)
    n_runs = sum(len(v) for v in runs_by_acc.values())
    print(f"  done: {len(runs_by_acc):,} accessions resolved, {n_runs:,} runs, "
          f"{len(failures):,} failures")

    # ---- collapse runs to samples ----------------------------------------
    # Keyed by (study, sample), NOT by sample alone: an accession pair such as
    # PRJNA123 + its secondary SRP456 resolves to the same samples, and a paper
    # citing both must see both studies. Deduping globally silently dropped the
    # second accession and could push a paper to "nothing resolved".
    print("\ncollapsing runs to samples")
    samples = {}          # (study, sample_accession) -> merged record
    study_samples = defaultdict(list)
    for acc, rows in runs_by_acc.items():
        for r in rows:
            s = r.get("sample_accession", "").strip()
            if not s:
                continue
            key = (acc, s)
            if key not in samples:
                samples[key] = {"source_study": acc, "sample_accession": s}
                study_samples[acc].append(s)
            rec = samples[key]
            for k, v in r.items():
                if k == "sample_accession" or not v:
                    continue
                prev = rec.get(k)
                if not prev:
                    rec[k] = v
                elif v not in prev.split("; "):   # lossless multi-run aggregation
                    rec[k] = prev + "; " + v
    n_uniq = len({s for _, s in samples})
    print(f"  {len(samples):,} study-sample rows "
          f"({n_uniq:,} distinct samples) across {len(study_samples):,} studies")

    # ---- Stage B ----------------------------------------------------------
    want = []
    for acc, sids in study_samples.items():
        want.extend(evenly_spaced(sids, MAX_SAMPLES_PER_STUDY))
    have = set()
    if os.path.exists(ATTR_CACHE):
        with open(ATTR_CACHE, encoding="utf8") as fh:
            for line in fh:
                try:
                    have.add(json.loads(line)["sample_accession"])
                except Exception:  # noqa: BLE001 - tolerate a truncated tail
                    pass
    todo = [s for s in want if s not in have]
    print(f"\nStage B: sample XML attributes")
    print(f"  targeted {len(want):,} samples (cap {MAX_SAMPLES_PER_STUDY}/study); "
          f"{len(have):,} cached, {len(todo):,} to fetch")

    batches = [todo[i:i + XML_BATCH] for i in range(0, len(todo), XML_BATCH)]
    fetched = 0
    with open(ATTR_CACHE, "a", encoding="utf8") as cache_fh:
        with ThreadPoolExecutor(WORKERS) as ex:
            futs = {ex.submit(fetch_attrs, b): b for b in batches}
            for i, f in enumerate(as_completed(futs), 1):
                recs, err = f.result()
                if err:
                    failures.append((",".join(futs[f][:3]) + "...", "sample_xml", err))
                for r in recs:
                    cache_fh.write(json.dumps(r) + "\n")
                    fetched += 1
                cache_fh.flush()
                if i % 50 == 0:
                    print(f"  batch {i:,}/{len(batches):,}  samples={fetched:,}", flush=True)
    print(f"  fetched {fetched:,} sample attribute records")

    # ---- merge and write --------------------------------------------------
    attrs = {}
    with open(ATTR_CACHE, encoding="utf8") as fh:
        for line in fh:
            try:
                r = json.loads(line)
                attrs[r["sample_accession"]] = r
            except Exception:  # noqa: BLE001
                pass

    import csv
    csv.field_size_limit(sys.maxsize)
    cols = (["source_study", "sample_accession", "custom_attributes", "sample_title_xml"]
            + [f for f in fields if f != "sample_accession"])
    n_with_attrs = 0
    with open(OUT, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for (_study, s), rec in samples.items():
            a = attrs.get(s)
            row = dict(rec)
            if a:
                row["custom_attributes"] = json.dumps(a["custom_attributes"])
                row["sample_title_xml"] = a["sample_title_xml"]
                n_with_attrs += 1
            else:
                row["custom_attributes"] = ""
                row["sample_title_xml"] = ""
            w.writerow(row)

    with open(FAIL, "w", encoding="utf8") as fh:
        fh.write("accession\tstage\treason\n")
        for a, st, r in failures:
            fh.write(f"{a}\t{st}\t{r}\n")

    print(f"\nwrote {OUT}")
    print(f"  {len(samples):,} samples; {n_with_attrs:,} carry XML attributes")
    print(f"wrote {FAIL}  ({len(failures):,} failures)")
    rate = 100 * len(runs_by_acc) / len(accs)
    print(f"\nretrieval rate: {len(runs_by_acc):,}/{len(accs):,} accessions ({rate:.1f}%)")


if __name__ == "__main__":
    main()
