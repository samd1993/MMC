#!/usr/bin/env python3
"""Pull ENA run + sample metadata (including fastq URLs) for a list of accessions.

    python3 scripts/fetch_ena_metadata.py data/accessions.tsv --tier 1,2

Two stages, both cached under cache/ so the run is resumable:

  A. portal filereport / read_run per study accession -> one row per run, with
     fastq_ftp, fastq_bytes and md5 so the download step needs no extra calls.
  B. sample XML in batches -> SAMPLE_ATTRIBUTES. This is the ONLY place the
     author's custom columns live -- `host_disease`, `gynecologic_disord`,
     `sample_cog_status` and friends -- and it is where the disease label in
     data/disease_field_map.tsv will be found.

Writes ena_runs.tsv and ena_sample_attributes.tsv into the output directory.
"""
import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed

PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
XML = "https://www.ebi.ac.uk/ena/browser/api/xml/"
RETRIES, WORKERS, XML_BATCH = 3, 8, 50

FIELDS = [
    "run_accession", "experiment_accession", "experiment_title", "experiment_alias",
    "study_accession", "secondary_study_accession", "study_title",
    "sample_accession", "secondary_sample_accession", "sample_alias",
    "sample_title", "sample_description", "run_alias", "library_name",
    "library_strategy", "library_source", "library_selection", "library_layout",
    "instrument_platform", "instrument_model", "read_count", "base_count",
    "scientific_name", "tax_id", "host_scientific_name", "host_sex",
    "host_body_site", "host_phenotype", "host_status", "age", "sex",
    "disease", "isolate", "isolation_source", "collection_date", "country",
    "fastq_ftp", "fastq_bytes", "fastq_md5", "submitted_ftp",
]


def get(url, tries=RETRIES):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode("utf8", "replace")
        except (urllib.error.URLError, TimeoutError) as e:
            if i == tries - 1:
                raise
            time.sleep(2 ** i)
    raise RuntimeError("unreachable")


def valid_fields():
    """An unknown field name 400s the whole request, so intersect with live."""
    try:
        txt = get(PORTAL.replace("filereport", "returnFields")
                  + "?result=read_run&format=tsv")
        live = {l.split("\t")[0].strip() for l in txt.splitlines()[1:] if l.strip()}
        keep = [f for f in FIELDS if f in live]
        if dropped := [f for f in FIELDS if f not in live]:
            print(f"  note: ENA no longer serves {', '.join(dropped)}")
        return keep or FIELDS
    except Exception:                                    # noqa: BLE001
        return FIELDS


def fetch_runs(acc, fields, cache):
    path = os.path.join(cache, "runs", f"{acc}.tsv")
    if os.path.exists(path):
        return acc, open(path, encoding="utf8").read(), None
    url = PORTAL + "?" + urllib.parse.urlencode(
        {"accession": acc, "result": "read_run", "fields": ",".join(fields),
         "format": "tsv", "limit": "0"})
    try:
        txt = get(url)
    except urllib.error.HTTPError as e:
        return acc, None, f"HTTP {e.code}"
    except Exception as e:                               # noqa: BLE001
        return acc, None, str(e)[:80]
    if not txt.strip() or len(txt.strip().splitlines()) < 2:
        return acc, None, "no runs returned"
    with open(path, "w", encoding="utf8") as fh:
        fh.write(txt)
    return acc, txt, None


def fetch_attrs(batch):
    """SAMPLE_ATTRIBUTES for up to XML_BATCH sample accessions."""
    try:
        txt = get(XML + ",".join(batch))
        root = ET.fromstring(txt)
    except Exception as e:                               # noqa: BLE001
        return [], str(e)[:80]
    out = []
    for s in root.iter("SAMPLE"):
        acc = s.get("accession") or ""
        attrs = {}
        for a in s.iter("SAMPLE_ATTRIBUTE"):
            tag = (a.findtext("TAG") or "").strip()
            val = (a.findtext("VALUE") or "").strip()
            if tag and val:
                attrs[tag] = val
        title = s.findtext("TITLE") or ""
        out.append({"sample_accession": acc, "sample_title_xml": title,
                    "attributes": attrs})
    return out, None


def read_accessions(path, tier_filter=""):
    """accessions.tsv (accession/record_id/tier) or a bare one-per-line list."""
    tiers = {t.strip() for t in tier_filter.split(",") if t.strip()}
    lines = [l.rstrip("\n") for l in open(path, encoding="utf8") if l.strip()]
    if not lines:
        return []
    head = lines[0].split("\t")
    if "accession" not in head:
        if tiers:
            print("  note: --tier ignored, this file has no tier column")
        return sorted({l.strip() for l in lines})
    ai, ti = head.index("accession"), (head.index("tier") if "tier" in head else None)
    out = set()
    for l in lines[1:]:
        f = l.split("\t")
        if len(f) <= ai:
            continue
        if tiers and ti is not None and (f[ti] if len(f) > ti else "") not in tiers:
            continue
        out.add(f[ai].strip())
    return sorted(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("accessions",
                    help="data/accessions.tsv, or any file with one accession "
                         "per line")
    ap.add_argument("--tier", default="",
                    help="comma-separated tiers to keep, e.g. 1,2 "
                         "(only meaningful for accessions.tsv)")
    ap.add_argument("-o", "--outdir", default="ena")
    ap.add_argument("--cache", default="cache")
    ap.add_argument("--skip-attributes", action="store_true")
    args = ap.parse_args()

    os.makedirs(os.path.join(args.cache, "runs"), exist_ok=True)
    os.makedirs(args.outdir, exist_ok=True)
    accs = read_accessions(args.accessions, args.tier)
    fields = valid_fields()
    print(f"Stage A: {len(accs)} accessions, {len(fields)} fields")

    rows, failures, header = [], [], None
    with ThreadPoolExecutor(WORKERS) as ex:
        futs = {ex.submit(fetch_runs, a, fields, args.cache): a for a in accs}
        for i, f in enumerate(as_completed(futs), 1):
            acc, txt, err = f.result()          # 3-tuple. Unpacking it as 2
            if err:                             # silently loses every accession.
                failures.append((acc, err))
                continue
            lines = txt.rstrip("\n").split("\n")
            header = header or lines[0].split("\t")
            for ln in lines[1:]:
                rows.append([acc] + ln.split("\t"))
            if i % 50 == 0:
                print(f"  {i:,}/{len(accs):,}", flush=True)

    run_path = os.path.join(args.outdir, "ena_runs.tsv")
    with open(run_path, "w", encoding="utf8") as fh:
        fh.write("\t".join(["queried_accession"] + (header or [])) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"  {len(rows):,} runs -> {run_path}   ({len(failures)} accessions failed)")
    if failures:
        with open(os.path.join(args.outdir, "failures.tsv"), "w", encoding="utf8") as fh:
            fh.write("accession\treason\n")
            for a, e in failures:
                fh.write(f"{a}\t{e}\n")
    if args.skip_attributes:
        return

    si = (header or []).index("sample_accession") + 1 if header else None
    samples = sorted({r[si] for r in rows if si and len(r) > si and r[si]})
    print(f"\nStage B: SAMPLE_ATTRIBUTES for {len(samples):,} samples")
    cache_f = os.path.join(args.cache, "sample_attrs.jsonl")
    have = set()
    if os.path.exists(cache_f):
        for line in open(cache_f, encoding="utf8"):
            try:
                have.add(json.loads(line)["sample_accession"])
            except Exception:                            # noqa: BLE001
                pass
    todo = [s for s in samples if s not in have]
    batches = [todo[i:i + XML_BATCH] for i in range(0, len(todo), XML_BATCH)]
    print(f"  {len(have):,} cached, {len(todo):,} to fetch")
    with open(cache_f, "a", encoding="utf8") as ch, ThreadPoolExecutor(WORKERS) as ex:
        futs = {ex.submit(fetch_attrs, b): b for b in batches}
        for i, f in enumerate(as_completed(futs), 1):
            recs, err = f.result()
            if err:
                print(f"  batch failed: {err}")
            for r in recs:
                ch.write(json.dumps(r) + "\n")
            ch.flush()
            if i % 25 == 0:
                print(f"  batch {i:,}/{len(batches):,}", flush=True)

    keys, recs = set(), []
    for line in open(cache_f, encoding="utf8"):
        try:
            r = json.loads(line)
        except Exception:                                # noqa: BLE001
            continue
        if r["sample_accession"] in set(samples):
            recs.append(r)
            keys |= set(r["attributes"])
    cols = ["sample_accession", "sample_title_xml"] + sorted(keys)
    attr_path = os.path.join(args.outdir, "ena_sample_attributes.tsv")
    with open(attr_path, "w", encoding="utf8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in recs:
            row = {"sample_accession": r["sample_accession"],
                   "sample_title_xml": r["sample_title_xml"], **r["attributes"]}
            fh.write("\t".join(str(row.get(c, "")).replace("\t", " ")
                               for c in cols) + "\n")
    print(f"  {len(recs):,} samples, {len(keys):,} distinct attribute keys "
          f"-> {attr_path}")


if __name__ == "__main__":
    main()
