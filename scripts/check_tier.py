#!/usr/bin/env python3
"""Which BHAM tier does a sequence deposit satisfy?

    python3 scripts/check_tier.py PRJEB10878
    python3 scripts/check_tier.py PRJNA389280 ERP012177 --json
    python3 scripts/check_tier.py --table my_sample_sheet.tsv

Give it an INSDC study accession (PRJ*/ERP/SRP/DRP/DRA...) and it pulls the
run and sample records from ENA, or give it the per-sample metadata sheet you
are about to submit (one row per sample, one column per attribute) and it reads
that instead. Either way it reports the tier the deposit earns, which keys
supplied the evidence, and what is missing for the next tier up.

The tiers are the ones used in the paper (Degregori et al., iMeta):

  Tier 1  per-sample host disease AND host age AND host sex
  Tier 2  per-sample host disease, from an attribute or from informative
          sample names (case_07 / control_12)
  Tier 3  valid accession whose sample IDs tell samples apart, no biology
  Tier 4  no accession, or samples that cannot be told apart

This is the paper's AUTOMATED pass, with the same rules and vocabularies (see
docs/TIERS.md). In the paper three curators then reviewed every Tier 1 and
Tier 2 call by hand and overruled the pipeline where they disagreed, so for a
study in data/studies.tsv the curated tier there takes precedence over this
script. Treat the output as a pre-submission check, not a verdict.

Two rules do most of the work and are worth knowing before you read a result:

  * A category only counts when its values DIFFER across samples. "Crohn's
    disease" written on every row is a study-level label copied down, and
    proves nothing about which sample came from whom.
  * In a disease field, "healthy" / "control" / "none" is information (it marks
    the control arm). "not applicable" / "missing" / blank is not.

Standard library only. Python 3.9+.
"""
import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict

csv.field_size_limit(sys.maxsize)
HERE = os.path.dirname(os.path.abspath(__file__))
VOCAB = os.path.join(HERE, "vocab")

PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
XML = "https://www.ebi.ac.uk/ena/browser/api/xml/"
XML_BATCH = 50

# ------------------------------------------------------------------ key rules
# Standard ENA portal columns surveyed as pseudo-keys, as in the paper.
STD_COLS = [
    "sample_alias", "sample_title", "sample_title_xml", "sample_description",
    "experiment_title", "experiment_alias", "run_alias", "library_name",
    "age", "sex", "host_sex", "submitted_host_sex", "dev_stage",
    "disease", "host_status", "host_phenotype", "host_body_site", "host_genotype",
    "host", "isolation_source", "host_scientific_name",
    "country", "location", "lat", "lon",
    "environment_biome", "environmental_medium", "sampling_site",
]
RUN_FIELDS = ["run_accession", "sample_accession", "study_accession"] + [
    c for c in STD_COLS if c != "sample_title_xml"]
ID_KEYS = ["sample_alias", "sample_title", "sample_title_xml", "sample_description"]
NAME_FIELDS = ["sample_alias", "sample_title", "sample_title_xml",
               "sample_description", "library_name", "experiment_title",
               "experiment_alias", "run_alias"]

DISEASE_KW = ["disease", "diagnos", "phenotype", "health", "condition", "clinical",
              "illness", "cancer", "disorder", "patholog", "host_status", "host_disease"]
AGE_RE = re.compile(r"(?:^|[^a-z])age(?:[^a-z]|$)")
SEX_RE = re.compile(r"(?:^|[^a-z])(?:sex|gender)(?:[^a-z]|$)")
DISGRP_RE = re.compile(r"(?:^|[^a-z])(?:group|status)(?:[^a-z]|$)")

ABX_DRUGS = (r"amoxicillin|azithromycin|ciprofloxacin|metronidazole|vancomycin|"
             r"doxycycline|clindamycin|cephalosporin|ceftriaxone|cefazolin|penicillin|"
             r"rifaximin|rifampin|gentamicin|levofloxacin|moxifloxacin|trimethoprim|"
             r"sulfamethoxazole|nitrofurantoin|meropenem|piperacillin|tazobactam|"
             r"erythromycin|clarithromycin|tetracycline|minocycline|linezolid|"
             r"colistin|imipenem|ertapenem|ampicillin|cefepime|aztreonam|bactrim|"
             r"cayston|tobramycin|neomycin|streptomycin|kanamycin|polymyxin")
ABX_VALUE_RX = re.compile(
    r"\bantibiotic|\bantimicrobial|\bantibacterial|\babx\b|" + ABX_DRUGS, re.I)
YESNO_RX = re.compile(r"^\s*(y|n|yes|no|true|false|0|1|pos|neg|positive|negative)\s*$", re.I)

NULLISH = {"", "na", "n/a", "not applicable", "not collected", "not provided",
           "missing", "none", "null", "unknown", "not available", "unspecified",
           "not determined", "restricted access", "-"}
MISSING_TOKENS = {"", "-", "na", "n/a", "nan", "null", "missing", "unknown",
                  "not available", "not provided", "not collected", "not reported",
                  "not determined", "not applicable", "unspecified",
                  "restricted access", "not specified", "no data"}
NEGATIVE_TOKENS = {"none", "no", "healthy", "normal", "negative", "control",
                   "controls", "absent", "nil", "false", "0", "non",
                   "disease-free", "disease free", "nondiseased", "non-diseased",
                   "no disease", "unaffected", "healthy control", "nonsmoker",
                   "never", "no antibiotics", "not treated", "untreated"}
PRESENCE_ABSENCE_CATS = {"disease", "antibiotic"}
CATS = ["disease", "age", "sex", "antibiotic", "geography"]


def _vocab(name, cols):
    out = [set() for _ in cols]
    path = os.path.join(VOCAB, name)
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf8") as fh:
        for row in csv.DictReader(fh):
            for i, c in enumerate(cols):
                v = (row.get(c) or "").strip().lower()
                if v and "summary" not in v and "included keys" not in v:
                    out[i].add(v)
    return out


AGE_INC, AGE_EXC, SEX_INC, SEX_EXC = _vocab(
    "age_sex_keys.csv",
    ["Age Included", "Age Excluded", "Sex/Gender Included", "Sex/Gender Excluded"])
ABX_INC, ABX_EXC, ABX_INSPECT = _vocab(
    "antibiotic_keys.csv",
    ["Antibiotic Included", "Antibiotic Excluded", "Antibiotic Value-Inspect"])
GEO_INC, GEO_EXC = _vocab("geography_keys.csv",
                          ["Geography Included", "Geography Excluded"])


def value_kind(v):
    t = str(v).strip().lower()
    if t in MISSING_TOKENS:
        return "missing"
    if t in NEGATIVE_TOKENS:
        return "negative"
    return "positive"


def classify(key, values):
    """Categories an attribute key belongs to, given the values it takes."""
    k = key.strip().lower()
    cats = []
    is_age = (k in AGE_INC) or (AGE_RE.search(k) and k not in AGE_EXC)
    is_sex = (k in SEX_INC) or (k == "sex") or (SEX_RE.search(k) and k not in SEX_EXC)
    if is_age:
        cats.append("age")
    if is_sex:
        cats.append("sex")
    if k in ABX_EXC:
        pass
    elif k in ABX_INC or re.search(r"antibiotic|antimicrob|antibacter|\babx\b", k):
        cats.append("antibiotic")
    elif k in ABX_INSPECT or re.search(r"treatment|medicat|drug|therap|regimen|"
                                       r"interven|exposure|prophyla", k):
        if any(ABX_VALUE_RX.search(v) for v in values if not YESNO_RX.match(v)):
            cats.append("antibiotic")
    if k in GEO_EXC:
        pass
    elif k in GEO_INC or re.search(r"geo_loc|geographic|lat_lon|latitude|longitude|"
                                   r"\bcountry\b|\bnation\b|\bcontinent\b", k):
        cats.append("geography")
    if not is_age and not is_sex:
        if any(w in k for w in DISEASE_KW):
            cats.append("disease")
        elif DISGRP_RE.search(k):
            cats.append("disease?")
    return cats


def _pa(cats):
    return any(c.rstrip("?") in PRESENCE_ABSENCE_CATS for c in cats)


def informative(values, cats=()):
    if _pa(cats):
        return {v for v in values if value_kind(v) != "missing"}
    return {v for v in values if v.strip().lower() not in NULLISH}


def contrast(values, cats=()):
    """Does this key express a real per-sample contrast?"""
    info = informative(values, cats)
    if not _pa(cats):
        return len(info) > 1
    pos = {v for v in info if value_kind(v) == "positive"}
    neg = {v for v in info if value_kind(v) == "negative"}
    return len(pos) >= 2 or (len(pos) >= 1 and len(neg) >= 1)


# ----------------------------------------------------- informative sample names
B, E = r"(?<![A-Za-z])", r"(?![A-Za-z])"
HEALTHY_RE = re.compile("|".join([
    B + r"control(?:s|ed|led)?" + E, B + r"ctrl" + E, B + r"ctl" + E,
    B + r"healthy" + E, B + r"normal" + E, B + r"donor" + E,
    B + r"unaffected" + E, B + r"negative" + E,
    B + r"non[\s_\-]?(?:asthma|allerg|diabet|smoker|cancer|disease|tumor|tumour|ad|aad|ibd|cd|uc|ms|crc|crs|als|nafld|hpv|hiv|ppt|msm)" + E,
    B + r"(?:hc|nc|hv)[_\s\-]?\d", B + r"(?:hc|nc|hv)\d",
    B + r"placebo" + E, B + r"baseline" + E, B + r"vehicle" + E, B + r"sham" + E,
    B + r"(?:con|cont|ctl|ctr)[_\s\-]?\d", B + r"(?:con|cont|ctl|ctr)\d",
    B + r"ova" + E,
]), re.I)
DISEASE_RE = re.compile("|".join([
    B + r"cases?" + E, B + r"patient(?:s)?" + E,
    B + r"tumou?r(?:s|ous)?" + E, B + r"cancer(?:s|ous)?" + E,
    B + r"carcinoma" + E, B + r"malignant" + E, B + r"neoplasm" + E,
    B + r"diseased?" + E, B + r"infected" + E, B + r"infection" + E,
    B + r"positive" + E, B + r"severe" + E, B + r"moderate" + E, B + r"mild" + E,
    B + r"lesional" + E, B + r"lesion" + E, B + r"symptomatic" + E,
    B + r"affected" + E, B + r"sick" + E, B + r"diagnosed" + E,
    B + r"(?:pt|pts)[_\s\-]?\d", B + r"inflam(?:ed|mation|matory)" + E,
    B + r"obes(?:e|ity)" + E, B + r"overweight" + E, B + r"smoker" + E,
    B + r"diabet(?:ic|es)" + E, B + r"nafld" + E, B + r"nash" + E,
    B + r"melanoma" + E, B + r"eczema" + E, B + r"crohn" + E, B + r"celiac" + E,
    B + r"lam" + E, B + r"ppt" + E, B + r"caox" + E, B + r"caries" + E,
    B + r"(?:cd|uc|ibd|ad|als|ms|crc|copd|hpv|hiv|sars|crs|t1d|t2d|nash|aad|asd|tb|eae|mog|pcos|dcis|bc)[_\s\-]?\d",
    B + r"(?:cd|uc|ibd|ad|als|ms|crc|copd|hpv|hiv|sars|crs|t1d|t2d|nash|aad|asd|tb|eae|mog|pcos|dcis)\d",
    B + r"glaucoma" + E, B + r"ischemia" + E, B + r"enterocolitis" + E,
    B + r"volvulus" + E, B + r"atresia" + E, B + r"perforation" + E,
    B + r"transmitter" + E, B + r"pcos" + E, B + r"preterm" + E, B + r"dry eye" + E,
    B + r"ascus" + E, B + r"lsil" + E, B + r"hsil" + E,
]), re.I)


def _sig(t):
    return (frozenset(m.group(0).lower() for m in HEALTHY_RE.finditer(t)),
            frozenset(m.group(0).lower() for m in DISEASE_RE.finditer(t)))


def names_evidence(rows):
    """Disease contrast hidden in sample names / free-text attribute values.

    A contrast needs more than one distinct token signature across samples AND
    at least one disease token. One signature on every row is a drag-down.
    """
    sig_names, sig_attrs, n_dis, tokens = set(), set(), 0, set()
    for r in rows:
        name_text = " ".join(str(r["std"].get(f) or "") for f in NAME_FIELDS)
        attr_text = " ".join(str(v) for v in r["attrs"].values())
        hs, ds = _sig(name_text + " " + attr_text)
        if ds:
            n_dis += 1
            tokens |= ds
        nh, nd = _sig(name_text)
        if nh or nd:
            sig_names.add((nh, nd))
        if r["attrs"] and (hs or ds):
            sig_attrs.add((hs, ds))
    varies = len(sig_names) > 1 or len(sig_attrs) > 1
    return varies and n_dis > 0, sorted(tokens)[:6]


# ------------------------------------------------------------------ assessment
def assess(rows):
    """rows: [{"std": {portal column: value}, "attrs": {custom key: value}}]."""
    keys = defaultdict(list)
    for r in rows:
        for c in STD_COLS:
            v = str(r["std"].get(c) or "").strip()
            if v:
                keys[c].append(v)
        for k, v in r["attrs"].items():
            if str(v).strip():
                keys[k].append(str(v).strip())

    best = {c: ("", 0, 0) for c in CATS}
    examples, id_varies = {}, False
    for key, values in sorted(keys.items()):
        cats = classify(key, values)
        distinct = informative(values, cats)
        varies = contrast(values, cats)
        filled_n = len([v for v in values
                        if v.strip().lower() not in NULLISH
                        or (_pa(cats) and value_kind(v) == "negative")])
        if key in ID_KEYS and varies:
            id_varies = True
        for c in cats:
            c = c.rstrip("?")
            if c not in best:
                continue
            cur = best[c]
            if (varies, filled_n) > (cur[1] > 1, cur[2]):
                best[c] = (key, 2 if varies else len(distinct), filled_n)
                examples[c] = sorted(distinct)[:4]

    present = {}
    for c in CATS:
        k, d, f = best[c]
        present[c] = (bool(k) and f > 0) if c == "geography" else (bool(k) and d > 1)
    from_names, tokens = names_evidence(rows)
    disease = present["disease"] or from_names

    if not rows:
        tier = "4"
    elif disease and present["age"] and present["sex"]:
        tier = "1"
    elif disease:
        tier = "2"
    elif id_varies:
        tier = "3"
    else:
        tier = "4"

    missing = []
    if tier == "4":
        missing.append("sample identifiers (alias/title) that differ between samples"
                       if rows else "a public INSDC deposit with samples")
    if tier in ("4", "3"):
        missing.append("a per-sample host disease field (e.g. host_disease) whose "
                       "values differ across samples, controls included")
    if tier in ("4", "3", "2"):
        for c in ("age", "sex"):
            if not present[c]:
                missing.append(f"per-sample host {c} (e.g. host_{c})")

    return {
        "tier": tier,
        "n_rows": len(rows),
        "n_samples": len({r["std"].get("sample_accession") or id(r) for r in rows}),
        "evidence": {c: {"key": best[c][0], "present": present[c],
                         "examples": examples.get(c, [])} for c in CATS},
        "disease_from_names": from_names,
        "name_tokens": tokens,
        "sample_ids_differ": id_varies,
        "missing_for_next_tier": missing,
    }


# ------------------------------------------------------------------ inputs
def _get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode("utf8", "replace")
        except urllib.error.HTTPError:
            raise
        except (urllib.error.URLError, TimeoutError):
            if i == tries - 1:
                raise
            time.sleep(2 ** i)


def _portal_fields():
    try:
        txt = _get(PORTAL.replace("filereport", "returnFields")
                   + "?result=read_run&format=tsv")
        live = {l.split("\t")[0].strip() for l in txt.splitlines()[1:] if l.strip()}
        return [f for f in RUN_FIELDS if f in live] or RUN_FIELDS
    except Exception:                                    # noqa: BLE001
        return RUN_FIELDS


def from_ena(acc, max_samples=0, quiet=False):
    """Run rows for an accession, each joined to its sample's attributes."""
    url = PORTAL + "?" + urllib.parse.urlencode(
        {"accession": acc, "result": "read_run", "fields": ",".join(_portal_fields()),
         "format": "tsv", "limit": "0"})
    try:
        txt = _get(url)
    except urllib.error.HTTPError as e:
        if e.code in (400, 404):
            return [], f"ENA has no public runs for {acc}"
        raise
    lines = [l for l in (txt or "").splitlines() if l.strip()]
    if len(lines) < 2:
        return [], f"ENA has no public runs for {acc}"
    head = lines[0].split("\t")
    runs = [dict(zip(head, l.split("\t"))) for l in lines[1:]]
    samples = sorted({r.get("sample_accession", "") for r in runs} - {""})
    if max_samples and len(samples) > max_samples:
        samples = samples[:max_samples]
        keep = set(samples)
        runs = [r for r in runs if r.get("sample_accession") in keep]
    attrs, titles = {}, {}
    for i in range(0, len(samples), XML_BATCH):
        batch = samples[i:i + XML_BATCH]
        if not quiet:
            print(f"  {acc}: sample records {i + len(batch):,}/{len(samples):,}",
                  file=sys.stderr, flush=True)
        try:
            root = ET.fromstring(_get(XML + ",".join(batch)))
        except Exception as e:                           # noqa: BLE001
            print(f"  warning: sample batch failed ({str(e)[:60]})", file=sys.stderr)
            continue
        for s in root.iter("SAMPLE"):
            a = s.get("accession") or ""
            titles[a] = s.findtext("TITLE") or ""
            attrs[a] = {(x.findtext("TAG") or "").strip(): (x.findtext("VALUE") or "").strip()
                        for x in s.iter("SAMPLE_ATTRIBUTE") if (x.findtext("TAG") or "").strip()}
    rows = []
    for r in runs:
        sa = r.get("sample_accession", "")
        std = dict(r, sample_title_xml=titles.get(sa, ""))
        rows.append({"std": std, "attrs": attrs.get(sa, {})})
    return rows, None


def from_table(path):
    """A submitter's sheet: one row per sample, one column per attribute."""
    delim = "," if path.lower().endswith(".csv") else "\t"
    with open(path, encoding="utf8", newline="") as fh:
        data = list(csv.DictReader(fh, delimiter=delim))
    std_set = set(STD_COLS) | {"sample_accession"}
    rows = []
    for d in data:
        d = {(k or "").strip(): (v or "").strip() for k, v in d.items()}
        # the usual sample-name column goes in as the alias
        for alt in ("sample_name", "sample_id", "#sampleid", "sampleid"):
            if alt in {k.lower() for k in d} and "sample_alias" not in d:
                d["sample_alias"] = next(v for k, v in d.items() if k.lower() == alt)
        rows.append({"std": {k: v for k, v in d.items() if k in std_set},
                     "attrs": {k: v for k, v in d.items() if k not in std_set}})
    return rows


# ------------------------------------------------------------------ report
LABEL = {"1": "Tier 1 - disease, age and sex per sample",
         "2": "Tier 2 - disease per sample",
         "3": "Tier 3 - distinguishable samples, no biology",
         "4": "Tier 4 - not reusable"}


def curated_tier(acc):
    """The paper's hand-curated tier(s) for an accession, if it is in data/."""
    path = os.path.join(os.path.dirname(HERE), "data", "accessions.tsv")
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf8") as fh:
        return [(r["record_id"], r["tier"]) for r in csv.DictReader(fh, delimiter="\t")
                if r.get("accession", "").upper() == acc.upper()]


def show(name, res, curated):
    print(f"\n{name}\n  {LABEL[res['tier']]}"
          f"   ({res['n_samples']:,} samples, {res['n_rows']:,} records)")
    for c in ("disease", "age", "sex"):
        e = res["evidence"][c]
        mark = "yes" if e["present"] else "no "
        ex = ", ".join(e["examples"])[:70]
        print(f"  {c:8s} {mark}  {e['key'] or '-'}" + (f"  [{ex}]" if e["key"] else ""))
    if res["disease_from_names"]:
        print(f"  disease  yes  from sample names  [{', '.join(res['name_tokens'])}]")
    print(f"  sample IDs differ: {'yes' if res['sample_ids_differ'] else 'no'}")
    for m in res["missing_for_next_tier"]:
        print(f"  to move up: {m}")
    for rid, t in curated:
        print(f"  curated tier in the paper: Tier {t} ({rid})")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accessions", nargs="*", help="INSDC study/project accessions")
    ap.add_argument("--table", help="a per-sample metadata sheet (.tsv or .csv) "
                                    "to check before submission")
    ap.add_argument("--max-samples", type=int, default=0,
                    help="read at most this many samples per accession (0 = all)")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    a = ap.parse_args()
    if not a.accessions and not a.table:
        ap.error("give an accession or --table")

    out = {}
    if a.table:
        out[a.table] = (assess(from_table(a.table)), [])
    for acc in a.accessions:
        rows, err = from_ena(acc.strip(), a.max_samples, quiet=a.json)
        res = assess(rows)
        if err:
            res["note"] = err
        out[acc] = (res, curated_tier(acc))

    if a.json:
        print(json.dumps({k: dict(v[0], curated=[{"record_id": r, "tier": t}
                                                 for r, t in v[1]])
                          for k, v in out.items()}, indent=2))
        return
    for name, (res, cur) in out.items():
        show(name, res, cur)
        if res.get("note"):
            print(f"  note: {res['note']}")


if __name__ == "__main__":
    main()
