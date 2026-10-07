"""Per-study survey of ENA attribute keys, classified into five categories.

Extends the three-category classifier in
`MMC/scripts/enaharm/scripts/survey_study_keys.py` (age / sex / disease) with
**antibiotic** and **geography**.

Two behaviours matter and are preserved from the MMC1 version:

  * Every key is surveyed -- custom_attributes keys AND standard portal columns
    AND alias/title/description -- because per-sample values frequently live in
    keys whose name carries no obvious token (e.g. an antibiotic arm encoded in
    `sample_alias`).
  * Generic keys (`treatment`, `medication`, `drug`, ...) are **value-inspected**
    rather than key-matched, so `treatment: amoxicillin` counts as antibiotic
    metadata while `treatment: placebo-controlled diet` does not.

Outputs:
  MMC2_key_survey.tsv        one row per (study, key) with fill rate + examples
  MMC2_study_column_map.csv  one row per study: chosen key per category
"""

import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import MMC, ROOT  # noqa: E402

csv.field_size_limit(sys.maxsize)

SAMPLES = os.path.join(ROOT, "MMC2_ena_samples.csv")
AGE_SEX_VOCAB = os.path.join(MMC, "claude dump", "age_sex_metadata_keys_all_studies.csv")
ABX_VOCAB = os.path.join(ROOT, "antibiotic_metadata_keys.csv")
GEO_VOCAB = os.path.join(ROOT, "geography_metadata_keys.csv")
SURVEY = os.path.join(ROOT, "MMC2_key_survey.tsv")
COLMAP = os.path.join(ROOT, "MMC2_study_column_map.csv")

# Standard portal columns worth surveying as pseudo-keys, mirroring STD_COLS in
# the MMC1 script plus the geo fields the portal exposes directly.
STD_COLS = [
    "sample_alias", "sample_title", "sample_title_xml", "sample_description",
    "experiment_title", "experiment_alias", "run_alias", "library_name",
    "age", "sex", "host_sex", "submitted_host_sex", "dev_stage",
    "disease", "host_status", "host_phenotype", "host_body_site", "host_genotype",
    "host", "isolation_source", "host_scientific_name",
    "country", "location", "lat", "lon",
    "environment_biome", "environmental_medium", "sampling_site",
]

DISEASE_KW = ["disease", "diagnos", "phenotype", "health", "condition", "clinical",
              "illness", "cancer", "disorder", "patholog", "host_status", "host_disease"]
AGE_RE = re.compile(r"(?:^|[^a-z])age(?:[^a-z]|$)")
SEX_RE = re.compile(r"(?:^|[^a-z])(?:sex|gender)(?:[^a-z]|$)")
DISGRP_RE = re.compile(r"(?:^|[^a-z])(?:group|status)(?:[^a-z]|$)")

# Antibiotic evidence inside a VALUE (used for generic keys such as `treatment`).
ABX_DRUGS = (r"amoxicillin|azithromycin|ciprofloxacin|metronidazole|vancomycin|"
             r"doxycycline|clindamycin|cephalosporin|ceftriaxone|cefazolin|penicillin|"
             r"rifaximin|rifampin|gentamicin|levofloxacin|moxifloxacin|trimethoprim|"
             r"sulfamethoxazole|nitrofurantoin|meropenem|piperacillin|tazobactam|"
             r"erythromycin|clarithromycin|tetracycline|minocycline|linezolid|"
             r"colistin|imipenem|ertapenem|ampicillin|cefepime|aztreonam|bactrim|"
             r"cayston|tobramycin|neomycin|streptomycin|kanamycin|polymyxin")
ABX_VALUE_RX = re.compile(
    r"\bantibiotic|\bantimicrobial|\bantibacterial|\babx\b|" + ABX_DRUGS, re.I
)
# Values that merely say yes/no are only meaningful when the KEY is antibiotic-
# specific; they carry no antibiotic evidence on their own.
YESNO_RX = re.compile(r"^\s*(y|n|yes|no|true|false|0|1|pos|neg|positive|negative)\s*$", re.I)

NULLISH = {"", "na", "n/a", "not applicable", "not collected", "not provided",
           "missing", "none", "null", "unknown", "not available", "unspecified",
           "not determined", "restricted access", "-"}

# Genuinely absent data. Distinguished from NEGATIVE below because in a
# disease or exposure column the two mean opposite things.
MISSING_TOKENS = {"", "-", "na", "n/a", "nan", "null", "missing", "unknown",
                  "not available", "not provided", "not collected", "not reported",
                  "not determined", "not applicable", "unspecified",
                  "restricted access", "not specified", "no data"}

# "No disease" / "no exposure". In a disease field `None` means the subject is
# a healthy control -- real information, not a blank. Treating it as null
# collapsed a genuine case/control contrast: PRJEB1415 records
# `disease status` as DOID:9256 on the cancer pool and None on the healthy
# pool, which scored as a single value and hid the contrast entirely.
NEGATIVE_TOKENS = {"none", "no", "healthy", "normal", "negative", "control",
                   "controls", "absent", "nil", "false", "0", "non",
                   "disease-free", "disease free", "nondiseased", "non-diseased",
                   "no disease", "unaffected", "healthy control", "nonsmoker",
                   "never", "no antibiotics", "not treated", "untreated"}

# Categories where a negative value is informative. Age, sex and geography have
# no meaningful "negative" -- a blank there is simply missing.
PRESENCE_ABSENCE_CATS = {"disease", "antibiotic"}


def value_kind(v):
    t = str(v).strip().lower()
    if t in MISSING_TOKENS:
        return "missing"
    if t in NEGATIVE_TOKENS:
        return "negative"
    return "positive"


def load_vocab(path, cols):
    out = [set() for _ in cols]
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf8") as fh:
        for row in csv.DictReader(fh):
            for i, c in enumerate(cols):
                v = (row.get(c) or "").strip().lower()
                if v and "summary" not in v and "included keys" not in v:
                    out[i].add(v)
    return out


AGE_INC, AGE_EXC, SEX_INC, SEX_EXC = load_vocab(
    AGE_SEX_VOCAB,
    ["Age Included", "Age Excluded", "Sex/Gender Included", "Sex/Gender Excluded"],
)
ABX_INC, ABX_EXC, ABX_INSPECT = load_vocab(
    ABX_VOCAB, ["Antibiotic Included", "Antibiotic Excluded", "Antibiotic Value-Inspect"]
)
GEO_INC, GEO_EXC = load_vocab(GEO_VOCAB, ["Geography Included", "Geography Excluded"])


def classify(key, values):
    """Categories for one attribute key, given the values it takes in a study."""
    k = key.strip().lower()
    cats = []

    is_age = (k in AGE_INC) or (AGE_RE.search(k) and k not in AGE_EXC)
    is_sex = (k in SEX_INC) or (k == "sex") or (SEX_RE.search(k) and k not in SEX_EXC)
    if is_age:
        cats.append("age")
    if is_sex:
        cats.append("sex")

    # --- antibiotic: key match, or value evidence under a generic key --------
    if k in ABX_EXC:
        pass
    elif k in ABX_INC or re.search(r"antibiotic|antimicrob|antibacter|\babx\b", k):
        cats.append("antibiotic")
    elif k in ABX_INSPECT or re.search(r"treatment|medicat|drug|therap|regimen|"
                                       r"interven|exposure|prophyla", k):
        # generic key -- only counts if a value actually names antibiotics
        if any(ABX_VALUE_RX.search(v) for v in values if not YESNO_RX.match(v)):
            cats.append("antibiotic")

    # --- geography ------------------------------------------------------------
    if k in GEO_EXC:
        pass
    elif k in GEO_INC or re.search(r"geo_loc|geographic|lat_lon|latitude|longitude|"
                                   r"\bcountry\b|\bnation\b|\bcontinent\b", k):
        cats.append("geography")

    # --- disease (never if already age/sex) -----------------------------------
    if not is_age and not is_sex:
        if any(w in k for w in DISEASE_KW):
            cats.append("disease")
        elif DISGRP_RE.search(k):
            cats.append("disease?")
    return cats


def informative(values, cats=()):
    """Distinct values that carry information.

    For disease and antibiotic keys a negative value ("None", "healthy", "no")
    is informative -- it says the subject is a control / was not exposed. For
    every other category a blank-like token is simply missing.
    """
    if any(c.rstrip("?") in PRESENCE_ABSENCE_CATS for c in cats):
        return {v for v in values if value_kind(v) != "missing"}
    return {v for v in values if v.strip().lower() not in NULLISH}


def contrast(values, cats=()):
    """Does this key express a real per-sample contrast?

    For presence/absence categories a contrast needs at least one POSITIVE
    value (a named disease, a drug) plus something to contrast it against --
    either a second distinct positive or an explicit negative. Variation
    between only negatives and missings is not a contrast: American Gut's
    `covid_chronic_conditions_*` fields take {false, Unknown, Not applicable},
    which says nobody was recorded as having the condition, not that samples
    differ by disease status.
    """
    info = informative(values, cats)
    if not any(c.rstrip("?") in PRESENCE_ABSENCE_CATS for c in cats):
        return len(info) > 1
    pos = {v for v in info if value_kind(v) == "positive"}
    neg = {v for v in info if value_kind(v) == "negative"}
    return len(pos) >= 2 or (len(pos) >= 1 and len(neg) >= 1)


def main():
    if not os.path.exists(SAMPLES):
        sys.exit(f"missing {SAMPLES} -- run fetch_ena_mmc2.py first")

    # study -> key -> list of values
    per_study = defaultdict(lambda: defaultdict(list))
    n_samples = Counter()

    with open(SAMPLES, encoding="utf8") as fh:
        for row in csv.DictReader(fh):
            study = row["source_study"]
            n_samples[study] += 1
            for c in STD_COLS:
                v = (row.get(c) or "").strip()
                if v:
                    per_study[study][c].append(v)
            raw = row.get("custom_attributes") or ""
            if raw:
                try:
                    for k, v in json.loads(raw).items():
                        if str(v).strip():
                            per_study[study][k].append(str(v).strip())
                except Exception:  # noqa: BLE001 - tolerate a malformed cell
                    pass

    survey_rows, colmap_rows = [], []
    CATS = ["age", "sex", "disease", "antibiotic", "geography"]

    for study, keys in sorted(per_study.items()):
        total = n_samples[study]
        best = {c: ("", 0, 0) for c in CATS}   # cat -> (key, distinct, filled)
        for key, values in sorted(keys.items()):
            cats = classify(key, values)
            distinct = informative(values, cats)
            varies = contrast(values, cats)
            filled = len(informative(values, cats))
            filled_n = len([v for v in values
                            if v.strip().lower() not in NULLISH
                            or (any(c.rstrip("?") in PRESENCE_ABSENCE_CATS for c in cats)
                                and value_kind(v) == "negative")])
            survey_rows.append({
                "study": study,
                "key": key,
                "n_samples": total,
                "n_filled": filled_n,
                "fill_rate": round(filled_n / total, 3) if total else 0,
                "n_distinct": len(distinct),
                "varies": varies,
                "categories": ";".join(cats),
                "examples": " | ".join(sorted(distinct)[:4])[:300],
            })
            for c in cats:
                c = c.rstrip("?")
                if c not in best:
                    continue
                # prefer a key that expresses a contrast, then higher fill
                score = (varies, filled_n)
                cur = best[c]
                if score > (cur[1] > 1, cur[2]):
                    best[c] = (key, 2 if varies else len(distinct), filled_n)

        row = {"study": study, "n_samples": total}
        for c in CATS:
            k, d, f = best[c]
            row[f"{c}_key"] = k
            row[f"{c}_distinct"] = d
            row[f"{c}_filled"] = f
            row[f"{c}_varies"] = bool(k) and d > 1
            row[f"{c}_any"] = bool(k) and f > 0
            # Presence semantics differ by what the field describes:
            #   age / sex / disease / antibiotic are per-SUBJECT attributes, so a
            #   single value repeated on every sample is a study-level drag-down,
            #   not per-sample evidence -> require variation.
            #   geography describes the SAMPLING SITE. A single-site study
            #   legitimately reports one country for every sample; demanding
            #   variation there would score complete metadata as missing.
            row[f"{c}_present"] = row[f"{c}_any"] if c == "geography" else row[f"{c}_varies"]
        colmap_rows.append(row)

    with open(SURVEY, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(survey_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(survey_rows)
    with open(COLMAP, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(colmap_rows[0].keys()))
        w.writeheader()
        w.writerows(colmap_rows)

    print(f"wrote {SURVEY}   ({len(survey_rows):,} study-key rows)")
    print(f"wrote {COLMAP}   ({len(colmap_rows):,} studies)")
    print(f"\nstudies surveyed: {len(colmap_rows):,}")
    print(f"{'category':12s} {'key found':>10s} {'populated':>10s} {'varies':>8s} "
          f"{'-> present':>11s}  rule")
    for c in CATS:
        any_key = sum(1 for r in colmap_rows if r[f"{c}_key"])
        filled = sum(1 for r in colmap_rows if r[f"{c}_any"])
        varies = sum(1 for r in colmap_rows if r[f"{c}_varies"])
        n = sum(1 for r in colmap_rows if r[f"{c}_present"])
        rule = "populated" if c == "geography" else "varies (drag-down rule)"
        print(f"{c:12s} {any_key:10,d} {filled:10,d} {varies:8,d} "
              f"{n:7,d} ({100 * n / len(colmap_rows):4.1f}%)  {rule}")

    cc = Counter()
    for r in survey_rows:
        for c in (r["categories"].split(";") if r["categories"] else []):
            cc[c] += 1
    print(f"\nclassified keys: {dict(cc)}")

    # Sanity: MMC1's collapsed-attribute-string bug, where a whole
    # "key:value key:value ..." run ended up as a single pseudo-key. The
    # signature is a colon inside the key, not mere length -- some studies
    # genuinely use very long descriptive attribute names.
    bad = [r for r in survey_rows if ":" in r["key"] and len(r["key"]) > 40]
    print(f"\ncollapsed key:value parse artefacts: {len(bad)}")
    for r in bad[:5]:
        print(f"   {r['study']}: {r['key'][:100]}")
    longk = [r for r in survey_rows if len(r["key"]) > 60]
    print(f"long but well-formed keys (>60 chars, no colon): "
          f"{len([r for r in longk if ':' not in r['key']])}")


if __name__ == "__main__":
    main()
