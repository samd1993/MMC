#!/usr/bin/env python3
"""Label every fetched sample healthy vs disease from its study's disease field.

    python3 scripts/label_samples.py --ena ena -o sample_labels.tsv

Reads data/disease_field_map.tsv (which column a human reviewer found the
per-sample disease in, per study), then for each sample pulls that column out of
the ENA record and classifies the value.

The label is only ever as good as the field. Three things to know:

  * `disease_field_is_free_text = yes` means the disease lives in a sample name,
    alias, title or description, or in a coded column -- it needs parsing or a
    lookup, not a dictionary hit. Those come back as `needs_parsing`.
  * `healthy_cohort = yes` means the whole study is healthy participants. Every
    sample is a control regardless of what the field says.
  * a sample whose value is missing gets `unknown`, NOT healthy. Absence of a
    disease label is not evidence of health, and conflating the two is the
    single most common way this dataset gets misused.
"""
import argparse
import csv
import os
import re
import sys

csv.field_size_limit(10 ** 9)

# Same vocabulary the paper used, so labels here agree with the published tiers.
CONTROL = re.compile(
    r"^\s*(healthy|health|control|controls|hc|ctrl|normal|nc|non[- ]?diseased|"
    r"negative|neg|uninfected|unaffected|no[nt][- ]?(case|disease|ibd|cancer)|"
    r"reference|baseline healthy)\b", re.I)
# Questionnaire cohorts (American Gut / Microsetta) put the CONDITION in the
# field name and only yes/no in the value.
SURVEY_NEG = re.compile(r"i do not have this condition|^\s*(no|never|false)\s*$", re.I)
SURVEY_POS = re.compile(r"diagnosed by|^\s*(yes|true)\s*$", re.I)
SURVEY_UNSURE = re.compile(r"self[- ]?diagnosed|unspecified|not sure", re.I)
# Two different things that both look like "empty". Keeping them apart is the
# whole point: MISSING means nobody recorded anything; ABSENT means somebody
# recorded that the condition is not present, which is a healthy control.
# The published harmonized file folds ABSENT into MISSING and then fills the
# blank with the study's disease -- 1,678 samples in 19 studies are labelled
# with a disease their own record denies. Do not repeat that here.
MISSING = {"", "na", "n/a", "not collected", "not provided", "missing",
           "unknown", "nan", "-", "null", "not available", "unspecified",
           "not specified", "not determined"}
ABSENT = {"none", "no", "not applicable", "absent", "not present", "nil",
          "negative", "n", "0", "false"}
NULLISH = MISSING

AGE_FIELDS = ["host_age", "age", "host age", "age_years", "host_age_years"]
SEX_FIELDS = ["host_sex", "sex", "gender", "host sex", "host_gender"]


def norm(k):
    return re.sub(r"[^a-z0-9]", "", str(k).lower())


def pick(row, names):
    """First non-empty value among `names`, matched loosely on the key."""
    idx = {norm(k): k for k in row}
    for n in names:
        k = idx.get(norm(n))
        if k and str(row[k]).strip() and str(row[k]).strip().lower() not in MISSING:
            # note MISSING, not MISSING|ABSENT -- "none" is a real answer
            return str(row[k]).strip(), k
    return "", ""


def classify(value, field, free_text):
    v = value.strip().lower()
    if not v or v in MISSING:
        return "unknown", "no_value"
    if v in ABSENT:
        return "healthy", "explicit_absence"
    if SURVEY_NEG.search(value):
        return "healthy", "survey_response_negative"
    if SURVEY_UNSURE.search(value):
        return "unknown", "survey_response_unsure"
    if SURVEY_POS.search(value):
        # the condition is the FIELD NAME, e.g. skin_condition = "Diagnosed by..."
        return "disease", f"survey_response_positive:{field}"
    if CONTROL.match(value):
        return "healthy", "control_value"
    if free_text:
        return "needs_parsing", "free_text_or_coded_field"
    return "disease", "sample_value"


def read_tsv(path):
    with open(path, encoding="utf8", errors="replace", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def main():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser()
    ap.add_argument("--ena", default="ena", help="output dir of fetch_ena_metadata.py")
    ap.add_argument("--map", default=os.path.join(here, "data", "disease_field_map.tsv"))
    ap.add_argument("-o", "--out", default="sample_labels.tsv")
    args = ap.parse_args()

    runs_p = os.path.join(args.ena, "ena_runs.tsv")
    attr_p = os.path.join(args.ena, "ena_sample_attributes.tsv")
    if not os.path.exists(runs_p):
        sys.exit(f"missing {runs_p} -- run fetch_ena_metadata.py first")

    # study accession -> the studies that cite it, and their disease field
    by_acc = {}
    for r in read_tsv(args.map):
        for a in [x.strip() for x in r["accession_codes"].split(";") if x.strip()]:
            by_acc.setdefault(a, []).append(r)

    attrs = {}
    if os.path.exists(attr_p):
        for r in read_tsv(attr_p):
            attrs[r["sample_accession"]] = r
    else:
        print(f"note: no {attr_p}; author-named columns will be unavailable")

    rows, seen = [], set()
    for run in read_tsv(runs_p):
        s = run.get("sample_accession", "")
        qa = run.get("queried_accession", "")
        if not s or (qa, s) in seen:
            continue
        seen.add((qa, s))
        merged = {**run, **attrs.get(s, {})}
        for st in by_acc.get(qa, [{}]) or [{}]:
            fields = [f.strip() for f in str(st.get("disease_field", "")).split(";")
                      if f.strip()]
            free = str(st.get("disease_field_is_free_text", "")).strip() == "yes"
            val, used = pick(merged, fields) if fields else ("", "")
            if not val:                       # fall back to the standard columns
                val, used = pick(merged, ["host_disease", "disease",
                                          "host_phenotype", "host_status"])
            label, how = classify(val, used, free)
            if str(st.get("healthy_cohort", "")).strip() == "yes":
                label, how = "healthy", "healthy_cohort_study"
            age, _ = pick(merged, AGE_FIELDS)
            sex, _ = pick(merged, SEX_FIELDS)
            rows.append({
                "record_id": st.get("record_id", ""),
                "tier": st.get("tier", ""),
                "study_accession": qa,
                "sample_accession": s,
                "run_accession": run.get("run_accession", ""),
                "expected_disease_field": st.get("disease_field", ""),
                "field_used": used,
                "raw_value": val,
                "label": label,
                "label_method": how,
                "host_age": age,
                "host_sex": sex,
                "fastq_ftp": run.get("fastq_ftp", ""),
                "fastq_bytes": run.get("fastq_bytes", ""),
            })

    cols = list(rows[0]) if rows else []
    with open(args.out, "w", encoding="utf8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    print(f"{n:,} sample rows -> {args.out}")
    from collections import Counter
    for k, v in Counter(r["label"] for r in rows).most_common():
        print(f"  {k:<14} {v:>8,}  ({v / n:5.1%})")
    print("\nby method:")
    for k, v in Counter(r["label_method"] for r in rows).most_common(8):
        print(f"  {k:<34} {v:>8,}")


if __name__ == "__main__":
    main()
