#!/usr/bin/env python3
"""One per-sample metadata table for EVERY tier, not just 1 and 2.

Supersedes harmonize_t12_samples.py, which it also regenerates as a strict
subset so the paper's supplementary file and the collaborator's master can
never disagree.

What changes against the old T1/T2-only harmonization:

  * Tiers 3 and 4 are included wherever the deposit resolves. They get their
    record's real values and NOTHING inferred -- no study-MeSH fallback, because
    "this study is about colorectal cancer" says nothing about a Tier 3 sample.
    A Tier 3/4 sample that does turn out to carry a usable disease value is
    marked `unexpected_annotation`, which is a tiering error worth chasing.

  * "the condition is absent" is no longer read as "nothing was recorded".
    NULLISH used to swallow the literal value `none`, after which the blank was
    filled from the study's MeSH heading -- so 1,678 samples across 19 studies
    carried a disease their own record denies (PRJEB31817 alone: 468 samples
    answering `none` to `gastrointestinal tract disorder`, counted as cases).
    ABSENT values now resolve to controls and are never overwritten.

Writes MMC2_sample_master.tsv (all tiers), MMC2_sample_master_by_study.tsv, and
regenerates MMC2_harmonized_T1T2_samples.tsv / _by_study.tsv.
"""
import csv
import json
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402
from mondo_map import NULLISH, MondoMapper  # noqa: E402

csv.field_size_limit(10 ** 9)
MASTER = os.path.join(ROOT, "MMC2_master_reviewed.tsv")
STUDY = os.path.join(ROOT, "MMC2_study_metadata.tsv")
ANALYSIS = os.path.join(ROOT, "MMC2_analysis_set.tsv")
RUNS = os.path.join(ROOT, "ena_cache", "runs")
ATTRS = os.path.join(ROOT, "ena_cache", "sample_attrs.jsonl")
INDEX = os.path.join(ROOT, "ontology", "mondo_index.json")
OUT = os.path.join(ROOT, "MMC2_sample_master.tsv")
OUT_STUDY = os.path.join(ROOT, "MMC2_sample_master_by_study.tsv")
OUT_T12 = os.path.join(ROOT, "MMC2_harmonized_T1T2_samples.tsv")
OUT_T12_STUDY = os.path.join(ROOT, "MMC2_harmonized_T1T2_by_study.tsv")

SURVEY_NEG = re.compile(r"i do not have this condition|^no$|^never$|"
                        r"^not diagnosed", re.I)
SURVEY_POS = re.compile(r"diagnosed by a medical professional|self[- ]diagnosed|"
                        r"^diagnosed\b|^yes\b", re.I)
SURVEY_UNSURE = re.compile(r"^unspecified$|^unsure$|^don'?t know$|^na$", re.I)

# NULLISH conflates three different facts. Split them:
#   MISSING        nobody wrote anything down          -> study-level fallback
#   UNINFORMATIVE  somebody wrote "we don't know"      -> no label, and no
#                  fallback either: the record HAS spoken and it said nothing
#   ABSENT         somebody wrote that the condition is not there -> control
# "not applicable" belongs in UNINFORMATIVE, not ABSENT: 3,611 of its 4,868 uses
# are American Gut's generic blank marker (Lynn, on MMC2-04547: "says 'not
# applicable' for sex and age of every sample"), not a negative answer.
ABSENT = {"none", "no", "absent", "not present", "nil", "negative", "false",
          "no disease", "non", "0"}
UNINFORMATIVE = {"not applicable", "unknown", "not collected", "not provided",
                 "missing", "na", "n/a", "not available", "unspecified", "nan",
                 "-", "null"}
MISSING = {""} | {v for v in NULLISH if v not in ABSENT | UNINFORMATIVE}
# A negatively-phrased field inverts the reading: `is_healthy = no` is a CASE,
# not a control, and `is_healthy = yes` is a control.
HEALTH_FIELD = re.compile(r"^\s*(is[_ ]?)?health(y|_?status|_?state)?\s*$", re.I)

AGE_FIELDS = ["age", "host_age", "host age", "age_years", "host_age_years"]
SEX_FIELDS = ["host_sex", "sex", "gender", "host sex", "submitted_host_sex"]
DISEASE_FALLBACK = ["host_disease", "disease", "host_phenotype", "host_status",
                    "host_health_state", "host_disease_stat"]


def _scan(rec, cust, names, skip):
    for n in names:
        v = rec.get(n)
        if v and str(v).strip().lower() not in skip:
            return str(v).strip(), n
        for k, vv in (cust or {}).items():
            if (k.strip().lower() == n.lower() and vv
                    and str(vv).strip().lower() not in skip):
                return str(vv).strip(), k
    return "", ""


def first_value(rec, cust, names):
    """First INFORMATIVE value among `names`; failing that, the uninformative
    token that was actually recorded, so the caller can tell "we don't know"
    apart from "nobody said"."""
    hit = _scan(rec, cust, names, MISSING | UNINFORMATIVE)
    return hit if hit[0] else _scan(rec, cust, names, MISSING)


def resolve(raw, used, fields, mm, s_mid, s_src, tier):
    """-> (mondo_id, case_control, mapping_method). Nothing is inferred for
    tiers 3 and 4: their defining property is that the record does not say."""
    infer_study = tier in ("1", "2")
    low = raw.strip().lower()
    if raw and low in UNINFORMATIVE:
        return None, "", "uninformative_value"
    if raw and low in ABSENT:
        if HEALTH_FIELD.match(str(used)):     # is_healthy = no -> a case
            mid, hw = (s_mid, "study_mesh") if (s_mid and infer_study) else (None, "")
            return mid, "case", ("health_field_negated_" + hw if mid
                                 else "health_field_negated")
        return None, "control", "explicit_absence"
    if raw and HEALTH_FIELD.match(str(used)) and re.match(r"^(yes|true|1)$", low, re.I):
        return None, "control", "health_field_affirmed"
    if raw and (SURVEY_NEG.search(raw) or SURVEY_POS.search(raw)
                or SURVEY_UNSURE.search(raw)):
        if SURVEY_NEG.search(raw):
            return None, "control", "survey_response_negative"
        if SURVEY_UNSURE.search(raw):
            return None, "", "survey_response_unspecified"
        mid, hw = mm.lookup(used)                  # the FIELD NAME is the condition
        if not mid and s_mid and infer_study:
            mid, hw = s_mid, "study_mesh"
        return mid, "case", (f"survey_response_positive_field_{hw}" if mid
                             else "survey_response_positive_unmapped")
    if raw:
        if mm.is_control(raw):
            return None, "control", "control_value"
        mid, hw = mm.lookup(raw)
        if mid:
            return mid, "case", f"sample_value_{hw}"
        if s_mid and infer_study:
            return s_mid, "case", "study_mesh_heading"
        return None, "", "unmapped"
    if s_mid and infer_study:
        return s_mid, "", "study_mesh_heading_no_sample_value"
    return None, "", ("no_disease_value" if infer_study
                      else f"tier{tier}_no_annotation")


def main():
    mm = MondoMapper(INDEX)
    m = pd.read_csv(MASTER, sep="\t", low_memory=False, keep_default_na=False,
                    dtype=str)
    sm = pd.read_csv(STUDY, sep="\t", low_memory=False, keep_default_na=False,
                     dtype=str)[["record_id", "accession_codes"]]
    m = m.drop(columns=[c for c in ("accession_codes",) if c in m.columns]) \
         .merge(sm, on="record_id", how="left")
    m = m[m.final_tier.isin(["1", "2", "3", "4"]) & m.remove_flag.ne("yes")].copy()

    an = pd.read_csv(ANALYSIS, sep="\t", low_memory=False, keep_default_na=False,
                     dtype=str).set_index("record_id")
    study_term = {}
    for r in m.itertuples():
        mid, src = None, ""
        if r.record_id in an.index:
            for h in str(an.loc[r.record_id, "mesh_disease_headings"]).split(";"):
                h = h.strip()
                if not h:
                    continue
                mid, _ = mm.lookup(h)
                if mid:
                    src = h
                    break
        study_term[r.record_id] = (mid, src)

    attrs = {}
    if os.path.exists(ATTRS):
        with open(ATTRS, encoding="utf8", errors="replace") as fh:
            for line in fh:
                try:
                    j = json.loads(line)
                except ValueError:
                    continue
                attrs[j.get("sample_accession", "")] = j.get("custom_attributes") or {}

    rows, per_study = [], []
    for r in m.itertuples():
        tier = r.final_tier
        fields = [f.strip() for f in str(r.disease_field).split(";") if f.strip()]
        accs = [a.strip() for a in str(r.accession_codes).split(";") if a.strip()]
        s_mid, s_src = study_term.get(r.record_id, (None, ""))
        seen = set()
        n_s = n_dis = n_ctrl = n_map = n_cust = n_unexp = 0
        for acc in accs:
            path = os.path.join(RUNS, f"{acc}.tsv")
            if not os.path.exists(path):
                continue
            with open(path, encoding="utf8", errors="replace") as fh:
                for rec in csv.DictReader(fh, delimiter="\t"):
                    s = rec.get("sample_accession", "")
                    if not s or s in seen:
                        continue
                    seen.add(s)
                    n_s += 1
                    cust = attrs.get(s)
                    if cust:
                        n_cust += 1
                    raw, used = first_value(rec, cust, fields + DISEASE_FALLBACK)
                    age, _ = first_value(rec, cust, AGE_FIELDS)
                    sex, _ = first_value(rec, cust, SEX_FIELDS)
                    mid, cc, how = resolve(raw, used, fields, mm, s_mid, s_src, tier)
                    if raw:
                        n_dis += 1
                    if cc == "control":
                        n_ctrl += 1
                    if mid:
                        n_map += 1
                    # a tier 3/4 sample that carries a real, mappable disease
                    # value contradicts its tier -- surface it, do not use it
                    unexpected = ("yes" if tier in ("3", "4")
                                  and how.startswith("sample_value_") else "")
                    if unexpected:
                        n_unexp += 1
                    d = mm.describe(mid) if mid else {}
                    rows.append(dict(
                        record_id=r.record_id, study_accession=acc,
                        sample_accession=s, tier=tier,
                        raw_disease_value=raw, disease_source_field=used,
                        case_control=cc, mondo_id=mid or "",
                        mondo_label=d.get("mondo_label", ""),
                        doid=d.get("doid", ""), mesh=d.get("mesh", ""),
                        host_age=age, host_sex=sex, mapping_method=how,
                        study_mesh_heading=s_src,
                        custom_attrs_available="yes" if cust else "no",
                        name_derived=r.name_derived_disease_flag,
                        in_analysis=r.in_analysis,
                        unexpected_annotation=unexpected,
                        run_accession=rec.get("run_accession", ""),
                        library_strategy=rec.get("library_strategy", ""),
                        instrument_platform=rec.get("instrument_platform", "")))
        if n_s:
            per_study.append(dict(
                record_id=r.record_id, tier=tier, in_analysis=r.in_analysis,
                n_samples=n_s,
                n_with_disease_value=n_dis, n_controls=n_ctrl,
                n_mondo_mapped=n_map, n_custom_attrs=n_cust,
                n_unexpected_annotation=n_unexp,
                study_mesh_heading=s_src, disease_field=r.disease_field))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    ps = pd.DataFrame(per_study)
    ps.to_csv(OUT_STUDY, sep="\t", index=False)
    t12_cols = [c for c in out.columns
                if c not in ("unexpected_annotation", "run_accession",
                             "library_strategy", "instrument_platform")]
    out[out.tier.isin(["1", "2"])][t12_cols].to_csv(OUT_T12, sep="\t", index=False)
    ps[ps.tier.isin(["1", "2"])].drop(columns=["n_unexpected_annotation"]) \
        .to_csv(OUT_T12_STUDY, sep="\t", index=False)

    n = len(out)
    print(f"wrote {os.path.basename(OUT)}  "
          f"({n:,} sample rows, {out.sample_accession.nunique():,} distinct "
          f"samples, {out.record_id.nunique():,} studies)")
    print(out.groupby("tier").agg(studies=("record_id", "nunique"),
                                  rows=("tier", "size"),
                                  with_value=("raw_disease_value",
                                              lambda s: int((s != "").sum())),
                                  controls=("case_control",
                                            lambda s: int((s == "control").sum())),
                                  mondo=("mondo_id",
                                         lambda s: int((s != "").sum()))).to_string())
    print(f"\nexplicit_absence (was silently a case): "
          f"{int((out.mapping_method == 'explicit_absence').sum()):,}")
    print(f"tier 3/4 samples with a real disease value (check the tier): "
          f"{int((out.unexpected_annotation == 'yes').sum()):,} "
          f"in {out[out.unexpected_annotation == 'yes'].record_id.nunique()} studies")
    print("\nmapping_method:")
    for k, v in out.mapping_method.value_counts().items():
        print(f"  {k:42s} {v:>7,}")


if __name__ == "__main__":
    main()
