"""Harmonized per-sample metadata for every Tier 1 and Tier 2 study.

The community resource the paper promises: one row per deposited sample, with
disease mapped onto MONDO (carrying DOID and MeSH cross-references), plus host
age and sex, drawn straight from the repository record.

Where each value comes from:

  portal fields   host_disease, disease, host_phenotype, host_status, age, sex,
                  host_sex, sample_title, sample_alias ... complete for every
                  sample, from the uncapped runs cache
  custom keys     author-named columns, from sample XML. PARTIAL at present: the
                  ENA XML API returned HTTP 500 corpus-wide during the backfill,
                  so 122,883 of 171,144 samples have no custom attributes yet.
                  `custom_attrs_available` marks which. Rerunning after the
                  outage fills them; the fetch is resumable.

Disease resolution, in order:
  1. the sample's own value in the study's recorded disease field, mapped to
     MONDO directly;
  2. failing that, if the value reads as a control, the sample is marked control
     and no disease term is assigned;
  3. failing that -- the common case for coded values like "CD_8" -- the study's
     disease comes from its MeSH disease headings, and the sample keeps the raw
     code so a reader can see what was decoded.

That third route is why `mapping_method` matters: a MONDO id assigned from a
study-level MeSH heading is a weaker claim than one read off the sample, and the
column says which happened.

Writes MMC2_harmonized_T1T2_samples.tsv and a per-study coverage summary.
"""
import csv
import json
import os
import re
import sys
from collections import defaultdict

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402
from mondo_map import NULLISH, MondoMapper  # noqa: E402

csv.field_size_limit(10 ** 9)
MASTER = os.path.join(ROOT, "MMC2_master_reviewed.tsv")
ANALYSIS = os.path.join(ROOT, "MMC2_analysis_set.tsv")
RUNS = os.path.join(ROOT, "ena_cache", "runs")
ATTRS = os.path.join(ROOT, "ena_cache", "sample_attrs.jsonl")
INDEX = os.path.join(ROOT, "ontology", "mondo_index.json")
OUT = os.path.join(ROOT, "MMC2_harmonized_T1T2_samples.tsv")
OUT_STUDY = os.path.join(ROOT, "MMC2_harmonized_T1T2_by_study.tsv")

# Cohort questionnaires (American Gut / Microsetta and similar) inverted the
# usual layout: the FIELD NAME carries the condition ("skin_condition",
# "ibd_diagnosis") and the VALUE carries only whether the subject reports it.
# 27,376 samples in this corpus answer "I do not have this condition", which is a
# control statement, not an unmappable disease name.
SURVEY_NEG = re.compile(r"i do not have this condition|^no$|^never$|"
                        r"^not diagnosed", re.I)
SURVEY_POS = re.compile(r"diagnosed by a medical professional|self[- ]diagnosed|"
                        r"^diagnosed\b|^yes\b", re.I)
SURVEY_UNSURE = re.compile(r"^unspecified$|^unsure$|^don'?t know$|^na$", re.I)

AGE_FIELDS = ["age", "host_age", "host age", "age_years", "host_age_years"]
SEX_FIELDS = ["host_sex", "sex", "gender", "host sex", "submitted_host_sex"]
DISEASE_FALLBACK = ["host_disease", "disease", "host_phenotype", "host_status",
                    "host_health_state", "host_disease_stat"]


def first_value(rec, cust, names):
    for n in names:
        v = rec.get(n)
        if v and str(v).strip().lower() not in NULLISH:
            return str(v).strip(), n
        for k, vv in (cust or {}).items():
            if k.strip().lower() == n.lower() and vv and str(vv).strip().lower() not in NULLISH:
                return str(vv).strip(), k
    return "", ""


def main():
    mm = MondoMapper(INDEX)
    m = pd.read_csv(MASTER, sep="\t", low_memory=False, keep_default_na=False, dtype=str)
    t12 = m[m.final_tier.isin(["1", "2"])].copy()
    an = pd.read_csv(ANALYSIS, sep="\t", low_memory=False,
                     keep_default_na=False, dtype=str).set_index("record_id")

    # study-level disease from MeSH headings, used for coded samples
    study_term = {}
    for r in t12.itertuples():
        mid, how, src = None, None, ""
        if r.record_id in an.index:
            for h in str(an.loc[r.record_id, "mesh_disease_headings"]).split(";"):
                h = h.strip()
                if not h:
                    continue
                mid, how = mm.lookup(h)
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
    for r in t12.itertuples():
        fields = [f.strip() for f in r.disease_field.split(";") if f.strip()]
        accs = [a.strip() for a in r.accession_codes.split(";") if a.strip()]
        s_mid, s_src = study_term.get(r.record_id, (None, ""))
        seen = set()
        n_s = n_dis = n_ctrl = n_map = n_cust = 0
        for acc in accs:
            p = os.path.join(RUNS, f"{acc}.tsv")
            if not os.path.exists(p):
                continue
            with open(p, encoding="utf8", errors="replace") as fh:
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

                    mid, how = (None, None)
                    cc = ""
                    if raw and (SURVEY_NEG.search(raw) or SURVEY_POS.search(raw)
                                or SURVEY_UNSURE.search(raw)):
                        # questionnaire layout: condition is the field name
                        n_dis += 1
                        if SURVEY_NEG.search(raw):
                            cc, how = "control", "survey_response_negative"
                            n_ctrl += 1
                        elif SURVEY_UNSURE.search(raw):
                            cc, how = "", "survey_response_unspecified"
                        else:
                            mid, hw = mm.lookup(used)      # map the FIELD NAME
                            if not mid and s_mid:
                                mid, hw = s_mid, "study_mesh"
                            cc = "case"
                            how = (f"survey_response_positive_field_{hw}" if mid
                                   else "survey_response_positive_unmapped")
                    elif raw:
                        n_dis += 1
                        if mm.is_control(raw):
                            cc, how = "control", "control_value"
                            n_ctrl += 1
                        else:
                            mid, hw = mm.lookup(raw)
                            if mid:
                                cc, how = "case", f"sample_value_{hw}"
                            elif s_mid:
                                mid, cc, how = s_mid, "case", "study_mesh_heading"
                            else:
                                how = "unmapped"
                    elif s_mid:
                        mid, cc, how = s_mid, "", "study_mesh_heading_no_sample_value"
                    else:
                        how = "no_disease_value"
                    if mid:
                        n_map += 1
                    d = mm.describe(mid) if mid else {}
                    rows.append(dict(
                        record_id=r.record_id, study_accession=acc,
                        sample_accession=s, tier=r.final_tier,
                        raw_disease_value=raw, disease_source_field=used,
                        case_control=cc,
                        mondo_id=mid or "", mondo_label=d.get("mondo_label", ""),
                        doid=d.get("doid", ""), mesh=d.get("mesh", ""),
                        host_age=age, host_sex=sex,
                        mapping_method=how or "",
                        study_mesh_heading=s_src,
                        custom_attrs_available="yes" if cust else "no",
                        name_derived=r.name_derived_disease_flag))
        per_study.append(dict(record_id=r.record_id, tier=r.final_tier,
                              n_samples=n_s, n_with_disease_value=n_dis,
                              n_controls=n_ctrl, n_mondo_mapped=n_map,
                              n_custom_attrs=n_cust,
                              study_mesh_heading=s_src,
                              disease_field=r.disease_field))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    ps = pd.DataFrame(per_study)
    ps.to_csv(OUT_STUDY, sep="\t", index=False)

    n = len(out)
    print(f"wrote {OUT}  ({n:,} samples from {out.record_id.nunique()} studies)")
    print(f"wrote {OUT_STUDY}\n")
    print(f"  samples with a disease value in the record : "
          f"{int((out.raw_disease_value != '').sum()):,} ({100*(out.raw_disease_value!='').mean():.1f}%)")
    print(f"  mapped to a MONDO term                     : "
          f"{int((out.mondo_id != '').sum()):,} ({100*(out.mondo_id!='').mean():.1f}%)")
    print(f"  identified as controls                     : "
          f"{int((out.case_control == 'control').sum()):,}")
    print(f"  custom attributes available                : "
          f"{int((out.custom_attrs_available == 'yes').sum()):,} "
          f"({100*(out.custom_attrs_available=='yes').mean():.1f}%)  <- ENA XML outage")
    print("\n  mapping_method:")
    for k, v in out.mapping_method.value_counts().items():
        print(f"    {k:38s} {v:>7,}")
    print(f"\n  distinct MONDO terms: {out[out.mondo_id!=''].mondo_id.nunique()}")
    print("  most common:")
    for k, v in out[out.mondo_id != ""].mondo_label.value_counts().head(8).items():
        print(f"    {k[:44]:46s} {v:>7,}")


if __name__ == "__main__":
    main()
