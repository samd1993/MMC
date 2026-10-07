"""Assign per-study metadata flags and reusability tiers for MMC2.

Tier definitions are MMC1's, unchanged
(`MMC/scripts/enaharm/references/tier_definitions.md`):

  Tier 1  per-sample disease + age + sex
  Tier 2  per-sample disease
  Tier 3  valid accession, sample IDs differentiate samples, but no biology
  Tier 4  no accession, or an accession whose samples carry no differentiating IDs

Antibiotics and geography are reported alongside but are NOT tier inputs.

The drag-down rule is enforced upstream in survey_keys_mmc2.py: a category only
counts as present when its key takes more than one distinct value across the
study's samples. A single value repeated on every row is the submitter copying a
study-level label down, not per-sample evidence.

Writes MMC2_study_metadata.tsv and MMC2_tier_review.tsv.
"""

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402

ANALYSIS = os.path.join(ROOT, "MMC2_analysis_set.tsv")
COLMAP = os.path.join(ROOT, "MMC2_study_column_map.csv")
SURVEY = os.path.join(ROOT, "MMC2_key_survey.tsv")
ABXTEXT = os.path.join(ROOT, "MMC2_antibiotic_text.tsv")
OUT = os.path.join(ROOT, "MMC2_study_metadata.tsv")
REVIEW = os.path.join(ROOT, "MMC2_tier_review.tsv")

# Keys whose variation shows the submitter gave samples distinguishable names.
ID_KEYS = ["sample_alias", "sample_title", "sample_title_xml", "sample_description"]


def main():
    an = pd.read_csv(ANALYSIS, sep="\t", low_memory=False)
    cm = pd.read_csv(COLMAP, low_memory=False)
    sv = pd.read_csv(SURVEY, sep="\t", low_memory=False)
    ab = pd.read_csv(ABXTEXT, sep="\t", low_memory=False)

    # --- do sample IDs differentiate samples? -----------------------------
    idvar = (sv[sv["key"].isin(ID_KEYS)]
             .groupby("study")["varies"]
             .any()
             .rename("ids_differentiate"))
    cm = cm.merge(idvar, left_on="study", right_index=True, how="left")
    cm["ids_differentiate"] = cm["ids_differentiate"].fillna(False)

    # --- map studies (accessions) back to papers --------------------------
    # a paper can cite several accessions; it gets the best evidence among them
    rows = []
    for _, r in an.iterrows():
        accs = [a for a in str(r.get("insdc_codes") or "").split("; ") if a]
        sub = cm[cm["study"].isin(accs)]
        rec = {
            "record_id": r["record_id"],
            "doi": r["doi"],
            "year": r["year"],
            "journal": r["journal"],
            "title": str(r["title"])[:200],
            "metadata_source": r["metadata_source"],
            "accession_codes": r["accession_codes"],
            "n_insdc_accessions": len(accs),
            "n_accessions_resolved": len(sub),
            "n_samples": int(sub["n_samples"].sum()) if len(sub) else 0,
        }
        for cat in ("disease", "age", "sex", "antibiotic", "geography"):
            rec[f"{cat}_present"] = bool(sub[f"{cat}_present"].any()) if len(sub) else False
            rec[f"{cat}_key"] = "; ".join(sorted({k for k in sub[f"{cat}_key"].dropna() if k})) if len(sub) else ""
        # Tier 2 is "per-sample disease via a metadata field OR an informative
        # sample name", so disease uses the rescued evidence, not the key alone.
        rec["disease_from_names"] = bool(sub["disease_from_names"].any()) if len(sub) else False
        rec["disease_evidence"] = rec["disease_present"] or rec["disease_from_names"]
        rec["ids_differentiate"] = bool(sub["ids_differentiate"].any()) if len(sub) else False
        rows.append(rec)
    df = pd.DataFrame(rows)

    # --- tier assignment ---------------------------------------------------
    def tier(r):
        if r["metadata_source"] == "non_insdc":
            return ""                       # excluded from the tier denominator
        if r["metadata_source"] == "none":
            return "4"                      # no accession -> not reusable
        if r["n_accessions_resolved"] == 0:
            return "4"                      # accession cited but nothing resolved
        if r["disease_evidence"] and r["age_present"] and r["sex_present"]:
            return "1"
        if r["disease_evidence"]:
            return "2"
        if r["ids_differentiate"]:
            return "3"
        return "4"

    df["tier"] = df.apply(tier, axis=1)
    df["tier_auto"] = df["tier"]            # preserved; hand review edits `tier`
    df["tier_reviewed"] = ""
    TIER_LABEL = {
        "1": "Tier 1: Repository biological annotation",
        "2": "Tier 2: Informative sample names",
        "3": "Tier 3: Unique sample IDs, no biology",
        "4": "Tier 4: Not reusable",
        "": "Not assessable (non-INSDC repository)",
    }
    df["reusability"] = df["tier"].map(TIER_LABEL)

    # --- antibiotics: three-state labeling --------------------------------
    ab = ab.set_index("record_id")
    df = df.join(
        ab[["abx_involves_study", "abx_universal_candidate"]], on="record_id"
    )
    df["abx_involves_study"] = df["abx_involves_study"].fillna(False)
    df["abx_universal_candidate"] = df["abx_universal_candidate"].fillna(False)

    def abx_state(r):
        if not r["abx_involves_study"]:
            return "not_applicable"
        if r["antibiotic_present"]:
            return "yes_per_sample"
        if r["abx_universal_candidate"]:
            return "all_participants"       # candidate; confirmed by hand review
        return "no"

    df["abx_labeling"] = df.apply(abx_state, axis=1)

    df.to_csv(OUT, sep="\t", index=False)

    # --- hand-review worksheet: every T1/T2 + the all_participants set -----
    need = df[df["tier"].isin(["1", "2"]) | df["abx_labeling"].eq("all_participants")].copy()
    need["review_reason"] = [
        "tier 1/2 candidate" if t in ("1", "2") else "all_participants candidate"
        for t in need["tier"]
    ]
    need["reviewer_tier"] = ""
    need["reviewer_abx"] = ""
    need["reviewer_notes"] = ""
    cols = ["record_id", "doi", "year", "title", "review_reason", "tier_auto",
            "disease_evidence", "disease_from_names", "age_present", "sex_present", "disease_key",
            "age_key", "sex_key", "antibiotic_present", "antibiotic_key",
            "abx_labeling", "n_samples", "accession_codes",
            "reviewer_tier", "reviewer_abx", "reviewer_notes"]
    need[cols].to_csv(REVIEW, sep="\t", index=False)

    # --- report -------------------------------------------------------------
    print(f"wrote {OUT}   ({len(df):,} papers)")
    print(f"wrote {REVIEW}  ({len(need):,} rows needing hand review)")

    assessable = df[df["tier"] != ""]
    print(f"\n=== TIERS (denominator {len(assessable):,}; "
          f"{int((df['tier'] == '').sum()):,} non-INSDC excluded) ===")
    for t in ("1", "2", "3", "4"):
        n = int((df["tier"] == t).sum())
        print(f"  Tier {t}: {n:5,d} ({100 * n / len(assessable):5.1f}%)")

    print("\n=== METADATA PRESENCE (INSDC studies only) ===")
    ins = df[df["metadata_source"] == "insdc"]
    for cat, col in [("disease", "disease_evidence"), ("age", "age_present"),
                     ("sex", "sex_present"), ("antibiotic", "antibiotic_present"),
                     ("geography", "geography_present")]:
        n = int(ins[col].sum())
        print(f"  {cat:11s} {n:5,d} / {len(ins):,} ({100 * n / len(ins):5.1f}%)")
    key_only = int(ins["disease_present"].sum())
    resc = int((ins["disease_from_names"] & ~ins["disease_present"]).sum())
    print(f"    (disease: {key_only:,} via attribute key + {resc:,} rescued from "
          f"informative sample names)")

    print("\n=== ANTIBIOTIC LABELING ===")
    vc = df["abx_labeling"].value_counts()
    inv = int(df["abx_involves_study"].sum())
    for k in ("yes_per_sample", "all_participants", "no", "not_applicable"):
        n = int(vc.get(k, 0))
        pct = f" ({100 * n / inv:5.1f}% of abx studies)" if k != "not_applicable" and inv else ""
        print(f"  {k:17s} {n:5,d}{pct}")
    y = int(vc.get("yes_per_sample", 0))
    print(f"\n  -> of {inv:,} studies involving antibiotics, {y:,} "
          f"({100 * y / inv:.1f}%) label exposure per sample" if inv else "")


if __name__ == "__main__":
    main()
