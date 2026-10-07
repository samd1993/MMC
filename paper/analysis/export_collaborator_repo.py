#!/usr/bin/env python3
"""Export the data pack for the sample-labelling repo.

The collaborator's job is: pull the reads for these studies, and label every
sample healthy vs disease from the metadata field the reviewers identified.
Everything he needs to do that -- and nothing about the manuscript.
"""
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.join(os.path.dirname(ROOT), "mmc-sample-labels")
p = lambda f: os.path.join(ROOT, f)
out = lambda f: os.path.join(REPO, "data", f)

m = pd.read_csv(p("MMC2_master_reviewed.tsv"), sep="\t", low_memory=False,
                keep_default_na=False, dtype=str)
sm = pd.read_csv(p("MMC2_study_metadata.tsv"), sep="\t", low_memory=False)
at = pd.read_csv(p("MMC2_study_attributes.tsv"), sep="\t", low_memory=False)

m = m[m["remove_flag"].ne("yes") & m["final_tier"].ne("")]
s = (m.merge(sm[["record_id", "doi", "n_samples", "accession_codes",
                 "n_accessions_resolved", "metadata_source"]],
             on="record_id", how="left", suffixes=("", "_p"))
      .merge(at[["record_id", "body_site_group", "sequencing_type", "country",
                 "continent"]], on="record_id", how="left", suffixes=("", "_a")))

studies = pd.DataFrame({
    "record_id": s["record_id"],
    "tier": s["final_tier"],
    "year": s["year"],
    "journal": s["journal"],
    "title": s["title"],
    "doi": s["doi"],
    "accession_codes": s["accession_codes_p"],
    "n_accessions_resolved": s["n_accessions_resolved"],
    "n_samples_in_repository": s["n_samples_p"],
    # the field the reviewers found the per-sample disease in -- the single
    # most useful column here, and the reason this repo exists
    "disease_field": s["disease_metadata_field"],
    "disease_field_source": s["disease_metadata_field_source"],
    "disease_field_is_free_text": s["name_derived_field"],
    "disease_group": s["disease_group_final"],
    "healthy_cohort": s["healthy_flag"],
    "body_site": s["body_site_group_a"],
    "sequencing_type": s["sequencing_type_a"],
    "country": s["country"],
    "tier_source": s["final_tier_source"],
    "needs_curator_check": s["adjudication_flag"],
    # Healthy-cohort studies are excluded from the PAPER's analyses but are kept
    # here on purpose: they are the cleanest source of controls for a labelling
    # task, which is the opposite of a reason to drop them.
    "in_paper_analysis": s["in_analysis"],
}).sort_values(["tier", "record_id"])
studies.to_csv(out("studies.tsv"), sep="\t", index=False)

# One per-sample table for every tier. Tier is a column, not a filename.
samples = pd.read_csv(p("MMC2_sample_master.tsv"), sep="\t", low_memory=False,
                      keep_default_na=False, dtype=str)
samples = samples[samples["record_id"].isin(set(studies["record_id"]))]
samples.to_csv(out("samples.tsv.gz"), sep="\t", index=False, compression="gzip")

# One accession list for every tier, with the tier alongside so a caller can
# filter instead of picking a different file.
INSDC = {"PRJ", "ERP", "SRP", "DRP", "ERS", "SRS", "DRS", "ERR", "SRR", "DRR",
         "SRA", "DRA", "CRA", "CNP"}
rows = []
for r in studies.itertuples(index=False):
    for a in [x.strip() for x in str(r.accession_codes).split(";") if x.strip()]:
        # DOIs and GEO series need a different route and would 404 the ENA API
        if a[:3].upper() in INSDC:
            rows.append({"accession": a, "record_id": r.record_id, "tier": r.tier})
acc = (pd.DataFrame(rows).drop_duplicates()
       .sort_values(["tier", "accession"]))
acc.to_csv(out("accessions.tsv"), sep="\t", index=False)
print("accessions.tsv:", len(acc), "rows,",
      acc["accession"].nunique(), "distinct accessions")
print(acc.groupby("tier")["accession"].nunique().to_string())

dfm = studies[studies["disease_field"].ne("")][
    ["record_id", "tier", "accession_codes", "disease_field",
     "disease_field_source", "disease_field_is_free_text", "disease_group",
     "healthy_cohort"]]
dfm.to_csv(out("disease_field_map.tsv"), sep="\t", index=False)

print(f"studies.tsv: {len(studies)} rows")
print(studies["tier"].value_counts().sort_index().to_string())
print(f"samples.tsv.gz: {len(samples):,} rows, "
      f"{samples['sample_accession'].nunique():,} distinct samples")
print(f"disease_field_map.tsv: {len(dfm)} rows")
