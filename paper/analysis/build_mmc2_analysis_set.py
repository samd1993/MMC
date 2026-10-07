"""Restrict the eligible set to the MMC2 phase-2 analysis set.

Applies the three phase-2 scope changes:
  1. scraped accessions only (`accession_scraped`; the union is not used)
  2. publication years 2012-2025, matching MMC1's window
  3. classify each record by whether its accessions are retrievable from INSDC

Writes MMC2_analysis_set.tsv.
"""

import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT, acc_set  # noqa: E402

ELIGIBLE = os.path.join(ROOT, "MMC2_eligible_v1.csv")
OUT = os.path.join(ROOT, "MMC2_analysis_set.tsv")

YEAR_MIN, YEAR_MAX = 2012, 2025

# Accessions the ENA portal / browser API can resolve. Project, study, sample,
# experiment and run accessions across the three INSDC nodes (NCBI/EBI/DDBJ).
INSDC_RX = re.compile(r"^(PRJ(NA|EB|DB)\d+|[SED]R[PRSXA]\d+)$")

# Everything else, labelled so the report can name the repository rather than
# lumping it into "no metadata".
NON_INSDC_REPOS = [
    (re.compile(r"^(CRA|PRJCA|HRA|OMIX)\d+", re.I), "NGDC (China National Center for Bioinformation)"),
    (re.compile(r"^CNP\d+", re.I), "CNGB Nucleotide Sequence Archive"),
    (re.compile(r"^GSE\d+", re.I), "NCBI GEO"),
    (re.compile(r"^(PHS|phs)\d+", re.I), "dbGaP (controlled access)"),
    (re.compile(r"^EGA[DS]\d+", re.I), "EGA (controlled access)"),
    (re.compile(r"^ZENODO", re.I), "Zenodo"),
    (re.compile(r"^OEP\d+", re.I), "OMIX/NODE"),
    (re.compile(r"^E-[A-Z]+-\d+", re.I), "ArrayExpress/BioStudies"),
    (re.compile(r"^SUB\d+", re.I), "NCBI submission stub (not a public accession)"),
    (re.compile(r"^SAM[NED]", re.I), "BioSample (sample-level, no study context)"),
    (re.compile(r"^\d*\.?\d+/", re.I), "DOI-style identifier"),
]


def repo_of(code):
    for rx, name in NON_INSDC_REPOS:
        if rx.match(code):
            return name
    return "other / unrecognized"


def main():
    df = pd.read_csv(ELIGIBLE, low_memory=False)
    n_in = len(df)

    year = pd.to_numeric(df["year"], errors="coerce")
    df = df[year.between(YEAR_MIN, YEAR_MAX)].copy()
    df["year"] = year[year.between(YEAR_MIN, YEAR_MAX)].astype(int)
    print(f"eligible {n_in:,} -> {len(df):,} in {YEAR_MIN}-{YEAR_MAX}")

    # --- corrections from the DOI audit -----------------------------------
    # 179 records carry a DOI that resolves to a different paper (57.6% of the
    # cervicovaginal arm, whose RIS import took DOIs from reference lists). The
    # original accession scrape keyed on DOI, so those accessions can belong to
    # the wrong article -- one "vaginal" record held 43 codes from the SIAMCAT
    # toolbox paper, while the real article states its data are not public.
    # rescrape_accessions.py re-derives them from the correct PMID.
    resc_path = os.path.join(ROOT, "MMC2_rescraped_accessions.tsv")
    n_fixed = 0
    if os.path.exists(resc_path):
        rs = pd.read_csv(resc_path, sep="\t", low_memory=False)
        fix_acc = dict(zip(rs["record_id"], rs["rescraped_accessions"].fillna("")))
        fix_doi = {k: v for k, v in zip(rs["record_id"], rs["correct_doi"].fillna(""))
                   if isinstance(v, str) and v.strip()}
        changed = set(rs.loc[~rs["unchanged"], "record_id"])
        df["doi_corrected"] = df["record_id"].isin(fix_doi)
        df["doi_as_recorded"] = df["doi"]
        df["doi"] = [fix_doi.get(rid, d) for rid, d in zip(df["record_id"], df["doi"])]
        df["accession_scraped"] = [
            fix_acc[rid] if rid in fix_acc else a
            for rid, a in zip(df["record_id"], df["accession_scraped"])
        ]
        df["accession_rescraped"] = df["record_id"].isin(changed)
        n_fixed = len(fix_acc)
        print(f"applied corrections to {n_fixed:,} records "
              f"({len(changed):,} with changed accessions, {len(fix_doi):,} with a corrected DOI)")
    else:
        df["doi_corrected"] = False
        df["doi_as_recorded"] = df["doi"]
        df["accession_rescraped"] = False
        print("no corrections file -- run verify_dois.py --all then rescrape_accessions.py")

    # --- manual curation --------------------------------------------------
    # Human decisions on the reuse/accession flags: `remove` drops a reused
    # meta-analysis outright, `trim` replaces the accession list with just the
    # datasets the authors generated, `keep` confirms the record as-is.
    cur_path = os.path.join(ROOT, "MMC2_curation.tsv")
    df["curation"] = ""
    df["curation_note"] = ""
    if os.path.exists(cur_path):
        cur = pd.read_csv(cur_path, sep="\t", low_memory=False).fillna("")
        dec = dict(zip(cur["record_id"], cur["decision"]))
        keep_acc = dict(zip(cur["record_id"], cur["keep_accessions"]))
        note = dict(zip(cur["record_id"], cur["note"]))
        df["curation"] = df["record_id"].map(dec).fillna("")
        df["curation_note"] = df["record_id"].map(note).fillna("")
        df["accession_scraped"] = [
            keep_acc[rid] if dec.get(rid) == "trim" else a
            for rid, a in zip(df["record_id"], df["accession_scraped"])
        ]
        n_rm = int((df["curation"] == "remove").sum())
        n_tr = int((df["curation"] == "trim").sum())
        df = df[df["curation"] != "remove"].copy()
        print(f"curation: removed {n_rm} reused meta-analyses, "
              f"trimmed {n_tr} accession lists, {len(df):,} papers remain")

    # --- collapse duplicate records ---------------------------------------
    # Phase 1 flagged duplicate DOIs but deliberately left them in place. They
    # must go before any rate is computed: the same paper counted twice inflates
    # every tier and presence statistic, and shows up here as an accession
    # "shared" between two records that are really one study.
    # Require BOTH the same normalized DOI and the same normalized title -- a
    # title alone collides (two distinct papers share the opening 70 characters
    # of "A Predictive Model Based on the Gut Microbiota Improves...").
    def _nt(s):
        if pd.isna(s):
            return None
        return re.sub(r"[^a-z0-9]", "", str(s).lower())[:70] or None

    df["_dk"] = [
        (str(d).strip().lower(), _nt(t)) if pd.notna(d) and _nt(t) else None
        for d, t in zip(df["doi"], df["title"])
    ]
    before = len(df)
    # keep the richest record in each group: most accessions, then most fields
    df["_rank"] = (df["accession_scraped"].map(lambda a: len(acc_set(a))) * 100
                   + df.notna().sum(axis=1))
    dup_mask = df["_dk"].notna() & df["_dk"].duplicated(keep=False)
    n_groups = df.loc[dup_mask, "_dk"].nunique()
    df = (df.sort_values("_rank", ascending=False)
            .drop_duplicates(subset="_dk", keep="first")
            .sort_index())
    df = df.drop(columns=["_dk", "_rank"])
    if before != len(df):
        print(f"deduplication: collapsed {n_groups} duplicate record groups, "
              f"removed {before - len(df)} rows -> {len(df):,} papers")

    # --- manual additions -------------------------------------------------
    # A cohort can be cited by papers in the corpus while its own originating
    # publication was never retrieved by any search arm. American Gut is the
    # case in point: two of the largest studies here are built on it, but the
    # 2018 mSystems paper is absent. Adding it once attributes the deposit to
    # its originator instead of to papers that reused it.
    add_path = os.path.join(ROOT, "MMC2_manual_additions.tsv")
    if os.path.exists(add_path):
        add = pd.read_csv(add_path, sep="\t", low_memory=False)
        for c in df.columns:
            if c not in add.columns:
                add[c] = pd.NA
        add = add[df.columns]
        df = pd.concat([df, add], ignore_index=True)
        print(f"manual additions: +{len(add)} record(s) -> {len(df):,} papers")

    # --- canonicalize accession aliases -----------------------------------
    # PRJEB11419 and ERP012803 are the same ENA study under primary and
    # secondary accessions. Left as distinct strings they hide shared cohorts:
    # 31 such pairs both appear in this corpus, American Gut among them.
    alias_path = os.path.join(ROOT, "accession_aliases.tsv")
    alias = {}
    if os.path.exists(alias_path):
        al = pd.read_csv(alias_path, sep="\t")
        alias = dict(zip(al["secondary"].str.upper(), al["primary"].str.upper()))

    def canon(cell):
        return "; ".join(sorted({alias.get(c, c) for c in acc_set(cell)}))

    if alias:
        before_codes = {c for cell in df["accession_scraped"] for c in acc_set(cell)}
        df["accession_scraped"] = df["accession_scraped"].map(canon)
        after_codes = {c for cell in df["accession_scraped"] for c in acc_set(cell)}
        print(f"alias canonicalization: {len(before_codes):,} -> {len(after_codes):,} "
              f"distinct codes ({len(before_codes) - len(after_codes)} secondaries folded)")

    # --- scraped accessions only ------------------------------------------
    codes = df["accession_scraped"].map(acc_set)
    df["accession_codes"] = codes.map(lambda s: "; ".join(sorted(s)))
    df["n_accession_codes"] = codes.map(len)

    insdc = codes.map(lambda s: sorted(c for c in s if INSDC_RX.match(c)))
    other = codes.map(lambda s: sorted(c for c in s if not INSDC_RX.match(c)))
    df["insdc_codes"] = insdc.map("; ".join)
    df["non_insdc_codes"] = other.map("; ".join)

    df["metadata_source"] = [
        "insdc" if i else "non_insdc" if o else "none"
        for i, o in zip(insdc, other)
    ]
    df["non_insdc_repos"] = other.map(
        lambda cs: "; ".join(sorted({repo_of(c) for c in cs}))
    )

    # Tier 4 by definition: no accession at all. Everything else is decided
    # downstream once the ENA harvest lands.
    df["tier_prelim"] = ["4" if m == "none" else "" for m in df["metadata_source"]]
    df["tier_assessable"] = df["metadata_source"].isin(["insdc", "none"])

    # --- collapse companion papers onto one row per cohort ----------------
    # Papers sharing an IDENTICAL accession set are publications off the same
    # deposit -- three HCHS/SOL papers on ERP117287, two SCAPIS papers on
    # PRJEB51353. Counting each separately multiplies one deposition event
    # across the tier statistics. Keep the earliest publication as the
    # cohort's representative and record the rest.
    df["cohort_key"] = ["; ".join(sorted(s)) if s else "" for s in codes]
    has = df["cohort_key"] != ""
    grp = df[has].groupby("cohort_key")["record_id"].apply(list).to_dict()
    multi = {k: v for k, v in grp.items() if len(v) > 1}
    keep_ids, collapsed_into = set(), {}
    for key, recs in multi.items():
        sub = df[df["record_id"].isin(recs)].sort_values(["year", "record_id"])
        rep = sub.iloc[0]["record_id"]
        keep_ids.add(rep)
        collapsed_into[rep] = [r for r in sub["record_id"] if r != rep]
    drop = {r for recs in multi.values() for r in recs} - keep_ids
    df["collapsed_companions"] = df["record_id"].map(
        lambda r: "; ".join(collapsed_into.get(r, []))
    ).fillna("")
    df["n_companion_papers"] = df["collapsed_companions"].map(
        lambda s: len(s.split("; ")) if s else 0
    )
    before_c = len(df)
    df = df[~df["record_id"].isin(drop)].copy()
    codes = df["accession_scraped"].map(acc_set)
    if drop:
        print(f"cohort collapse: {len(multi)} shared-accession groups -> "
              f"removed {before_c - len(df)} companion papers, {len(df):,} remain")

    df.to_csv(OUT, sep="\t", index=False)

    ms = df["metadata_source"].value_counts()
    all_insdc = sorted({c for cs in insdc for c in cs})
    all_other = sorted({c for cs in other for c in cs})

    print(f"\nwrote {OUT}  ({len(df):,} rows)")
    print("\nmetadata_source:")
    for k in ("insdc", "non_insdc", "none"):
        n = int(ms.get(k, 0))
        print(f"  {k:10s} {n:6,d} ({100 * n / len(df):5.1f}%)")
    print(f"\ndistinct codes: {len(all_insdc):,} INSDC / {len(all_other):,} non-INSDC")
    print("\nnon-INSDC repositories (excluded from the tier denominator):")
    rc = {}
    for c in all_other:
        rc[repo_of(c)] = rc.get(repo_of(c), 0) + 1
    for k, v in sorted(rc.items(), key=lambda x: -x[1]):
        print(f"  {v:4d}  {k}")

    # the accession list the harvest will consume
    with open(os.path.join(ROOT, "MMC2_insdc_accessions.txt"), "w") as fh:
        fh.write("\n".join(all_insdc) + "\n")
    print(f"\nwrote MMC2_insdc_accessions.txt ({len(all_insdc):,} accessions to harvest)")


if __name__ == "__main__":
    main()
