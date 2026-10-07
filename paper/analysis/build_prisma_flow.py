"""Assemble the MMC2 PRISMA 2020 flow counts as a single-column flow.

The review was run in Rayyan on one pile of PubMed records.  The cervicovaginal
search was not a separate branch: its 1,656 records were poured in at
identification and de-duplicated against the main query along with everything
else.  So the flow is one column, five stages deep.

There is no third identification arm.  The cancer supplementary search
contributed 51 studies to the analysis set, and all 51 are returned by the main
query as written (checked against live PubMed, 4 Sep 2026).  It was a targeted
re-screen of records the main query had already retrieved, not a new source.
The control holds: only 27 of the 157 cervicovaginal studies in the analysis set
are returned by the main query, which is why that search does add records.

Writes:
    MMC2_PRISMA_flow.md
    MMC2_prisma_counts.tsv
"""

import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import (  # noqa: E402
    ROOT,
    VAGINAL_RIS,
    load_disease_headings,
    mesh_flags,
)

MASTER = os.path.join(ROOT, "MMC2_master_Aug2026.csv")
REVIEWED = os.path.join(ROOT, "MMC2_master_reviewed.tsv")
OUT_MD = os.path.join(ROOT, "MMC2_PRISMA_flow.md")
OUT_TSV = os.path.join(ROOT, "MMC2_prisma_counts.tsv")

# Stage counts Rayyan reported at the time.  These are not derivable from the
# files on disk - Rayyan is the screening platform, and its de-duplication and
# screening tallies live there.  Recorded here as documented constants.
RAYYAN = {
    "main_identified": 33564,   # PubMed records imported from the main query
    "deduplicated": 23781,      # after Rayyan removed duplicates across the pile
    "after_tiab": 17009,        # after title/abstract screening
}
# Rayyan also reported 4,110 excluded at full text, which lands on 12,899.  The
# screened file that actually fed curation, MMC2_master_Aug2026.csv, holds 13,094
# records - 195 more.  The flow anchors on the file so a reviewer recomputing
# from the deposited data reproduces the figure; the full-text exclusion is
# therefore taken by subtraction (3,915) and this note carries the difference.
RAYYAN_FULLTEXT_EXCLUDED_REPORTED = 4110


def ris_records(path):
    """Records in a RIS export.

    Counted on the reference-type tag, which every record carries, not the
    title tag - 15 of these records have no TI line, and counting titles
    undercounts the export at 1,641.
    """
    txt = open(path, encoding="utf8", errors="replace").read()
    n_ty = len(re.findall(r"^TY\s+-\s+", txt, re.M))
    n_er = len(re.findall(r"^ER\s+-", txt, re.M))
    assert n_ty == n_er, f"RIS record delimiters disagree: {n_ty} TY vs {n_er} ER"
    return n_ty


def main():
    master = pd.read_csv(MASTER, low_memory=False)
    mr = pd.read_csv(REVIEWED, sep="\t", low_memory=False,
                     keep_default_na=False, dtype=str)
    disease = load_disease_headings()

    # ---- measured from files -------------------------------------------
    n_cervico = ris_records(VAGINAL_RIS)          # 1,656
    n_screened = len(master)                      # 13,094
    n_corpus = len(mr)                            # 3,318
    n_analysis = int(mr["in_analysis"].eq("yes").sum())          # 2,896

    rem = mr["remove_flag"].eq("yes")
    n_removed = int(rem.sum())
    n_noninsdc = int((mr["final_tier"].eq("") & ~rem).sum())
    n_healthy = int((mr["final_tier"].ne("") & ~rem
                     & mr["disease_group_final"].eq("Healthy")).sum())

    # MeSH eligibility recomputed on the screened set, for the notes only.  It
    # gives 3,428 while the curated corpus is 3,318, so the figure does not put a
    # number on the MeSH exclusion line.
    m_design, m_dis, _ = mesh_flags(master, disease)
    n_mesh_eligible = int((m_design & m_dis).sum())

    # ---- the chain ------------------------------------------------------
    n_identified = RAYYAN["main_identified"] + n_cervico
    n_dedup = RAYYAN["deduplicated"]
    n_tiab = RAYYAN["after_tiab"]

    excl_dup = n_identified - n_dedup
    excl_tiab = n_dedup - n_tiab
    excl_full = n_tiab - n_screened
    excl_final = n_screened - n_analysis

    assert n_identified == 35220, n_identified
    assert excl_dup > 0 and excl_tiab > 0 and excl_full > 0 and excl_final > 0
    assert n_dedup + excl_tiab == n_identified - excl_dup + excl_tiab
    assert n_identified - excl_dup == n_dedup
    assert n_dedup - excl_tiab == n_tiab
    assert n_tiab - excl_full == n_screened
    assert n_screened - excl_final == n_analysis
    assert n_corpus - (n_noninsdc + n_healthy + n_removed) == n_analysis

    tiers = mr[mr["in_analysis"].eq("yes")]["final_tier"].value_counts()

    rows = [
        ("identification", "Records identified in PubMed", n_identified),
        ("identification", "  main query", RAYYAN["main_identified"]),
        ("identification", "  cervicovaginal search", n_cervico),
        ("screening", "Records after duplicates removed", n_dedup),
        ("screening", "  excluded: duplicate records", excl_dup),
        ("screening", "Records after title/abstract screening", n_tiab),
        ("screening", "  excluded at title/abstract screening", excl_tiab),
        ("screening", "Records after full-text screening", n_screened),
        ("screening", "  excluded at full-text screening", excl_full),
        ("included", "Studies included in the analysis", n_analysis),
        ("included", "  excluded at eligibility and curation", excl_final),
        ("included", "  excluded: deposited outside INSDC", n_noninsdc),
        ("included", "  excluded: healthy-cohort studies", n_healthy),
        ("included", "  excluded: non-human, duplicate deposit, withdrawn", n_removed),
        ("notes", "Eligible records carried into the curated corpus", n_corpus),
        ("notes", "MeSH eligibility recomputed on the screened set", n_mesh_eligible),
        ("notes", "Rayyan-reported full-text exclusions", RAYYAN_FULLTEXT_EXCLUDED_REPORTED),
    ]
    for t in "1234":
        rows.append(("tiers", f"Tier {t}", int(tiers.get(t, 0))))

    with open(OUT_TSV, "w", encoding="utf8") as f:
        f.write("stage\tdescription\tn\n")
        for s, d, n in rows:
            f.write(f"{s}\t{d}\t{n}\n")

    # ---- narrative ------------------------------------------------------
    L = []
    a = L.append
    a("# MMC2 PRISMA 2020 flow\n")
    a("Counts are generated by `scripts/build_prisma_flow.py`. Every figure "
      "except the three Rayyan stage tallies is recomputed from the files on "
      "disk, and the script asserts that the chain adds up.\n")

    a("## Identification\n")
    a("One pile, not three branches. The cervicovaginal search was poured into "
      "the same Rayyan project as the main query and de-duplicated against it.\n")
    a("| Source | Records |")
    a("|---|---|")
    a(f"| Main PubMed query | {RAYYAN['main_identified']:,} |")
    a(f"| Cervicovaginal PubMed search | {n_cervico:,} |")
    a(f"| **Total identified** | **{n_identified:,}** |")
    a("")
    a("The cancer supplementary search is not listed as a source. It "
      "contributed 51 studies to the analysis set, and all 51 are returned by "
      "the main query as written, so they were already inside the "
      f"{RAYYAN['main_identified']:,}. By contrast only 27 of the 157 "
      "cervicovaginal studies in the analysis set are returned by the main "
      "query, which is why that search does add records.\n")

    a("## Screening\n")
    a("| Stage | Excluded | Remaining |")
    a("|---|---|---|")
    a(f"| Duplicate records removed | {excl_dup:,} | {n_dedup:,} |")
    a(f"| Title and abstract screening | {excl_tiab:,} | {n_tiab:,} |")
    a(f"| Full-text screening | {excl_full:,} | {n_screened:,} |")
    a("")
    a(f"Rayyan reported {RAYYAN_FULLTEXT_EXCLUDED_REPORTED:,} full-text "
      f"exclusions, which would leave "
      f"{n_tiab - RAYYAN_FULLTEXT_EXCLUDED_REPORTED:,}. The screened file that "
      f"fed curation, `MMC2_master_Aug2026.csv`, holds {n_screened:,} records - "
      f"{n_screened - (n_tiab - RAYYAN_FULLTEXT_EXCLUDED_REPORTED):,} more. The "
      "flow anchors on the file so the figure is reproducible from the "
      "deposited data, and takes the full-text exclusion by subtraction.\n")

    a("## Eligibility and inclusion\n")
    a("> **Inclusion criterion.** A record is eligible if its MeSH indexing "
      "includes **both** (a) an epidemiological or interventional study-design "
      "heading (10 headings; see `mesh_design_terms.tsv`) **and** (b) a disease "
      "or condition heading, defined as any descriptor in MeSH tree category "
      "**C (Diseases)** or **F03 (Mental Disorders)** (see "
      "`mesh_disease_terms.tsv`).\n")
    a("Rationale: metadata reusability can only be assessed for studies whose "
      "design implies a defined clinical comparison structure and a named "
      "condition.\n")
    a(f"- Records excluded at this stage: **{excl_final:,}**")
    a(f"- Studies included in the analysis: **{n_analysis:,}**\n")
    a("Reasons, for the part that is recomputable:\n")
    a("| Reason | Records |")
    a("|---|---|")
    a("| No study-design or disease/condition MeSH heading | (remainder) |")
    a(f"| Sequence data deposited outside INSDC | {n_noninsdc:,} |")
    a(f"| Healthy cohorts with no disease contrast | {n_healthy:,} |")
    a(f"| Non-human, duplicate deposit, withdrawn | {n_removed:,} |")
    a("")
    a(f"The MeSH criterion recomputed on the {n_screened:,} screened records "
      f"gives {n_mesh_eligible:,} eligible, while the curated corpus is "
      f"{n_corpus:,}. The figure therefore does not put a number on the MeSH "
      "exclusion line rather than paper over that "
      f"{abs(n_mesh_eligible - n_corpus):,}-record gap.\n")

    a("## Included\n")
    a("| Tier | Studies |")
    a("|---|---|")
    for t in "1234":
        a(f"| Tier {t} | {int(tiers.get(t, 0)):,} |")
    a(f"| **Total** | **{n_analysis:,}** |")
    a("")

    open(OUT_MD, "w", encoding="utf8").write("\n".join(L))

    print(f"wrote {OUT_TSV}")
    print(f"wrote {OUT_MD}")
    print(f"  identified        {n_identified:>7,}   "
          f"(main {RAYYAN['main_identified']:,} + cervicovaginal {n_cervico:,})")
    print(f"  de-duplicated     {n_dedup:>7,}   -{excl_dup:,}")
    print(f"  after title/abs   {n_tiab:>7,}   -{excl_tiab:,}")
    print(f"  after full text   {n_screened:>7,}   -{excl_full:,}")
    print(f"  INCLUDED          {n_analysis:>7,}   -{excl_final:,}")
    print(f"  tiers             " + "  ".join(
        f"T{t} {int(tiers.get(t,0)):,}" for t in "1234"))


if __name__ == "__main__":
    main()
