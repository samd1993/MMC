"""Build the MMC2 master sheet: merge scraped accession codes into the Rayyan output.

Keeps every original column (including the old `Accession Code`) and adds a
merged column plus provenance flags. Writes:
    MMC2_master_Aug2026.csv
    MMC2_accession_merge_report.md
"""

import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import (  # noqa: E402
    CANCER_CSV,
    EXTRA_CSV,
    MOHAK,
    ORIGINAL_QUERY,
    RAYYAN_OUT,
    ROOT,
    SCRAPED_COL,
    VAGINAL_RIS,
    acc_set,
    fmt_acc,
    norm_doi,
    norm_pmid,
    norm_title,
)

OUT_CSV = os.path.join(ROOT, "MMC2_master_Aug2026.csv")
OUT_MD = os.path.join(ROOT, "MMC2_accession_merge_report.md")


def load_keys(path, title_col="title", doi_col="doi", pmid_col=None):
    """Return (titles, dois, pmids) key sets from a CSV source arm."""
    df = pd.read_csv(path, low_memory=False)
    cols = {c.lower(): c for c in df.columns}
    t = cols.get(title_col, None)
    d = cols.get(doi_col, None)
    p = pmid_col or next((c for c in df.columns if "pubmed" in c.lower() or c.lower() == "pmid"), None)
    return (
        set(df[t].map(norm_title).dropna()) if t else set(),
        set(df[d].map(norm_doi).dropna()) if d else set(),
        set(df[p].map(norm_pmid).dropna()) if p else set(),
        len(df),
    )


def load_ris_keys(path):
    """Parse titles / DOIs / accession numbers out of an RIS export."""
    txt = open(path, encoding="utf8", errors="replace").read()
    titles = {norm_title(m.group(1)) for m in re.finditer(r"^TI\s+-\s+(.+)$", txt, re.M)}
    dois = {norm_doi(m.group(1)) for m in re.finditer(r"^DO\s+-\s+(.+)$", txt, re.M)}
    pmids = {norm_pmid(m.group(1)) for m in re.finditer(r"^(?:AN|ID)\s+-\s+(\d{6,9})", txt, re.M)}
    n = len(re.findall(r"^TI\s+-\s+", txt, re.M))
    return titles - {None}, dois - {None}, pmids - {None}, n


def main():
    out = pd.read_csv(RAYYAN_OUT, low_memory=False)
    mo = pd.read_csv(MOHAK, low_memory=False)
    n_rows = len(out)

    out["_doi"] = out["doi"].map(norm_doi)
    out["_title"] = out["title"].map(norm_title)
    out["_pmid"] = out["pubmed_id"].map(norm_pmid)
    mo["_doi"] = mo["doi [scraped sheet]"].map(norm_doi)

    # --- integrity: the two Rayyan sheets must describe the same corpus -----
    so, sm = set(out["_doi"].dropna()), set(mo["_doi"].dropna())
    assert not (so ^ sm), f"DOI sets diverge: {len(so - sm)} out-only, {len(sm - so)} mohak-only"

    # --- accession sets ----------------------------------------------------
    out["_old"] = out["Accession Code"].map(acc_set)
    mo["_scr"] = mo[SCRAPED_COL].map(acc_set)

    by_doi = mo.dropna(subset=["_doi"]).groupby("_doi")["_scr"].agg(lambda s: frozenset().union(*s))
    cert = mo.dropna(subset=["_doi"]).groupby("_doi")["certainty"].agg(
        lambda s: s.dropna().iloc[0] if s.notna().any() else None
    )
    flag = mo.dropna(subset=["_doi"]).groupby("_doi")["flag_detail"].agg(
        lambda s: s.dropna().iloc[0] if s.notna().any() else None
    )

    empty = frozenset()
    out["_scr"] = out["_doi"].map(by_doi).map(lambda v: v if isinstance(v, frozenset) else empty)
    out["scrape_certainty"] = out["_doi"].map(cert)
    out["scrape_flag_detail"] = out["_doi"].map(flag)

    out["_merged"] = [a | b for a, b in zip(out["_old"], out["_scr"])]
    out["accession_scraped"] = out["_scr"].map(fmt_acc)
    out["accession_merged"] = out["_merged"].map(fmt_acc)

    has_old = out["_old"].map(bool)
    has_scr = out["_scr"].map(bool)
    out["accession_source"] = [
        "both" if o and s else "old_only" if o else "scraped_only" if s else "neither"
        for o, s in zip(has_old, has_scr)
    ]
    out["accession_conflict"] = [
        bool(o) and bool(s) and not (o & s) for o, s in zip(out["_old"], out["_scr"])
    ]

    out["accession_status"] = [
        "present"
        if m
        else "absent_verified"
        if c == "verified-empty"
        else "unknown_not_assessed"
        if c == "needs-review"
        else "unknown_no_scrape_record"
        for m, c in zip(out["_merged"].map(bool), out["scrape_certainty"])
    ]

    # --- source arm --------------------------------------------------------
    oq = pd.read_csv(ORIGINAL_QUERY, low_memory=False)
    oq_t = set(oq["title"].map(norm_title).dropna())
    oq_d = set(oq["doi"].map(norm_doi).dropna())
    oq_p = set(oq["pubmed_id"].map(norm_pmid).dropna())

    vag_t, vag_d, vag_p, n_vag = load_ris_keys(VAGINAL_RIS)
    can_t, can_d, can_p, n_can = load_keys(CANCER_CSV)
    ext_t, ext_d, ext_p, n_ext = load_keys(EXTRA_CSV)

    def hits(t, d, p):
        return out["_title"].isin(t) | out["_doi"].isin(d) | out["_pmid"].isin(p)

    in_main, in_vag = hits(oq_t, oq_d, oq_p), hits(vag_t, vag_d, vag_p)
    in_cancer = hits(can_t, can_d, can_p) | hits(ext_t, ext_d, ext_p)

    out["source_arm"] = [
        "main_query" if m else "cervicovaginal" if v else "cancer" if c else "unresolved"
        for m, v, c in zip(in_main, in_vag, in_cancer)
    ]
    out["in_main_query"] = in_main

    # --- ids / duplicate groups -------------------------------------------
    out["record_id"] = [f"MMC2-{i + 1:05d}" for i in range(n_rows)]
    dup_key = out["_doi"].where(out["_doi"].notna())
    counts = dup_key.value_counts()
    dupes = set(counts[counts > 1].index)
    codes = {d: f"DUP-{i + 1:03d}" for i, d in enumerate(sorted(dupes))}
    out["dup_group"] = [codes.get(d, "") if pd.notna(d) else "" for d in dup_key]

    drop = [c for c in out.columns if c.startswith("_")]
    out.drop(columns=drop).to_csv(OUT_CSV, index=False)

    # --- report ------------------------------------------------------------
    src = out["accession_source"].value_counts()
    n_old, n_scr, n_uni = int(has_old.sum()), int(has_scr.sum()), int(out["_merged"].map(bool).sum())

    # DOI-level view (rows collapse to unique works)
    d = pd.DataFrame(
        {
            "old": out.dropna(subset=["_doi"]).groupby("_doi")["_old"].agg(lambda s: frozenset().union(*s)),
            "new": out.dropna(subset=["_doi"]).groupby("_doi")["_scr"].agg(lambda s: frozenset().union(*s)),
        }
    )
    both = d[(d.old.map(bool)) & (d.new.map(bool))]
    identical = int((both.old == both.new).sum())

    pmc = out["pmc_id"].notna()
    present = out["accession_status"].eq("present")

    L = []
    a = L.append
    a("# MMC2 accession merge — overlap and gap report\n")
    a(f"Source rows: **{n_rows:,}** ({out['_doi'].nunique():,} unique DOIs, "
      f"{int(out['_doi'].isna().sum())} without a DOI)\n")
    a("## Coverage\n")
    a("| Source | rows with a real code |")
    a("|---|---|")
    a(f"| `Accession Code` (existing) | {n_old:,} ({100 * n_old / n_rows:.1f}%) |")
    a(f"| `Accession Code [scraped sheet]` | {n_scr:,} ({100 * n_scr / n_rows:.1f}%) |")
    a(f"| **Union (`accession_merged`)** | **{n_uni:,} ({100 * n_uni / n_rows:.1f}%)** |\n")
    a("## Provenance of the merged column\n")
    a("| `accession_source` | rows |")
    a("|---|---|")
    for k in ("both", "scraped_only", "old_only", "neither"):
        a(f"| {k} | {int(src.get(k, 0)):,} |")
    a("")
    a(f"Among the {len(both):,} DOIs where **both** sheets carry codes, "
      f"**{identical:,} ({100 * identical / len(both):.1f}%) are identical sets**. "
      f"{int(out['accession_conflict'].sum())} rows have zero overlap and are flagged "
      "`accession_conflict` for manual review.\n")
    a("## Assessment status\n")
    a("`accession_status` separates *verified absent* from *never assessed* so the "
      "un-scraped papers are not silently counted as non-depositing.\n")
    a("| status | rows |")
    a("|---|---|")
    for k, v in out["accession_status"].value_counts().items():
        a(f"| {k} | {int(v):,} |")
    a("")
    a("| scrape certainty | rows |")
    a("|---|---|")
    for k, v in out["scrape_certainty"].value_counts(dropna=False).items():
        a(f"| {k} | {int(v):,} |")
    a("")
    a("## Open-access confound\n")
    a(f"Accession present: **{100 * present[pmc].mean():.1f}%** of PMC-available papers "
      f"vs **{100 * present[~pmc].mean():.1f}%** of non-PMC papers. Deposition rates are "
      "therefore entangled with full-text retrievability and must be reported as such.\n")
    a("## Gap by year\n")
    y = pd.to_numeric(out["year"], errors="coerce")
    g = out.groupby(y.astype("Int64"))["accession_source"].value_counts().unstack(fill_value=0)
    g["total"] = g.sum(axis=1)
    g = g[g.total >= 50]
    a("| year | total | both | old_only | scraped_only | neither | % gained by scrape |")
    a("|---|---|---|---|---|---|---|")
    for yr, r in g.iterrows():
        a(f"| {yr} | {int(r.total)} | {int(r.get('both', 0))} | {int(r.get('old_only', 0))} | "
          f"{int(r.get('scraped_only', 0))} | {int(r.get('neither', 0))} | "
          f"{100 * r.get('scraped_only', 0) / r.total:.1f}% |")
    a("")
    a("## Source arms\n")
    a("| arm | rows | identification source |")
    a("|---|---|---|")
    lbl = {
        "main_query": f"`original_query.csv` (n={len(oq):,})",
        "cervicovaginal": f"`cervicovaginal_microbiome_full.ris` (n={n_vag:,})",
        "cancer": f"`cancer_articles.csv` (n={n_can:,}) + `extra_articles.csv` (n={n_ext:,})",
        "unresolved": "not traced to any identification file",
    }
    for k, v in out["source_arm"].value_counts().items():
        a(f"| {k} | {int(v):,} | {lbl.get(k, '')} |")
    a("")
    unres = out[out["source_arm"].eq("unresolved")]
    if len(unres):
        a(f"### The {len(unres)} unresolved records\n")
        for t in unres["title"].dropna():
            a(f"- {str(t)[:150]}")
        a("")
    a("## Conflicts requiring manual review\n")
    cf = out[out["accession_conflict"]]
    a(f"{len(cf)} rows. Many are cosmetic (whitespace, or a BioProject/SRA alias pair).\n")
    a("| DOI | old | scraped |")
    a("|---|---|---|")
    for _, r in cf.iterrows():
        a(f"| {r['doi']} | {fmt_acc(r['_old'])} | {fmt_acc(r['_scr'])} |")
    a("")

    open(OUT_MD, "w").write("\n".join(L))

    print(f"wrote {OUT_CSV}  ({n_rows:,} rows x {len(out.columns) - len(drop)} cols)")
    print(f"wrote {OUT_MD}")
    print(f"  old={n_old}  scraped={n_scr}  union={n_uni}")
    print(f"  {dict(src)}")
    print(f"  arms: {dict(out['source_arm'].value_counts())}")


if __name__ == "__main__":
    main()
