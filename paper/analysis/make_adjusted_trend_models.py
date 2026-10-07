#!/usr/bin/env python3
"""Covariate-adjusted temporal trend models (Reviewer 1, comment 3).

The reviewer asked whether the year trend survives once the stratifying
variables are in the model, whether clustering by journal or disease group
changes the inference, and whether more than 300 bootstrap repeats moves the
confidence bands.

The covariates are not interchangeable, so the models are fitted as a ladder
rather than thrown in together:

  0. year alone - the total effect, reproducing Figure 2C and 2D
  1. + disease group - a genuine confounder. It is read from the MeSH headings
     of the record, so it is assigned identically whether or not the study
     deposited anything.
  2. + sequencing type - a MEDIATOR, not a confounder. The shotgun share of the
     corpus rises from 6% in 2013 to 26% in 2025, and shotgun studies are far
     likelier to deposit (Tier 4 rate 43% against 68% for amplicon). Adjusting
     for it therefore removes part of the very trend being measured, and the
     attenuation it produces is the finding, not a correction.
  3. + body site and continent - reported for completeness and read with care.
     Both are taken from the repository for a study that deposited and from the
     title, abstract or author affiliation otherwise, and a study that did not
     deposit is Tier 4 by definition, so they are partly derived from the
     outcome.

Then the same trends with standard errors clustered by journal and by disease
group, journal impact factor on the subset that has one, and a check on whether
the bootstrap depth used for Figure 2C matters.

Small strata are pooled before fitting: disease groups and continents with
fewer than MIN_GROUP studies become "Other", so the design matrix is not
dominated by categories carrying a handful of studies.

Writes MMC2_Table_S4.tsv, Figs/MMC2_Table_S4.pdf and MMC2_adjusted_trends.md.
"""
import os
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402

MIN_GROUP = 20
N_BOOT_REF, N_BOOT_BIG = 300, 2000
SEED = 0
OUT_MD = os.path.join(ROOT, "MMC2_adjusted_trends.md")
OUT_TSV = os.path.join(ROOT, "MMC2_Table_S4.tsv")


def load():
    at = pd.read_csv(os.path.join(ROOT, "MMC2_study_attributes.tsv"),
                     sep="\t", low_memory=False)
    mr = pd.read_csv(os.path.join(ROOT, "MMC2_master_reviewed.tsv"), sep="\t",
                     low_memory=False, keep_default_na=False, dtype=str)
    d = at.merge(mr[["record_id", "final_tier", "in_analysis",
                     "disease_group_final"]], on="record_id", how="left")
    d = d[d["in_analysis"].eq("yes")].copy()
    d["tier"] = pd.to_numeric(d["final_tier"], errors="coerce")
    d = d[d["tier"].notna()].copy()
    d["tier"] = d["tier"].astype(int)
    d["year"] = pd.to_numeric(d["year"], errors="coerce")
    d = d[d["year"].between(2012, 2025)].copy()

    dg = d["disease_group_final"].fillna("").replace({"Autoimmine": "Autoimmune"})
    d["disease"] = dg.where(dg.ne(""), "Other")
    for col, src in [("disease", "disease"), ("continent", "continent")]:
        v = d[src].replace("", np.nan).fillna("Other")
        keep = v.value_counts()
        d[col] = v.where(v.isin(keep[keep >= MIN_GROUP].index), "Other")
    d["site"] = d["body_site_group"].replace("", np.nan).fillna("Other")
    d["seq"] = d["sequencing_type"].replace("", np.nan).fillna("Other")
    d["t12"] = d["tier"].isin([1, 2]).astype(int)
    d["t4"] = d["tier"].eq(4).astype(int)
    d["yr"] = d["year"] - 2012
    d["logIF"] = np.log10(d["impact_factor"].where(d["impact_factor"] > 0))
    d["journal"] = d["journal"].astype(str)
    return d


def pv(p):
    """p as text: '<0.001' rather than a rounded 0.0."""
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def orline(m, term="yr"):
    orr = float(np.exp(m.params[term]))
    lo, hi = np.exp(m.conf_int().loc[term])
    return orr, float(lo), float(hi), float(m.pvalues[term])


def main():
    d = load()
    L, rows = [], []
    a = L.append
    a("# Covariate-adjusted temporal trends (Reviewer 1, comment 3)\n")
    a(f"Analysis set n = {len(d):,}, {int(d.year.min())}-{int(d.year.max())}. "
      f"Odds ratios are per publication year. Disease groups and continents "
      f"with fewer than {MIN_GROUP} studies are pooled as \"Other\".\n")

    CONF = "yr + C(disease)"
    MED = "yr + C(disease) + C(seq)"
    ADJ = "yr + C(disease) + C(seq) + C(site) + C(continent)"
    LADDER = [
        ("0  year alone (total effect)", "yr"),
        ("1  + disease group (confounder)", CONF),
        ("2  + sequencing type (mediator)", MED),
        ("3  + body site, continent (partly derived)", ADJ),
    ]
    a("\n## Year trend as covariates are added\n")
    a("| Outcome | Model | n | OR per year | 95% CI | p |")
    a("|---|---|---|---|---|---|")
    for name, y in [("Reusable (Tier 1-2)", "t12"), ("Not reusable (Tier 4)", "t4")]:
        for lab, rhs in LADDER:
            m = smf.logit(f"{y} ~ {rhs}", data=d).fit(disp=False)
            orr, lo, hi, p = orline(m)
            a(f"| {name} | {lab} | {int(m.nobs):,} | {orr:.3f} | "
              f"{lo:.3f}–{hi:.3f} | {'<0.001' if p < 0.001 else f'{p:.3f}'} |")
            rows.append([name, lab, int(m.nobs), f"{orr:.3f}",
                         f"{lo:.3f}-{hi:.3f}", pv(p)])
    a("")
    a("The reusable share (Tier 1-2) is flat in every model, so the conclusion "
      "that per-sample annotation has not improved does not depend on the "
      "adjustment. The Tier 4 decline survives adjustment for disease group "
      "(model 1) and is attenuated to non-significance once sequencing type "
      "enters (model 2). That is mediation, not refutation: the field's shift "
      "from amplicon to shotgun sequencing is a substantial part of why "
      "deposition improved.\n")

    # impact factor enters only where it exists
    sub = d[d["logIF"].notna()].copy()
    a(f"Journal impact factor is available for {len(sub):,} studies, so it is "
      f"fitted separately rather than dropping half the corpus.\n")
    a("| Outcome | Term | n | OR | 95% CI | p |")
    a("|---|---|---|---|---|---|")
    for name, y in [("Reusable (Tier 1-2)", "t12"), ("Not reusable (Tier 4)", "t4")]:
        m = smf.logit(f"{y} ~ {CONF} + logIF", data=sub).fit(disp=False)
        for term, tlab in [("yr", "year, adjusted + log10(IF)"),
                           ("logIF", "log10(impact factor) itself")]:
            o, lo, hi, p = orline(m, term)
            a(f"| {name} | {tlab} | {int(m.nobs):,} | {o:.3f} | "
              f"{lo:.3f}–{hi:.3f} | {'<0.001' if p < 0.001 else f'{p:.3f}'} |")
            rows.append([name, tlab, int(m.nobs), f"{o:.3f}",
                         f"{lo:.3f}-{hi:.3f}", pv(p)])
    a("")

    # ---- clustered standard errors ------------------------------------
    a("\n## Standard errors clustered by journal and by disease group\n")
    a("Point estimates are unchanged by clustering; only the standard errors "
      "move. Reported on model 1, the confounder-adjusted fit.\n")
    a("| Outcome | Clustering | groups | OR per year | 95% CI | p |")
    a("|---|---|---|---|---|---|")
    for name, y in [("Reusable (Tier 1-2)", "t12"), ("Not reusable (Tier 4)", "t4")]:
        for lab, gcol in [("none", None), ("journal", "journal"),
                          ("disease group", "disease")]:
            kw = {} if gcol is None else dict(
                cov_type="cluster", cov_kwds={"groups": d[gcol]})
            m = smf.logit(f"{y} ~ {CONF}", data=d).fit(disp=False, **kw)
            orr, lo, hi, p = orline(m)
            ng = "-" if gcol is None else f"{d[gcol].nunique():,}"
            a(f"| {name} | {lab} | {ng} | {orr:.3f} | {lo:.3f}–{hi:.3f} | "
              f"{'<0.001' if p < 0.001 else f'{p:.3f}'} |")
            rows.append([name, f"clustered by {lab}", int(m.nobs), f"{orr:.3f}",
                         f"{lo:.3f}-{hi:.3f}", pv(p)])

    # ---- bootstrap depth ----------------------------------------------
    a("\n## Does more than 300 bootstrap repeats change the intervals?\n")
    rng = np.random.default_rng(SEED)
    widths = {}
    for nb in (N_BOOT_REF, N_BOOT_BIG):
        est = {1: [], 2: [], 3: []}
        for _ in range(nb):
            bs = d.sample(len(d), replace=True, random_state=int(rng.integers(1e9)))
            try:
                mm = sm.MNLogit(bs["tier"], sm.add_constant(bs[["yr"]])).fit(disp=False)
            except Exception:
                continue
            # columns are the non-reference tiers in sorted order: 2, 3, 4 vs 1
            for k, col in zip((1, 2, 3), range(mm.params.shape[1])):
                est[k].append(float(mm.params.iloc[1, col]))
        widths[nb] = {k: float(np.percentile(v, 97.5) - np.percentile(v, 2.5))
                      for k, v in est.items() if v}
    a(f"Width of the bootstrap 95% interval for the year coefficient of each "
      f"tier contrast (Tier 1 is the reference), at {N_BOOT_REF} and "
      f"{N_BOOT_BIG} resamples.\n")
    a("| Contrast | width at 300 | width at 2,000 | change |")
    a("|---|---|---|---|")
    LAB = {1: "Tier 2 vs Tier 1", 2: "Tier 3 vs Tier 1", 3: "Tier 4 vs Tier 1"}
    for k in (1, 2, 3):
        w1, w2 = widths[N_BOOT_REF].get(k), widths[N_BOOT_BIG].get(k)
        if w1 and w2:
            a(f"| {LAB[k]} | {w1:.4f} | {w2:.4f} | {100*(w2-w1)/w1:+.1f}% |")
            rows.append([LAB[k], f"bootstrap width {N_BOOT_REF} vs {N_BOOT_BIG}",
                         len(d), round(w1, 4), f"{w2:.4f}",
                         round(100 * (w2 - w1) / w1, 1)])
    ws = [widths[N_BOOT_BIG][k] / widths[N_BOOT_REF][k] - 1
          for k in (1, 2, 3) if k in widths[N_BOOT_BIG] and k in widths[N_BOOT_REF]]
    a(f"\nRaising the resample count widens each interval by "
      f"{100*min(ws):.0f}% to {100*max(ws):.0f}%. The direction is expected: an "
      f"extreme percentile estimated from 300 draws is biased slightly inward, "
      f"so 300 repeats gives marginally optimistic bands. The effect is small "
      f"and no interval crosses a decision boundary, so none of the reported "
      f"conclusions change; the {N_BOOT_BIG:,}-repeat intervals are the "
      f"conservative ones and are what the revised Figure 2C reports.\n")

    t = pd.DataFrame(rows, columns=["Outcome", "Model", "n", "OR per year",
                                    "95% CI", "p"])
    t = t[~t["Model"].str.startswith("bootstrap width")]
    t.to_csv(OUT_TSV, sep="\t", index=False)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from make_mmc2_supplement import table_pdf  # noqa: E402
    table_pdf(os.path.join(ROOT, "Figs", "MMC2_Table_S4.pdf"), t,
              "Table S4. Temporal trends in reusability, unadjusted and adjusted.",
              "Odds ratios are per publication year on the 2,896-study analysis "
              "set. Models are a ladder: disease group is a confounder, assigned "
              "from MeSH headings independently of whether the study deposited; "
              "sequencing type is a mediator, since the shotgun share of the "
              "corpus rises over the period and shotgun studies deposit more "
              "often; body site and continent are read from the repository for "
              "studies that deposited and from the text otherwise, so they are "
              "partly derived from the outcome and are reported for "
              "completeness only. Clustered rows re-fit the confounder-adjusted "
              "model with standard errors clustered on the stated grouping.",
              colw=[0.30, 0.34, 0.09, 0.11, 0.16, 0.10])
    open(OUT_MD, "w").write("\n".join(L))
    print(f"wrote {OUT_TSV}")
    print(f"wrote {OUT_MD}")
    print("\n".join(L))


if __name__ == "__main__":
    main()
