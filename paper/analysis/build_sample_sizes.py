#!/usr/bin/env python3
"""Per-study sample sizes under Sam's 2026-08-31 counting rules.

The old count was every sample under every accession a paper cited, which
credited a study with datasets it reused rather than generated -- American Gut
alone was 18% of all Tier 1+2 samples. The rules now:

  Tiers 1-3  the study's OWN deposit: the FIRST accession listed, which is what
             `repository_link` pointed the reviewer at in all 177
             multi-accession review rows. Named exceptions below override it.
             (Sam, 2 Sep: one rule for every tier. Tier 3's all-accession total
             is still computed, as `n_samples_all_accessions`, because that is
             the size of the curation job still to do -- a workload figure, not
             a reusability figure.)
  Tier 4     no deposit by definition -> the mean of Tiers 1-3, over studies of
             >=10 samples so pooled or broken deposits do not drag it.

Studies whose disease category is "Healthy" are excluded from every total
(`in_analysis` in the master); they are still written out, flagged, so the
PRISMA flow can account for them.

Writes MMC2_study_sample_sizes.tsv: record_id, tier, n_samples, rule.
"""
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402

p = lambda f: os.path.join(ROOT, f)
OUT = p("MMC2_study_sample_sizes.tsv")

# --- Sam's per-study calls, 2026-08-31 ------------------------------------
# Both are now applied in the master by apply_adjudication.py; kept here as
# no-ops so the rule set stays readable in one place.
DROP, RETIER = {}, {}
# a study's own new data, where it is not simply the first accession
OWN_ACCESSIONS = {
    "MMC2-04315": ["PRJNA595346", "PRJNA726955"],
    "MMC2-00496": ["PRJNA731589"],
    "MMC2-03978": ["PRJNA701982"],
}
# counts Sam supplied directly; these beat anything derived from ENA
FIXED_N = {
    "MMC2-12782": 146,   # only PRJNA737052 is new, 146 sequenced
    # Tier 2 top-20
    "MMC2-03989": 96,    "MMC2-02710": 850,   "MMC2-03717": 121,
    "MMC2-00205": 1232,  "MMC2-04458": 1702,  "MMC2-02003": 377,
    "MMC2-03131": 500,   "MMC2-04397": 143,   "MMC2-03919": 195,
    "MMC2-04344": 679,   "MMC2-04364": 280,   "MMC2-02724": 102,
}
MIN_FOR_MEAN = 10
FILL_POOL_TIERS = ("1", "2", "3")
# Sam, 2 Sep: MMC2-03277 sits in Tier 4 because its own data is blocked and
# every accession it cites is public reuse, so it must not be credited with
# those samples -- it takes the Tier 4 fill like any other study with no deposit.
NO_OWN_DEPOSIT = {"MMC2-03277"}


def reported_n(s):
    v = [int(x) for x in re.findall(r"\bn\s*=\s*(\d{1,6})", str(s))]
    return max(v) if v else np.nan


def main():
    m = pd.read_csv(p("MMC2_master_reviewed.tsv"), sep="\t", low_memory=False,
                    keep_default_na=False, dtype=str)
    st = pd.read_csv(p("MMC2_study_metadata.tsv"), sep="\t", low_memory=False,
                     keep_default_na=False, dtype=str)[["record_id", "accession_codes"]]
    an = pd.read_csv(p("MMC2_analysis_set.tsv"), sep="\t", low_memory=False,
                     keep_default_na=False, dtype=str)[["record_id", "abstract"]]
    sm = pd.read_csv(p("MMC2_sample_master.tsv"), sep="\t", low_memory=False,
                     keep_default_na=False, dtype=str,
                     usecols=["record_id", "study_accession", "sample_accession"])
    per_acc = sm.groupby(["record_id", "study_accession"])["sample_accession"].nunique()

    d = m.merge(st, on="record_id", how="left", suffixes=("", "_p")) \
         .merge(an, on="record_id", how="left")
    d = d[d.final_tier.isin(["1", "2", "3", "4"]) & d.remove_flag.ne("yes")].copy()
    for rid, (t, _why) in RETIER.items():
        d.loc[d.record_id.eq(rid), "final_tier"] = t
    d["abstract_n"] = d["abstract"].map(reported_n)

    def accs(r):
        col = r.get("accession_codes_p") or r.get("accession_codes") or ""
        return [a.strip() for a in str(col).split(";") if a.strip()]

    def count(rid, a):
        return int(per_acc.get((rid, a), 0))

    rows = []
    for r in d.to_dict("records"):
        rid, tier, A = r["record_id"], r["final_tier"], accs(r)
        n, rule = np.nan, ""
        if rid in FIXED_N:
            n, rule = FIXED_N[rid], "sam_supplied"
        elif rid in OWN_ACCESSIONS:
            n = sum(count(rid, a) for a in OWN_ACCESSIONS[rid])
            rule = "own_accessions:" + "+".join(OWN_ACCESSIONS[rid])
        elif rid in NO_OWN_DEPOSIT:
            n, rule = np.nan, "no_own_deposit"
        elif tier in ("1", "2", "3"):
            # the first accession that actually resolves -- a leading GEO series
            # or Zenodo DOI has no ENA samples and would otherwise score zero
            for a in A:
                if count(rid, a):
                    n, rule = count(rid, a), f"first_accession:{a}"
                    break
            if pd.isna(n):
                n, rule = (r["abstract_n"], "abstract_n") if pd.notna(r["abstract_n"]) \
                    else (np.nan, "no_count")
        else:                                    # tier 4
            tot = sum(count(rid, a) for a in A)
            n, rule = (tot, "all_accessions") if tot else (np.nan, "fill_tier1_3_mean")
        rows.append(dict(record_id=rid, tier=tier, n_samples=n, rule=rule,
                         n_accessions=len(A), abstract_n=r["abstract_n"],
                         n_samples_all_accessions=sum(count(rid, a) for a in A),
                         in_analysis=r["in_analysis"],
                         disease_group=r["disease_group_final"]))

    out = pd.DataFrame(rows)
    ana = out[out.in_analysis.eq("yes")]        # Healthy never feeds the mean
    t3 = ana[ana.tier.isin(FILL_POOL_TIERS) & (ana.n_samples >= MIN_FOR_MEAN)]
    fill = float(t3.n_samples.mean())
    # Tier 4 gets the Tier 3 mean. A Tier 1/2/3 study whose accession simply did
    # not resolve is a different problem -- it is filled from its OWN tier, so a
    # fetch failure cannot inflate Tier 1 to look like Tier 3.
    own = out[out.n_samples.notna()].groupby("tier")["n_samples"].mean()
    miss = out.n_samples.isna()
    out.loc[miss & out.tier.eq("4"), "n_samples"] = fill
    out.loc[miss & out.tier.eq("4"), "rule"] = "fill_tier1_3_mean"
    for t_ in ("1", "2", "3"):
        sel = miss & out.tier.eq(t_)
        out.loc[sel, "n_samples"] = own[t_]
        out.loc[sel, "rule"] = f"tier{t_}_own_mean_fill(unresolved accession)"
    out["n_samples"] = out["n_samples"].round().astype(int)
    out.to_csv(OUT, sep="\t", index=False)

    print(f"Tier 1-3 mean over {len(t3):,} analysis-set studies with "
          f">={MIN_FOR_MEAN} samples = {fill:.1f} samples/study -> Tier 4 fill")
    print(f"excluded as Healthy: {int(out.in_analysis.ne('yes').sum())}; "
          f"analysis set {int(out.in_analysis.eq('yes').sum()):,} studies")
    ana = out[out.in_analysis.eq("yes")]
    g = ana.groupby("tier").agg(studies=("record_id", "size"),
                                samples=("n_samples", "sum"),
                                median=("n_samples", "median"))
    g["pct_samples"] = (g.samples / g.samples.sum() * 100).round(1)
    print(g.to_string())
    T = g.samples.sum()
    r = g.samples.get("1", 0) + g.samples.get("2", 0)
    print(f"\ntotal {T:,}   T1+T2 {r:,} = {r/T:.1%} reusable   "
          f"NOT reusable {T-r:,} = {1-r/T:.1%}")
    print(f"\nrule counts:\n{ana.rule.str.split(':').str[0].value_counts().to_string()}")
    t3all = int(ana.loc[ana.tier.eq("3"), "n_samples_all_accessions"].sum())
    print(f"\nTier 3 workload figure (all accessions, for 'samples still to "
          f"check'): {t3all:,}  vs {int(ana.loc[ana.tier.eq('3'),'n_samples'].sum()):,} counted")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
