#!/usr/bin/env python3
"""Journal impact factor against reusability, as a plain scatter.

Two renders from one builder so the main figure and the supplement cannot
drift apart:

  Figs/Fig2B_if_scatter.{pdf,png}   Tier 1-2 and Tier 4      -> main Figure 2B
  Figs/MMC2_Fig_S3.{pdf,png}        Tier 1-2, Tier 3, Tier 4 -> supplement

Every point is a journal, all the same size; the line is the study-level
logistic fit on log10(impact factor). Restricted to journals carrying a
Clarivate JIF and at least ten studies.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402
from scipy.special import expit  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402

FIGS = os.path.join(ROOT, "Figs")
INK, MUTED, GRID = "#1a1a1a", "#4d4d4d", "#d9d9d9"
DOT, LINE = "#9a9a9a", "#000000"
MIN_STUDIES = 10

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 9,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.edgecolor": GRID, "figure.dpi": 120,
})

PANELS = {
    "t12": ("Tier 1–2", lambda t: t.isin([1, 2])),
    "t3": ("Tier 3", lambda t: t.eq(3)),
    "t4": ("Tier 4", lambda t: t.eq(4)),
}


def data():
    at = pd.read_csv(os.path.join(ROOT, "MMC2_study_attributes.tsv"),
                     sep="\t", low_memory=False)
    mr = pd.read_csv(os.path.join(ROOT, "MMC2_master_reviewed.tsv"), sep="\t",
                     low_memory=False, keep_default_na=False, dtype=str)
    d = at.merge(mr[["record_id", "final_tier", "in_analysis"]],
                 on="record_id", how="left")
    d = d[d["in_analysis"].eq("yes")].copy()
    d["tier"] = pd.to_numeric(d["final_tier"], errors="coerce")
    d = d[d["tier"].notna() & d["impact_factor"].notna()
          & (d["impact_factor"] > 0)].copy()
    d["tier"] = d["tier"].astype(int)
    vc = d["journal"].value_counts()
    d = d[d["journal"].isin(vc[vc >= MIN_STUDIES].index)].copy()
    d["logIF"] = np.log10(d["impact_factor"])
    for k, (_, fn) in PANELS.items():
        d[k] = fn(d["tier"]).astype(int)
    return d


def panel(ax, d, key, first):
    label, _ = PANELS[key]
    m = smf.logit(f"{key} ~ logIF", data=d).fit(disp=False)
    orr = float(np.exp(m.params["logIF"]))
    lo, hi = np.exp(m.conf_int().loc["logIF"])
    p = float(m.pvalues["logIF"])

    per_j = d.groupby("journal").agg(x=("logIF", "first"), y=(key, "mean"),
                                     n=(key, "size"))
    ax.scatter(10 ** per_j.x, 100 * per_j.y, s=34, color=DOT,
               edgecolor="none", alpha=0.85, zorder=3)
    grid = np.linspace(d.logIF.min(), d.logIF.max(), 200)
    ax.plot(10 ** grid, 100 * expit(m.params["Intercept"]
                                    + m.params["logIF"] * grid),
            color=LINE, lw=1.8, zorder=4)
    ax.set_xscale("log")
    ax.set_xlabel("Journal impact factor")
    if first:
        ax.set_ylabel("Percentage of studies (%)")
    ax.set_title(label, fontsize=11, fontweight="bold", pad=14)

    ticks = [t for t in (2, 3, 5, 10, 20, 50)
             if d.impact_factor.min() * 0.85 <= t <= d.impact_factor.max() * 1.15]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())

    # put the stats in whichever corner holds the fewest journals, so the box
    # never lands on the data
    xr = (per_j.x - per_j.x.min()) / max(per_j.x.max() - per_j.x.min(), 1e-9)
    yr = per_j.y
    corners = {("left", "top"): ((xr < .45) & (yr > .62)).sum(),
               ("right", "top"): ((xr > .55) & (yr > .62)).sum(),
               ("left", "bottom"): ((xr < .45) & (yr < .38)).sum(),
               ("right", "bottom"): ((xr > .55) & (yr < .38)).sum()}
    (ha, va) = min(corners, key=corners.get)
    ptxt = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
    ax.annotate(f"OR = {orr:.2f} [{lo:.2f}–{hi:.2f}]\n{ptxt}",
                xy=(0.035 if ha == "left" else 0.965,
                    0.97 if va == "top" else 0.03),
                xycoords="axes fraction", ha=ha, va=va, fontsize=8,
                color=MUTED,
                bbox=dict(boxstyle="square,pad=0.28", facecolor="white",
                          edgecolor="none", alpha=0.85))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.yaxis.grid(True, color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    return orr, lo, hi, p, int(m.nobs)


def render(d, keys, stem, width):
    fig, axs = plt.subplots(1, len(keys), figsize=(width, 3.5), sharey=True)
    axs = np.atleast_1d(axs)
    out = []
    for i, k in enumerate(keys):
        out.append((k,) + panel(axs[i], d, k, i == 0))
    axs[0].set_ylim(0, 100)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.86, bottom=0.17,
                        wspace=0.10)
    os.makedirs(FIGS, exist_ok=True)
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"{stem}.{e}"), dpi=300,
                    facecolor="white")
    plt.close(fig)
    return out


def main():
    d = data()
    print(f"journals: {d.journal.nunique()}   studies: {len(d):,}\n")
    for stem, keys, w in [("Fig2B_if_scatter", ["t12", "t4"], 7.4),
                          ("MMC2_Fig_S3", ["t12", "t3", "t4"], 10.4)]:
        print(f"{stem}")
        for k, orr, lo, hi, p, n in render(d, keys, stem, w):
            print(f"   {PANELS[k][0]:9s} OR={orr:.2f} [{lo:.2f}-{hi:.2f}] "
                  f"p={p:.3f}  n={n:,}")
    print(f"\nwrote {FIGS}/Fig2B_if_scatter.pdf and {FIGS}/MMC2_Fig_S3.pdf")


if __name__ == "__main__":
    main()
