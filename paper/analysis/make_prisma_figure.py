#!/usr/bin/env python3
"""PRISMA 2020 flow diagram for MMC2.

Single column, five stages.  The cervicovaginal search is not drawn as a branch:
its records entered the same Rayyan pile at identification and were
de-duplicated against the main query there, so it belongs inside the
identification box.

Every count is read from MMC2_prisma_counts.tsv (build_prisma_flow.py) and the
master, so the figure cannot drift from the text.

Writes Figs/MMC2_Fig_PRISMA.pdf / .png and the same as MMC2_Fig_S7.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
import pandas as pd                                    # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402

p = lambda f: os.path.join(ROOT, f)
FIGS = p("Figs")
MAIN, EXCL = "#eaf1f7", "#f6f0e8"      # flow boxes vs exclusion boxes
EDGE, RULE = "#2c3e50", "#8a9aa8"
BAND = "#dfe6ec"

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42,
                     "font.family": "sans-serif",
                     "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"]})


def counts():
    c = pd.read_csv(p("MMC2_prisma_counts.tsv"), sep="\t")
    d = dict(zip(c["description"], c["n"]))
    return {
        "identified": d["Records identified in PubMed"],
        "main": d["  main query"],
        "cervico": d["  cervicovaginal search"],
        "dedup": d["Records after duplicates removed"],
        "excl_dup": d["  excluded: duplicate records"],
        "tiab": d["Records after title/abstract screening"],
        "excl_tiab": d["  excluded at title/abstract screening"],
        "fulltext": d["Records after full-text screening"],
        "excl_full": d["  excluded at full-text screening"],
        "analysis": d["Studies included in the analysis"],
        "excl_final": d["  excluded at eligibility and curation"],
        "noninsdc": d["  excluded: deposited outside INSDC"],
        "healthy": d["  excluded: healthy-cohort studies"],
        "removed": d["  excluded: non-human, duplicate deposit, withdrawn"],
        "tiers": {t: d[f"Tier {t}"] for t in "1234"},
    }


LINE = 0.0130          # vertical space one line of box text needs
PAD = 0.0130           # padding above and below the text block


def box(ax, x, y, w, text, fill=MAIN, weight="normal", size=8.6):
    """Draw a box sized to its own content, so text can never overflow."""
    h = len(text.split("\n")) * LINE + PAD
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.003,rounding_size=0.008",
                                facecolor=fill, edgecolor=EDGE, linewidth=0.9))
    ax.text(x, y, text, ha="center", va="center", fontsize=size,
            fontweight=weight, linespacing=1.42, color="#1b2733")
    return h


def vline(ax, x, y0, y1, arrow=True):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1),
                                 arrowstyle="-|>" if arrow else "-",
                                 mutation_scale=11, linewidth=0.9, color=EDGE,
                                 shrinkA=0, shrinkB=0))


def hline(ax, x0, x1, y, arrow=True):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y),
                                 arrowstyle="-|>" if arrow else "-",
                                 mutation_scale=11, linewidth=0.9, color=EDGE,
                                 shrinkA=0, shrinkB=0))


def main():
    k = counts()
    f = lambda n: f"{n:,}"
    t = k["tiers"]
    fig, ax = plt.subplots(figsize=(8.9, 8.6))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    LX, RX = 0.400, 0.795
    LW, RW = 0.40, 0.345

    # ---- identification --------------------------------------------------
    y = 0.965
    h1 = box(ax, LX, y, 0.52,
             f"Records identified in PubMed\nn = {f(k['identified'])}",
             weight="bold", size=9.6)
    ax.text(LX, y - h1 / 2 - 0.010,
            f"main query {f(k['main'])}  ·  cervicovaginal search {f(k['cervico'])}",
            ha="center", va="top", fontsize=7.8, color="#4a5b6b")
    top = y + h1 / 2
    prev_bot = y - h1 / 2 - 0.016
    bands = [(prev_bot, top, "Identification")]

    ROWS = [
        (f"Records after duplicates removed\nn = {f(k['dedup'])}",
         f"Duplicate records removed\nn = {f(k['excl_dup'])}",
         "Screening", False),
        (f"Records after title and abstract\nscreening    n = {f(k['tiab'])}",
         f"Records excluded at title and\nabstract screening  n = {f(k['excl_tiab'])}",
         "Screening", False),
        (f"Records after full-text screening\nn = {f(k['fulltext'])}",
         f"Records excluded at full-text\nscreening  n = {f(k['excl_full'])}",
         "Screening", False),
        (f"Studies included in the analysis\nn = {f(k['analysis'])}",
         f"Records excluded  n = {f(k['excl_final'])}\n"
         f"no epidemiological or interventional\nstudy-design MeSH heading, or no\n"
         f"disease or condition MeSH heading;\n"
         f"sequence data outside INSDC ({f(k['noninsdc'])});\n"
         f"healthy cohorts with no disease\ncontrast ({f(k['healthy'])}); non-human, duplicate\n"
         f"deposit, withdrawn ({f(k['removed'])})",
         "Eligibility", True),
    ]

    # Left-column boxes are spaced evenly on their own bottoms, so a tall
    # exclusion box on the right never opens a gap in the flow.
    GAP = 0.052
    stage_top, prev_stage = prev_bot, None
    floor = prev_bot
    for left, right, stage, bold in ROWS:
        hl = len(left.split("\n")) * LINE + PAD
        hr = len(right.split("\n")) * LINE + PAD
        yc = prev_bot - GAP - hl / 2
        vline(ax, LX, prev_bot, yc + hl / 2)
        box(ax, LX, yc, LW, left, weight="bold" if bold else "normal",
            size=9.6 if bold else 8.6)
        box(ax, RX, yc, RW, right, fill=EXCL, size=7.6 if hr < 0.09 else 7.1)
        hline(ax, LX + LW / 2, RX - RW / 2, yc)
        if prev_stage and stage != prev_stage:
            bands.append((prev_bot, stage_top, prev_stage))
            stage_top = yc + max(hl, hr) / 2
        prev_stage = stage
        prev_bot = yc - hl / 2
        floor = min(floor, yc - hr / 2)

    y6 = min(prev_bot - GAP, floor - 0.030) - 0.0
    h6 = len(f"x\nx".split("\n")) * LINE + PAD
    y6 = y6 - h6 / 2
    h6 = box(ax, LX, y6, LW,
             f"Tier 1  {f(t['1'])}        Tier 2  {f(t['2'])}\n"
             f"Tier 3  {f(t['3'])}        Tier 4  {f(t['4'])}", size=8.8)
    vline(ax, LX, prev_bot, y6 + h6 / 2)
    bands.append((y6 - h6 / 2, stage_top, "Included"))

    for y0, y1, label in bands:
        ax.add_patch(FancyBboxPatch((0.022, y0), 0.044, y1 - y0,
                                    boxstyle="round,pad=0,rounding_size=0.006",
                                    facecolor=BAND, edgecolor="none"))
        ax.text(0.044, (y0 + y1) / 2, label, rotation=90, ha="center",
                va="center", fontsize=9.2, fontweight="bold", color="#33475b")

    ax.set_ylim(min(y6 - h6 / 2, floor) - 0.02, 1.0)
    os.makedirs(FIGS, exist_ok=True)
    for ext in ("pdf", "png"):
        for stem in ("MMC2_Fig_PRISMA", "MMC2_Fig_S7"):
            fig.savefig(os.path.join(FIGS, f"{stem}.{ext}"), dpi=300,
                        bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote Figs/MMC2_Fig_PRISMA.pdf / .png")
    for lab, key in [("identified", "identified"), ("de-duplicated", "dedup"),
                     ("after title/abs", "tiab"), ("after full text", "fulltext"),
                     ("INCLUDED", "analysis")]:
        print(f"  {lab:<16} {k[key]:>7,}")


if __name__ == "__main__":
    main()
