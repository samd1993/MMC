#!/usr/bin/env python3
"""Figure 1 redesigns: the BHAM tier schematic, rebuilt with the final numbers.

Numbers are read from the master and the sample-size table so they cannot
drift. Six designs, one per page, Nature/Science-style: white ground, restrained
palette, thin rules, drawn (not pasted) icons, captions only.

Writes MMC2_Figure1_redesigns.pdf and Figs/Figure1_v*.png.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                        # noqa: E402
import numpy as np                                      # noqa: E402
import pandas as pd                                     # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages    # noqa: E402
from matplotlib.patches import (Circle, FancyBboxPatch, Rectangle, Polygon,   # noqa: E402
                                Ellipse, PathPatch, FancyArrowPatch)
from matplotlib.path import Path                        # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402

p = lambda f: os.path.join(ROOT, f)
FIGS = p("Figs")
OUT = p("MMC2_Figure1_redesigns.pdf")

# Nature-ish palette: muted, print-safe, colour-blind distinguishable
C = {1: "#3B6FA0", 2: "#D9A441", 3: "#E07B39", 4: "#B93A32"}
CL = {1: "#DCE7F2", 2: "#F7ECD3", 3: "#F9E1D2", 4: "#F3D8D5"}   # light fills
INK, MUTED, RULE, GREY = "#1F2933", "#5F6B76", "#C9D1D8", "#B8C0C8"
OK, NO = "#2E7D5B", "#B93A32"

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42,
                     "font.family": "sans-serif",
                     "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial",
                                         "DejaVu Sans"]})

TITLE = {1: "Sample-level metadata\nin repository",
         2: "Disease recoverable\nper sample",
         3: "Sequence data accessible,\nno per-sample biology",
         4: "No reusable\nsequence data"}
DESC = {1: "Per-sample disease, host age and host sex\nare all recorded in the repository.",
        2: "Per-sample disease is recoverable from a\nmetadata field or an informative sample name.",
        3: "A valid accession with unique sample IDs,\nbut no per-sample biological annotation.",
        4: "No accession, an inaccessible one, or\nsamples that carry no differentiating IDs."}
VERDICT = {1: "Reusable as deposited", 2: "Reusable for disease",
           3: "Reads only", 4: "Not reusable"}


def numbers():
    m = pd.read_csv(p("MMC2_master_reviewed.tsv"), sep="\t", low_memory=False,
                    keep_default_na=False, dtype=str)
    a = m[m.in_analysis.eq("yes")]
    ss = pd.read_csv(p("MMC2_study_sample_sizes.tsv"), sep="\t", dtype={"tier": str})
    ss = ss[ss.in_analysis.eq("yes")]
    st = {t: int((a.final_tier == str(t)).sum()) for t in (1, 2, 3, 4)}
    sa = {t: int(ss.loc[ss.tier == str(t), "n_samples"].sum()) for t in (1, 2, 3, 4)}
    return st, sa, sum(st.values()), sum(sa.values())


# ------------------------------------------------------------------ icons --
def icon_disease(ax, x, y, r, col, on=True):
    c = col if on else GREY
    ax.add_patch(Circle((x, y), r * 0.55, facecolor=c, edgecolor="none"))
    for k in range(8):
        a = k * np.pi / 4
        ax.plot([x + r * 0.62 * np.cos(a), x + r * 0.95 * np.cos(a)],
                [y + r * 0.62 * np.sin(a), y + r * 0.95 * np.sin(a)],
                color=c, lw=1.6, solid_capstyle="round")


def icon_age(ax, x, y, r, col, on=True):
    c = col if on else GREY
    w, h = r * 1.7, r * 1.5
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0,rounding_size=%.3f" % (r * 0.18),
                                facecolor="white", edgecolor=c, lw=1.6))
    ax.add_patch(Rectangle((x - w / 2, y + h / 2 - r * 0.42), w, r * 0.42,
                           facecolor=c, edgecolor="none"))
    for i in range(3):
        for j in range(2):
            ax.add_patch(Rectangle((x - w / 2 + r * 0.28 + i * r * 0.5,
                                    y - h / 2 + r * 0.2 + j * r * 0.42),
                                   r * 0.28, r * 0.22, facecolor=c, edgecolor="none"))


def icon_sex(ax, x, y, r, col, on=True):
    """Venus and Mars drawn as strokes, so no font dependency."""
    c = col if on else GREY
    lw = 1.6
    # venus: circle with a cross below, offset left
    vx, vy = x - r * 0.42, y + r * 0.12
    ax.add_patch(Circle((vx, vy), r * 0.36, facecolor="none", edgecolor=c, lw=lw))
    ax.plot([vx, vx], [vy - r * 0.36, vy - r * 0.95], color=c, lw=lw)
    ax.plot([vx - r * 0.25, vx + r * 0.25], [vy - r * 0.7, vy - r * 0.7], color=c, lw=lw)
    # mars: circle with an arrow to the upper right, offset right
    mx, my = x + r * 0.42, y - r * 0.12
    ax.add_patch(Circle((mx, my), r * 0.36, facecolor="none", edgecolor=c, lw=lw))
    ex, ey = mx + r * 0.85, my + r * 0.85
    ax.plot([mx + r * 0.26, ex], [my + r * 0.26, ey], color=c, lw=lw)
    ax.plot([ex - r * 0.36, ex], [ey, ey], color=c, lw=lw)
    ax.plot([ex, ex], [ey - r * 0.36, ey], color=c, lw=lw)


def icon_db(ax, x, y, r, col, on=True):
    c = col if on else GREY
    w, h = r * 1.5, r * 1.5
    ax.add_patch(Rectangle((x - w / 2, y - h / 2 + r * 0.25), w, h - r * 0.5,
                           facecolor="white", edgecolor=c, lw=1.6))
    for dy in (-h / 2 + r * 0.25, y - y + h / 2 - r * 0.25):
        ax.add_patch(Ellipse((x, y + dy), w, r * 0.5, facecolor="white",
                             edgecolor=c, lw=1.6))
    ax.add_patch(Ellipse((x, y + h / 2 - r * 0.25), w, r * 0.5, facecolor=c,
                         edgecolor=c, lw=1.6))


def icon_id(ax, x, y, r, col, on=True):
    c = col if on else GREY
    w, h = r * 1.7, r * 1.1
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0,rounding_size=%.3f" % (r * 0.15),
                                facecolor="white", edgecolor=c, lw=1.6))
    xs = np.linspace(x - w / 2 + r * 0.25, x + w / 2 - r * 0.25, 9)
    for i, xx in enumerate(xs):
        ax.plot([xx, xx], [y - h / 2 + r * 0.22, y + h / 2 - r * 0.22], color=c,
                lw=1.0 if i % 3 else 2.2)


def mark(ax, x, y, ok, r=0.011, fill=None):
    """Tick or cross drawn as strokes inside a disc."""
    ax.add_patch(Circle((x, y), r, facecolor=fill or (OK if ok else NO), edgecolor="none"))
    lw = max(0.9, r * 130)
    if ok:
        ax.plot([x - r * 0.5, x - r * 0.1, x + r * 0.55], [y, y - r * 0.42, y + r * 0.42],
                color="white", lw=lw, solid_capstyle="round", solid_joinstyle="round")
    else:
        ax.plot([x - r * 0.45, x + r * 0.45], [y - r * 0.45, y + r * 0.45], color="white", lw=lw, solid_capstyle="round")
        ax.plot([x - r * 0.45, x + r * 0.45], [y + r * 0.45, y - r * 0.45], color="white", lw=lw, solid_capstyle="round")


def caption(fig, txt):
    fig.text(0.04, 0.035, txt, fontsize=8.6, color=MUTED, va="bottom", ha="left",
             linespacing=1.45)


def frame(fig, title):
    fig.text(0.04, 0.985, title, fontsize=12.5, fontweight="bold", color=INK, va="top")


# --------------------------------------------------------------- design 1 --
def v1_columns(st, sa, N, S):
    fig = plt.figure(figsize=(12.6, 7.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    frame(fig, "V1  Refined columns  —  the current layout, rebuilt")
    # reusability gradient rule
    xs = np.linspace(0.06, 0.94, 300)
    for i in range(299):
        t = i / 298
        col = np.array(matplotlib.colors.to_rgb(C[1])) * (1 - t) + np.array(matplotlib.colors.to_rgb(C[4])) * t
        ax.plot(xs[i:i + 2], [0.905, 0.905], color=col, lw=3, solid_capstyle="butt")
    ax.text(0.06, 0.925, "MORE REUSABLE", fontsize=9, color=C[1], fontweight="bold", va="bottom")
    ax.text(0.94, 0.925, "LESS REUSABLE", fontsize=9, color=C[4], fontweight="bold", va="bottom", ha="right")
    W, G, x0, top, bot = 0.205, 0.022, 0.06, 0.875, 0.10
    for i, t in enumerate((1, 2, 3, 4)):
        x = x0 + i * (W + G)
        ax.add_patch(FancyBboxPatch((x, bot), W, top - bot, boxstyle="round,pad=0,rounding_size=0.012",
                                    facecolor="white", edgecolor=RULE, lw=1.0))
        ax.add_patch(Rectangle((x, top - 0.115), W, 0.115, facecolor=C[t], edgecolor="none"))
        ax.text(x + W / 2, top - 0.035, f"TIER {t}", ha="center", va="center", fontsize=14,
                fontweight="bold", color="white")
        ax.text(x + W / 2, top - 0.083, TITLE[t], ha="center", va="center", fontsize=9.3,
                color="white", linespacing=1.25)
        # icon row
        iy = top - 0.19; r = 0.0185
        need = {1: (1, 1, 1), 2: (1, 0, 0), 3: (0, 0, 0), 4: (0, 0, 0)}[t]
        for k, (fn, lab) in enumerate([(icon_disease, "Disease"), (icon_age, "Age"), (icon_sex, "Sex")]):
            ix = x + W * (0.22 + 0.28 * k)
            ax.add_patch(Circle((ix, iy), r * 1.55, facecolor=CL[t] if need[k] else "#F3F5F7",
                                edgecolor=C[t] if need[k] else GREY, lw=1.2,
                                linestyle="-" if need[k] else (0, (2, 2))))
            fn(ax, ix, iy, r * 0.75, C[t], on=bool(need[k]))
            ax.text(ix, iy - r * 2.05, lab, ha="center", va="top", fontsize=7.6, color=INK)
            if not need[k]:
                mark(ax, ix + r * 1.15, iy - r * 1.15, False, r=0.0085)
        # mock record
        ry0 = top - 0.30
        ax.text(x + W / 2, ry0 + 0.008, "Repository sample record", ha="center", va="bottom",
                fontsize=7.8, color=C[t], fontweight="bold")
        rows = {1: [("CD_treated_01", "34", "F", "Crohn's"), ("CD_placebo_02", "41", "M", "Crohn's"), ("Control_03", "29", "F", "Healthy")],
                2: [("SRS123456", "—", "—", "Crohn's"), ("SRS123457", "—", "—", "Crohn's"), ("SRS123458", "—", "—", "Healthy")],
                3: [("SRS0001", "—", "—", "—"), ("SRS0002", "—", "—", "—"), ("SRS0003", "—", "—", "—")],
                4: None}[t]
        if rows:
            cw = [0.42, 0.16, 0.16, 0.26]; rh = 0.036
            cx = x + W * 0.05; tw = W * 0.90
            heads = ["Sample", "Age", "Sex", "Disease"]
            ax.add_patch(Rectangle((cx, ry0 - rh), tw, rh, facecolor="#F3F5F7", edgecolor=RULE, lw=0.6))
            xx = cx
            for h_, w_ in zip(heads, cw):
                ax.text(xx + tw * w_ / 2, ry0 - rh / 2, h_, ha="center", va="center", fontsize=6.8,
                        color=INK, fontweight="bold"); xx += tw * w_
            for j, row in enumerate(rows):
                yy = ry0 - rh * (j + 2)
                ax.add_patch(Rectangle((cx, yy), tw, rh, facecolor="white", edgecolor=RULE, lw=0.6))
                xx = cx
                for k, (v, w_) in enumerate(zip(row, cw)):
                    if v == "—":
                        mark(ax, xx + tw * w_ / 2, yy + rh / 2, False, r=0.0075)
                    else:
                        ax.text(xx + tw * w_ / 2, yy + rh / 2, v, ha="center", va="center",
                                fontsize=6.6, color=INK if k else MUTED)
                    xx += tw * w_
        else:
            # tier 4: the broken chain
            yy = ry0 - 0.08
            for k, (fn, lab) in enumerate([(icon_id, "No sample IDs"), (icon_db, "No / invalid\naccession")]):
                ix = x + W * (0.3 + 0.4 * k)
                fn(ax, ix, yy, 0.02, GREY, on=False)
                mark(ax, ix + 0.024, yy + 0.02, False, r=0.0085)
                ax.text(ix, yy - 0.042, lab, ha="center", va="top", fontsize=7.4, color=MUTED, linespacing=1.2)
        # verdict
        ax.text(x + W / 2, bot + 0.135, VERDICT[t], ha="center", va="center", fontsize=9.6,
                fontweight="bold", color=C[t])
        # footer stats
        ax.add_patch(Rectangle((x, bot), W, 0.095, facecolor=CL[t], edgecolor="none"))
        ax.text(x + W * 0.07, bot + 0.062, f"{100*st[t]/N:.1f}%", ha="left", va="center",
                fontsize=17, fontweight="bold", color=C[t])
        ax.text(x + W * 0.07, bot + 0.026, f"{st[t]:,} studies", ha="left", va="center",
                fontsize=8.2, color=INK)
        ax.text(x + W * 0.93, bot + 0.062, f"{100*sa[t]/S:.1f}%", ha="right", va="center",
                fontsize=11, fontweight="bold", color=MUTED)
        ax.text(x + W * 0.93, bot + 0.026, f"{sa[t]/1000:,.0f}k samples", ha="right", va="center",
                fontsize=8.2, color=MUTED)
    caption(fig, f"BHAM reusability tiers. Per tier, the metadata a repository record must carry, what the record looks like, and the share of the "
                 f"{N:,} analysed studies (left) and {S/1000:,.0f}k samples (right) that land there. "
                 f"Dashed icons are absent fields. Tiers 1–2 are reusable for disease analysis; Tiers 3–4 ({100*(st[3]+st[4])/N:.1f}% of studies, "
                 f"{100*(sa[3]+sa[4])/S:.1f}% of samples) are not.")
    return fig


# --------------------------------------------------------------- design 2 --
def v2_matrix(st, sa, N, S):
    fig = plt.figure(figsize=(12.6, 7.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    frame(fig, "V2  Requirement matrix  —  what each tier guarantees, at a glance")
    reqs = ["Public\naccession", "Unique\nsample IDs", "Per-sample\ndisease", "Per-sample\nage", "Per-sample\nsex"]
    icons = [icon_db, icon_id, icon_disease, icon_age, icon_sex]
    have = {1: (1, 1, 1, 1, 1), 2: (1, 1, 1, 0, 0), 3: (1, 1, 0, 0, 0), 4: (0, 0, 0, 0, 0)}
    x0, colw = 0.29, 0.076; rows_y = [0.70, 0.565, 0.43, 0.295]; rh = 0.105
    # header icons
    for k, (fn, lab) in enumerate(zip(icons, reqs)):
        cx = x0 + colw * (k + 0.5)
        fn(ax, cx, 0.85, 0.016, INK)
        ax.text(cx, 0.812, lab, ha="center", va="top", fontsize=8, color=INK, linespacing=1.2)
    bx, bw = x0 + colw * 5 + 0.035, 0.085
    bx2 = bx + bw + 0.085
    ax.text(bx + bw / 2, 0.815, "Share of studies", ha="center", va="center", fontsize=8.6, fontweight="bold", color=INK)
    ax.text(bx2 + bw / 2, 0.815, "Share of samples", ha="center", va="center", fontsize=8.6, fontweight="bold", color=INK)
    ax.plot([0.05, 0.97], [0.768, 0.768], color=INK, lw=0.9)
    for t, y in zip((1, 2, 3, 4), rows_y):
        ax.add_patch(Rectangle((0.05, y - rh / 2), 0.92, rh, facecolor=CL[t] if t in (1, 2) else "#FAFBFC", edgecolor="none"))
        ax.add_patch(Rectangle((0.05, y - rh / 2), 0.012, rh, facecolor=C[t], edgecolor="none"))
        ax.text(0.075, y + 0.022, f"Tier {t}", fontsize=13, fontweight="bold", color=C[t], va="center")
        ax.text(0.075, y - 0.014, TITLE[t].replace("\n", " "), fontsize=8.3, color=MUTED, va="center")
        for k in range(5):
            cx = x0 + colw * (k + 0.5)
            if have[t][k]:
                mark(ax, cx, y, True, r=0.017, fill=C[t])
            else:
                ax.add_patch(Circle((cx, y), 0.017, facecolor="white", edgecolor=GREY, lw=1.1, linestyle=(0, (2, 2))))
        # study bar and sample bar, each with its own label; verdict lives under the tier title
        ax.add_patch(Rectangle((bx, y - 0.014), bw, 0.028, facecolor="#EEF1F4", edgecolor="none"))
        ax.add_patch(Rectangle((bx, y - 0.014), bw * st[t] / N, 0.028, facecolor=C[t], edgecolor="none"))
        ax.text(bx + bw + 0.007, y, f"{100*st[t]/N:.1f}%  ({st[t]:,})", va="center", fontsize=8.2, color=INK)
        ax.add_patch(Rectangle((bx2, y - 0.014), bw, 0.028, facecolor="#EEF1F4", edgecolor="none"))
        ax.add_patch(Rectangle((bx2, y - 0.014), bw * sa[t] / S, 0.028, facecolor=C[t], alpha=0.55, edgecolor="none"))
        ax.text(bx2 + bw + 0.007, y, f"{100*sa[t]/S:.1f}%  ({sa[t]/1000:,.0f}k)", va="center", fontsize=8.2, color=INK)
        ax.text(0.075, y - 0.041, VERDICT[t], fontsize=8.3, va="center", fontweight="bold",
                color=OK if t in (1, 2) else (MUTED if t == 3 else NO))
    ax.plot([0.05, 0.97], [rows_y[-1] - rh / 2 - 0.01, rows_y[-1] - rh / 2 - 0.01], color=INK, lw=0.9)
    # bracket
    ax.annotate("", xy=(0.03, rows_y[0] + rh / 2), xytext=(0.03, rows_y[1] - rh / 2),
                arrowprops=dict(arrowstyle="-", color=OK, lw=2.2))
    ax.text(0.022, (rows_y[0] + rows_y[1]) / 2, f"reusable\n{100*(st[1]+st[2])/N:.1f}% of studies",
            rotation=90, ha="right", va="center", fontsize=7.8, color=OK, linespacing=1.2)
    ax.annotate("", xy=(0.03, rows_y[2] + rh / 2), xytext=(0.03, rows_y[3] - rh / 2),
                arrowprops=dict(arrowstyle="-", color=NO, lw=2.2))
    ax.text(0.022, (rows_y[2] + rows_y[3]) / 2, f"not reusable\n{100*(sa[3]+sa[4])/S:.1f}% of samples",
            rotation=90, ha="right", va="center", fontsize=7.8, color=NO, linespacing=1.2)
    caption(fig, "Each row is a tier; each column a requirement checked in the repository record. A filled circle means the tier guarantees "
                 "that element for every sample. Bars give the tier's share of the 2,896 analysed studies (solid) and of the 752k deposited "
                 "samples (faded). The five requirements nest: every Tier 1 study also satisfies Tiers 2 and 3.")
    return fig


# --------------------------------------------------------------- design 3 --
def v3_funnel(st, sa, N, S):
    fig = plt.figure(figsize=(12.6, 7.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    frame(fig, "V3  The funnel  —  where the corpus is lost")
    stages = [("All analysed studies", N, S, INK, ""),
              ("Deposited with a public accession", st[1] + st[2] + st[3], sa[1] + sa[2] + sa[3], C[3],
               f"Tier 4 falls out: {st[4]:,} studies, {sa[4]/1000:,.0f}k samples — {DESC[4].replace(chr(10),' ')}"),
              ("… with per-sample disease", st[1] + st[2], sa[1] + sa[2], C[2],
               f"Tier 3 falls out: {st[3]:,} studies, {sa[3]/1000:,.0f}k samples — reads are there, biology is not"),
              ("… with disease, age and sex", st[1], sa[1], C[1],
               f"Tier 2 falls out: {st[2]:,} studies, {sa[2]/1000:,.0f}k samples — disease only")]
    ys = [0.78, 0.60, 0.42, 0.24]; H = 0.13; cx = 0.55; maxw = 0.50
    for i, (lab, n, s, col, lost) in enumerate(stages):
        w = maxw * n / N; y = ys[i]
        nxt = stages[i + 1][1] / N * maxw if i < 3 else None
        ax.add_patch(Rectangle((cx - w / 2, y - H / 2), w, H, facecolor=col, edgecolor="white", lw=1.5, alpha=0.92))
        if nxt is not None:
            ax.add_patch(Polygon([(cx - w / 2, y - H / 2), (cx + w / 2, y - H / 2),
                                  (cx + nxt / 2, ys[i + 1] + H / 2), (cx - nxt / 2, ys[i + 1] + H / 2)],
                                 closed=True, facecolor=col, alpha=0.18, edgecolor="none"))
        if w > 0.16:
            ax.text(cx, y + 0.02, f"{n:,} studies", ha="center", va="center", fontsize=13, fontweight="bold", color="white")
            ax.text(cx, y - 0.028, f"{s/1000:,.0f}k samples  ·  {100*n/N:.1f}% of studies", ha="center", va="center",
                    fontsize=8.6, color="white")
        else:   # too narrow to hold text: put it beside the band
            ax.text(cx + w / 2 + 0.012, y + 0.016, f"{n:,} studies", ha="left", va="center", fontsize=11, fontweight="bold", color=col)
            ax.text(cx + w / 2 + 0.012, y - 0.02, f"{s/1000:,.0f}k samples  ·  {100*n/N:.1f}%", ha="left", va="center", fontsize=8.4, color=MUTED)
        ax.text(0.27, y, lab, ha="right", va="center", fontsize=10.5, color=INK, fontweight="bold" if i == 0 else "normal")
        if lost:
            yy = (y - H / 2 + ys[i] - 0.18 + H / 2) / 2 if False else y - H / 2 - 0.025
            ax.text(0.27, yy, lost.split(" — ")[0], ha="right", va="center", fontsize=8.2, color=NO, style="italic")
            ax.text(0.27, yy - 0.024, lost.split(" — ")[1], ha="right", va="center", fontsize=7.6, color=MUTED, style="italic")
    ax.text(0.04, 0.078, f"Reusable for a disease-associated analysis: Tiers 1 + 2 = {st[1]+st[2]:,} studies ({100*(st[1]+st[2])/N:.1f}%), "
                         f"{(sa[1]+sa[2])/1000:,.0f}k samples ({100*(sa[1]+sa[2])/S:.1f}%).",
            fontsize=9.6, color=OK, fontweight="bold")
    ax.text(0.04, 0.050, f"Not reusable: {100*(sa[3]+sa[4])/S:.1f}% of all deposited samples.", fontsize=9.6, color=NO, fontweight="bold")
    fig.text(0.04, 0.035, "", fontsize=1)
    caption2 = lambda f, t: f.text(0.04, 0.012, t, fontsize=8.6, color=MUTED, va="bottom", ha="left", linespacing=1.45)
    caption2(fig, "Each band is the set of studies that survives one more requirement; widths are proportional to study count. What drops out at each "
                 "step is a BHAM tier. Tier 4 sample counts are imputed at the Tier 1–3 mean, so the top band's sample total is an estimate.")
    return fig


# --------------------------------------------------------------- design 4 --
def v4_cards(st, sa, N, S):
    fig = plt.figure(figsize=(12.6, 7.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    frame(fig, "V4  The record itself  —  four deposits, four fates")
    W, G, x0 = 0.215, 0.018, 0.045
    heads = ["sample_id", "host_disease", "age", "sex"]
    data = {1: [("CRC_case_01", "CRC", "58", "F"), ("CRC_case_02", "CRC", "61", "M"),
                ("HC_03", "control", "55", "F"), ("HC_04", "control", "49", "M")],
            2: [("SAMN0812", "adenoma", "", ""), ("SAMN0813", "adenoma", "", ""),
                ("SAMN0814", "control", "", ""), ("SAMN0815", "CRC", "", "")],
            3: [("SAMN2201", "", "", ""), ("SAMN2202", "", "", ""), ("SAMN2203", "", "", ""), ("SAMN2204", "", "", "")],
            4: [("sample", "", "", ""), ("sample", "", "", ""), ("sample", "", "", ""), ("sample", "", "", "")]}
    for i, t in enumerate((1, 2, 3, 4)):
        x = x0 + i * (W + G); top = 0.88; base = 0.30
        # card with shadow
        ax.add_patch(FancyBboxPatch((x + 0.004, base - 0.004), W, top - base, boxstyle="round,pad=0,rounding_size=0.01",
                                    facecolor="#E9EDF1", edgecolor="none"))
        ax.add_patch(FancyBboxPatch((x, base), W, top - base, boxstyle="round,pad=0,rounding_size=0.01",
                                    facecolor="white", edgecolor=RULE, lw=0.9))
        ax.add_patch(Rectangle((x, top - 0.012), W, 0.012, facecolor=C[t], edgecolor="none"))
        ax.text(x + 0.012, top - 0.045, f"Tier {t}", fontsize=13, fontweight="bold", color=C[t], va="center")
        ax.text(x + W - 0.012, top - 0.045, "BioProject PRJNA…" if t < 4 else "accession: —", fontsize=7.6,
                color=MUTED if t < 4 else NO, va="center", ha="right", family="monospace")
        # table
        ty = top - 0.085; rh = 0.05; cw = [0.36, 0.30, 0.15, 0.19]; tw = W - 0.024; cx = x + 0.012
        ax.add_patch(Rectangle((cx, ty - rh * 0.8), tw, rh * 0.8, facecolor="#F3F5F7", edgecolor=RULE, lw=0.5))
        xx = cx
        for h_, w_ in zip(heads, cw):
            ax.text(xx + 0.004, ty - rh * 0.4, h_, fontsize=6.6, color=INK, va="center", family="monospace"); xx += tw * w_
        for j, row in enumerate(data[t]):
            yy = ty - rh * 0.8 - rh * (j + 1)
            ax.add_patch(Rectangle((cx, yy), tw, rh, facecolor="white" if j % 2 else "#FBFCFD", edgecolor=RULE, lw=0.5))
            xx = cx
            for k, (v, w_) in enumerate(zip(row, cw)):
                if v:
                    col = INK
                    if t == 4 and k == 0: col = NO
                    ax.text(xx + 0.004, yy + rh / 2, v, fontsize=6.6, color=col, va="center", family="monospace")
                else:
                    ax.add_patch(Rectangle((xx + 0.004, yy + rh * 0.3), tw * w_ - 0.008, rh * 0.4, facecolor="#F0F2F4", edgecolor="none"))
                xx += tw * w_
        # what a machine can do with it
        vy = 0.455
        ax.plot([x + 0.012, x + W - 0.012], [vy + 0.05, vy + 0.05], color=RULE, lw=0.7)
        can = {1: [("disease", True), ("age", True), ("sex", True)], 2: [("disease", True), ("age", False), ("sex", False)],
               3: [("disease", False), ("age", False), ("sex", False)], 4: [("reads", False), ("identify samples", False), ("disease", False)]}[t]
        for k, (lab, ok) in enumerate(can):
            yy = vy + 0.02 - k * 0.03
            mark(ax, x + 0.03, yy, ok, r=0.009)
            ax.text(x + 0.05, yy, lab, fontsize=8, color=INK if ok else MUTED, va="center")
        ax.text(x + W / 2, 0.335, VERDICT[t], ha="center", fontsize=9.6, fontweight="bold",
                color=OK if t in (1, 2) else (MUTED if t == 3 else NO))
        # stats under the card
        ax.text(x + 0.012, 0.235, f"{100*st[t]/N:.1f}%", fontsize=18, fontweight="bold", color=C[t], va="center")
        ax.text(x + 0.012, 0.20, f"of studies  ({st[t]:,})", fontsize=8, color=MUTED, va="center")
        ax.text(x + W - 0.012, 0.235, f"{100*sa[t]/S:.1f}%", fontsize=12, fontweight="bold", color=MUTED, va="center", ha="right")
        ax.text(x + W - 0.012, 0.20, f"of samples  ({sa[t]/1000:,.0f}k)", fontsize=8, color=MUTED, va="center", ha="right")
    caption(fig, "The tier is nothing more than what the repository record looks like. Grey cells are fields a pipeline finds empty. The ticks below each card "
                 "are what an automated reader can recover for every sample without opening the paper. Tier 4 records either do not exist or cannot be "
                 "told apart.")
    return fig


# --------------------------------------------------------------- design 5 --
def v5_ladder(st, sa, N, S):
    fig = plt.figure(figsize=(12.6, 7.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    frame(fig, "V5  The staircase  —  each step adds one guarantee")
    steps = [(4, 0.06, 0.16), (3, 0.28, 0.30), (2, 0.50, 0.44), (1, 0.72, 0.58)]
    W = 0.215
    for t, x, h in steps:
        ax.add_patch(Rectangle((x, 0.16), W, h, facecolor=C[t], edgecolor="none", alpha=0.95))
        ax.add_patch(Rectangle((x, 0.16 + h - 0.012), W, 0.012, facecolor="white", alpha=0.25, edgecolor="none"))
        ax.text(x + 0.012, 0.16 + h - 0.035, f"Tier {t}", fontsize=14, fontweight="bold", color="white", va="center")
        ax.text(x + 0.012, 0.16 + h - 0.075, VERDICT[t], fontsize=8.6, color="white", va="center")
        # stats inside the step foot
        ax.text(x + W / 2, 0.215, f"{100*st[t]/N:.1f}%  ·  {st[t]:,} studies", ha="center", fontsize=8.8, color="white", fontweight="bold")
        ax.text(x + W / 2, 0.185, f"{100*sa[t]/S:.1f}% of samples", ha="center", fontsize=8, color="white", alpha=0.9)
        # what the step adds (above)
        add = {4: ("nothing usable", icon_db, False), 3: ("a public accession and\nunique sample IDs", icon_id, True),
               2: ("per-sample disease", icon_disease, True), 1: ("per-sample age and sex", icon_age, True)}[t]
        iy = 0.16 + h + 0.07
        add[1](ax, x + 0.03, iy, 0.017, C[t], on=add[2])
        tx = x + 0.06
        if t == 1:
            icon_sex(ax, x + 0.072, iy, 0.017, C[t]); tx = x + 0.105
        ax.text(tx, iy, ("+ " if t < 4 else "") + add[0], fontsize=8.6, color=INK, va="center", linespacing=1.2)
    # baseline & arrow
    ax.plot([0.04, 0.96], [0.16, 0.16], color=INK, lw=1.0)
    ax.annotate("", xy=(0.96, 0.10), xytext=(0.04, 0.10), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2, mutation_scale=14))
    ax.text(0.5, 0.075, r"increasing per-sample annotation  $\rightarrow$  increasing reusability by automated tools", ha="center",
            fontsize=8.8, color=MUTED)
    ax.text(0.04, 0.905, f"{N:,} studies  ·  {S/1000:,.0f}k samples", ha="left", fontsize=9, color=MUTED)
    ax.text(0.04, 0.875, f"{100*(sa[3]+sa[4])/S:.1f}% of samples never reach Tier 2", ha="left", fontsize=9.2, color=NO, fontweight="bold")
    caption(fig, "Read left to right: each step adds one guarantee to the repository record. A study stands on the highest step whose guarantee its "
                 "deposit meets for every sample. Step heights are schematic; the percentages are the measured shares of studies and samples.")
    return fig


# --------------------------------------------------------------- design 6 --
def v6_rings(st, sa, N, S):
    fig = plt.figure(figsize=(12.6, 7.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    frame(fig, "V6  Nested rings  —  tiers as nested requirements, sized by the data")
    # nested: a Tier 1 study satisfies everything -> innermost disc; areas by cumulative study share
    cum = [(1, st[1]), (2, st[1] + st[2]), (3, st[1] + st[2] + st[3]), (4, N)]
    axr = fig.add_axes([0.04, 0.12, 0.52, 0.80]); axr.set_aspect("equal"); axr.axis("off")
    axr.set_xlim(-1.05, 1.05); axr.set_ylim(-1.05, 1.05)
    R = 1.0
    for t, n in reversed(cum):
        axr.add_patch(Circle((0, 0), R * np.sqrt(n / N), facecolor=C[t], edgecolor="white", lw=1.8, alpha=0.96))
    # anchor points in ring middles, mapped back to figure coordinates for the leaders
    def fig_xy(r, ang):
        px, py = r * np.cos(ang), r * np.sin(ang)
        return fig.transFigure.inverted().transform(axr.transData.transform((px, py)))
    cx, cy = 0.30, 0.52
    # labels with leaders
    lab = {1: (f"Tier 1 — disease + age + sex", st[1], sa[1]), 2: (f"Tier 2 — disease", st[2], sa[2]),
           3: (f"Tier 3 — accession, IDs only", st[3], sa[3]), 4: (f"Tier 4 — no usable deposit", st[4], sa[4])}
    ys = {1: 0.50, 2: 0.615, 3: 0.73, 4: 0.845}
    for t in (1, 2, 3, 4):
        r_out = R * np.sqrt(cum[t - 1][1] / N)
        r_in = R * np.sqrt(cum[t - 2][1] / N) if t > 1 else 0
        rr = (r_out + r_in) / 2 if t > 1 else 0
        px, py = fig_xy(rr, np.pi / 4)
        ax.plot([px, 0.62], [py, ys[t]], color=MUTED, lw=0.7)
        ax.add_patch(Ellipse((px, py), 0.006, 0.006 * 12.6 / 7.4, facecolor=INK, edgecolor="none"))
        ax.add_patch(Rectangle((0.63, ys[t] - 0.02), 0.012, 0.04, facecolor=C[t], edgecolor="none"))
        ax.text(0.655, ys[t] + 0.012, lab[t][0], fontsize=10, fontweight="bold", color=INK, va="center")
        ax.text(0.655, ys[t] - 0.016, f"{lab[t][1]:,} studies ({100*lab[t][1]/N:.1f}%)   ·   {lab[t][2]/1000:,.0f}k samples ({100*lab[t][2]/S:.1f}%)",
                fontsize=8.4, color=MUTED, va="center")
    ax.text(0.655, 0.30, f"reusable  =  Tier 1 + Tier 2", fontsize=10.5, fontweight="bold", color=OK)
    ax.text(0.655, 0.265, f"{st[1]+st[2]:,} studies ({100*(st[1]+st[2])/N:.1f}%)  ·  {(sa[1]+sa[2])/1000:,.0f}k samples ({100*(sa[1]+sa[2])/S:.1f}%)",
            fontsize=8.6, color=OK)
    ax.text(0.655, 0.21, f"not reusable  =  Tier 3 + Tier 4", fontsize=10.5, fontweight="bold", color=NO)
    ax.text(0.655, 0.175, f"{st[3]+st[4]:,} studies ({100*(st[3]+st[4])/N:.1f}%)  ·  {(sa[3]+sa[4])/1000:,.0f}k samples ({100*(sa[3]+sa[4])/S:.1f}%)",
            fontsize=8.6, color=NO)
    caption(fig, "Tiers nest: every requirement of an outer ring is also met by the rings inside it. Ring area is proportional to the number of studies "
                 "that meet at least that tier's requirements, so the innermost disc is the fully annotated core and the outer band is the deposit-less majority.")
    return fig


def main():
    st, sa, N, S = numbers()
    os.makedirs(FIGS, exist_ok=True)
    designs = [v1_columns, v2_matrix, v3_funnel, v4_cards, v5_ladder, v6_rings]
    with PdfPages(OUT) as pdf:
        for i, fn in enumerate(designs, 1):
            fig = fn(st, sa, N, S)
            pdf.savefig(fig, bbox_inches="tight", facecolor="white")
            fig.savefig(os.path.join(FIGS, f"Figure1_v{i}.png"), dpi=170, bbox_inches="tight", facecolor="white")
            plt.close(fig)
    print(f"wrote {OUT}  ({len(designs)} designs)")
    print("studies:", st, "=", N); print("samples:", sa, "=", S)


if __name__ == "__main__":
    main()
