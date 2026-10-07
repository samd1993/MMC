"""MMC2 figures: metadata presence, antibiotic labeling, tier distribution.

Colour assignment follows the data's job rather than MMC1's per-series hues:
the x-axis already carries the metadata dimension, so colour encodes the
*dataset* (MMC1 vs MMC2). That keeps the categorical palette at two hues and
avoids MMC1's own red/brown pair, which is only ΔE 2.4 apart under protanopia
and is separable in Figure 2E solely because of marker shape.

Palettes were checked with the dataviz validator (all six checks pass):
  MMC1 vs MMC2        #377eb8 / #d73027
  antibiotic 3-state  #238b45 / #d94801 / #6a51a3

Every bar is directly labelled, so identity never rests on colour alone.

Writes Figs/MMC2_Fig2C_presence.{pdf,png}, Figs/MMC2_antibiotic_labeling.{pdf,png},
Figs/MMC2_tier_distribution.{pdf,png} and MMC2_figure_source_data.tsv.
"""

import io
import os
import sys
import tarfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import MMC, ROOT  # noqa: E402

STUDY = os.path.join(ROOT, "MMC2_study_metadata.tsv")
FIGS = os.path.join(ROOT, "Figs")
SRC = os.path.join(ROOT, "MMC2_figure_source_data.tsv")

C_MMC1, C_MMC2 = "#377eb8", "#d73027"
C_ABX = {"yes_per_sample": "#238b45", "all_participants": "#6a51a3", "no": "#d94801"}
INK, MUTED, GRID = "#1a1a1a", "#4d4d4d", "#d9d9d9"

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 10, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.edgecolor": GRID, "axes.linewidth": 0.8,
})


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, lw=0.7)
    ax.xaxis.grid(False)


def label_bars(ax, bars, vals, fmt="{:.0f}%"):
    for b, v in zip(bars, vals):
        ax.annotate(fmt.format(v), (b.get_x() + b.get_width() / 2, b.get_height()),
                    ha="center", va="bottom", fontsize=8.5, color=INK,
                    xytext=(0, 2), textcoords="offset points")


def load_mmc1():
    tf = tarfile.open(os.path.join(MMC, "MMC_resubmission_handoff.tar.gz"))
    return pd.read_csv(
        io.BytesIO(tf.extractfile("MMC_resubmission_handoff/MMC1_study_data_final.tsv").read()),
        sep="\t", low_memory=False, keep_default_na=False, na_values=[],
    )


def main():
    os.makedirs(FIGS, exist_ok=True)
    df = pd.read_csv(STUDY, sep="\t", low_memory=False)
    m1 = load_mmc1()

    # MMC1 reference rates (its 2,046 studies)
    m1_acc = 100 * (~m1["AccessionCode"].astype(str).str.strip().str.upper()
                    .isin(["", "N/A", "NA"])).mean()
    m1_dis = 100 * m1["Tier"].astype(int).isin([1, 2]).mean()
    m1_age = 100 * m1["Age"].eq("Yes").mean()
    m1_sex = 100 * m1["Sex"].eq("Yes").mean()

    # MMC2 rates over the assessable denominator (non-INSDC excluded).
    # `tier` round-trips through TSV as float with NaN for the excluded rows, so
    # test for null -- comparing str(nan) to "" silently keeps all 175 of them.
    df["tier"] = df["tier"].map(lambda v: "" if pd.isna(v) else str(int(float(v))))
    assess = df[df["tier"] != ""]
    ins = df[df["metadata_source"] == "insdc"]
    m2_acc = 100 * assess["metadata_source"].eq("insdc").mean()
    m2 = {
        "Accession": m2_acc,
        "Disease": 100 * assess["disease_evidence"].mean(),
        "Age": 100 * assess["age_present"].mean(),
        "Sex": 100 * assess["sex_present"].mean(),
        "Antibiotics": 100 * assess["antibiotic_present"].mean(),
        "Geography": 100 * assess["geography_present"].mean(),
    }
    m1v = {"Accession": m1_acc, "Disease": m1_dis, "Age": m1_age, "Sex": m1_sex,
           "Antibiotics": None, "Geography": None}

    # ---------------- Figure 2C-style presence panel ----------------------
    cats = list(m2)
    x = range(len(cats))
    w = 0.38
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    fig.subplots_adjust(left=0.09, right=0.98, top=0.88, bottom=0.14)

    b1 = ax.bar([i - w / 2 for i in x], [m1v[c] or 0 for c in cats], w,
                color=C_MMC1, label=f"MMC1 (n={len(m1):,})",
                edgecolor="white", linewidth=1.4)
    b2 = ax.bar([i + w / 2 for i in x], [m2[c] for c in cats], w,
                color=C_MMC2, label=f"MMC2 (n={len(assess):,})",
                edgecolor="white", linewidth=1.4)
    label_bars(ax, [b for b, c in zip(b1, cats) if m1v[c] is not None],
               [m1v[c] for c in cats if m1v[c] is not None])
    label_bars(ax, b2, [m2[c] for c in cats])
    for i, c in enumerate(cats):
        if m1v[c] is None:
            ax.annotate("not\nmeasured", (i - w / 2, 1.5), ha="center", va="bottom",
                        fontsize=7.5, color=MUTED, style="italic")

    ax.set_xticks(list(x))
    ax.set_xticklabels(cats)
    ax.set_ylabel("Studies with metadata present (%)")
    ax.set_ylim(0, max(max(m2.values()), m1_acc) * 1.22)
    ax.set_title("Per-sample metadata present in the repository",
                 loc="left", fontsize=12, color=INK, pad=12)
    ax.legend(frameon=False, loc="upper right", fontsize=9)
    style(ax)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"MMC2_Fig2C_presence.{ext}"), dpi=300)
    plt.close(fig)

    # ---------------- antibiotic labeling ---------------------------------
    inv = df[df["abx_involves_study"]]
    order = ["yes_per_sample", "all_participants", "no"]
    lbl = {"yes_per_sample": "Labelled per sample",
           "all_participants": "All participants exposed\n(labelling uninformative)",
           "no": "Not labelled"}
    counts = [int((inv["abx_labeling"] == k).sum()) for k in order]
    pcts = [100 * c / len(inv) if len(inv) else 0 for c in counts]

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    fig.subplots_adjust(left=0.30, right=0.97, top=0.84, bottom=0.13)
    bars = ax.barh([lbl[k] for k in order], pcts,
                   color=[C_ABX[k] for k in order], height=0.6,
                   edgecolor="white", linewidth=1.4)
    for b, c, p in zip(bars, counts, pcts):
        ax.annotate(f"{c:,}  ({p:.1f}%)", (b.get_width(), b.get_y() + b.get_height() / 2),
                    va="center", ha="left", fontsize=9, color=INK,
                    xytext=(4, 0), textcoords="offset points")
    ax.invert_yaxis()
    ax.set_xlabel("Share of studies involving antibiotics (%)")
    ax.set_xlim(0, max(pcts) * 1.28 if any(pcts) else 1)
    ax.set_title(f"Antibiotic exposure labelling\n{len(inv):,} of {len(df):,} studies "
                 f"involve antibiotics", loc="left", fontsize=12, color=INK, pad=12)
    style(ax)
    ax.xaxis.grid(True, color=GRID, lw=0.7)
    ax.yaxis.grid(False)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"MMC2_antibiotic_labeling.{ext}"), dpi=300)
    plt.close(fig)

    # ---------------- tier distribution -----------------------------------
    t1 = m1["Tier"].astype(int).value_counts(normalize=True).mul(100)
    t2 = assess["tier"].astype(str).value_counts(normalize=True).mul(100)
    tiers = ["1", "2", "3", "4"]
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.13)
    xa = range(len(tiers))
    b1 = ax.bar([i - w / 2 for i in xa], [t1.get(int(t), 0) for t in tiers], w,
                color=C_MMC1, label=f"MMC1 (n={len(m1):,})", edgecolor="white", linewidth=1.4)
    b2 = ax.bar([i + w / 2 for i in xa], [t2.get(t, 0) for t in tiers], w,
                color=C_MMC2, label=f"MMC2 (n={len(assess):,})", edgecolor="white", linewidth=1.4)
    label_bars(ax, b1, [t1.get(int(t), 0) for t in tiers], "{:.1f}%")
    label_bars(ax, b2, [t2.get(t, 0) for t in tiers], "{:.1f}%")
    ax.set_xticks(list(xa))
    ax.set_xticklabels([f"Tier {t}" for t in tiers])
    ax.set_ylabel("Studies (%)")
    ax.set_title("Reusability tier distribution", loc="left", fontsize=12, color=INK, pad=12)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    style(ax)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"MMC2_tier_distribution.{ext}"), dpi=300)
    plt.close(fig)

    # ---------------- source data ------------------------------------------
    src = [{"panel": "presence", "category": c, "MMC1_pct": m1v[c],
            "MMC2_pct": round(m2[c], 2)} for c in cats]
    src += [{"panel": "antibiotic", "category": lbl[k].replace("\n", " "),
             "MMC1_pct": None, "MMC2_pct": round(p, 2), "n": c}
            for k, c, p in zip(order, counts, pcts)]
    src += [{"panel": "tier", "category": f"Tier {t}",
             "MMC1_pct": round(t1.get(int(t), 0), 2),
             "MMC2_pct": round(t2.get(t, 0), 2)} for t in tiers]
    pd.DataFrame(src).to_csv(SRC, sep="\t", index=False)

    print(f"wrote 3 figures to {FIGS}/")
    print(f"wrote {SRC}")
    print(f"\npresence (MMC2, n={len(assess):,} assessable):")
    for c in cats:
        ref = f"   MMC1 {m1v[c]:.1f}%" if m1v[c] is not None else "   MMC1 n/a"
        print(f"  {c:12s} {m2[c]:5.1f}%{ref}")
    print(f"\nantibiotic labeling among {len(inv):,} antibiotic studies:")
    for k, c, p in zip(order, counts, pcts):
        print(f"  {k:17s} {c:4,d} ({p:5.1f}%)")


if __name__ == "__main__":
    main()
