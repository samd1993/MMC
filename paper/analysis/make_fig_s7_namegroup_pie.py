"""Figure S7: how much of Tier 1 and Tier 2 rests on a name or free-text field.

Four wedges: each of Tier 1 and Tier 2 split by whether the per-sample disease
label sits in a structured repository metadata column, or only in a name or
free-text field (sample_alias, sample_title, sample_description,
experiment_title/alias, run_alias, library_name, study_title, sample name/ID).

Both halves satisfy the tier definition as written. The distinction matters
because an automated pipeline reading the archive cannot decode "HC_24" against
"CD_8" without a human deciding what the prefixes mean, so the name-derived half
is reusable in principle and not in practice.

Tiers 1-3 are complete manual censuses, so this is measured, not sampled.

Writes Figs/MMC2_Fig_S6.{pdf,png}.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402

MASTER = os.path.join(ROOT, "MMC2_master_reviewed.tsv")
FIGS = os.path.join(ROOT, "Figs")
# tier colour, then a desaturated version of it for the name-derived wedge
COLORS = {("1", False): "#91bfdb", ("1", True): "#cfe3ee",
          ("2", False): "#fee090", ("2", True): "#fef3d0"}

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42,
                     "font.family": "sans-serif",
                     "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"]})


def main():
    m = pd.read_csv(MASTER, sep="\t", low_memory=False,
                    keep_default_na=False, dtype=str)
    m = m[m.final_tier.isin(["1", "2"]) & m.in_analysis.eq("yes")].copy()
    m["nd"] = m.name_derived_disease_flag == "yes"

    order = [("1", False), ("1", True), ("2", True), ("2", False)]
    vals, labels, cols = [], [], []
    for t, nd in order:
        n = int(((m.final_tier == t) & (m.nd == nd)).sum())
        vals.append(n)
        labels.append(f"Tier {t}\n{'name / free text' if nd else 'metadata column'}\n"
                      f"n = {n}")
        cols.append(COLORS[(t, nd)])
    total = sum(vals)

    fig, ax = plt.subplots(figsize=(7.4, 6.4))
    wedges, _txt, auto = ax.pie(
        vals, colors=cols, startangle=90, counterclock=False,
        autopct=lambda p: f"{p:.0f}%", pctdistance=0.72,
        wedgeprops=dict(linewidth=0), textprops=dict(fontsize=11))
    for a in auto:
        a.set_fontweight("bold")
    ax.legend(wedges, labels, loc="center left", bbox_to_anchor=(1.0, 0.5),
              frameon=False, fontsize=10, labelspacing=1.1)
    ax.set_title(
        f"Reusable studies (Tier 1 and Tier 2), n = {total:,}\n"
        f"{100 * m.nd.mean():.0f}% carry the disease label only in a name or "
        f"free-text field",
        fontsize=11.5, loc="center", pad=16)
    ax.axis("equal")
    fig.subplots_adjust(left=0.02, right=0.62, top=0.86, bottom=0.04)

    os.makedirs(FIGS, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"MMC2_Fig_S6.{ext}"), dpi=300)
    plt.close(fig)
    print(f"wrote {FIGS}/MMC2_Fig_S6.pdf / .png")
    for (t, nd), v in zip(order, vals):
        print(f"  Tier {t} {'name/free-text' if nd else 'metadata column'}: "
              f"{v:>4}  ({100 * v / total:4.1f}% of reusable)")
    print(f"  total reusable {total:,}; name-derived overall "
          f"{100 * m.nd.mean():.1f}%")


if __name__ == "__main__":
    main()
