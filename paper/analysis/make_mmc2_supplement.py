"""MMC2 supplement — Tables S1-S3 and Figures S1-S5, mirroring MMC1's.

  Table S1  study counts by tier x body site x sequencing type
  Table S2  per-journal tier rates for journals with >=10 studies, + impact factor
  Table S3  study counts by disease group
  Fig S1    study counts by country (log-scaled bar; MMC1 used a choropleth,
            which needs a GeoJSON download -- a ranked bar carries the same
            information without a network dependency)
  Fig S2    tier composition, (A) top-30 journals by impact factor,
            (B) continents
  Fig S3    journal impact factor vs reusability, (A) Tier 1-2 and (B) Tier 4
  Fig S4    chi-square page: standardized residuals of tier against
            (A) continent, (B) disease group, (C) body site, (D) hemisphere

Fig S3 uses `impact_factor` (Clarivate) ONLY. The mixed
`impact_factor_combined` column shifts the Tier 1-2 odds ratio from 1.03 to 1.51
and toward significance -- a metric artifact, since `jif_like` is a rank proxy
running ~1.2x higher. Fig S2 uses the combined column, where the metric only
orders journals on an axis and the extra coverage is a real gain.

Writes Figs/MMC2_Table_S1..S3.pdf, Figs/MMC2_Fig_S1..S5.{pdf,png},
(S5 is the tier-over-time panel, formerly S6),
MMC2_supplement_stats.md and the underlying TSVs.
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from scipy.special import expit  # noqa: E402
from scipy.stats import chi2_contingency  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402

ATTR = os.path.join(ROOT, "MMC2_study_attributes.tsv")
STUDY = os.path.join(ROOT, "MMC2_study_metadata.tsv")
FIGS = os.path.join(ROOT, "Figs")
STATS = os.path.join(ROOT, "MMC2_supplement_stats.md")

TIER_COLORS = {1: "#91bfdb", 2: "#fee090", 3: "#fc8d59", 4: "#d73027"}
TIERS = [1, 2, 3, 4]
TLAB = {1: "Tier 1", 2: "Tier 2", 3: "Tier 3", 4: "Tier 4"}
SITES = ["Gut", "Multiple", "Pulmonary", "Vaginal", "Oral", "Nasal", "Skin", "Other"]
SEQ = ["Amplicon", "Shotgun", "Both"]
# Chi-square subgroups are kept only where the group has at least this many
# studies; Sam, 4 Sep. Below it the expected counts are too thin to test.
MIN_GROUP = 20
NORTH = {"North America", "Europe", "Asia", "Oceania", "Russia"}
SOUTH = {"Southeast Asia", "Middle East", "South America", "Africa"}
INK, MUTED, GRID = "#1a1a1a", "#4d4d4d", "#d9d9d9"
# diverging ramp for standardized residuals: two hues, neutral grey midpoint
DIVERGE = "RdBu_r"

# Helvetica throughout; titles are dropped from the plots because the numbered
# captions in the Supplementary Information already carry them.
plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 9,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.edgecolor": GRID, "figure.dpi": 120,
})
L = []


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_axisbelow(True)


def table_pdf(path, df, title, note, colw=None):
    h = max(3.0, 0.28 * len(df) + 2.2)
    fig, ax = plt.subplots(figsize=(11, h))
    ax.axis("off")
    t = ax.table(cellText=df.astype(str).values, colLabels=df.columns,
                 cellLoc="center", loc="upper center",
                 colWidths=colw or [1 / len(df.columns)] * len(df.columns))
    t.auto_set_font_size(False)
    t.set_fontsize(7.5)
    t.scale(1, 1.25)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor("#e8e8e8")
        if r == 0:
            cell.set_facecolor("#f2f2f2")
            cell.set_text_props(weight="bold")
    ax.set_title(title, loc="left", fontsize=11, pad=18)
    fig.text(0.06, 0.045, note, fontsize=7.5, color=MUTED, wrap=True, va="top")
    fig.subplots_adjust(top=0.90, bottom=0.14)
    with PdfPages(path) as pp:
        pp.savefig(fig)
    plt.close(fig)


def pval(p, lo=0.001):
    """Bare p-value: '<0.001' or '0.021'. Never '0.000'."""
    return f"<{lo}" if p < lo else f"{p:.3f}"


def pfmt(p, lo=0.001):
    """p-values print as 'p<0.001', never as 'p=<0.001' or 'p=0.000'."""
    return f"p<{lo}" if p < lo else f"p={p:.3f}"


def simulated_p(ct, n_sim=20000, seed=0):
    """Monte-Carlo p for a contingency table, for use when expected counts are
    thin enough that the chi-square approximation is not trustworthy."""
    rng = np.random.default_rng(seed)
    obs = np.asarray(ct, dtype=float)
    stat = chi2_contingency(obs)[0]
    rows, cols = obs.sum(1), obs.sum(0)
    n = obs.sum()
    pool = np.repeat(np.arange(len(cols)), cols.astype(int))
    hits = 0
    for _ in range(n_sim):
        perm = rng.permutation(pool)
        tab = np.zeros_like(obs)
        start = 0
        for i, r in enumerate(rows.astype(int)):
            idx, cnt = np.unique(perm[start:start + r], return_counts=True)
            tab[i, idx] = cnt
            start += r
        keep = tab.sum(0) > 0
        if chi2_contingency(tab[:, keep])[0] >= stat - 1e-9:
            hits += 1
    return (hits + 1) / (n_sim + 1)


def residual_heatmap(ax, ct, title):
    chi2, p, dof, exp = chi2_contingency(ct)
    res = (ct.values - exp) / np.sqrt(exp)
    vmax = max(3.0, np.abs(res).max())
    im = ax.imshow(res, cmap=DIVERGE, vmin=-vmax, vmax=vmax, aspect="equal")
    ax.set_anchor("N")
    ax.set_xticks(range(ct.shape[1]))
    ax.set_xticklabels([TLAB.get(c, c) for c in ct.columns], fontsize=9)
    ax.set_yticks(range(ct.shape[0]))
    ax.set_yticklabels(ct.index, fontsize=9)
    for i in range(res.shape[0]):
        for j in range(res.shape[1]):
            v = res[i, j]
            ax.text(j, i, f"{v:+.1f}", ha="center", va="center", fontsize=8,
                    color="white" if abs(v) > vmax * 0.55 else INK)
    head = f"{title}\n" if title else ""
    ax.set_title(f"{head}χ²={chi2:.1f}, df={dof}, {pfmt(p)}",
                 loc="left", fontsize=10.5, pad=8)
    return im, chi2, p, dof, exp.min()


def main():
    os.makedirs(FIGS, exist_ok=True)
    at = pd.read_csv(ATTR, sep="\t", low_memory=False)
    sm_ = pd.read_csv(STUDY, sep="\t", low_memory=False)
    at = at.merge(sm_[["record_id", "metadata_source"]], on="record_id",
                  how="left", suffixes=("", "_s"))
    # Tier comes from the consolidated master so the supplement agrees with
    # Figure 2; the pipeline column is superseded wherever review moved a study.
    mr = pd.read_csv(os.path.join(ROOT, "MMC2_master_reviewed.tsv"), sep="\t",
                     low_memory=False, keep_default_na=False, dtype=str)
    at = at.merge(mr[["record_id", "final_tier", "in_analysis",
                      "disease_group_final"]], on="record_id", how="left")
    # MMC2_study_attributes.tsv carries a pre-review disease_group that
    # disagrees with the reviewed master for 29 studies, 13 of them still
    # labelled "Healthy" after review moved them to a real disease group. The
    # reviewed column wins. "Autoimmine" is a typo for "Autoimmune" in one
    # record and is normalised here rather than left to split the group.
    dgf = at["disease_group_final"].fillna("").replace(
        {"Autoimmine": "Autoimmune"})
    at["disease_group"] = dgf.where(dgf.ne(""), at["disease_group"])
    # Sam, 2 Sep: Healthy-category studies are out of every analysis, as are the
    # non-INSDC deposits and the withdrawn records
    at = at[at["in_analysis"].eq("yes")].copy()
    at["tier"] = pd.to_numeric(at["final_tier"], errors="coerce")
    d = at[at["tier"].notna()].copy()
    d["tier"] = d["tier"].astype(int)
    N = len(d)
    L.append(f"# MMC2 supplement — statistics\n\nAssessable studies: **{N:,}** "
             f"({int(at['tier'].isna().sum()):,} non-INSDC excluded).\n")

    # ---------------- Table S1 --------------------------------------------
    rows = []
    for s in SITES:
        g = d[d["body_site_group"] == s]
        if not len(g):
            continue
        rows.append({"Body site": s, "Sequencing": "— all —", "n": len(g),
                     **{TLAB[t]: f"{int((g.tier == t).sum())} "
                                 f"({100 * (g.tier == t).mean():.0f}%)" for t in TIERS}})
        for q in SEQ:
            gg = g[g["sequencing_type"] == q]
            if not len(gg):
                continue
            rows.append({"Body site": "", "Sequencing": q, "n": len(gg),
                         **{TLAB[t]: f"{int((gg.tier == t).sum())} "
                                     f"({100 * (gg.tier == t).mean():.0f}%)"
                            for t in TIERS}})
    t1 = pd.DataFrame(rows)
    t1.to_csv(os.path.join(ROOT, "MMC2_Table_S1.tsv"), sep="\t", index=False)
    table_pdf(os.path.join(FIGS, "MMC2_Table_S1.pdf"), t1,
              "Table S1. Study counts by reusability tier, body site and sequencing type.",
              f"Human-microbiome studies (n={N:,}) cross-tabulated by tier, body site "
              f"(8 collapsed groups) and sequencing approach. Percentages are within row. "
              f"Sequencing: {', '.join(f'{q} {100 * d.sequencing_type.eq(q).mean():.0f}%' for q in SEQ)}.")

    # ---------------- Table S2 --------------------------------------------
    vc = d["journal"].value_counts()
    big = vc[vc >= 10]
    rows = []
    for j in big.index:
        g = d[d["journal"] == j]
        rows.append({
            "Journal": str(j)[:52], "n": len(g),
            "Tier 1-2 %": f"{100 * g.tier.isin([1, 2]).mean():.1f}",
            "Tier 4 %": f"{100 * g.tier.eq(4).mean():.1f}",
            "IF": ("" if g["impact_factor_combined"].isna().all()
                   else f"{g['impact_factor_combined'].dropna().iloc[0]:.1f}"),
            "IF source": (g["impact_factor_source"].dropna().iloc[0]
                          if g["impact_factor_source"].notna().any() else ""),
        })
    t2 = pd.DataFrame(rows).sort_values("n", ascending=False)
    t2.to_csv(os.path.join(ROOT, "MMC2_Table_S2.tsv"), sep="\t", index=False)
    table_pdf(os.path.join(FIGS, "MMC2_Table_S2.pdf"), t2,
              f"Table S2. Reusability by journal ({len(big)} journals with >=10 studies).",
              "IF is Clarivate JIF where available, otherwise an OpenAlex-derived rank "
              "proxy (mean 2024 citations to 2022-23 research articles); the 'IF source' "
              "column records which. The proxy is NOT interchangeable with JIF and is not "
              "used for the Fig S3 regression.",
              colw=[0.42, 0.08, 0.14, 0.12, 0.10, 0.14])

    # ---------------- Table S3 --------------------------------------------
    dg = d["disease_group"].value_counts()
    t3 = pd.DataFrame([{"Disease group": k, "n": int(v),
                        "% of corpus": f"{100 * v / N:.1f}",
                        "Tier 1-2 %": f"{100 * d.loc[d.disease_group == k, 'tier'].isin([1, 2]).mean():.1f}"}
                       for k, v in dg.items()])
    t3.to_csv(os.path.join(ROOT, "MMC2_Table_S3.tsv"), sep="\t", index=False)
    table_pdf(os.path.join(FIGS, "MMC2_Table_S3.pdf"), t3,
              "Table S3. Study counts by disease group.",
              f"Disease groups assigned from MeSH disease headings. Counts sum to "
              f"{int(dg.sum()):,}.", colw=[0.46, 0.14, 0.20, 0.20])

    # ---------------- Fig S1: countries -----------------------------------
    cc = d["country"].replace("", np.nan).dropna().value_counts()
    top = cc.head(30)[::-1]
    fig, ax = plt.subplots(figsize=(8.2, 9))
    ax.barh(top.index, top.values, color="#377eb8", edgecolor="white", linewidth=0.8)
    for y, v in enumerate(top.values):
        ax.annotate(f"{v:,}", (v, y), va="center", ha="left", fontsize=8,
                    xytext=(3, 0), textcoords="offset points")
    ax.set_xscale("log")
    ax.set_xlabel("Studies (log scale)")
    style(ax)
    ax.xaxis.grid(True, color=GRID, lw=0.7)
    fig.subplots_adjust(left=0.28, right=0.95, top=0.94, bottom=0.07)
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"MMC2_Fig_S1.{e}"), dpi=300)
    plt.close(fig)
    two = cc.head(2)
    L.append(f"\n## Fig S1 — geography of study output\n\n{cc.size} countries. "
             f"Top two ({', '.join(two.index)}) account for "
             f"**{100 * two.sum() / cc.sum():.0f}%** of studies with a country.\n")

    # ---------------- Fig S2: tier composition, journals and continents ---
    # Panel B was Figure S5A. Its legend used to sit inside the axes on
    # "lower right" and covered the Tier 3/4 bars; both panels now share one
    # legend above the figure.
    jj = t2.copy()
    jj["IFv"] = pd.to_numeric(jj["IF"], errors="coerce")
    jj = jj.dropna(subset=["IFv"]).nlargest(30, "IFv").sort_values("IFv")
    ins_geo = d[d["metadata_source"] == "insdc"].copy()
    gc = ins_geo[ins_geo["continent"].replace("", np.nan).notna()].copy()
    corder = gc["continent"].value_counts().index.tolist()[::-1]

    fig, axs = plt.subplots(2, 1, figsize=(9.5, 11.6),
                            gridspec_kw={"height_ratios": [30, len(corder) + 2]})
    left = np.zeros(len(jj))
    for t in TIERS:
        vals = []
        for j in jj["Journal"]:
            g = d[d["journal"].astype(str).str.slice(0, 52) == j]
            vals.append(100 * g.tier.eq(t).mean() if len(g) else 0)
        axs[0].barh(range(len(jj)), vals, left=left, color=TIER_COLORS[t],
                    edgecolor="white", linewidth=0.7, label=TLAB[t])
        left += np.array(vals)
    axs[0].set_yticks(range(len(jj)))
    axs[0].set_yticklabels([f"{r.Journal[:40]}  ({r.IFv:.1f})" for r in jj.itertuples()],
                           fontsize=8)
    axs[0].set_xlabel("Proportion of studies (%)")
    axs[0].set_xlim(0, 100)
    axs[0].set_title("A", loc="left", fontsize=12, fontweight="bold", pad=6)
    style(axs[0])
    axs[0].xaxis.grid(True, color=GRID, lw=0.7)

    left = np.zeros(len(corder))
    for t in TIERS:
        vals = [100 * gc.loc[gc.continent == c, "tier"].eq(t).mean() for c in corder]
        axs[1].barh(range(len(corder)), vals, left=left, color=TIER_COLORS[t],
                    edgecolor="white", linewidth=0.7)
        left += np.array(vals)
    axs[1].set_yticks(range(len(corder)))
    axs[1].set_yticklabels(
        [f"{c}  (n={int((gc.continent == c).sum()):,})" for c in corder], fontsize=9)
    axs[1].set_xlabel("Proportion of studies (%)")
    axs[1].set_xlim(0, 100)
    axs[1].set_title("B", loc="left", fontsize=12, fontweight="bold", pad=6)
    style(axs[1])
    axs[1].xaxis.grid(True, color=GRID, lw=0.7)

    h, lb = axs[0].get_legend_handles_labels()
    fig.legend(h, lb, frameon=False, fontsize=9, ncol=4,
               loc="upper center", bbox_to_anchor=(0.55, 0.995))
    fig.subplots_adjust(left=0.40, right=0.97, top=0.955, bottom=0.05, hspace=0.16)
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"MMC2_Fig_S2.{e}"), dpi=300)
    plt.close(fig)
    L.append(f"\n## Fig S2 \u2014 tier composition by journal and continent\n\n"
             f"(A) the 30 journals with the highest impact factor. "
             f"(B) continents, accession-bearing studies only (n={len(gc):,}), "
             f"moved here from the geography figure.\n")

    # ---------------- Fig S3: IF vs reusability (Clarivate only) ----------
    cl = d[d["impact_factor"].notna() & (d["impact_factor"] > 0)].copy()
    vcl = cl["journal"].value_counts()
    cl = cl[cl["journal"].isin(vcl[vcl >= 10].index)].copy()
    cl["logIF"] = np.log10(cl["impact_factor"])
    cl["t12"] = cl.tier.isin([1, 2]).astype(int)
    cl["t4"] = cl.tier.eq(4).astype(int)
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2))
    L.append(f"\n## Fig S3 — impact factor vs reusability\n\n"
             f"Clarivate JIF only, {cl.journal.nunique()} journals with >=10 studies, "
             f"n={len(cl):,}.\n")
    for k_, (ax, col, name, mmc1) in enumerate(
            [(axs[0], "t12", "P(Tier 1-2)", "OR=0.99, p=0.96"),
             (axs[1], "t4", "P(Tier 4)", "OR=0.57, p=0.001")]):
        ax.text(-0.13, 1.14, "AB"[k_], transform=ax.transAxes, fontsize=12,
                fontweight="bold", va="top")
        m = smf.logit(f"{col} ~ logIF", data=cl).fit(disp=False)
        orr, p = np.exp(m.params["logIF"]), m.pvalues["logIF"]
        lo, hi = np.exp(m.conf_int().loc["logIF"])
        grid = np.linspace(cl.logIF.min(), cl.logIF.max(), 100)
        pred = expit(m.params["Intercept"] + m.params["logIF"] * grid)
        per_j = cl.groupby("journal").agg(x=("logIF", "first"), y=(col, "mean"),
                                         n=(col, "size"))
        ax.scatter(10 ** per_j.x, 100 * per_j.y, s=per_j.n * 1.6, color="#377eb8",
                   edgecolor="black", linewidth=0.5, alpha=0.75, zorder=3)
        ax.plot(10 ** grid, 100 * pred, color="#d73027", lw=2.4, zorder=2)
        ax.set_xscale("log")
        ax.set_xlabel("Journal impact factor (log scale)")
        ax.set_ylabel(f"{name} (%)")
        ax.set_title(f"{name}\nOR={orr:.2f} [{lo:.2f}–{hi:.2f}] per log10(IF), "
                     f"p={'<0.001' if p < 0.001 else f'{p:.3f}'}",
                     loc="left", fontsize=11, pad=10)
        style(ax)
        ax.yaxis.grid(True, color=GRID, lw=0.7)
        L.append(f"- {name}: OR={orr:.2f} [{lo:.2f}–{hi:.2f}], "
                 f"p={p:.3f}  (MMC1: {mmc1})\n")
    fig.subplots_adjust(left=0.07, right=0.98, top=0.85, bottom=0.16, wspace=0.25)
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"MMC2_Fig_S3.{e}"), dpi=300)
    plt.close(fig)

    # ---------------- Fig S4: the chi-square page -------------------------
    # Every test on this page is restricted to accession-bearing (INSDC)
    # studies. On the full corpus these tests are CONFOUNDED BY CONSTRUCTION:
    # a paper with an accession gets its body site and country from ENA, while
    # a paper without one gets them from the title, abstract or affiliation and
    # is Tier 4 by definition. The same variable drives both axes, manufacturing
    # chi2=202.5 for body site where the uniform-source stratum gives 74.7.
    ins = d[d["metadata_source"] == "insdc"].copy()

    def keep_big(frame, col):
        """Groups with at least MIN_GROUP studies, largest first."""
        vc = frame[col].replace("", np.nan).dropna().value_counts()
        return vc[vc >= MIN_GROUP].index.tolist()

    cont_keep = keep_big(ins, "continent")
    dis_keep = keep_big(ins, "disease_group")
    gcont = ins[ins["continent"].isin(cont_keep)]
    gdis = ins[ins["disease_group"].isin(dis_keep)]
    gsite = ins[ins["body_site_group"].notna()]
    ghemi = ins.assign(hemi=ins["continent"].map(
        lambda x: "Global North" if x in NORTH
        else ("Global South" if x in SOUTH else None))).dropna(subset=["hemi"])

    site_ct = pd.crosstab(pd.Categorical(gsite["body_site_group"], SITES), gsite["tier"])
    site_ct = site_ct.loc[[x for x in SITES if x in site_ct.index]]

    PANELS = [
        ("A", f"Continent  (\u2265{MIN_GROUP} studies, n={len(gcont):,})",
         pd.crosstab(gcont["continent"], gcont["tier"]).loc[cont_keep], "continent"),
        ("B", f"Disease group  (\u2265{MIN_GROUP} studies, n={len(gdis):,})",
         pd.crosstab(gdis["disease_group"], gdis["tier"]).loc[dis_keep], "disease group"),
        ("C", f"Body site  (n={len(gsite):,})", site_ct, "body site"),
        ("D", f"Hemisphere  (n={len(ghemi):,})",
         pd.crosstab(ghemi["hemi"], ghemi["tier"]), "hemisphere"),
    ]

    fig, axs = plt.subplots(2, 2, figsize=(11.2, 10.6))
    axs = axs.ravel()
    stats = []
    for ax, (letter, title, ct, name) in zip(axs, PANELS):
        im, c2, pv, dof, minexp = residual_heatmap(ax, ct, title)
        ax.text(-0.34, 1.13, letter, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="top")
        # The chi-square approximation needs every expected count above 5; below
        # that a permutation p is reported alongside rather than instead, so the
        # reader can see both.
        sim = simulated_p(ct) if minexp < 5 else None
        tail = ("" if sim is None else
                f"\nsimulated {pfmt(sim)}  (min expected {minexp:.1f})")
        ax.set_title(f"{title}\nχ²={c2:.1f}, df={dof}, {pfmt(pv)}{tail}",
                     loc="left", fontsize=10.2, pad=8)
        stats.append((name, len(ct), c2, dof, pv, minexp, sim))
    fig.subplots_adjust(left=0.155, right=0.875, top=0.90, bottom=0.05,
                        wspace=0.60, hspace=0.28)
    cax = fig.add_axes([0.915, 0.34, 0.016, 0.32])
    fig.colorbar(im, cax=cax, label="standardized residual")
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"MMC2_Fig_S4.{e}"), dpi=300)
    plt.close(fig)

    L.append(f"\n## Fig S4 \u2014 chi-square tests of reusability tier\n\n"
             f"All four tests are restricted to accession-bearing studies "
             f"(n={len(ins):,}) so that the grouping variable is assigned from "
             f"repository metadata throughout. Groups with fewer than "
             f"{MIN_GROUP} studies are dropped.\n\n"
             f"| Grouping | Groups | \u03c7\u00b2 | df | p | min expected | simulated p |\n"
             f"|---|---|---|---|---|---|---|\n")
    for name, k, c2, dof, pv, minexp, sim in stats:
        L.append(f"| {name} | {k} | {c2:.1f} | {dof} | "
                 f"{pval(pv)} | {minexp:.1f} | "
                 f"{'-' if sim is None else pval(sim)} |\n")
    # size of the assignment-source confound, recomputed every run
    def _full(frame, col, order=None):
        f = frame[frame[col].replace("", np.nan).notna()]
        idx = pd.Categorical(f[col], order) if order else f[col]
        c = pd.crosstab(idx, f["tier"])
        if order:
            c = c.loc[[x for x in order if x in c.index]]
        return chi2_contingency(c)[0]

    full_site = _full(d, "body_site_group", SITES)
    here_site = [x for x in stats if x[0] == "body site"][0][2]
    ins_mult = 100 * ins["body_site_group"].eq("Multiple").mean()
    oth_mult = 100 * d.loc[d["metadata_source"] != "insdc",
                           "body_site_group"].eq("Multiple").mean()
    L.append(f"\nWhere an expected count falls below 5 the chi-square "
             f"approximation is not reliable, so a Monte-Carlo p (20,000 "
             f"permutations of the table) is reported beside it. "
             f"(MMC1 body site: \u03c7\u00b2=24.4, df=21, p=0.28 \u2014 homogeneous.)\n\n"
             f"Why the restriction matters: run on the full corpus the body-site "
             f"test gives \u03c7\u00b2={full_site:.1f} against {here_site:.1f} here. "
             f"Body site is read from the repository for accession-bearing studies "
             f"({ins_mult:.1f}% 'Multiple') and from the title or abstract for the "
             f"rest ({oth_mult:.1f}% 'Multiple'), and the rest are Tier 4 by "
             f"definition, so the same variable drives both axes. The full-corpus "
             f"values must not be reported.\n")

    reu_n = 100 * ghemi.loc[ghemi.hemi == "Global North", "tier"].isin([1, 2]).mean()
    reu_s = 100 * ghemi.loc[ghemi.hemi == "Global South", "tier"].isin([1, 2]).mean()
    L.append(f"\nReusable (Tier 1-2): Global North {reu_n:.1f}% "
             f"(n={int((ghemi.hemi == 'Global North').sum()):,}) vs Global South "
             f"{reu_s:.1f}% (n={int((ghemi.hemi == 'Global South').sum()):,}).\n")

    # ---------------- Fig S5: tier composition over time -------------------
    # ---------------- Fig S6: tier composition over time -------------------
    # Was panel B of Figure 2; moved here so the main figure can carry the
    # impact-factor result instead. Same 3-year centred rolling mean as before,
    # with the tier legend the panel version inherited from panel A.
    dy = d.copy()
    dy["Year"] = pd.to_numeric(dy["year"], errors="coerce")
    dy = dy[dy["Year"].between(2012, 2025)]
    cnt = dy.groupby(["Year", "tier"]).size().unstack(fill_value=0)
    for t in TIERS:
        if t not in cnt:
            cnt[t] = 0
    cnt = cnt[[4, 3, 2, 1]]
    pctc = (cnt.div(cnt.sum(axis=1), axis=0) * 100
            ).rolling(window=3, center=True, min_periods=1).mean()
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.stackplot(pctc.index, [pctc[t] for t in [4, 3, 2, 1]],
                 colors=[TIER_COLORS[t] for t in [4, 3, 2, 1]],
                 labels=[TLAB[t] for t in [4, 3, 2, 1]], alpha=0.9)
    ax.set_xlim(pctc.index.min(), pctc.index.max())
    ax.set_ylim(0, 100)
    ax.set_xlabel("Year")
    ax.set_ylabel("Proportion (%)")
    h, lb = ax.get_legend_handles_labels()
    ax.legend(h[::-1], lb[::-1], frameon=False, fontsize=8, ncol=4,
              loc="lower center", bbox_to_anchor=(0.5, 1.0))
    style(ax)
    fig.subplots_adjust(left=0.11, right=0.97, top=0.90, bottom=0.12)
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"MMC2_Fig_S5.{e}"), dpi=300)
    plt.close(fig)
    first, last = pctc.index.min(), pctc.index.max()
    L.append(f"\n## Fig S5 — tier composition over time\n\n"
             f"3-year centred rolling mean of the within-year tier split, "
             f"{int(first)}-{int(last)} (n={len(dy):,}). Tier 4 goes from "
             f"{pctc.loc[first, 4]:.0f}% to {pctc.loc[last, 4]:.0f}%, Tier 1-2 "
             f"from {pctc.loc[first, 1] + pctc.loc[first, 2]:.0f}% to "
             f"{pctc.loc[last, 1] + pctc.loc[last, 2]:.0f}%. The first and last "
             f"points average over fewer years (min_periods=1), so read the "
             f"interior of the series, not its endpoints. Moved out of Figure 2.\n")

    open(STATS, "w").write("".join(L))
    print(f"wrote 3 tables + 5 figures to {FIGS}/")
    print(f"wrote {STATS}")
    print("".join(L))


if __name__ == "__main__":
    main()
