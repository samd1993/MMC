"""Figure 2 — six panels, built on the original Figure2.py layout.

Three changes against that layout, all requested by the user:
  * the smoothed tier-composition panel (old B) moves to the supplement as
    Figure S6;
  * the metadata present/absent count bars (old C) become an INSET inside the
    headroom of panel A rather than a panel of their own;
  * the journal impact-factor result takes the top-right space the two of them
    used to occupy.

  A  stacked study counts by year and tier
     inset: metadata present / absent, raw counts (light grey = Yes, dark
     grey = No)
  B  journal impact factor against reusability            (promoted from Fig S3)
  C  Tiers 1-3 over time: multinomial-logit trend + bootstrap CI, bubbles sized
     by study count
  D  Tier 4 and missing-accession trends (circles / squares) plus per-sample
     age, sex and disease annotation (triangles)
  E  waffle of SAMPLES by tier, one block = 5,000 samples
  F  Sankey of accession -> manuscript metadata -> sample-ID match

Panel F is built from the validation dataset (the 2,046-study curated corpus),
which is the only source here with manuscript-level metadata. Every other panel
uses the current corpus.

Inset denominators: accession, age, sex, disease and geography are scored over
all assessable studies; ANTIBIOTICS is scored only over studies whose text shows
antibiotic exposure, since "labelled per sample" is meaningless for a study that
never involved antibiotics. The inset is too small to carry that note, so the
dagger is explained in the Figure 2 caption (Figure2_caption.md).

Panel B uses the Clarivate JIF only (`impact_factor`), never the mixed column —
the OpenAlex value is a rank proxy and must not enter a regression.

Panel E sample sizes follow the earlier method: take the per-study sample count,
then fill studies that have none with the mean of their own tier. Here the count
is the deposited ENA sample count where an accession exists, otherwise the
largest "n = " reported in the abstract; Tier 4 studies have no deposit by
definition, so 82% of them are filled from the mean deposit size of studies
that DID deposit -- not from their own abstracts, which run 2-3x low.

Writes Figs/Figure2_revised.{pdf,png} and figure2_revised_source_data.tsv.
"""

import io
import os
import re
import sys
import tarfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.gridspec as gridspec  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from matplotlib.path import Path as MplPath  # noqa: E402
from scipy.special import expit  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import MMC, ROOT  # noqa: E402

MASTER = os.path.join(ROOT, "MMC2_master_reviewed.tsv")
ATTR = os.path.join(ROOT, "MMC2_study_attributes.tsv")
STUDY = os.path.join(ROOT, "MMC2_study_metadata.tsv")
ANALYSIS = os.path.join(ROOT, "MMC2_analysis_set.tsv")
FIGS = os.path.join(ROOT, "Figs")
SRC = os.path.join(ROOT, "figure2_revised_source_data.tsv")

TIER_COLORS = {1: "#91bfdb", 2: "#fee090", 3: "#fc8d59", 4: "#d73027"}
TIER_LABELS = {1: "Tier 1", 2: "Tier 2", 3: "Tier 3", 4: "Tier 4"}
STACK_ORDER = [4, 3, 2, 1]
PRED_COL = {1: 0, 2: 1, 3: 2, 4: 3}
GREY_YES, GREY_NO = "#d9d9d9", "#636363"
# the inset is small enough that the full category names collide; the caption
# carries the long form
INSET_LABELS = {"Accession\nCode": "Accession", "Age": "Age", "Sex": "Sex",
                "Disease": "Disease", "Antibiotics†": "Abx†",
                "Geography": "Geography"}
SCALE = 3.0
N_BOOT = 2000   # raised from 200 for Reviewer 1; see MMC2_adjusted_trends.md
rows = []

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
})


def style(ax):
    ax.tick_params(axis="both", labelsize=11)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# Tier 4 sample-size imputation. "all" = every deposited study;
# "single_accession" = one-deposit studies only, excluding American Gut, which
# strips the meta-analysis inflation out of the average.
T4_POOL = "all"
MIN_POOL_SAMPLES = 2


def panel_label(ax, letter, dx=-0.12, dy=1.08):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=20,
            fontweight="bold", va="top", ha="left")


def reported_n(s):
    if pd.isna(s):
        return np.nan
    m = [int(x) for x in re.findall(r"\bn\s*=\s*(\d{1,5})", str(s))]
    return max(m) if m else np.nan


def load():
    at = pd.read_csv(ATTR, sep="\t", low_memory=False)
    st = pd.read_csv(STUDY, sep="\t", low_memory=False)
    an = pd.read_csv(ANALYSIS, sep="\t", low_memory=False)[["record_id", "abstract"]]
    cols = ["record_id", "disease_evidence", "age_present", "sex_present",
            "antibiotic_present", "geography_present", "abx_involves_study",
            "metadata_source", "n_samples", "accession_codes"]
    d = (at.merge(st[[c for c in cols if c in st.columns]], on="record_id",
                  how="left", suffixes=("", "_s"))
           .merge(an, on="record_id", how="left"))
    # Tier comes from MMC2_master_reviewed.tsv: the manual validation moved 124
    # studies, so the pipeline column is superseded everywhere it disagrees.
    mr = pd.read_csv(MASTER, sep="\t", low_memory=False,
                     keep_default_na=False, dtype=str)
    d = d.merge(mr[["record_id", "final_tier", "name_derived_disease_flag",
                    "healthy_flag", "wrong_accession_flag", "in_analysis",
                    "orphan_accession_recovered"]],
                on="record_id", how="left")
    # Sam, 2 Sep: Healthy-category studies are out of every analysis
    d = d[d["in_analysis"].eq("yes")].copy()
    d["tier"] = pd.to_numeric(d["final_tier"], errors="coerce")
    d = d[d["tier"].notna()].copy()
    d["tier"] = d["tier"].astype(int)
    d["Year"] = pd.to_numeric(d["year"], errors="coerce").astype(int)
    d = d[d["Year"].between(2012, 2025)].copy()
    # six studies gained an accession that reviewers wrote into sheet 1 and no
    # extraction pass had picked up, so deposition is scored after that fix
    d["AccPresent"] = (d["metadata_source"].eq("insdc")
                       | d["orphan_accession_recovered"].eq("yes")).astype(int)
    return d


def main():
    os.makedirs(FIGS, exist_ok=True)
    df = load()
    N = len(df)
    counts = df.groupby(["Year", "tier"]).size().unstack(fill_value=0)
    for t in (1, 2, 3, 4):
        if t not in counts:
            counts[t] = 0
    counts = counts[STACK_ORDER]
    pct = counts.div(counts.sum(axis=1), axis=0) * 100
    pct_smooth = pct.rolling(window=3, center=True, min_periods=1).mean()
    yearly_totals = df.groupby("Year").size()
    yearly_props = (df.groupby(["Year", "tier"]).size().unstack(fill_value=0)
                    .div(yearly_totals, axis=0))
    year_grid = np.linspace(df["Year"].min(), df["Year"].max(), 200)
    X_grid = sm.add_constant(year_grid)

    # ---- multinomial model for panel D ----------------------------------
    cat = pd.Categorical(df["tier"], categories=[1, 2, 3, 4])
    Xc, yc = sm.add_constant(df["Year"]).values, cat.codes
    mnl = sm.MNLogit(yc, Xc).fit(disp=False, maxiter=100)
    mnl_pred = mnl.predict(X_grid)
    # RRR per year for Tiers 1-3 against Tier 4 as the reference level.
    # MNLogit params are indexed (exog, outcome): row 1 is Year, column k is the
    # k-th non-base outcome.
    cat4 = pd.Categorical(df["tier"], categories=[4, 1, 2, 3])
    m4 = sm.MNLogit(cat4.codes, Xc).fit(disp=False, maxiter=100)
    par = np.asarray(m4.params)
    pv = np.asarray(m4.pvalues)
    mnl_stat = {t: (float(np.exp(par[1, k])), float(pv[1, k]))
                for k, t in enumerate((1, 2, 3))}

    rng = np.random.RandomState(42)
    boot = np.full((N_BOOT, len(year_grid), 4), np.nan)
    for b in range(N_BOOT):
        idx = rng.choice(N, N, replace=True)
        try:
            boot[b] = sm.MNLogit(yc[idx], Xc[idx]).fit(disp=False, maxiter=50).predict(X_grid)
        except Exception:  # noqa: BLE001
            pass
    ci_lo = np.nanpercentile(boot, 2.5, axis=0)
    ci_hi = np.nanpercentile(boot, 97.5, axis=0)

    # ---- logistic trends for panel E -------------------------------------
    df["is_t4"] = (df["tier"] == 4).astype(int)
    mt4 = smf.logit("is_t4 ~ Year", data=df).fit(disp=False)
    t4_or, t4_p = np.exp(mt4.params["Year"]), mt4.pvalues["Year"]
    lo4 = mt4.params["Intercept"] + mt4.params["Year"] * year_grid
    se4 = np.sqrt(np.sum(X_grid @ mt4.cov_params().values * X_grid, axis=1))
    t4_line, t4_lo, t4_hi = expit(lo4), expit(lo4 - 1.96 * se4), expit(lo4 + 1.96 * se4)
    t4_obs = df.groupby("Year")["is_t4"].mean()

    ma = smf.logit("AccPresent ~ Year", data=df).fit(disp=False)
    acc_or_no, acc_p = np.exp(-ma.params["Year"]), ma.pvalues["Year"]
    loa = ma.params["Intercept"] + ma.params["Year"] * year_grid
    sea = np.sqrt(np.sum(X_grid @ ma.cov_params().values * X_grid, axis=1))
    acc_line = 1 - expit(loa)
    acc_lo, acc_hi = 1 - expit(loa + 1.96 * sea), 1 - expit(loa - 1.96 * sea)
    acc_obs = 1 - df.groupby("Year")["AccPresent"].mean()

    # per-sample biological annotation trends shown alongside the two
    # reusability series; geography and antibiotics are deliberately excluded
    def logistic_trend(col):
        m = smf.logit(f"{col} ~ Year", data=df).fit(disp=False)
        lin = m.params["Intercept"] + m.params["Year"] * year_grid
        se = np.sqrt(np.sum(X_grid @ m.cov_params().values * X_grid, axis=1))
        return dict(line=expit(lin), lo=expit(lin - 1.96 * se),
                    hi=expit(lin + 1.96 * se), OR=float(np.exp(m.params["Year"])),
                    p=float(m.pvalues["Year"]), obs=df.groupby("Year")[col].mean())

    df["age_p"] = df["age_present"].astype(int)
    df["sex_p"] = df["sex_present"].astype(int)
    # after review, "this study has per-sample disease" IS the Tier 1/2 test
    df["dis_p"] = df["tier"].isin([1, 2]).astype(int)
    BIO = {"dis": ("dis_p", "#a65628", "Disease"),
           "age": ("age_p", "#4daf4a", "Age"),
           "sex": ("sex_p", "#984ea3", "Sex")}
    bio = {k: logistic_trend(c) for k, (c, _, _) in BIO.items()}

    # ---- panel F sample sizes --------------------------------------------
    df["reported_n"] = df["abstract"].map(reported_n)
    df["sample_size"] = np.where(df["n_samples"].fillna(0) > 0, df["n_samples"],
                                 df["reported_n"])
    # Tier 4 is DEFINED by having no usable deposit, so the only size it can
    # state is a text scrape of its abstract -- which runs 2-3x below the same
    # study's deposit would (median ENA/abstract ratio 3.45x in T1, 2.22x in T2,
    # 2.59x in T3). Filling Tier 4 from its own abstracts therefore made it the
    # SMALLEST tier per study (141) when the hand-curated MMC1 pass found it the
    # LARGEST (365). Impute it instead from what a deposited study actually
    # holds. Studies of 1-2 samples are deposit errors, not one-sample studies,
    # and are excluded from the average.
    # Sizes come from build_sample_sizes.py, which is the single place the
    # 2026-08-31 counting rules live: Tier 1/2 count the study's OWN deposit
    # (the first accession, which is the one the reviewer was shown), Tier 3
    # counts every accession, Tier 4 is filled at the Tier 3 mean.
    ss = pd.read_csv(os.path.join(ROOT, "MMC2_study_sample_sizes.tsv"), sep="\t")
    ss = ss[ss["in_analysis"].eq("yes")]
    df = df.merge(ss[["record_id", "n_samples", "rule"]], on="record_id",
                  how="left", suffixes=("", "_rule"))
    df["sample_size_filled"] = df["n_samples_rule"]
    df["sample_size"] = np.where(df["rule"].fillna("").str.contains("fill"),
                                 np.nan, df["n_samples_rule"])
    t4_fill = float(ss.loc[ss["rule"].eq("fill_tier1_3_mean"), "n_samples"].iloc[0])
    print(f"  sample sizes from MMC2_study_sample_sizes.tsv "
          f"(Tier 4 fill {t4_fill:.0f}/study)")
    tier_samples = [int(round(df.loc[df["tier"] == t, "sample_size_filled"].sum()))
                    for t in (1, 2, 3, 4)]
    total_samples = sum(tier_samples)
    # Observed-only totals, reported alongside. The tier-mean fill is 37% of the
    # corpus total and 81% of Tier 4, because Tier 4 means "no usable deposit" so
    # only 366 of its 1,969 studies state a size anywhere. Both numbers go into
    # the source data so neither is quoted without the other.
    obs = df[df["sample_size"].notna()]
    tier_obs = [int(round(obs.loc[obs["tier"] == t, "sample_size"].sum()))
                for t in (1, 2, 3, 4)]
    total_obs = sum(tier_obs)
    n_imputed = int(df["sample_size"].isna().sum())
    n_filled = int(df["sample_size"].isna().sum())

    # ---- panel G: validation dataset --------------------------------------
    tf = tarfile.open(os.path.join(MMC, "MMC_resubmission_handoff.tar.gz"))
    v9 = pd.read_csv(io.BytesIO(tf.extractfile(
        "MMC_resubmission_handoff/MMC1_data_merged_v9_May23_26.tsv").read()),
        sep="\t", low_memory=False)
    v9["Year"] = pd.to_numeric(v9["Year"], errors="coerce")
    g = v9[v9["Year"].between(2012, 2030) & v9["Tier_new"].isin([1, 2, 3, 4])].copy()
    g["AccessionCode"] = g["AccessionCode"].fillna("").astype(str).str.strip()
    g["AccPres"] = ~g["AccessionCode"].isin(["", "N/A", "n/a", "N/a"])
    g["MetadataPresent"] = g["MetadataPresent"].fillna("No").astype(str).str.strip()
    g["SampleIDmatches"] = g["SampleIDmatches"].fillna("No").astype(str).str.strip()
    g["SIDYes"] = g["SampleIDmatches"].isin(["Yes", "Yes in sample link"])
    g_total = len(g)
    acc_yes = int(g["AccPres"].sum())
    acc_no = g_total - acc_yes
    ga = g[g["AccPres"]]
    meta_yes = int((ga["MetadataPresent"] == "Yes").sum())
    meta_no = len(ga) - meta_yes
    gmy = ga[ga["MetadataPresent"] == "Yes"]
    sid_yes = int(gmy["SIDYes"].sum())
    sid_no = len(gmy) - sid_yes
    gna = g[~g["AccPres"]]
    na_meta_yes = int((gna["MetadataPresent"] == "Yes").sum())
    na_meta_no = len(gna) - na_meta_yes

    # =================== FIGURE ==========================================
    # Top row is now A (with its inset) on the left and the impact-factor panel
    # on the right, where the composition panel and the count bars used to sit.
    fig = plt.figure(figsize=(16.5, 16.2))
    gs = gridspec.GridSpec(3, 6, height_ratios=[1, 1.18, 0.95],
                           hspace=0.30, wspace=0.80)
    axA = fig.add_subplot(gs[0, 0:3])
    axIF = fig.add_subplot(gs[0, 3:6])
    axD = fig.add_subplot(gs[1, 0:3])
    axE = fig.add_subplot(gs[1, 3:6])
    axF = fig.add_subplot(gs[2, 0:3])
    axG = fig.add_subplot(gs[2, 3:6])

    # ---- A ---------------------------------------------------------------
    bottom = np.zeros(len(counts))
    xs = counts.index.astype(int)
    for t in STACK_ORDER:
        axA.bar(xs, counts[t], bottom=bottom, color=TIER_COLORS[t],
                linewidth=0, width=0.75)
        bottom += counts[t].values
    axA.set_xticks(list(xs))
    axA.set_xticklabels([int(y) for y in xs], rotation=45, ha="right")
    axA.set_xlabel("Year", fontsize=14, fontweight="bold")
    axA.set_ylabel("Number of Studies", fontsize=14, fontweight="bold")
    # The bar chart keeps its own scale -- no added headroom. The inset sits in
    # the white space the rising bars leave open, above the early years.
    # laid along the top so the whole upper-left is free for the inset
    axA.legend(handles=[Patch(facecolor=TIER_COLORS[t], linewidth=0,
                              label=TIER_LABELS[t]) for t in (1, 2, 3, 4)],
               fontsize=11, loc="lower center", bbox_to_anchor=(0.5, 1.005),
               ncol=4, frameon=False, handlelength=1.0, handleheight=1.0,
               columnspacing=1.6, borderaxespad=0.0)
    panel_label(axA, "A")
    style(axA)
    for y, v in zip(xs, bottom):
        rows.append({"panel": "A", "x": int(y), "series": "total", "value": int(v)})

    # ---- A inset: metadata present / absent ------------------------------
    # Was its own panel C. Shrunk into A's headroom; the dagger denominator note
    # is carried by the figure caption instead of the panel.
    # Sam, 3 Sep: more room, centre-left of the panel
    axC = axA.inset_axes([0.105, 0.500, 0.480, 0.430], zorder=6)
    inv = df[df["abx_involves_study"].fillna(False).astype(bool)]
    specs = [("Accession\nCode", int(df["AccPresent"].sum()), N),
             ("Age", int(df["age_present"].sum()), N),
             ("Sex", int(df["sex_present"].sum()), N),
             # after review, per-sample disease IS the Tier 1/2 criterion
             ("Disease", int(df["tier"].isin([1, 2]).sum()), N),
             ("Antibiotics†", int(inv["antibiotic_present"].sum()), len(inv)),
             ("Geography", int(df["geography_present"].sum()), N)]
    yes = [s[1] for s in specs]
    no = [s[2] - s[1] for s in specs]
    x = np.arange(len(specs))
    w = 0.40
    by = axC.bar(x - w / 2, yes, w, color=GREY_YES, linewidth=0, label="Yes")
    bn = axC.bar(x + w / 2, no, w, color=GREY_NO, linewidth=0, label="No")
    # Counts above the bars are replaced by a real y axis -- the numbers are
    # in figure2_revised_source_data.tsv for anyone who needs them exactly.
    axC.set_xticks(x)
    axC.set_xticklabels([INSET_LABELS[s[0]] for s in specs], fontsize=8.5,
                        fontweight="bold", rotation=90, ha="center",
                        va="top")
    top = max(max(yes), max(no))
    step = 500 if top <= 2600 else 1000
    ticks = list(range(0, int(top) + step, step))
    axC.set_yticks(ticks)
    axC.set_yticklabels([f"{t:,}" for t in ticks], fontsize=8.0)
    axC.set_xlim(-0.65, len(specs) - 0.35)
    axC.set_ylim(0, top * 1.10)
    axC.set_ylabel("Studies", fontsize=9.5, fontweight="bold", labelpad=1)
    axC.legend(fontsize=8.0, loc="upper right", frameon=False, borderpad=0.1,
               handlelength=0.9, handleheight=0.9, borderaxespad=0.2)
    axC.patch.set_facecolor("white")
    axC.patch.set_alpha(1.0)
    for s in ("top", "right"):
        axC.spines[s].set_visible(False)
    axC.spines["bottom"].set_linewidth(0.7)
    axC.spines["left"].set_linewidth(0.7)
    axC.tick_params(axis="x", labelsize=8.5, length=0, pad=2.0)
    axC.tick_params(axis="y", labelsize=8.0, length=2.5, pad=1.5)
    axC.grid(axis="y", color="#d9d9d9", linewidth=0.5, alpha=0.8)
    axC.set_axisbelow(True)
    for lab, y_, tot in specs:
        rows.append({"panel": "A inset", "x": lab.replace("\n", " "),
                     "series": "Yes", "value": y_})
        rows.append({"panel": "A inset", "x": lab.replace("\n", " "),
                     "series": "No", "value": tot - y_})

    # ---- B: journal impact factor vs reusability -------------------------
    # Promoted out of the supplement into the space the composition panel and
    # the count bars vacated. Clarivate JIF only.
    jr = df[df["impact_factor"].notna() & (df["impact_factor"] > 0)].copy()
    jvc = jr["journal"].value_counts()
    jr = jr[jr["journal"].isin(jvc[jvc >= 10].index)].copy()
    jr["logIF"] = np.log10(jr["impact_factor"])
    jr["reusable"] = jr["tier"].isin([1, 2]).astype(int)
    mj = smf.logit("reusable ~ logIF", data=jr).fit(disp=False)
    j_or, j_p = float(np.exp(mj.params["logIF"])), float(mj.pvalues["logIF"])
    j_lo, j_hi = np.exp(mj.conf_int().loc["logIF"])
    grid_if = np.linspace(jr["logIF"].min(), jr["logIF"].max(), 120)
    pred_if = expit(mj.params["Intercept"] + mj.params["logIF"] * grid_if)
    per_j = jr.groupby("journal").agg(x=("logIF", "first"),
                                      y=("reusable", "mean"),
                                      n=("reusable", "size"))
    axIF.scatter(10 ** per_j.x, 100 * per_j.y, s=per_j.n * 2.4, color="#377eb8",
                 linewidth=0, alpha=0.75, zorder=3)
    axIF.plot(10 ** grid_if, 100 * pred_if, color="#d73027", lw=2.6, zorder=2)
    axIF.set_xscale("log")
    axIF.set_xlabel("Journal Impact Factor (log scale)", fontsize=14,
                    fontweight="bold")
    axIF.set_ylabel("Reusable studies,\nTier 1 or 2 (%)", fontsize=14,
                    fontweight="bold")
    axIF.legend([plt.Line2D([], [], color="#d73027", lw=2.6)],
                [f"Odds ratio = {j_or:.2f} [{j_lo:.2f}–{j_hi:.2f}] "
                 f"per log10(Impact Factor),\n"
                 f"p = {'<0.001' if j_p < 0.001 else f'{j_p:.3f}'};  "
                 f"{jr.journal.nunique()} journals, {len(jr):,} studies"],
                frameon=False, fontsize=10, loc="upper left")
    axIF.set_ylim(-3, 100 * per_j.y.max() * 1.30)   # headroom for the legend
    axIF.set_title("Marker area = studies per journal", loc="left",
                   fontsize=10, color="#4d4d4d", pad=8)
    panel_label(axIF, "B")
    style(axIF)
    rows.append({"panel": "B", "x": "P(Tier 1-2)", "series": "OR",
                 "value": round(j_or, 3)})
    rows.append({"panel": "B", "x": "P(Tier 1-2)", "series": "p", "value": j_p})
    rows.append({"panel": "B", "x": "n_journals", "series": "count",
                 "value": int(jr.journal.nunique())})

    # ---- C: tier 1-3 trends ----------------------------------------------
    for t in (1, 2, 3):
        c = TIER_COLORS[t]
        axD.fill_between(year_grid, ci_lo[:, PRED_COL[t]], ci_hi[:, PRED_COL[t]],
                         color=c, alpha=0.2, zorder=1)
        axD.plot(year_grid, mnl_pred[:, PRED_COL[t]], color=c, linewidth=2.5, zorder=2)
        obs = yearly_props[t] if t in yearly_props.columns else pd.Series(dtype=float)
        axD.scatter(obs.index, obs.values,
                    s=yearly_totals.reindex(obs.index) * SCALE, color=c,
                    linewidth=0, alpha=0.8, zorder=3)
        rows.append({"panel": "C", "x": TIER_LABELS[t], "series": "RRR",
                     "value": round(mnl_stat[t][0], 3)})
        rows.append({"panel": "C", "x": TIER_LABELS[t], "series": "p",
                     "value": float(mnl_stat[t][1])})
    d_handles = [plt.scatter([], [], s=100, color=TIER_COLORS[t], linewidth=0)
                 for t in (1, 2, 3)]
    d_labels = [f"{TIER_LABELS[t]}   RRR={mnl_stat[t][0]:.3f}, p={mnl_stat[t][1]:.3f}"
                for t in (1, 2, 3)]
    size_handles = [axD.scatter([], [], s=n * SCALE, facecolor="#cfcfcf",
                                linewidth=0) for n in (50, 150, 400)]
    size_leg = axD.legend(size_handles, [f"n={n}" for n in (50, 150, 400)],
                          fontsize=9, loc="upper right", title="Study count",
                          frameon=False, handlelength=1.5, handleheight=1.5,
                          labelspacing=1.2)
    axD.add_artist(size_leg)
    axD.legend(handles=d_handles, labels=d_labels, fontsize=9.5, loc="upper left",
               frameon=False, handlelength=1.0, handleheight=1.0)
    ymax = float(np.nanmax([yearly_props[t].max() for t in (1, 2, 3)
                            if t in yearly_props.columns]))
    axD.set_ylim(-0.01, max(0.25, ymax * 1.25))
    axD.set_xlabel("Year", fontsize=14, fontweight="bold")
    axD.set_ylabel("Proportion", fontsize=14, fontweight="bold")
    panel_label(axD, "C")
    style(axD)

    # ---- D: tier 4 / missing-accession trends ----------------------------
    axE.fill_between(year_grid, t4_lo, t4_hi, color="#d73027", alpha=0.15, zorder=1)
    axE.plot(year_grid, t4_line, color="#d73027", linewidth=2.5, zorder=2)
    axE.scatter(t4_obs.index, t4_obs.values, s=yearly_totals * SCALE, color="#d73027",
                linewidth=0, alpha=0.8, zorder=3)
    axE.fill_between(year_grid, acc_lo, acc_hi, color="#d73027", alpha=0.1, zorder=1)
    axE.plot(year_grid, acc_line, color="#d73027", linewidth=2.5, linestyle="--", zorder=2)
    axE.scatter(acc_obs.index, acc_obs.values, s=yearly_totals * SCALE, color="#d73027",
                linewidth=0, alpha=0.55, zorder=3, marker="s")
    for k, (col, color, _name) in BIO.items():
        b = bio[k]
        axE.fill_between(year_grid, b["lo"], b["hi"], color=color, alpha=0.15, zorder=1)
        axE.plot(year_grid, b["line"], color=color, linewidth=2.4, zorder=2)
        axE.scatter(b["obs"].index, b["obs"].values,
                    s=yearly_totals.reindex(b["obs"].index) * SCALE, color=color,
                    linewidth=0, alpha=0.8, zorder=3, marker="^")
        rows.append({"panel": "D", "x": _name, "series": "OR",
                     "value": round(b["OR"], 3)})
        rows.append({"panel": "D", "x": _name, "series": "p", "value": b["p"]})

    def pf(pv):
        return "p<0.001" if pv < 0.001 else f"p={pv:.3f}"

    handles = [plt.Line2D([], [], color="#d73027", marker="o", ls="-", lw=2.2,
                          markersize=9, markeredgewidth=0),
               plt.Line2D([], [], color="#d73027", marker="s", ls="--", lw=2.2,
                          markersize=9, markeredgewidth=0, alpha=0.7)]
    labels = [f"Tier 4 (not reusable)   OR={t4_or:.3f}, {pf(t4_p)}",
              f"No accession code   OR={acc_or_no:.3f}, {pf(acc_p)}"]
    for k, (col, color, name) in BIO.items():
        handles.append(plt.Line2D([], [], color=color, marker="^", ls="-", lw=2.2,
                                  markersize=9, markeredgewidth=0))
        labels.append(f"{name}   OR={bio[k]['OR']:.3f}, {pf(bio[k]['p'])}")
    axE.legend(handles, labels, fontsize=9.5, loc="upper right", frameon=False,
               handlelength=1.6, labelspacing=0.55)
    axE.set_ylim(-0.02, 1.0)
    axE.set_xlabel("Year", fontsize=14, fontweight="bold")
    axE.set_ylabel("Proportion", fontsize=14, fontweight="bold")
    panel_label(axE, "D")
    style(axE)
    rows.append({"panel": "D", "x": "Tier 4", "series": "OR", "value": round(t4_or, 3)})
    rows.append({"panel": "D", "x": "Tier 4", "series": "p", "value": float(t4_p)})
    rows.append({"panel": "D", "x": "No accession", "series": "OR",
                 "value": round(acc_or_no, 3)})
    rows.append({"panel": "D", "x": "No accession", "series": "p", "value": float(acc_p)})

    # ---- E: waffle -------------------------------------------------------
    BLOCK_SAMPLES, ROW_H, BLOCK = 5000, 6, 0.9
    tier_blocks = [int(round(n / BLOCK_SAMPLES)) for n in tier_samples]
    total_cols = int(np.ceil(sum(tier_blocks) / ROW_H))
    bi = 0
    edges = []
    for n_blocks, t in zip(tier_blocks, (1, 2, 3, 4)):
        for _ in range(n_blocks):
            axF.add_patch(Rectangle((bi // ROW_H, bi % ROW_H), BLOCK, BLOCK,
                                    facecolor=TIER_COLORS[t], edgecolor="white",
                                    linewidth=0.3))
            bi += 1
        edges.append(bi / ROW_H)
    prev = 0
    for t, edge, n_samp in zip((1, 2, 3, 4), edges, tier_samples):
        axF.text((prev + edge) / 2, ROW_H + 0.4,
                 f"T{t}\n{n_samp / 1000:.0f}k ({n_samp / total_samples:.0%})",
                 ha="center", va="bottom", fontsize=10, fontweight="bold",
                 color=TIER_COLORS[t])
        prev = edge
        rows.append({"panel": "E", "x": f"Tier {t}", "series": "samples",
                     "value": n_samp})
        rows.append({"panel": "E", "x": f"Tier {t}", "series": "samples_observed",
                     "value": tier_obs[t - 1]})
    lost_start = edges[1]
    axF.annotate("", xy=(lost_start, -0.35), xytext=(edges[3], -0.35),
                 arrowprops=dict(arrowstyle="-", color="#4d4d4d", lw=1.0))
    axF.text((lost_start + edges[3]) / 2, -0.75,
             "Samples lost due to lack of repository metadata\nand/or sequencing data",
             ha="center", va="top", fontsize=10)
    axF.set_xlim(-0.4, total_cols + 0.4)
    axF.set_ylim(-2.6, ROW_H + 2.2)
    axF.set_aspect("equal")
    axF.set_xticks([])
    axF.set_yticks([])
    for s in axF.spines.values():
        s.set_visible(False)
    axF.text(total_cols, -2.35, f"= {BLOCK_SAMPLES // 1000}k samples   ·   "
             f"total ≈ {total_samples / 1000:.0f}k", ha="right", va="bottom",
             fontsize=9)
    axF.text(0, -2.35, f"{total_obs / 1000:.0f}k reported by {len(obs):,} studies; "
             f"{n_imputed:,} studies without a reported size are filled at "
             f"{t4_fill:.0f} samples (Tier 4) or their tier mean",
             ha="left", va="bottom", fontsize=7.5, color="#4d4d4d")
    axF.add_patch(Rectangle((total_cols - 0.0, -2.45), BLOCK, BLOCK,
                            facecolor="#808080", edgecolor="white", linewidth=0.3,
                            clip_on=False))
    panel_label(axF, "E", dy=1.02)

    # ---- F: Sankey (validation dataset) ----------------------------------
    axG.set_xlim(0, 13.6)
    axG.set_ylim(0, 10)
    axG.axis("off")
    c_yes, c_no, c_grey = "#2C7BB6", "#D7191C", "#999999"
    ALPHA = 0.25
    x_cols = [0.5, 3.3, 6.3, 9.2]
    nw, gap, total_h, y_base = 0.18, 0.30, 7.0, 1.6

    def sc(n):
        return max(n / g_total * total_h, 0.15)

    ay_h, an_h = sc(acc_yes), sc(acc_no)
    my_h, mn_h = sc(meta_yes), sc(meta_no)
    sy_h, sn_h = sc(sid_yes), sc(sid_no)
    na_y_h, na_n_h = sc(na_meta_yes), sc(na_meta_no)
    n_all = (x_cols[0], y_base, sc(g_total))
    n_ay = (x_cols[1], y_base + an_h + gap, ay_h)
    n_an = (x_cols[1], y_base, an_h)
    n_my = (x_cols[2], n_ay[1] + mn_h + gap, my_h)
    n_mn = (x_cols[2], n_ay[1], mn_h)
    n_sy = (x_cols[3], n_my[1] + sn_h + gap, sy_h)
    n_sn = (x_cols[3], n_my[1], sn_h)
    x_grey = (x_cols[1] + x_cols[2]) / 2
    n_nan = (x_grey, y_base, na_n_h)
    n_nay = (x_grey, y_base + na_n_h + gap * 0.5, na_y_h)

    def flow(x0, y0, h0b, h0t, x1, y1, h1, color):
        xm = (x0 + nw + x1) / 2
        sx, tx = x0 + nw, x1
        verts = [(sx, y0 + h0t), (xm, y0 + h0t), (xm, y1 + h1), (tx, y1 + h1),
                 (tx, y1), (xm, y1), (xm, y0 + h0b), (sx, y0 + h0b), (sx, y0 + h0t)]
        codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                 MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                 MplPath.CLOSEPOLY]
        axG.add_patch(mpatches.PathPatch(MplPath(verts, codes), facecolor=color,
                                         edgecolor="none", alpha=ALPHA))

    flow(n_all[0], n_all[1], an_h, an_h + ay_h, n_ay[0], n_ay[1], ay_h, c_yes)
    flow(n_all[0], n_all[1], 0, an_h, n_an[0], n_an[1], an_h, c_no)
    flow(n_ay[0], n_ay[1], mn_h, mn_h + my_h, n_my[0], n_my[1], my_h, c_yes)
    flow(n_ay[0], n_ay[1], 0, mn_h, n_mn[0], n_mn[1], mn_h, c_no)
    flow(n_my[0], n_my[1], sn_h, sn_h + sy_h, n_sy[0], n_sy[1], sy_h, c_yes)
    flow(n_my[0], n_my[1], 0, sn_h, n_sn[0], n_sn[1], sn_h, c_no)
    flow(n_an[0], n_an[1], na_n_h, na_n_h + na_y_h, n_nay[0], n_nay[1], na_y_h, c_grey)
    flow(n_an[0], n_an[1], 0, na_n_h, n_nan[0], n_nan[1], na_n_h, c_grey)

    for (x, y, h), col in [(n_all, "#555555"), (n_ay, c_yes), (n_an, c_no),
                           (n_my, c_yes), (n_mn, c_no), (n_sy, c_yes),
                           (n_sn, c_no), (n_nay, c_grey), (n_nan, c_grey)]:
        axG.add_patch(mpatches.FancyBboxPatch((x, y), nw, h,
                                              boxstyle="round,pad=0.02",
                                              facecolor=col, edgecolor="white",
                                              linewidth=0.5))
    def lab(node, text, side="right", fs=10):
        x, y, h = node
        if side == "right":
            axG.text(x + nw + 0.12, y + h / 2, text, va="center", ha="left", fontsize=fs)
        else:
            axG.text(x - 0.12, y + h / 2, text, va="center", ha="right", fontsize=fs)

    lab(n_all, f"All studies\n(n={g_total:,})", side="left")
    lab(n_ay, f"Accession present\n(n={acc_yes:,})")
    lab(n_an, f"Accession absent\n(n={acc_no:,})")
    lab(n_my, f"Manuscript metadata\npresent (n={meta_yes:,})")
    lab(n_mn, f"Metadata absent\n(n={meta_no:,})")
    lab(n_sy, f"Sample ID\nmatch (n={sid_yes:,})")
    lab(n_sn, f"Sample ID\nno match (n={sid_no:,})")
    lab(n_nay, f"Metadata present (n={na_meta_yes:,})", fs=9)
    lab(n_nan, f"Metadata absent (n={na_meta_no:,})", fs=9)
    panel_label(axG, "F", dx=-0.02, dy=1.02)
    for k, v in [("all", g_total), ("acc_yes", acc_yes), ("acc_no", acc_no),
                 ("meta_yes", meta_yes), ("meta_no", meta_no), ("sid_yes", sid_yes),
                 ("sid_no", sid_no), ("na_meta_yes", na_meta_yes),
                 ("na_meta_no", na_meta_no)]:
        rows.append({"panel": "F", "x": k, "series": "n", "value": int(v)})

    fig.subplots_adjust(left=0.055, right=0.985, top=0.965, bottom=0.035)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"Figure2_revised.{ext}"), dpi=300)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(SRC, sep="\t", index=False)

    print(f"wrote {FIGS}/Figure2_revised.pdf / .png")
    print(f"wrote {SRC}")
    print(f"\nA/B/D/E  n={N:,} studies, {counts.values.sum():,} in the year grid")
    print("C yes/no:")
    for lab_, y_, tot in specs:
        print(f"   {lab_.replace(chr(10), ' '):16s} Yes={y_:5d}  No={tot - y_:5d}  (of {tot:,})")
    print("D  RRR: " + " | ".join(f"T{t} {mnl_stat[t][0]:.3f} p={mnl_stat[t][1]:.3f}"
                                  for t in (1, 2, 3)))
    print(f"E  Tier4 OR={t4_or:.3f} p={t4_p:.3f} | no-accession OR={acc_or_no:.3f} p={acc_p:.3f}")
    print(f"F  samples per tier T1..T4: {tier_samples}  total={total_samples:,}")
    print(f"   observed only:            {tier_obs}  total={total_obs:,}")
    print(f"   {n_filled:,}/{N:,} studies had no per-study count and were filled "
          f"from the deposited-study mean ({t4_fill:.0f}) in Tier 4, "
          f"their tier mean elsewhere")
    print(f"G  validation dataset n={g_total:,}: accession {acc_yes}/{acc_no}, "
          f"metadata {meta_yes}/{meta_no}, sample-ID {sid_yes}/{sid_no}")


if __name__ == "__main__":
    main()
