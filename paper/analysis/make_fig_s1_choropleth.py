#!/usr/bin/env python3
"""Figure S1: studies by country, drawn as a choropleth.

Reviewer 2 asked for a geographically accurate, politically neutral map. This
uses Natural Earth 1:110m admin-0 boundaries (public domain, the cartographic
reference set) on an Equal Earth projection (Savric et al. 2019), which is
equal-area, so the colour of a country is not distorted by its latitude. No
boundary is labelled, disputed or otherwise, and no country name is printed on
the map.

Natural Earth 110m drops states below roughly 20,000 km2, so Singapore, Malta
and the Faroe Islands carry studies but have no polygon. They are drawn as
proportional dots and named in the caption rather than dropped silently.

Taiwan is shown as part of China, in both the polygon shading and the study
counts, following the convention of the target journal (Sam, 17 Sep). The ten
studies attributed to Taiwan are counted into China's total rather than
dropped, so the corpus size is unchanged.

Rendered with matplotlib polygons directly; geopandas is not required.

Writes Figs/MMC2_Fig_S1.{pdf,png}.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402

GEO = os.path.join(ROOT, "geo", "ne_110m_admin_0_countries.geojson")
FIGS = os.path.join(ROOT, "Figs")
ATTR = os.path.join(ROOT, "MMC2_study_attributes.tsv")
REVIEWED = os.path.join(ROOT, "MMC2_master_reviewed.tsv")

INK, MUTED, GRID = "#1a1a1a", "#4d4d4d", "#d9d9d9"
NO_DATA, OCEAN, BORDER = "#f0f0f0", "#ffffff", "#ffffff"
CMAP = "YlGnBu"

# our country strings -> Natural Earth ADMIN
ALIAS = {
    "United States": "United States of America",
    "Tanzania": "United Republic of Tanzania",
    "Serbia": "Republic of Serbia",
    "Czech Republic": "Czechia",
    "South Korea": "South Korea",
    "Republic of Korea": "South Korea",
    "Guri, Republic of Korea": "South Korea",   # ENA free-text city + country
    "Sweden [GAZ": "Sweden",                    # truncated GAZ ontology string
    "Bosnia and Herzegovina": "Bosnia and Herzegovina",
    "Ivory Coast": "Ivory Coast",
    "Democratic Republic of the Congo": "Democratic Republic of the Congo",
}
# territories shaded and counted as part of another state
MERGE_INTO = {"Taiwan": "China"}
# too small for the 110m polygon set; drawn as dots at (lon, lat)
SMALL = {"Singapore": (103.8, 1.35), "Malta": (14.4, 35.9),
         "Faroe Islands": (-6.9, 62.0), "Bahrain": (50.6, 26.1),
         "Hong Kong": (114.2, 22.3), "Luxembourg": (6.1, 49.8)}
# unresolvable free text, reported rather than mapped
DROP = {"Gaz"}
# Antarctica has no studies and Equal Earth smears it across the foot of the
# frame, so it is not drawn. Nothing is hidden: its count is zero.
SKIP_POLY = {"Antarctica"}

A1, A2, A3, A4 = 1.340264, -0.081106, 0.000893, 0.003796


def equal_earth(lon, lat):
    """Equal Earth projection (Savric, Patterson & Jenny 2019), radians in."""
    lon, lat = np.radians(lon), np.radians(lat)
    th = np.arcsin(np.sqrt(3) / 2 * np.sin(lat))
    t2 = th * th
    x = (2 * np.sqrt(3) * lon * np.cos(th)
         / (3 * (9 * A4 * t2**4 + 7 * A3 * t2**3 + 3 * A2 * t2 + A1)))
    y = th * (A1 + A2 * t2 + t2**3 * (A3 + A4 * t2))
    return x, y


def rings(geom):
    """Every outer ring of a Polygon / MultiPolygon, as (lon, lat) arrays."""
    if geom["type"] == "Polygon":
        polys = [geom["coordinates"]]
    elif geom["type"] == "MultiPolygon":
        polys = geom["coordinates"]
    else:
        return []
    out = []
    for poly in polys:
        if poly:
            a = np.asarray(poly[0], dtype=float)
            if a.ndim == 2 and len(a) >= 3:
                out.append(a)
    return out


def counts():
    at = pd.read_csv(ATTR, sep="\t", low_memory=False)
    mr = pd.read_csv(REVIEWED, sep="\t", low_memory=False,
                     keep_default_na=False, dtype=str)
    d = at.merge(mr[["record_id", "in_analysis"]], on="record_id", how="left")
    d = d[d["in_analysis"].eq("yes")]
    c = d["country"].replace("", np.nan).dropna()
    c = c[~c.isin(DROP)].map(lambda x: ALIAS.get(x, x))
    c = c.map(lambda x: MERGE_INTO.get(x, x))
    return c.value_counts(), len(d)


def main():
    vc, n_total = counts()
    gj = json.load(open(GEO, encoding="utf8"))

    verts, vals = [], []
    unmatched = dict(vc)
    for f in gj["features"]:
        name = f["properties"]["ADMIN"]
        if name in SKIP_POLY:
            continue
        name = MERGE_INTO.get(name, name)
        n = int(vc.get(name, 0))
        unmatched.pop(name, None)
        for ring in rings(f["geometry"]):
            x, y = equal_earth(ring[:, 0], ring[:, 1])
            verts.append(np.column_stack([x, y]))
            vals.append(n)
    vals = np.asarray(vals, dtype=float)

    fig, ax = plt.subplots(figsize=(11.6, 5.6))
    pos = vals > 0
    norm = LogNorm(vmin=1, vmax=max(vals.max(), 2))
    # zero-count countries first, then the data layer
    ax.add_collection(PolyCollection([v for v, p in zip(verts, pos) if not p],
                                     facecolors=NO_DATA, edgecolors=BORDER,
                                     linewidths=0.28, zorder=1))
    pc = PolyCollection([v for v, p in zip(verts, pos) if p],
                        array=vals[pos], cmap=CMAP, norm=norm,
                        edgecolors=BORDER, linewidths=0.28, zorder=2)
    ax.add_collection(pc)

    # states with studies but no polygon at 1:110m
    smalls = [(k, v) for k, v in unmatched.items() if k in SMALL]
    for name, n in smalls:
        lon, lat = SMALL[name]
        x, y = equal_earth(np.array([lon]), np.array([lat]))
        ax.scatter(x, y, s=26, facecolor=plt.get_cmap(CMAP)(norm(n)),
                   edgecolor="#404040", linewidth=0.6, zorder=4)

    allv = np.vstack(verts)
    ax.set_xlim(allv[:, 0].min() * 1.02, allv[:, 0].max() * 1.02)
    ax.set_ylim(allv[:, 1].min() * 1.06, allv[:, 1].max() * 1.06)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_facecolor(OCEAN)

    cax = fig.add_axes([0.33, 0.075, 0.34, 0.024])
    cb = fig.colorbar(pc, cax=cax, orientation="horizontal")
    cb.set_label("Studies per country (log scale)", fontsize=9.5, color=INK)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8.5, color=MUTED, labelcolor=MUTED)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.16)

    os.makedirs(FIGS, exist_ok=True)
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"MMC2_Fig_S1.{e}"), dpi=300,
                    facecolor="white")
    plt.close(fig)

    mapped = int(sum(v for k, v in vc.items() if k not in unmatched))
    top = vc.head(2)
    print(f"wrote {FIGS}/MMC2_Fig_S1.pdf / .png")
    print(f"  countries with >=1 study: {int((vc > 0).sum())}")
    for a, b in MERGE_INTO.items():
        print(f"  {a} shaded and counted as part of {b} "
              f"({b} total {int(vc.get(b, 0)):,})")
    print(f"  studies placed on the map: {mapped:,} of {int(vc.sum()):,} with a country")
    print(f"  top two ({', '.join(top.index)}) = {100 * top.sum() / vc.sum():.0f}%")
    if smalls:
        print(f"  drawn as dots (no 110m polygon): "
              f"{', '.join(f'{k} {v}' for k, v in smalls)}")
    still = {k: v for k, v in unmatched.items() if k not in SMALL}
    if still:
        print(f"  STILL UNMATCHED: {still}")


if __name__ == "__main__":
    main()
