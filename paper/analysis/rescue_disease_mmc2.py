"""Recover Tier 2 evidence from informative sample NAMES.

MMC1 defines Tier 2 as "per-sample disease present, via a metadata field **or an
informative sample name**". The key survey only inspects attribute *keys*, so a
study whose biology lives entirely in aliases like `case16` / `control44` scores
as Tier 3 -- which undercounts Tier 2.

This pass reads the *values* of the naming fields and looks for case/control and
disease/healthy contrasts, using the pattern lists from MMC1's
`MMC/scripts/enaharm/scripts/rescue_disease.py` verbatim so the two corpora are
scored by the same definition.

Evidence requires a CONTRAST, not a mention: either both a healthy-token and a
disease-token sample exist, or two or more distinct disease tokens appear. One
token repeated across every sample is a study-level label dragged down, and is
rejected the same way the key survey rejects a single-valued column.

Augments MMC2_study_column_map.csv in place with disease_from_names columns.
"""

import csv
import json
import os
import re
import sys
from collections import defaultdict

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mmc2_common import ROOT  # noqa: E402

csv.field_size_limit(sys.maxsize)

SAMPLES = os.path.join(ROOT, "MMC2_ena_samples.csv")
COLMAP = os.path.join(ROOT, "MMC2_study_column_map.csv")

NAME_FIELDS = ["sample_alias", "sample_title", "sample_title_xml",
               "sample_description", "library_name", "experiment_title",
               "experiment_alias", "run_alias"]

# --- verbatim from enaharm/scripts/rescue_disease.py ------------------------
# Letter-only boundaries so '_Tumor' matches ('\b' fails: '_' is a word char).
B = r"(?<![A-Za-z])"
E = r"(?![A-Za-z])"

HEALTHY_PATTS = [
    B + r"control(?:s|ed|led)?" + E, B + r"ctrl" + E, B + r"ctl" + E,
    B + r"healthy" + E, B + r"normal" + E, B + r"donor" + E,
    B + r"unaffected" + E, B + r"negative" + E,
    B + r"non[\s_\-]?(?:asthma|allerg|diabet|smoker|cancer|disease|tumor|tumour|ad|aad|ibd|cd|uc|ms|crc|crs|als|nafld|hpv|hiv|ppt|msm)" + E,
    B + r"(?:hc|nc|hv)[_\s\-]?\d", B + r"(?:hc|nc|hv)\d",
    B + r"placebo" + E, B + r"baseline" + E,
    B + r"vehicle" + E, B + r"sham" + E,
    B + r"(?:con|cont|ctl|ctr)[_\s\-]?\d",
    B + r"(?:con|cont|ctl|ctr)\d",
    B + r"ova" + E,
]
DISEASE_PATTS = [
    B + r"cases?" + E, B + r"patient(?:s)?" + E,
    B + r"tumou?r(?:s|ous)?" + E, B + r"cancer(?:s|ous)?" + E,
    B + r"carcinoma" + E, B + r"malignant" + E, B + r"neoplasm" + E,
    B + r"diseased?" + E, B + r"infected" + E, B + r"infection" + E,
    B + r"positive" + E,
    B + r"severe" + E, B + r"moderate" + E, B + r"mild" + E,
    B + r"lesional" + E, B + r"lesion" + E, B + r"symptomatic" + E,
    B + r"affected" + E, B + r"sick" + E, B + r"diagnosed" + E,
    B + r"(?:pt|pts)[_\s\-]?\d",
    B + r"inflam(?:ed|mation|matory)" + E,
    B + r"obes(?:e|ity)" + E, B + r"overweight" + E,
    B + r"smoker" + E,
    B + r"diabet(?:ic|es)" + E,
    B + r"nafld" + E, B + r"nash" + E,
    B + r"melanoma" + E, B + r"eczema" + E,
    B + r"crohn" + E, B + r"celiac" + E,
    B + r"lam" + E, B + r"ppt" + E,
    B + r"caox" + E, B + r"caries" + E,
    B + r"(?:cd|uc|ibd|ad|als|ms|crc|copd|hpv|hiv|sars|crs|t1d|t2d|nash|aad|asd|tb|eae|mog|pcos|dcis|bc)[_\s\-]?\d",
    B + r"(?:cd|uc|ibd|ad|als|ms|crc|copd|hpv|hiv|sars|crs|t1d|t2d|nash|aad|asd|tb|eae|mog|pcos|dcis)\d",
    B + r"glaucoma" + E, B + r"ischemia" + E, B + r"enterocolitis" + E,
    B + r"volvulus" + E, B + r"atresia" + E, B + r"perforation" + E,
    B + r"transmitter" + E, B + r"pcos" + E,
    B + r"preterm" + E, B + r"dry eye" + E,
    B + r"ascus" + E, B + r"lsil" + E, B + r"hsil" + E,
]
HEALTHY_RE = re.compile("|".join(HEALTHY_PATTS), re.I)
DISEASE_RE = re.compile("|".join(DISEASE_PATTS), re.I)
# ---------------------------------------------------------------------------


def main():
    n_healthy = defaultdict(int)
    n_disease = defaultdict(int)
    tokens = defaultdict(set)
    n_samples = defaultdict(int)
    examples = defaultdict(list)
    # Per-sample signature = the exact set of tokens that matched. If every
    # sample carries the SAME signature the text is a study-level description
    # dragged down onto every row, not per-sample evidence -- even when it
    # happens to contain both a healthy and a disease token (e.g. a single
    # abstract-like sample_description mentioning "inpatients" and "baseline"
    # repeated 246 times). Counting healthy and disease samples independently
    # misreads that as a contrast, so compare signatures instead.
    #
    # Two signature sets, because only a capped subset of samples was XML-
    # fetched. Mixing them lets attribute AVAILABILITY masquerade as biological
    # variation: a study with one dragged-down attribute value splits into two
    # signatures purely because the uncapped samples have no attributes at all.
    # Comparing name-only across all samples, and name+attributes only among
    # samples that actually have attributes, keeps each comparison like-for-like.
    sig_names = defaultdict(set)
    sig_attrs = defaultdict(set)

    with open(SAMPLES, encoding="utf8") as fh:
        for row in csv.DictReader(fh):
            study = row["source_study"]
            n_samples[study] += 1
            name_text = " ".join(str(row.get(f) or "") for f in NAME_FIELDS)
            raw = row.get("custom_attributes") or ""
            attr_text = ""
            if raw:
                try:
                    attr_text = " ".join(str(v) for v in json.loads(raw).values())
                except Exception:  # noqa: BLE001
                    pass
            text = name_text + " " + attr_text

            def sig(t):
                return (frozenset(m.group(0).lower() for m in HEALTHY_RE.finditer(t)),
                        frozenset(m.group(0).lower() for m in DISEASE_RE.finditer(t)))

            hs, ds = sig(text)
            if hs:
                n_healthy[study] += 1
            if ds:
                n_disease[study] += 1
                tokens[study] |= ds
            nh, nd = sig(name_text)
            if nh or nd:
                sig_names[study].add((nh, nd))
            if raw and (hs or ds):
                sig_attrs[study].add((hs, ds))
            if (hs or ds) and len(examples[study]) < 4:
                examples[study].append(
                    (str(row.get("sample_alias") or row.get("sample_title") or "")[:40])
                )

    cm = pd.read_csv(COLMAP, low_memory=False)
    cm["names_healthy_n"] = cm["study"].map(n_healthy).fillna(0).astype(int)
    cm["names_disease_n"] = cm["study"].map(n_disease).fillna(0).astype(int)
    cm["names_disease_tokens"] = cm["study"].map(
        lambda s: "; ".join(sorted(tokens.get(s, ()))[:6])
    ).fillna("")
    cm["names_examples"] = cm["study"].map(
        lambda s: " | ".join(examples.get(s, ()))
    ).fillna("")

    # A real contrast means samples differ from each other. Require:
    #   (a) more than one distinct token signature across the study, AND
    #   (b) at least one sample carrying a disease token.
    # (a) alone rejects the drag-down case; (b) alone rejects studies where the
    # only variation is in healthy-side wording.
    cm["sig_from_names"] = cm["study"].map(lambda s: len(sig_names.get(s, ()))).fillna(0).astype(int)
    cm["sig_from_attrs"] = cm["study"].map(lambda s: len(sig_attrs.get(s, ()))).fillna(0).astype(int)
    cm["names_signatures"] = cm[["sig_from_names", "sig_from_attrs"]].max(axis=1)
    varies = (cm["sig_from_names"] > 1) | (cm["sig_from_attrs"] > 1)
    has_disease = cm["names_disease_n"] > 0
    cm["disease_from_names"] = varies & has_disease
    cm["names_dragdown"] = (~varies) & has_disease   # reported, not counted

    # the flag harmonize_mmc2.py consumes
    cm["disease_evidence"] = cm["disease_present"] | cm["disease_from_names"]

    cm.to_csv(COLMAP, index=False)

    n = len(cm)
    print(f"updated {COLMAP}  ({n:,} studies)")
    print(f"\n  disease via attribute key      : {int(cm['disease_present'].sum()):5,d} "
          f"({100 * cm['disease_present'].mean():5.1f}%)")
    print(f"  disease via informative names  : {int(cm['disease_from_names'].sum()):5,d} "
          f"({100 * cm['disease_from_names'].mean():5.1f}%)")
    print(f"  -> combined disease evidence   : {int(cm['disease_evidence'].sum()):5,d} "
          f"({100 * cm['disease_evidence'].mean():5.1f}%)")
    gained = int((cm["disease_from_names"] & ~cm["disease_present"]).sum())
    print(f"     rescued (names only, no key): {gained:5,d}")
    print(f"\n  rejected as drag-down (disease token, one signature on every "
          f"sample): {int(cm['names_dragdown'].sum()):,}")
    print("\nexamples of rescued studies:")
    r = cm[cm["disease_from_names"] & ~cm["disease_present"]].head(6)
    for _, x in r.iterrows():
        print(f"   {x['study']}: tokens=[{x['names_disease_tokens']}]  "
              f"names={x['names_examples'][:70]}")


if __name__ == "__main__":
    main()
