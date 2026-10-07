"""Shared normalization + MeSH helpers for the MMC2 resubmission pipeline.

ROOT convention matches the MMC1 handoff: ROOT is the parent of scripts/.
"""

import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MMC = os.path.dirname(ROOT)
ONEDRIVE = os.path.dirname(MMC)

# ---------------------------------------------------------------- inputs
RAYYAN_OUT = os.path.join(MMC, "Rayyan_output_Aug6_2026.csv")
MOHAK = os.path.join(MMC, "Rayyan_mohak_accession_codes.csv")
ORIGINAL_QUERY = os.path.join(MMC, "original_query.csv")
VAGINAL_RIS = os.path.join(MMC, "claude dump", "cervicovaginal_microbiome_full.ris")
CANCER_CSV = os.path.join(
    ONEDRIVE, "Cancer Qiita", "Cancer_rayyan_output_2026-06-10_20-58-47", "cancer_articles.csv"
)
EXTRA_CSV = os.path.join(ONEDRIVE, "MMC2", "extra_add", "extra_articles.csv")
MTREES = os.path.join(ROOT, "mtrees2025.bin")

SCRAPED_COL = "Accession Code [scraped sheet] "  # trailing space is in the source header

# Placeholders that mean "no accession", not a code.
_NULL_ACC = {"", "-", "–", "—", "no accession", "none", "na", "n/a", "nan"}


# ---------------------------------------------------------------- normalizers
def norm_doi(s):
    """Lowercase, strip, drop a leading doi.org resolver prefix."""
    if pd.isna(s):
        return None
    s = re.sub(r"^https?://(dx\.)?doi\.org/", "", str(s).strip().lower())
    return s or None


def norm_title(s):
    """Aggressive title key for fuzzy cross-file matching (first 80 alnum chars)."""
    if pd.isna(s):
        return None
    return re.sub(r"[^a-z0-9]", "", str(s).lower())[:80] or None


def norm_pmid(s):
    if pd.isna(s):
        return None
    m = re.search(r"(\d{6,9})", str(s))
    return m.group(1) if m else None


def acc_set(s):
    """Parse an accession cell into a normalized frozenset of codes.

    Collapses internal whitespace so 'PRJNA 761547' == 'PRJNA761547', drops
    tokens that are too short or carry no digit (prose, 'the', 'available').

    DOI-style identifiers are kept whole. Splitting them on '/' turned
    '10.5281/zenodo.12345' into the bare registrant prefix '10.5281' plus
    'ZENODO.12345' -- and because every Zenodo deposit shares that prefix, 29
    unrelated papers looked like they shared an accession. Same for figshare
    (10.6084) and Dryad (10.5061).
    """
    if pd.isna(s):
        return frozenset()
    raw = str(s).strip()
    if raw.lower() in _NULL_ACC:
        return frozenset()
    keep = set()
    # peel off DOI-style identifiers first so the '/' split can't break them
    doi_like = re.findall(r"10\.\d{4,9}/[^\s;,|]+", raw)
    for d in doi_like:
        raw = raw.replace(d, " ")
        keep.add(d.strip(".").upper())
    for part in re.split(r"[;,/|]|\s+and\s+", raw):
        tok = re.sub(r"\s+", "", part.strip().strip(".")).upper()
        if not tok or len(tok) < 5 or not re.search(r"\d", tok):
            continue
        if re.fullmatch(r"10\.\d{4,9}", tok):   # bare DOI registrant prefix
            continue
        keep.add(tok)
    return frozenset(keep)


def fmt_acc(fs):
    return "; ".join(sorted(fs)) if fs else ""


# ---------------------------------------------------------------- MeSH
DESIGN_TERMS = [
    "Case-Control Studies",
    "Cohort Studies",
    "Prospective Studies",
    "Retrospective Studies",
    "Cross-Sectional Studies",
    "Longitudinal Studies",
    "Randomized Controlled Trials as Topic",
    "Double-Blind Method",
    "Controlled Clinical Trial",
    "Clinical Trials as Topic",
]

# Cross-check only. The frozen tree-C/F03 list is the primary criterion.
DISEASE_STEM_RX = (
    r"neoplasm|carcinom|diseases?(?:/|;|$)|syndrome|infections?(?:/|;|$)"
    r"|disorder|colitis|diabet|obesity|dysbiosis|inflammat"
)


def load_disease_headings(path=MTREES):
    """MeSH headings in tree category C (Diseases) or F03 (Mental Disorders).

    mtrees<year>.bin is one 'Heading;TreeNumber' pair per line.
    """
    heads = set()
    with open(path, encoding="utf8", errors="replace") as fh:
        for line in fh:
            if ";" not in line:
                continue
            heading, tree = line.rstrip("\n").rsplit(";", 1)
            tree = tree.strip()
            if tree.startswith("C") or tree.startswith("F03"):
                heads.add(heading.strip().lower())
    return heads


def split_mesh(cell):
    """Split a Rayyan keywords cell into bare MeSH main headings.

    'Humans;*Mouth/microbiology;Female' -> ['humans', 'mouth', 'female']
    Strips the '*' major-topic marker and any '/subheading' qualifiers.
    """
    if pd.isna(cell):
        return []
    out = []
    for part in str(cell).split(";"):
        head = part.strip().lstrip("*").split("/")[0].strip().lower()
        if head:
            out.append(head)
    return out


def mesh_flags(df, disease_headings):
    """Return (has_design, has_disease, has_mesh) boolean Series for a dataframe."""
    design_lc = {t.lower() for t in DESIGN_TERMS}
    heads = df["keywords"].map(split_mesh)
    has_design = heads.map(lambda hs: any(h in design_lc for h in hs))
    has_disease = heads.map(lambda hs: any(h in disease_headings for h in hs))
    has_mesh = heads.map(bool)
    return has_design, has_disease, has_mesh
