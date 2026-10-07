"""Map free-text disease strings onto MONDO, with DOID and MeSH cross-references.

MONDO is the cross-walk ontology: its terms carry xrefs to DOID, MeSH, NCIt,
OMIM and Orphanet, so one id interoperates with all of them -- including the
DOID codes that already appear in this corpus's own ENA records.

Matching runs through a cascade, because the strings come from three registers
that spell the same disease differently:

  MeSH headings     "colorectal neoplasms", "genital neoplasms, female"
  ENA free text     "Crohn's disease(tissue)", "CRC (Stage III/IV)"
  author codes      "HC", "CD", "case3"

so the cascade is: exact normalised -> singular/plural -> comma inversion
("genital neoplasms, female" -> "female genital neoplasms") -> adjective forms
(prostatic -> prostate) -> reversed noun phrase -> token-set.

NOTE the normaliser strips disease/disorder/syndrome/condition INCLUDING their
plurals. Stripping only the singular was a real bug: "inflammatory bowel
diseases" then failed to match MONDO's "inflammatory bowel disease", losing one
of the commonest conditions in the corpus.

Roughly 14% of MeSH disease headings have no MONDO equivalent at all -- MONDO
has no term labelled "bacteremia", and its only Helicobacter entry is
"Helicobacter infection, non-human animal". Those are ontology gaps, not
matching failures, and are reported as unmapped rather than forced.
"""
import json
import os
import re

_STRIP = r"\b(the|a|an|of|with|and|diseases?|disorders?|syndromes?|conditions?)\b"
_ADJ = {"prostatic": "prostate", "gastric": "stomach", "hepatic": "liver",
        "renal": "kidney", "pulmonary": "lung", "cardiac": "heart",
        "cerebral": "brain", "colonic": "colon", "oesophageal": "esophageal"}

# values that mean "no disease", not a disease name
CONTROL = re.compile(
    r"^\s*(healthy|health|control|controls|hc|ctrl|normal|nc|non[- ]?diseased|"
    r"negative|neg|uninfected|unaffected|no[nt][- ]?(case|disease|ibd|cancer)|"
    r"reference|baseline healthy)\b", re.I)
NULLISH = {"", "na", "n/a", "none", "not applicable", "not collected",
           "not provided", "missing", "unknown", "nan", "-", "null"}


class MondoMapper:
    def __init__(self, index_path):
        d = json.load(open(index_path))
        self.meta = d["meta"]
        self.index = dict(d["index"])
        # re-normalise the shipped keys so the stored index and the live
        # normaliser agree on plurals
        for k, v in list(d["index"].items()):
            nk = re.sub(r"\s+", " ", re.sub(_STRIP, " ", k)).strip()
            self.index.setdefault(nk, v)
        self.tokens = {}
        for k, v in self.index.items():
            self.tokens.setdefault(tuple(sorted(k.split())), v)

    @staticmethod
    def norm(s):
        s = str(s).lower().strip()
        s = re.sub(r"\(.*?\)", " ", s)
        s = re.sub(r"[^a-z0-9, ]+", " ", s)
        s = re.sub(_STRIP, " ", s)
        return re.sub(r"\s+", " ", s).strip()

    def _base(self, k):
        yield k
        w = k.split()
        if not w:
            return
        if w[-1].endswith("s") and len(w[-1]) > 3:
            yield " ".join(w[:-1] + [w[-1][:-1]])
        if w[-1].endswith("ies"):
            yield " ".join(w[:-1] + [w[-1][:-3] + "y"])
        sw = [_ADJ.get(x, x) for x in w]
        if sw != w:
            yield " ".join(sw)
            if sw[-1].endswith("s") and len(sw[-1]) > 3:
                yield " ".join(sw[:-1] + [sw[-1][:-1]])
        if len(w) > 1:
            yield " ".join(reversed(w))

    def _variants(self, s):
        k = self.norm(s)
        if "," in k:
            parts = [p.strip() for p in k.split(",")]
            yield from self._base(re.sub(r"\s+", " ", " ".join(reversed(parts))).strip())
        yield from self._base(re.sub(r"\s+", " ", k.replace(",", " ")).strip())

    def lookup(self, s):
        """(mondo_id, how) or (None, None)."""
        if not s or str(s).strip().lower() in NULLISH:
            return None, None
        seen = []
        for v in self._variants(s):
            if v and v not in seen:
                seen.append(v)
                if v in self.index:
                    return self.index[v], "exact"
        for v in seen:
            key = tuple(sorted(v.split()))
            if key in self.tokens:
                return self.tokens[key], "tokens"
        return None, None

    def describe(self, mondo_id):
        m = self.meta.get(mondo_id) or {}
        x = m.get("xref", {})
        return dict(mondo_label=m.get("label", ""),
                    doid=x.get("DOID", ""), mesh=x.get("MESH", ""),
                    ncit=x.get("NCIT", ""))

    @staticmethod
    def is_control(v):
        return bool(CONTROL.match(str(v).strip()))
