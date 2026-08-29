# Tier definitions

A tier describes what the **repository record** carries, not the quality of the
science. It is assigned per study and inherited by every sample in it.

**Tier 1** — per-sample disease, host age and host sex are all recorded in the
repository. 185 studies.

**Tier 2** — per-sample disease is recoverable, either from a metadata field or
from an informative sample name/alias/title. Age and sex may be absent.
380 studies.

**Tier 3** — a valid accession with unique sample identifiers, but no per-sample
biological annotation. The reads are downloadable and usable; the labels are
not there. 607 studies.

**Tier 4** — no accession at all, or an accession whose samples carry no
differentiating identifiers. 1,974 studies. Nothing to fetch.

## How a tier was decided

An automated pass over the ENA/SRA record proposed a tier for all 3,146
studies. Three curators then reviewed 1,322 of them against the live repository
record and the paper. Reviewers are the final authority; the automated tier is
superseded wherever they disagree.

The rules that turn a reviewer's verdict into a tier, applied 2026-08-27:

- A reviewer's **"no" is a statement about disease.** Host age and sex were
  independently correct ~99% of the time and are trusted.
- **"no" on a Tier 1 or 2 study → Tier 3**, unless the note is really naming
  *where* the disease lives (then the tier stands and the field is recorded), or
  objects only to age/sex (then Tier 1 → Tier 2, since Tier 1 requires both).
- **"no" on a Tier 3 study → Tier 1** if host age and sex are both present,
  otherwise **Tier 2**.
- A label that is **the same on every sample** ("all patients had COPD") is a
  study-level fact, not per-sample annotation, and never promotes a tier.
- On the Tier 4 sheet, the reviewer's own call stands with nothing inferred.

`tier_source` in `studies.tsv` records which of these fired for each study.
