# mmc-sample-labels

Reusable per-sample **healthy vs disease** labels for the human microbiome
literature, 2012–2025, plus the tooling to pull the reads they point at.

This is the data-engineering half of the Microbiome Metadata Crisis (MMC)
project. The paper asks *how much of the published human microbiome corpus can
actually be reused for a disease-associated analysis*; this repo is the answer
in a form you can run against. **3,146 studies** were screened and tiered by
hand, and for the 565 that carry per-sample biological annotation we recorded
**which metadata field the disease label lives in** — the column you need to
read to turn a pile of FASTQs into a labelled training set.

---

## The one idea to take away

There is no standard column for disease. Across these studies the per-sample
disease label lives in `host_disease`, or `gastrointestinal tract disorder`, or
`env_medium`, or `CASECTL2`, or the sample's own alias with no separate field at
all. Three curators read every study and wrote down where it actually is. That
is [`data/disease_field_map.tsv`](data/disease_field_map.tsv), and it is what
makes this corpus labellable at scale.

---

## Quick start

```bash
python3 scripts/fetch_ena_metadata.py data/accessions_tier1_tier2.txt -o ena
```
```bash
python3 scripts/label_samples.py --ena ena -o sample_labels.tsv
```
```bash
python3 scripts/download_fastq.py sample_labels.tsv --label disease,healthy --dry-run
```

Stage 1 takes a few hours for the full Tier 1+2 list (783 project accessions,
~187k samples) and is cached and resumable — kill it and rerun. Stage 3 always
reports the byte total before it moves anything; drop `--dry-run` when the
number looks sane. Python 3.9+, standard library only for the fetch and
download; `pandas` for the data files.

---

## What is in `data/`

| file | rows | what it is |
|---|---|---|
| `studies.tsv` | 3,146 | every screened study: tier, accessions, disease field, body site, sequencing type, country |
| `disease_field_map.tsv` | 569 | study → the metadata field carrying per-sample disease. **Start here.** |
| `samples_tier1_tier2.tsv.gz` | 186,673 | one row per deposited sample, disease mapped to MONDO with DOID/MeSH xrefs, plus host age and sex |
| `accessions_tier1_tier2.txt` | 783 | INSDC project accessions with per-sample disease annotation |
| `accessions_tier3.txt` | 671 | valid deposits with **no** per-sample annotation — reads are usable, labels are not |

### The tiers

| tier | studies | samples in repo | meaning |
|---|---|---|---|
| **1** | 185 | 95,187 | per-sample disease **and** host age **and** sex in the repository |
| **2** | 380 | 91,946 | per-sample disease recoverable, from a field or an informative sample name |
| **3** | 607 | 185,115 | valid accession, unique sample IDs, **no** per-sample biology |
| **4** | 1,974 | — | no accession, or an accession whose samples don't differentiate |

Tiers 1 and 2 are the labellable corpus. Tier 3 is where you go if you need more
reads and are willing to supply labels another way. Tier 4 has nothing to fetch.

Of the 565 Tier 1+2 studies, **236 have `disease_field_is_free_text = yes`** —
the label is inside a sample alias, title or description, or in a coded column
(`0=HC, 1=BPH`). Those need a per-study parser, not a dictionary lookup;
`label_samples.py` returns them as `needs_parsing` rather than guessing.
**130 studies are `healthy_cohort = yes`**: every sample is a control.

---

## What `label_samples.py` gives you

One row per sample per citing study:

```
record_id  tier  study_accession  sample_accession  run_accession
expected_disease_field  field_used  raw_value  label  label_method
host_age  host_sex  fastq_ftp  fastq_bytes
```

`label` is `disease` / `healthy` / `needs_parsing` / `unknown`, and
`label_method` says how it was decided, so you can keep or drop each route:

| method | means |
|---|---|
| `sample_value` | the study's disease field held a real value on this sample |
| `explicit_absence` | the field says the condition is **not** present → control |
| `control_value` | the value matched the control vocabulary (`healthy`, `HC`, `NC`, …) |
| `survey_response_negative` | questionnaire cohort; the subject answered no |
| `survey_response_positive:<field>` | questionnaire cohort; the **field name** is the condition |
| `healthy_cohort_study` | the whole study is healthy participants |
| `free_text_or_coded_field` | the label is in a name or a code — parse it yourself |
| `no_value` | nothing recorded. **Not** healthy. |

---

## Read this before you trust a label

Long version in [`docs/CAVEATS.md`](docs/CAVEATS.md). The four that will bite you:

1. **Missing ≠ healthy.** A sample with no disease value is `unknown`. Folding
   those into the control arm is the fastest way to a meaningless model.

2. **Half the MONDO terms in `samples_tier1_tier2.tsv` are study-level.**
   `mapping_method` is the column that tells you: only **8.8%** of rows are
   `sample_value_exact` — the disease read off the sample itself. **54.1%** are
   `study_mesh_heading*`, meaning the study's own MeSH disease heading was
   broadcast to every sample in it. That is a much weaker claim. Filter on
   `mapping_method` before you use `mondo_id` as ground truth.

3. **Samples are double-counted across studies.** 186,673 rows cover 174,812
   distinct sample accessions: 40 project accessions are cited by more than one
   paper, usually a re-analysis alongside the original. Deduplicate on
   `sample_accession`, not on row count, before computing anything.

4. **A study's sample count is the whole BioProject**, not what the paper
   analysed. Where both numbers exist, the repository holds a median 2–3.5×
   more samples than the paper's own stated *n*, because papers cite the
   projects they reused as well as the one they deposited.

---

## Provenance

Tiers come from a three-curator manual review (Isabella, Lynn, Sam) of 1,322
studies on top of an automated pass over all 3,146. `studies.tsv` carries
`tier_source` for every row — the rule that set it — and `needs_curator_check`
marks the **24 studies still awaiting adjudication**. The governing principle:
a reviewer's "no" is a statement about **disease**; host age and sex were
verified correct ~99% of the time and were trusted.

Upstream working directory is `MMC/MMC2_resubmission/` (not public). Regenerate
this data pack with `scripts/export_collaborator_repo.py` there.

Questions on the curation → Sam Degregori. Questions on a specific study → check
its `reviewer_notes` in the upstream master before assuming the field is wrong.
