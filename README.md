# MMC — the Microbiome Metadata Crisis

Data, code and tools for **Degregori et al., *iMeta* (in revision)**, which asks
how much of the published human microbiome literature can actually be reused
for a disease-associated analysis, and answers: very little.

A PRISMA search of PubMed returned **35,220 records**; screening left
**2,896 human microbiome studies**, each tiered by what its public sequence
deposit records about the people sampled (the BHAM tiers):

| tier | studies | what the deposit carries |
|---|---|---|
| **1** | 170 (5.9%) | per-sample host disease **and** age **and** sex |
| **2** | 360 (12.4%) | per-sample host disease, from a field or informative sample names |
| **3** | 539 (18.6%) | downloadable, distinguishable samples, **no** biology |
| **4** | 1,827 (63.1%) | no accession, or samples that cannot be told apart |

Of the 752,361 samples behind these studies, **645,997 (85.9%) cannot be reused**
for a disease analysis, because the metadata that would label them was never
deposited.

---

## Check a deposit before it is published

`scripts/check_tier.py` takes an accession and reports which tier the deposit
satisfies, which metadata keys supplied the evidence, and what is missing for
the next tier up. Authors, reviewers and editors can run it on a submission.

```bash
python3 scripts/check_tier.py PRJEB10878
```

Or check the per-sample sheet you are about to upload (one row per sample,
one column per attribute, `.tsv` or `.csv`):

```bash
python3 scripts/check_tier.py --table my_sample_metadata.tsv
```

`--json` gives machine-readable output. Standard library only, Python 3.9+.

It applies the paper's automated rules and vocabularies exactly: on the 1,448
ENA studies in the paper it reproduces every automated disease, age, sex and
sample-name call. In the paper, curators then reviewed the Tier 1 and 2
candidates by hand. Where a study appears in `data/`, the script also prints
its **curated** tier, which takes precedence. Two rules decide most outcomes
(details in [`docs/TIERS.md`](docs/TIERS.md)):

- A field only counts when its values **differ across samples**. "Crohn's
  disease" on every row is a study-level label, not per-sample metadata.
- In a disease field, `healthy` / `control` / `none` is information (it marks
  the control arm), but `not applicable` / `missing` is not.

---

## What is here

```
paper/                      the paper's supporting data and code
  Data_S1_studies.tsv           Data S1: curated study-level variables, 2,896 studies
  Data_S2_harmonized_samples.tsv.gz
                                Data S2: harmonized per-sample host disease, age, sex
                                for Tier 1 and 2 (124,298 samples, 491 studies)
  queries/                      full PubMed query strings and PRISMA counts
  analysis/                     the scripts behind the tiers, figures and statistics
data/                       per-sample healthy/disease labels, all screened studies
scripts/                    check_tier.py, plus fetch / label / download tooling
docs/                       tier definitions and caveats
```

### `paper/`

- **`queries/main_query.txt`** returned 33,564 records and
  **`queries/cervicovaginal_query.txt`** returned 1,656, for 35,220 combined.
  `prisma_counts.tsv` carries every count in the PRISMA flow (Figure S2).
- **Data S1** has one row per analyzed study: tier, year, journal, DOI,
  accessions, the field holding per-sample disease, disease group, body site,
  sequencing type, country and continent. For Tier 1 and 2 studies whose
  disease label sits in prose or a field that could not be pinned down,
  `disease_field` reads `free-text`.
- **Data S2** has one row per sample per citing study, with disease mapped to
  MONDO (DOID and MeSH cross-references), plus host age and sex.
  `mapping_method` says whether the disease came from the sample's own record
  or was inherited from the study's MeSH heading. Read
  [`docs/CAVEATS.md`](docs/CAVEATS.md) before treating it as ground truth.
- **`analysis/`** holds the pipeline as it ran for the paper: ENA harvest
  (`fetch_ena_mmc2.py`), key survey and tiering (`survey_keys_mmc2.py`,
  `rescue_disease_mmc2.py`, `harmonize_mmc2.py`), the sample master and
  harmonized file (`build_sample_master.py`, `harmonize_t12_samples.py`),
  sample sizes, PRISMA, figures, supplementary tables and the adjusted trend
  models. The scripts read from the authors' working directory and are
  published as a record of the method, not as a one-command rebuild.

### `data/` — reusable labels for every screened study

`data/` covers all **3,145 screened studies**, including the 249 the paper
excluded (healthy cohorts and non-INSDC deposits), because healthy cohorts are
the cleanest source of controls. Filter on `in_paper_analysis = yes` to get the
paper's 2,896. **Tier is a column, never a filename.**

| file | rows | what it is |
|---|---|---|
| `samples.tsv.gz` | 336,984 | **the master.** One row per sample per citing study, all tiers, disease mapped to MONDO with DOID/MeSH xrefs, plus host age and sex |
| `studies.tsv` | 3,145 | every screened study: tier, accessions, disease field, body site, sequencing type, country |
| `disease_field_map.tsv` | 568 | study → the metadata field carrying per-sample disease. **Start here.** |
| `accessions.tsv` | 1,567 | `accession` / `record_id` / `tier`, one list, filterable |

| tier | studies | with a resolvable deposit | sample rows |
|---|---|---|---|
| **1** | 183 | 174 | 57,051 |
| **2** | 380 | 349 | 91,941 |
| **3** | 607 | 554 | 184,572 |
| **4** | 1,975 | 15 | 3,420 |

There is no standard column for disease. Across these studies the per-sample
label lives in `host_disease`, or `gastrointestinal tract disorder`, or
`env_medium`, or `CASECTL2`, or the sample's own alias with no separate field
at all. Curators recorded where it actually is for each Tier 1 and 2 study:
[`data/disease_field_map.tsv`](data/disease_field_map.tsv).

Of the 563 Tier 1 and 2 studies, **211 have `disease_field_is_free_text = yes`**:
the label sits in a sample alias, title or description, or in a coded column
(`0=HC, 1=BPH`). Those need a per-study parser, and `label_samples.py` returns
them as `needs_parsing` rather than guessing. **49 are `healthy_cohort = yes`**,
meaning every sample is a control.

---

## Pull labelled reads

```bash
python3 scripts/fetch_ena_metadata.py data/accessions.tsv --tier 1,2 -o ena
```
```bash
python3 scripts/label_samples.py --ena ena --tier 1,2 -o sample_labels.tsv
```
```bash
python3 scripts/download_fastq.py sample_labels.tsv --label disease,healthy --dry-run
```

Drop `--tier` to take the whole corpus. Stage 1 covers 776 project accessions
for Tier 1 and 2, takes a few hours, and is cached and resumable. Stage 3
reports the byte total before it moves anything; drop `--dry-run` once the
number looks right. Fetch and download use only the standard library;
`pandas` is needed to read the data files.

`label_samples.py` writes one row per sample per citing study:

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
| `health_field_negated` | `is_healthy = no` → a **case**, not a control |
| `healthy_cohort_study` | the whole study is healthy participants |
| `free_text_or_coded_field` | the label is in a name or a code — parse it yourself |
| `uninformative_value` | the record says "not applicable"/"unknown". **Not** healthy. |
| `no_value` | nothing recorded. **Not** healthy. |

---

## Read this before you trust a label

Long version in [`docs/CAVEATS.md`](docs/CAVEATS.md). The four that will bite you:

1. **Missing ≠ healthy, and "not applicable" ≠ healthy.** Both come back
   `unknown`. Folding them into the control arm is the fastest way to a
   meaningless model.

2. **Most MONDO terms are study-level, not per-sample.** The `mapping_method`
   column tells you which. Over the Tier 1 and 2 rows, only **10.5%** are
   `sample_value_exact`, meaning the disease was read off the sample itself.
   **63%** are `study_mesh_heading*`, where the study's own MeSH heading was
   applied to every sample in it. Filter on `mapping_method` before you use
   `mondo_id` as ground truth.

3. **Samples are double-counted across studies.** 336,984 rows cover 314,167
   distinct sample accessions, because some project accessions are cited by
   more than one paper (usually a re-analysis alongside the original).
   Deduplicate on `sample_accession`, not on row count.

4. **A study's sample count is the whole BioProject**, not what the paper
   analysed, and it is often several times the paper's own stated *n*.

---

## Provenance and contact

Tiers come from an automated pass over the ENA/SRA record of every study,
followed by a three-curator manual review. 1,318 studies carry a reviewer-set
tier. `studies.tsv` records `tier_source`, the rule that set each tier, and
`needs_curator_check` marks the **24 studies still awaiting adjudication**.
A reviewer's "no" was read as a statement about **disease**: host age and sex
were verified correct about 99% of the time and were trusted.

This repository previously held the MMC1 study table and Figure 2 notebook.
Those are preserved in the `MMC1` legacy repository.

Questions → Sam Degregori (corresponding details in the paper).
