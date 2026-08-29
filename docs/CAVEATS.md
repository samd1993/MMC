# Caveats

Ordered by how much damage each one does if you miss it.

## 1. Three kinds of "empty", and only one of them is healthy

| `label_method` | what the record says | treat as |
|---|---|---|
| `no_value` | nothing at all was written | **unknown** |
| `uninformative_value` | `not applicable`, `unknown`, `not collected` | **unknown** |
| `explicit_absence` | `none`, `no`, `0`, `absent` | **healthy** |

Only the third is evidence of health. `not applicable` looks like a negative and
is not: 3,611 of its 4,868 uses are American Gut's generic blank marker (Lynn,
reviewing MMC2-04547: *"says 'not applicable' for sex and age of every sample"*).
Any pipeline that fills unknowns into the control arm trains on noise. If you
need more controls, take them from the 130 `healthy_cohort = yes` studies.

**Watch the field name, not just the value.** `is_healthy = no` means the subject
is *not* healthy — a case. The master handles this (`health_field_negated`), but
if you write your own parser against a coded column, check which way it points.

This is a change from the published supplementary file. `MMC2_harmonized_T1T2_samples.tsv`
as originally built folded absence into missing and then filled the blank from
the study's MeSH heading, so **~6,800 samples carried a disease their own record
denies** — PRJEB31817 alone contributed 468 samples answering `none` to
`gastrointestinal tract disorder`, counted as GI-disease cases. Rebuilt from
`build_sample_master.py`, Tier 1+2 cases fall 77,489 → 70,686 and controls rise
33,301 → 35,662. The data here is the corrected version.

## 2. `mapping_method` decides how much a MONDO term is worth

Over the 186,888 Tier 1+2 rows:

| mapping_method | share | strength |
|---|---|---|
| `sample_value_exact` | 8.4% | the sample's own value mapped to MONDO. Trust it. |
| `survey_response_positive_field_exact` | 3.5% | questionnaire; the field name is the condition. Good. |
| `control_value` + `survey_response_negative` + `explicit_absence` | 18.9% | a control statement. Good. |
| `study_mesh_heading` | 25.9% | sample had a value, it didn't map, so the **study's** disease was used |
| `study_mesh_heading_no_sample_value` | 26.3% | sample had **no** value; the study's disease was used anyway |
| `unmapped` / `no_disease_value` / `uninformative_value` | 16.8% | no term |

The `study_mesh_heading*` rows are study-level attributions, not per-sample
facts. For a supervised model, `sample_value_exact` plus the control classes is
the honest training set.

Tier 3 and 4 rows get **no** study-level fallback at all — their defining
property is that the record does not say. They come back
`tier3_no_annotation` / `tier4_no_annotation`.

## 3. Deduplicate on `sample_accession`

371,509 rows, 348,692 distinct samples. Project accessions cited by more than
one paper — PRJNA629344 under two studies, PRJNA763023 under three, PRJEB31817
under two — get a full copy of their samples per citing study. Grouping by
`record_id` double-counts; group by `sample_accession`.

## 4. Study sample counts are repository counts

`n_samples_in_repository` in `studies.tsv` is every sample under every accession
the paper cites, including datasets it reused rather than generated. Median
ratio of repository count to the paper's own stated *n*, where both are known:
**3.45× in Tier 1, 2.22× in Tier 2, 2.59× in Tier 3**. MMC2-03277 cites seven
BioProjects for a study whose abstract says *n* = 94; the repository count is
3,371. Do not read these as cohort sizes.

## 5. Free-text and coded disease fields need per-study work

236 of 565 Tier 1+2 studies have `disease_field_is_free_text = yes`. What is
actually in there:

- `sample_alias`: `FD` vs `L` = donor vs long-term care patient
- `host_disease`: `0` = healthy control, `1` = BPH
- `sample_cog_status`: opaque codes, curator noted "but codes"
- `experiment_title`: `CRC vs. healthy subject`

`label_samples.py` marks these `needs_parsing` and stops. The curator's note in
the upstream master usually says how to decode them.

## 6. 336 Tier 3/4 samples carry a real disease value

`unexpected_annotation = yes` in `samples.tsv.gz` — 336 samples across 9 studies
whose record holds a mappable disease value even though the study was tiered as
having none. Largest: MMC2-02990 (100, Crohn disease), MMC2-02917 (95, COPD),
MMC2-03343 (77, cystic fibrosis). These are candidate tiering errors, not a
usable label set. Flag them back to Sam rather than using them.

## 7. 24 studies are not finally adjudicated

`needs_curator_check` in `studies.tsv` is non-empty for 24 studies where the
reviewer was unsure, left an empty note, or where a study-level label was
recorded as if per-sample. Their tier may still move. Exclude them from anything
you intend to publish.

## 8. Sequencing type is mixed, and so is its provenance

Tier 1+2 is 414 amplicon, 101 shotgun, 50 both — not interchangeable inputs, so
filter `sequencing_type` before batching a processing run. Note that it and
`body_site` come from ENA where an accession resolves and from the paper's text
otherwise, so both are confounded with tier: fine for deciding what to download,
not for comparing tiers against each other.
