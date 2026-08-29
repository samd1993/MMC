# Caveats

Ordered by how much damage each one does if you miss it.

## 1. Missing is not healthy

`label = unknown` means nothing was recorded. It does **not** mean the subject
was well. Any pipeline that fills unknowns into the control arm will train on
noise. If you need controls and there aren't enough, take them from the 130
`healthy_cohort = yes` studies instead.

The published harmonized file gets this wrong in one direction and this repo
fixes it: `MMC2_harmonized_T1T2_samples.tsv` treats the literal value `none` as
missing, then fills the blank with the study's disease heading. **1,678 samples
across 19 studies carry a disease their own record denies** — largest offenders
MMC2-02971 and MMC2-09276 (468 each, the same 468 samples via a shared
accession), MMC2-02929 (164), MMC2-03092 (112). `label_samples.py` classifies
those as `explicit_absence` → healthy. If you use the harmonized file directly
rather than the labeller, exclude those 19 studies or re-read the raw field.

## 2. `mapping_method` decides how much a MONDO term is worth

In `samples_tier1_tier2.tsv.gz`:

| mapping_method | share | strength |
|---|---|---|
| `sample_value_exact` | 8.8% | the sample's own value mapped to MONDO. Trust it. |
| `survey_response_positive_field_exact` | 3.5% | questionnaire; the field name is the condition. Good. |
| `control_value` / `survey_response_negative` | 17.8% | a control statement. Good. |
| `study_mesh_heading` | 29.1% | sample had a value, it didn't map, so the **study's** disease was used |
| `study_mesh_heading_no_sample_value` | 25.0% | sample had **no** value; the study's disease was used anyway |
| `unmapped` / `no_disease_value` | 15.3% | no term |

The bottom three are study-level attributions, not per-sample facts. For a
supervised model, `sample_value_exact` + the control classes is the honest
training set; the `study_mesh_heading*` rows are at best weak labels.

## 3. Deduplicate on `sample_accession`

186,673 rows, 174,812 distinct samples. 40 project accessions are cited by more
than one paper — PRJNA629344 appears under two studies, PRJNA763023 under three
— and each citing study gets a full copy of the samples. Grouping by
`record_id` double-counts; group by `sample_accession`.

## 4. Study sample counts are repository counts

`n_samples_in_repository` is every sample under every accession the paper cites,
including datasets it reused rather than generated. Median ratio of repository
count to the paper's own stated *n*, where both are known: **3.45× in Tier 1,
2.22× in Tier 2, 2.59× in Tier 3**. MMC2-03277 cites seven BioProjects for a
study whose abstract says *n* = 94; the repository count is 3,371. Do not read
these as cohort sizes.

## 5. Free-text and coded disease fields need per-study work

236 of 565 Tier 1+2 studies have `disease_field_is_free_text = yes`. Examples of
what is actually in there:

- `sample_alias`: `FD` vs `L` = donor vs long-term care patient
- `host_disease`: `0` = healthy control, `1` = BPH
- `sample_cog_status`: opaque codes, curator noted "but codes"
- `experiment_title`: `CRC vs. healthy subject`

`label_samples.py` marks these `needs_parsing` and stops. The curator's note in
the upstream master usually says how to decode them — worth reading before
writing a parser.

## 6. 24 studies are not finally adjudicated

`needs_curator_check` in `studies.tsv` is non-empty for 24 studies where the
reviewer was unsure, left an empty note, or where a study-level label was
recorded as if per-sample. Their tier may still move. Exclude them from anything
you intend to publish, or ask before using them.

## 7. Sequencing type is mixed

Tier 1+2 is 414 amplicon, 101 shotgun, 50 both. Amplicon and shotgun are not
interchangeable inputs — filter `sequencing_type` before batching a processing
run.

## 8. `body_site` and `sequencing_type` are mixed provenance

Both come from ENA where an accession resolves and from the paper's text
otherwise, so they are confounded with tier. Fine for filtering what to
download; do not use them to compare tiers against each other.
