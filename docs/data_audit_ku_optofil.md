# KU-Optofil data audit (2026-09-19)

Evidence: `results/inspect_ku_optofil.txt`, `results/dups_ku_optofil.txt`, `results/folds_v1/`. Scripts: `scripts/01_*.py`, `scripts/02_*.py`.

## Verified facts
| Property | Proposal | Verified | Note |
|---|---|---|---|
| Cells | 31,489 | 31,488 jpg on disk (+1 Thumbs.db) = 31,484 in metadata + 4 `- Copy.jpg` duplicates | The proposal's 31,489 probably counts every file in the zip, Thumbs.db included. Report **31,484** labelled cells. |
| Classes | 13 | 13 | segmented_neutrophil 16,256; lymphocyte 5,762; giant_platelet 2,835; eosinophil 1,862; monocyte 1,381; reactive_lymphocyte 1,038; **blast 814**; basophil 546; metamyelocyte 262; erythroblast 243; band_neutrophil 199; platelet_cluster 188; myelocyte 98 |
| Patients | not stated | 275 identified + 1 pseudo-ID `unknown` (3,204 cells) | median 110 cells per identified patient (IQR 103–117), i.e. each patient is roughly one ~100-cell differential |
| Image size | not stated | 368×368 px (199/200 sampled) | single-cell crops |
| Spatial or field info | assumed absent | **absent**: no slide, field or coordinate columns | the filename index (`<patient>_<nnn>`) is not contiguous for 272/276 patients; what it means is unknown, so **do not** treat it as acquisition order |

## Problems found
1. **The folder split (`dataset/train|val|test`) is not patient-level.** 265 of 276 patients appear in more than one folder, so using the folders would leak patients between training and evaluation. Use the `split` column (0 patients span splits) or the folds below. Never use the folders.
2. **Blasts are concentrated in very few patients.** Only 12 IDs have any blast: 11 identified patients (489 blasts) and `unknown` (325).
   - The provided test split has only **3 blast patients (86 blasts)**, and validation has **none**.
   - Evaluating critical evidence on the provided split would therefore reuse blasts from three patients in every episode. **Methodological vulnerability:** the episodes would not be independent at the source level.
3. **`unknown` (3,204 cells, train split only)** has filename prefixes EO_, VLY_, BA_, BL_, ERB_, MY_, ERC_, MMY_, BNE_, … These cells can't be grouped by patient. *Suspicion, not verified:* the naming scheme resembles the Barcelona PBC datasets (Acevedo et al.), i.e. a second source. This needs checking against the data descriptor (Yarıkan et al., Sci. Data 2026). Until then they are **excluded**.
4. **Exact duplicates:** the only exact-content duplicates found are the 4 `- Copy` files plus their originals, which are ignored because they are not in the metadata. Near-duplicates (perceptually similar crops) have **not** been checked.

## Decision (Gate A): patient-grouped 5-fold cross-fitting
- Every identified patient goes to one fold. Blast patients are dealt round-robin by blast count; the rest are shuffled with seed 20260919.
- Each cell's embedding and classifier score come from the model trained on the other 4 folds. All 11 blast patients can then be evaluated with **no patient leakage**.

| fold | patients | cells | blasts | blast patients |
|---|---|---|---|---|
| 0 | 55 | 5,236 | 126 | 3 |
| 1 | 55 | 5,659 | 113 | 2 |
| 2 | 55 | 5,634 | 102 | 2 |
| 3 | 55 | 6,010 | 83 | 2 |
| 4 | 55 | 5,741 | 65 | 2 |

This deviates from the provided split and is logged in `DEVIATIONS.md`. The provided split is kept in `cells.csv` for a sensitivity analysis.

## Implications for the design
- **External validity (vulnerability):** critical evidence comes from 11 patients. The morphological diversity of blasts is bounded by those 11 sources. Episode-level CIs must therefore be complemented by a **patient-clustered** analysis: a cluster bootstrap over source patients, and results reported per source patient.
- **Composite fields are synthetic in spatial structure.** KU-Optofil provides no geometry at all, so spatial arrangement, clustering and stain drift in composites are *our* modelling choices. RQ4 cannot be answered on KU-Optofil; it needs RV-PBS or other whole-field data.
- **Part III:** real per-patient blast fractions exist (8 patients ≥ 20%, 9 ≥ 10%; the counts are small, ~100 cells). They can anchor realistic composite proportions near the WHO-HAEM5 (20%) and ICC 2022 (10%) boundaries.
- **Background composition:** a patient has only ~110 cells, so a specimen of N fields with m cells each must draw its background cells from several patients in the same fold. This mixes inter-patient heterogeneity into a "specimen". It is a design choice to pre-register (options: single-host with replacement, or a fold-level pool).
