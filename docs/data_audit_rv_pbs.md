# RV-PBS audit (2026-09-21) — consequences for RQ4 and Part III

Source: RV-PBS (Ramakrishna Vivekananda Peripheral Blood Smear), the dataset of Sen et al., *Expert Systems with Applications* (2024), doi 10.1016/j.eswa.2024.123660 (GitHub `Jimut123/RV-PBS`; a Figshare mirror exists). Leishman stain, 1000× oil immersion, photographed with an iPhone XR. Local copy: `data/rv_pbs/`.

## What is actually in it
| Property | Value |
|---|---|
| Annotated images | **583** (4032×3024 for 575; 3302×2476 for 8) |
| Unannotated images | 140 (`TEST_MIXED/`) |
| Annotation format | CVAT XML, **polygon instance masks** with full pixel coordinates |
| Annotated cell instances | **752** total, mean **1.3 per image** (median 1, max 9) |
| Class totals | blast 153, neutrophil 116, myelocyte 95, eosinophil 93, band 63, lymphocyte 62, promyelocyte 56, monocyte 54, metamyelocyte 32, basophil 28 |
| Images containing an annotated blast | 68 (all in `BLAST CELLS/`) |

## The decisive limitation: annotation is class-partial, not exhaustive
Each image sits in a class folder and **only cells of that folder's class are annotated** (mean distinct classes per image = 1.00, across all 583 images). A 4032×3024 field at 1000× contains many leukocytes, but on average only 1.3 are labelled, and never from more than one class.

**This is a methodological vulnerability for the proposal's plan, because:**
1. **The field-level blast proportion is unknown.** Part III's estimand — the blast fraction over a finite population of cells — cannot be computed on any RV-PBS field, and nor can a threshold decision be scored against ground truth.
2. **Counterfactual auditing of unexamined fields is impossible.** Unannotated cells are indistinguishable from absent cells, so "what was missed" cannot be recovered.
3. **Spatial analyses are biased.** Polygon coordinates do exist, so within-class positions are recoverable, but the annotated set is a class-conditioned, non-random subset of the cells actually present. Moran's I or any clustering statistic computed on it measures the annotation process as much as the biology. It does **not** support the RQ4 claim about spatial dependence of abnormal cells.
4. `TEST_MIXED` has no annotations at all, so it cannot serve as a mixed-composition evaluation set as published.

The proposal's framing of RV-PBS as "underpowered but external" understates the issue: the problem is not power, it is that the **estimand is not identified** on these data.

## What RV-PBS can still support
- **Acquisition / stain shift (RQ4, partially).** These are real Leishman-stained fields imaged with a phone camera, a genuinely different acquisition pipeline from KU-Optofil. Cropping the 752 annotated instances gives a labelled out-of-distribution set (including 153 blasts) on which the cross-fitted classifier's error can be measured. That directly tests whether classifier error is stable across sources, which is the assumption behind prediction-powered inference.
- **Morphology-dependent error.** The class list covers the immature-granulocyte continuum (promyelocyte, myelocyte, metamyelocyte, band), which is exactly where the KU-Optofil analysis found the blast confusions (`docs/separability_v1.md`).
- **A realistic picture of within-field geometry** for illustration and for setting the composite generator's field size, but not as ground truth.

## Options for RQ4 (to settle with the supervisor)
1. **Exhaustive re-annotation of a subset.** 20–40 fields fully annotated (every WBC, any class) would restore the estimand and make a genuine field-level population. The dataset authors explicitly invite haematologist contributions. This is the only option that makes Part III testable on real fields.
2. **Detector-derived pseudo-population.** Run a segmentation model to locate all WBCs, treat the output as the population, and verify a subset by expert. This is circular if the same classifier is being audited, and the circularity must be stated; it is defensible only if the verification sample is expert-labelled and the certificate is conditioned on the detector's recall.
3. **Narrow RQ4 to what is identified.** Test the assumptions that do not need exhaustive counts: classifier-error heterogeneity across sources (KU-Optofil vs RV-PBS vs PBC-DIB), error vs morphology class, error vs image quality within class, and error correlation across cells of the same patient. Report spatial dependence as **not testable** with available data and say why.

Recommended: option 3 as the guaranteed deliverable, with option 1 pursued in parallel if a haematologist can be found. Option 2 only as a clearly-labelled sensitivity analysis.
