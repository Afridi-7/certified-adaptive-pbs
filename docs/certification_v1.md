# Part III — certified threshold decisions (RQ3), exploratory run v1 (2026-09-23)

**Status: EXPLORATORY** (pre-registration amendment v2 for Part III is not yet written).
Code: `src/cadpbs/alpha.py`, `scripts/07_certify.py`. Results: `results/certify_theta010/`, `results/certify_theta020/`.

## Setup
- **Population**: one specimen = N = 100 fields × m = 5 cells = **M = 500 cells** (the upper end of the CLSI/NCCLS H20-A convention of 100–500 cells).
- **Estimand**: the blast proportion p over those 500 cells. **Assertion**: p > θ.
- **θ is stated per run and never merged**: θ = 0.10 (the ICC 2022 MDS/AML boundary) and θ = 0.20 (the WHO 5th edition 2022 boundary for most AML contexts) are run separately. This certifies a statement about a finite population of cells; it does not establish diagnostic validity.
- **Assorter**: a = b/(2θ) ∈ [0, 1/(2θ)], so mean(a) > 1/2 ⟺ p > θ.
- **Test**: ALPHA (Stark, *Ann. Appl. Stat.* 17(1), 2023, DOI 10.1214/22-AOAS1646), a non-negative supermartingale for a finite-population mean under sampling without replacement; certify when T ≥ 1/α, α = 0.05. Anytime-valid, so stopping may happen at any draw.
- **Unit and verification**: a field is sampled; an expert verifies **all** its cells; the observation is the field's mean assorter. The classifier is used only to *choose* fields, never to supply labels, so validity does not depend on classifier accuracy.
- 500 episodes per (θ, p), 5 folds × 100.

## Validation of the implementation
On synthetic populations at θ = 0.10, the rejection rate is 0.006 at p = 0.05, 0.016 at p = 0.08 and 0.023 at p = θ exactly — at or below the nominal α = 0.05 — and 1.000 for p ≥ 0.11. An earlier version of the code rejected 92% of the time at the exact null; the cause was treating an exhausted null mean as a contradiction rather than as a consistent outcome. Fixed and retested.

## Results (certification rate [95% CI], and median cells examined when certified)

**θ = 0.10**

| true p | uniform SRSWOR | cells | prioritised + HT correction | cells | prioritised, **no** correction | 100-cell count calls "above" |
|---|---|---|---|---|---|---|
| 0.05 | 0.000 [0.00,0.01] | – | 0.000 | – | 0.028 | 0.004 |
| 0.08 | 0.002 [0.00,0.01] | 25 | 0.000 | – | **0.550** | 0.144 |
| 0.10 (at θ) | **0.034** [0.02,0.05] | 110 | 0.002 | 20 | **0.996** | 0.412 |
| 0.12 | 1.000 | **395** | 0.020 | 87 | 1.000 | 0.708 |
| 0.15 | 1.000 | 160 | 0.222 | 70 | 1.000 | 0.914 |
| 0.20 | 1.000 | 62 | 0.972 | 60 | 1.000 | 0.998 |

**θ = 0.20**: same pattern. At the exact null, uniform 0.024, prioritised-uncorrected **0.906**. At p = 0.22 the uniform procedure needs a median of **440 of 500 cells**; at p = 0.30, 105 cells. The 100-cell count calls "above" 42.2% of the time when p = θ exactly.

## Findings
1. **Certification cost is governed by the margin |p − θ|, not by the blast count.** One fifth above the boundary (p = 0.12 vs θ = 0.10) costs a median of 395 of 500 cells; double the boundary costs 62. Near the boundary, certification approaches a census — the same phenomenon as Part II, now for a proportion rather than for absence.
2. **Adaptive selection without probability correction destroys validity, and the size of the failure is large.** Prioritising fields by classifier score and then testing as if the sample were uniform certifies a *false* assertion in **99.6%** of episodes at θ = 0.10 (90.6% at θ = 0.20) and in 55% of episodes when p is well below the boundary. This is the single clearest empirical illustration in the thesis of the difference between an informative sampler and a valid one.
3. **Correcting for the selection probabilities restores validity but costs most of the power.** Our Horvitz-Thompson-corrected variant holds its size (0.002 at the exact null) but certifies only 2% of episodes at p = 0.12 and 22% at p = 0.15, where uniform sampling certifies all of them. **Adaptivity did not help certification here.** Whether a better weighted procedure closes this gap is an open question to check against Shekhar et al. (UAI 2023, arXiv:2305.06884) and Spertus, Sridhar & Stark (arXiv:2409.06680); our variant is a construction validated empirically, not a proven procedure.
4. **The conventional fixed counts are fast and unguaranteed.** A 100-cell differential calls "above threshold" 41–42% of the time when the truth sits exactly on the boundary, 14.4% when p = 0.08 against θ = 0.10, and 70.8% when p = 0.12. This is Rümke's variability made concrete in decision terms: the convention answers quickly, with no control of the error rate near the boundary, which is where the decision matters.

## Threats to validity
- **Pseudo-replication of blasts.** A specimen at p = 0.20 needs 100 blast cells but the folds hold 65–126, so the higher proportions reuse blast cells (reuse factor 0.27–1.62 of the fold pool; with replacement above 1.0). Duplicate cells share their embedding and classifier score, which makes the specimen less variable than a real one. Mitigation for the confirmatory run: cap p per fold by the available blasts, or reduce M.
- **Our HT-weighted variant is not a proven procedure** (point 3).
- Composites have no spatial structure, so nothing here tests spatial dependence.
- Verification is simulated by revealing ground-truth labels: it assumes a perfect expert and no verification cost asymmetry.
