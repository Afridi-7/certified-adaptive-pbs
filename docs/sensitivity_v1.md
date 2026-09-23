# Sensitivity analyses for the kill-test result (exploratory, 2026-09-23)

Purpose: test whether the main Part I finding — that InfoDPP-PAC's advantage over matched-cost uniform sampling
disappears as fields hold more cells — is an artefact of two choices I had to make: the **field representation**
(mean pooling) and the **τ calibration point** (the paper's reported ~16.3% operating point).
Runs: `results/sens_pool_max`, `results/sens_pool_single`, `results/sens_tau_*` (500 episodes per cell, 5 folds).

## 1. Field representation (m = cells per field, k = critical fields, miss = premature-stop rate)

| pooling | m | k | cost | InfoDPP-PAC | uniform @ cost | difference | AUC(info gain) |
|---|---|---|---|---|---|---|---|
| mean (primary) | 10 | 1 | 33 | 0.805 | 0.798 | **+0.007** | 0.534 |
| max | 10 | 1 | 32 | 0.730 | 0.820 | −0.090 | 0.579 |
| single cell | 10 | 1 | 33 | 0.824 | 0.822 | +0.002 | 0.513 |
| mean | 20 | 1 | 33 | 0.817 | 0.833 | −0.017 | 0.522 |
| max | 20 | 1 | 33 | 0.764 | 0.814 | −0.050 | 0.567 |
| single cell | 20 | 1 | 33 | 0.814 | 0.804 | +0.010 | 0.512 |
| any | 1 | 1 | 33 | 0.344–0.397 | 0.850–0.855 | ≈ −0.5 | 0.75–0.80 |

**Conclusion: the finding is not an artefact of mean pooling.** Max pooling, which preserves the most signal of
the three, recovers only part of the advantage (difference −0.09 instead of +0.007 at m = 10) and still leaves a
**73% miss rate** on specimens with one critical field. Single-cell representation behaves like the mean. The
information-gain AUC stays in 0.51–0.58 under every representation at m ≥ 10, against 0.75–0.80 at m = 1.

## 2. τ calibration (m = 10, k = 1; the operating point sets the examination cost)

| target cost | realised cost | InfoDPP-PAC miss | uniform @ same cost | difference | 1 − n/N (Part II) |
|---|---|---|---|---|---|
| 5% | 10 | 0.922 | 0.940 | −0.018 | 0.950 |
| 10% | 20 | 0.866 | 0.892 | −0.026 | 0.900 |
| 16.3% (paper) | 33 | 0.805 | 0.798 | +0.007 | 0.835 |
| 30% | 59 | 0.610 | 0.690 | −0.080 | 0.705 |
| 50% | 100 | 0.384 | 0.480 | −0.096 | 0.500 |

**Conclusion: the finding is not an artefact of the calibration point.** Across a tenfold range of examination
cost the difference never exceeds 0.1, while the absolute miss rate is governed almost entirely by how much of
the specimen was examined. Uniform's miss rate tracks the finite-population prediction 1 − n/N at every operating
point (0.940 vs 0.950, 0.892 vs 0.900, 0.798 vs 0.835, 0.690 vs 0.705, 0.480 vs 0.500), an independent
confirmation of the Part II derivation at five different budgets (H1d).

## 3. Implication for the thesis
Adaptive selection buys, at most, a small constant reduction in miss rate at fixed cost; it does not change the
**rate** at which missing rare evidence decays with examination effort. Reducing the miss rate for a single
critical field from 0.80 to 0.05 requires examining ~95% of the specimen, whichever rule chooses the fields.
This is the empirical counterpart of the Part II bound and the reason Part III reformulates the question from
detection to certification of a proportion.

Caveat: still exploratory (`PREREG.md` is not yet frozen); episode counts here are 500 per cell rather than the
1,000 pre-registered for the confirmatory run.
