# Part III — certified threshold decisions (RQ3), exploratory run v2 (2026-09-25)

**Status: EXPLORATORY.** Supersedes `certification_v1.md` for adaptive sampling. The uniform-sampling results in v1
stand and are reproduced here.

## Correction to v1 (read first)
The v1 "prioritised + HT correction" design was **not a valid construction**. It fed importance-weighted draws
into ALPHA's null-mean update, which assumes unweighted draws. Its size happened to be conservative, but its
power numbers — and v1's conclusion that *"adaptivity did not help certification"* — are **withdrawn**. The
error was found by reading Shekhar et al. (below) and is logged in `DEVIATIONS.md`.

## Method now used for adaptive sampling
Shekhar, Xu, Lipton, Liang & Ramdas, *Risk-limiting Financial Audits via Weighted Sampling without Replacement*,
arXiv:2305.06884v1 (2023). [The proposal cites it as UAI 2023; the arXiv PDF does not state a venue, so **the
venue still needs verifying**.] Their construction (Sec. 2, Thm 1): observation Z_t = π(I_t) f(I_t)/q_t(I_t);
null updated with the **raw** values, μ_t(m) = m − Σ_{j<t} π(I_j) f(I_j); wealth W_t = W_{t−1}(1 + λ_t(Z_t − μ_t)).
This is a test martingale under **any** predictable sampling distribution q_t, which is exactly what adaptive
examination needs. We use it one-sided (λ_t ≥ 0) for H0: p ≤ θ, with their ApproxKelly bet, their logical bound,
their prop-MS sampling (q_t ∝ side information) and optionally their control variate (Thm 2).
Side information = the cross-fitted classifier's predicted blast content of each field; labels still come only
from expert verification of sampled fields. `src/cadpbs/weighted_betting.py`.

**Validity check on synthetic populations**, including side information built to steer *away* from blasts:
rejection at the exact null 0.003–0.013 across all designs (nominal 0.05).

## Results (M = 500 cells = 100 fields × 5; α = 0.05; 5 folds)
Certification rate at or below the boundary [95% CI]; median cells examined to certify above it.

**θ = 0.10 (ICC 2022 MDS/AML boundary), 300 episodes per p**

| p | ALPHA, uniform | betting, uniform | betting, classifier-driven | + control variate | classifier-driven, **uncorrected** | 100-cell count calls "above" |
|---|---|---|---|---|---|---|
| 0.08 | 0.000 | 0.003 | 0.003 | 0.000 | **0.563** | 0.127 |
| 0.10 (at θ) | 0.020 [0.01,0.04] | 0.013 [0.01,0.03] | 0.003 [0.00,0.02] | 0.000 | **0.997** | 0.403 |
| 0.12 | 400 cells | 420 | **220** | 220 | (80, invalid) | 0.707 |
| 0.15 | 180 | 252 | **75** | 100 | (50, invalid) | 0.917 |
| 0.20 | 60 | 105 | **40** | 60 | (35, invalid) | 0.997 |

**θ = 0.20 (WHO 5th ed. 2022 boundary), 250 episodes per p**

| p | ALPHA, uniform | betting, uniform | betting, classifier-driven | + control variate | uncorrected | 100-cell count |
|---|---|---|---|---|---|---|
| 0.18 | 0.004 | 0.004 | 0.000 | 0.000 | **0.260** | 0.244 |
| 0.20 (at θ) | 0.024 [0.01,0.05] | 0.016 | 0.000 | 0.000 | **0.872** | 0.436 |
| 0.22 | 435 | 455 | **335** | 335 | (172, invalid) | 0.652 |
| 0.25 | 255 | 340 | **165** | 200 | (112, invalid) | 0.896 |
| 0.30 | 105 | 130 | **60** | 85 | (70, invalid) | 0.988 |

## Findings (revised)
1. **Certification cost is governed by the margin |p − θ|.** Unchanged from v1: at one fifth above the boundary,
   uniform sampling needs 400–435 of 500 cells; at 1.5–2× the boundary, 60–105.
2. **Adaptive selection helps — provided the estimator accounts for it.** Classifier-driven sampling with the
   importance-weighted test stays valid (0.000–0.003 at the exact null) and cuts the cost of certification by
   **45% at p = 0.12 (400 → 220 cells)**, 58% at p = 0.15 and 35% at p = 0.25 (θ = 0.20). This **reverses**
   v1's finding 3.
3. **The same selection without the correction is invalid**: 99.7% false certification at the exact null
   (θ = 0.10), 87.2% (θ = 0.20), and 56% when p = 0.08, well below the boundary. Same samples, same classifier;
   the only difference is whether the analysis knows how the sample was drawn.
4. **The control variate added nothing here** (equal or worse cost than prop-MS alone): once the classifier drives
   the sampling, most of its information is already used, and the extra term adds variance. Shekhar et al.
   present control variates as helpful when side information is correlated with f; our result does not
   contradict that, it shows the gain is not additive with prop-MS in this setting.
5. **Near the boundary, adaptivity does not escape the census.** Even the best valid design needs 220 of 500 cells
   at p = 0.12 vs θ = 0.10, and 335 at p = 0.22 vs θ = 0.20.
6. **Conventional counts:** a 100-cell count calls "above" 40–44% of the time when p sits exactly on the boundary.

## Connection to Part I
Part I found that selecting by **information gain** does not find rare evidence once it is diluted. Part III finds
that selecting by **predicted relevance** does reduce certification cost — *if* the selection probabilities are
carried into the test. The contrast is the thesis in one line: the adaptive signal must be relevance, not novelty,
and the guarantee must come from the sampling design, not from the selector's confidence.

## Open threats
- Blast reuse at high p (see v1). The confirmatory run must cap p per fold.
- Side information is a cross-fitted classifier whose recall ranges 0.385–1.000 across patients (RQ4). How much
  of the 45% saving survives on the worst patients is not yet measured, and is the next analysis.
- Spertus, Sridhar & Stark (arXiv:2409.06680v5, 2026; submitted to *Ann. Appl. Stat.*) give union-of-intersection
  tests for **stratified** sequential sampling. They are the principled route for real smears with regional
  heterogeneity (RQ4), because validity then needs random sampling only *within* strata, not exchangeability
  across the whole smear. Not needed for composites, which have no regions.
