# Kill-tests 1 and 2 on composite specimens — exploratory run v1 (2026-09-21)

> **CORRECTION (2026-09-26).** The grid was run in two batches (folds 0–1, then 2–4) and the second batch
> overwrote the combined `episodes.csv`, so every number below originally came from **folds 2–4 only (12,000 of
> 20,000 episodes)**. The script now combines all per-fold files, and the tables below are recomputed from all
> **20,000 episodes (5 folds, 1,000 per cell)**. The qualitative conclusions are unchanged; the main numeric
> changes are in the KT2 table (e.g. m = 10, k = 1: 0.784 vs 0.811, difference −0.027, p = 0.15; previously
> 0.805 vs 0.798) and the KT1 AUCs (0.775 / 0.589 / 0.543 / 0.527).

**Status: EXPLORATORY.** `PREREG.md` is not yet frozen, so these results inform the pre-registration; they are
not the confirmatory analysis. The confirmatory run must use fresh seeds under the frozen protocol and be
reported separately.

Code: `src/cadpbs/infodpp.py` (faithful Eqs. 3–7 of arXiv:2608.23574), `src/cadpbs/composites.py`,
`scripts/05_killtests.py`, `scripts/06_killtest_analysis.py`. Results: `results/killtests_v1/`.

**Setup.** N = 200 fields per specimen; rarity k ∈ {1,2,5,10,20} critical fields; dilution m ∈ {1,5,10,20} cells
per field; mean-pooled field embeddings; background cells from the blast-free patients of the evaluation fold;
blasts from that fold's blast patients; 5 folds × 200 episodes per cell = **20,000 episodes**. τ is calibrated on
blast-free specimens from *other* folds to the paper's operating point (~16.3% of the budget); all other
InfoDPP-PAC hyperparameters are the paper's tuned values, except those the paper does not report (see DEVIATIONS).
Baselines are evaluated at the **same cost** as InfoDPP-PAC in the same episode (paired, exact McNemar).

## H1a — the stop point is content-blind
Median stop = **33 of 200 fields (16.5%) in every one of the 20 (m, k) cells**. Kruskal–Wallis of stop step against
rarity: p = 0.000 at m = 1 (median range 33–34), p = 0.247, 0.603, 0.935 at m = 5, 10, 20. At m = 1 the dependence
is statistically detectable with 3,000 episodes but is **one field wide**; under any dilution it disappears.
**H1a is supported**, and this is what the equations predict: Γ̂ (Eq. 6) contains no relevance term.

## KT1 — the three signals are separable, and the halting signal decays with dilution
AUC for critical vs non-critical fields, mean over 3,000 episodes per m:

| m | information gain (Eq. 5, before any selection) | teacher relevance (Eq. 3) | corr(info, relevance) |
|---|---|---|---|
| 1 | 0.775 | 0.977 | 0.085 |
| 5 | 0.589 | 0.941 | 0.025 |
| 10 | 0.543 | 0.919 | 0.015 |
| 20 | 0.527 | 0.894 | 0.010 |

Information gain and relevance are **nearly orthogonal** (r ≤ 0.08), and the information signal falls to chance
(0.52) at m = 20 while relevance stays high (0.89). **Kill-test 1 passes**: the signal the certificate is built on
carries almost no information about critical evidence once the observation unit holds more than one cell.

## KT2 — premature-stop rate at matched cost (selected cells; full table in `results/killtests_v1/summary_kt2.csv`)

| m | k | InfoDPP-PAC miss [95% CI] | uniform @ same cost | difference | McNemar b/c | p |
|---|---|---|---|---|---|---|
| 1 | 1 | 0.323 [0.295, 0.353] | 0.838 [0.814, 0.860] | −0.515 | 42/557 | 8e-116 |
| 1 | 5 | 0.008 [0.004, 0.016] | 0.421 [0.391, 0.452] | −0.413 | 5/418 | 1e-116 |
| 5 | 1 | 0.691 [0.662, 0.719] | 0.841 [0.817, 0.862] | −0.150 | 112/262 | 6e-15 |
| 10 | 1 | 0.784 [0.757, 0.808] | 0.811 [0.786, 0.834] | **−0.027** | 146/173 | 0.15 |
| 20 | 1 | 0.795 [0.769, 0.819] | 0.818 [0.793, 0.841] | −0.023 | 142/165 | 0.21 |
| 20 | 2 | 0.615 [0.584, 0.645] | 0.693 [0.664, 0.721] | −0.078 | 194/272 | 3e-04 |

**Findings.**
1. **The apparent advantage of information-saturation selection is an artefact of the observation unit.** At m = 1
   (a "field" is a single cell) InfoDPP-PAC misses far less than matched-cost uniform sampling, because a blast is
   then a genuine embedding-space outlier and the diversity term seeks outliers. At m ≥ 10 the advantage is
   reduced to a few points and no longer significant (m = 10, k = 1: 0.784 vs 0.811, p = 0.15; m = 20: 0.795 vs 0.818, p = 0.21).
2. **Neither method protects against rare evidence.** At m ≥ 10 and k = 1, both miss ~80% of positive specimens
   while examining 16.5% of the fields — which is exactly the finite-population prediction 1 − n/N = 0.835 of
   Part II. The empirical uniform miss rates (0.798–0.855 at k = 1) match that bound across all m, an independent
   check of the Part II derivation.
3. **Critical-evidence recall of InfoDPP-PAC collapses to the uniform rate under dilution**: 0.60 at m = 1, k = 1
   down to 0.18–0.21 at m = 20, against uniform's 0.165 ≈ n/N.
4. **Revised H1b.** The gap against matched-cost uniform is governed by **dilution**, not by rarity: it is large at
   m = 1 for every k, and absent at m ≥ 10 for k = 1. Rarity drives the *absolute* miss rate of every method. The
   proposal's original H1b ("the gap widens as evidence becomes rarer") is **not supported** in this form.
5. **Confidence stopping** is cheaper (median 8–41 fields) and misses more (0.78–0.88 at m ≥ 10, k = 1). It is a
   different failure: it stops when the classifier is confident, which is guaranteed early when nothing abnormal
   has yet been seen.

## Threats to validity of this run
- **Field representation.** Mean pooling is the primary specification and is the direct cause of dilution. Max and
  single-cell pooling are implemented (`pooling=` in `make_specimen`) and must be run as a sensitivity analysis.
- **τ is our calibration, not the paper's.** The paper reports no τ, δ, ρ, kernel or length-scale. Stop cost is
  therefore pinned to the paper's *reported operating point*, not to its actual configuration. A τ sweep is needed.
- **ρ = 1 by fiat** in ε_{k,n}. Since ε enters only through τ + ε and τ is calibrated to a target cost, the
  calibration absorbs ρ; but the nominal PAC guarantee is then not the paper's guarantee, and must not be claimed.
- **Eleven blast patients.** Fold-clustered ranges are reported; the patient-level bootstrap should be repeated here.
- **Composites have no spatial structure**, so "fixed_order" is an arbitrary order, not a raster scan.
