# Pre-registration v1 — Part I (RQ1) and Part II (RQ2)
Thesis: *Certified Adaptive Examination of Peripheral Blood Smears.* Author: Qudrat Ullah. Supervisor: Balázs Harangi.
Status: **DRAFT — to be frozen with `git tag prereg-v1` before the confirmatory run.** Part III (RQ3) and RQ4
will be added as amendment v2, written before those analyses are run.

## 0. What is already known, and why this document exists
An exploratory run (`docs/killtests_v1.md`, 12,000 episodes, seed family 20260919) has already been performed and
its results are reported as exploratory. This pre-registration fixes the protocol for a **confirmatory** run on a
disjoint seed family. Where an exploratory result motivates a hypothesis, that is stated explicitly. No
confirmatory claim will rest on the exploratory run.

## 1. Estimands and units
- **Specimen**: a finite population of N = 200 fields with exact ground truth.
- **Field**: a bag of m cells drawn from KU-Optofil; represented by the **mean** of its cells' embeddings (primary
  specification). Rationale: a field is observed as a whole and no cell-level detector is assumed. Sensitivity:
  max-pooling and single-cell (nearest-to-mean) representations.
- **Critical field**: a field containing at least one blast.
- **Rarity** k ∈ {1, 2, 5, 10, 20}: number of critical fields per specimen.
- **Dilution** m ∈ {1, 5, 10, 20}: cells per field.
- **Premature stop (miss)**: the examination stopped without any critical field having been examined.
- **Critical-evidence recall**: fraction of the k critical fields examined at the stop point.
- **Observation cost**: number of fields examined.

## 2. Hypotheses (confirmatory)
- **H1a (content-blind stopping).** The InfoDPP-PAC stop point does not depend on rarity. Operationalised as
  equivalence: the median stop count differs by **at most 2 fields (1% of N)** across all k, within each m.
  *Directional prediction from the equations*: Γ̂ (Eq. 6) has no relevance term, so any dependence is indirect.
- **H1b (revised; the exploratory H1b of the proposal is abandoned).** The miss-rate difference
  Δ(m, k) = miss(InfoDPP-PAC) − miss(uniform at matched cost) is governed by **dilution**, not rarity:
  Δ(1, k) < Δ(20, k) for every k, i.e. the advantage shrinks as fields hold more cells.
  *The proposal's original H1b* ("the gap widens as evidence becomes rarer") *is pre-registered as rejected on
  exploratory evidence and is not re-tested as a confirmatory hypothesis.*
- **H1c (no protection under dilution).** At m ≥ 10 and k = 1, InfoDPP-PAC and matched-cost uniform sampling are
  **equivalent within ±0.05** in miss rate.
- **H1d (Part II agreement, RQ2).** The empirical uniform miss rate at k = 1 equals the finite-population value
  1 − n/N within Monte-Carlo error, where n is the realised cost.
- **KT1 (signal separability).** The information gain and the teacher relevance are distinguishable as predictors
  of criticality: AUC(information gain) < AUC(relevance) at every m, and AUC(information gain) → 0.5 as m grows.

## 3. Design
Five patient-grouped folds (`scripts/02_make_folds.py`, 11 blast patients, no patient in two folds). Background
cells from the blast-free patients of the evaluation fold (primary; "single host patient" as sensitivity).
**200 episodes per (fold, m, k)** → 1,000 episodes per (m, k), 20,000 in total. Baselines evaluated at the
**same cost** as InfoDPP-PAC within the same episode (paired). Exhaustive inspection is the reference (miss = 0
by construction). Confidence stopping is reported at its own cost, not matched.

**Confirmatory seeds**: `--seed 20260401`, disjoint from the exploratory family 20260919. Any episode from the
exploratory run is excluded from confirmatory analysis.

## 4. Fixed analysis plan
| Hypothesis | Test | Decision rule |
|---|---|---|
| H1a | median stop per (m, k); Kruskal–Wallis reported for completeness | supported if the median range across k is ≤ 2 fields within each m; a significant p with a range ≤ 2 is reported as "significant, negligible" |
| H1b | Δ(1, k) and Δ(20, k) with fold-clustered bootstrap (2,000 resamples of folds) | supported if the 95% CI of Δ(20, k) − Δ(1, k) excludes 0 and is positive, for a majority of k |
| H1c | exact McNemar on paired miss indicators; CI for the difference of correlated proportions | equivalence declared if the 95% CI ⊆ [−0.05, +0.05] |
| H1d | empirical uniform miss at k = 1 vs 1 − n/N | agreement if \|empirical − predicted\| < 0.02 at every (m) |
| KT1 | per-episode AUC of information gain and of relevance | supported if mean AUC(info) < mean AUC(relevance) at every m and AUC(info) at m = 20 ∈ [0.45, 0.55] |

All proportions are reported with Wilson 95% intervals; all cross-fold quantities with a **bootstrap clustered on
fold** (folds are the unit of resampling because blast patients are nested in folds).

## 5. Fixed InfoDPP-PAC configuration
Paper values: λ1 = 1.75, λ2 = 0.25, α1 = 0.5, α2 = 1.5 (Sec. 6.2). Not reported in the paper and therefore fixed
here: RBF kernel, median-heuristic length-scale per specimen, σn² = 0.01, δ = 0.05, ρ = 1, m0 = max(5, 0.01n).
τ is **calibrated per (fold, m) on blast-free specimens from the other folds** to the paper's reported operating
point (median stop ≈ 16.3% of N). τ is never tuned on evaluation episodes or on any outcome.
**Claim restriction:** because ρ and τ are ours, the implementation reproduces the paper's *operating point*, not
its nominal PAC guarantee. No statement of the form "the certificate guarantees ..." will be made about our runs.

## 6. Pre-specified sensitivity analyses (reported whatever they show)
1. Field representation: max-pooling and single-cell, at all (m, k).
2. τ sweep: operating points 5%, 10%, 16.3%, 30%, 50% of N.
3. Background pool: single host patient instead of fold pool.
4. Patient-level (rather than fold-level) clustered bootstrap.

## 7. Exclusions and stopping rules
No episode is excluded after the fact. Episodes with k > number of blast cells available in the fold are not
generated. The analysis runs once, after the full grid completes; no interim looks and no optional stopping.

## 8. What this design cannot establish
It cannot establish clinical diagnostic validity; composites have no spatial structure, so spatial dependence and
stain drift are outside its scope; critical evidence comes from 11 patients, bounding external validity; results
are conditional on one frozen encoder (DINOv2 ViT-S/14) and on mean-pooled field representation.

---
Frozen: 2026-09-26 (git tag prereg-v1)
Supervisor informed: 2026-09-23 
