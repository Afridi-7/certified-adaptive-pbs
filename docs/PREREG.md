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

---
# Amendment 1 (v1.1) — 2026-09-26, BEFORE any confirmatory data were generated
Filed after a re-read of the frozen v1 against the exploratory data. No confirmatory episode (seed family
20260401) exists at the time of writing. The hypotheses' **content and margins are unchanged**; the amendments
fix three drafting errors in decision rules, set sample sizes from a power calculation, and fix a data-handling bug.

**A0. Correction to the exploratory record.** The exploratory kill-test run generated 20,000 episodes, but a
multi-batch bug meant the analysis used only folds 2–4 (12,000 episodes). The exploratory numbers cited in §0 are
superseded by the corrected 5-fold analysis in `docs/killtests_v1.md`. `scripts/05_killtests.py` now combines all
per-fold files; the confirmatory `episodes.csv` must contain all 5 folds, which the analysis checks.

**A1. H1d decision rule replaced.** v1 required |empirical − predicted| < 0.02. At 1,000 episodes the Monte-Carlo
standard error of the uniform miss rate is ≈ 0.012, so a formula that is *exactly* correct would fail that rule by
chance with substantial probability (it failed at m = 1 and m = 10 in the exploratory data; a re-run of those cells
with an unrelated seed gave z = −0.86 and +1.10, i.e. chance). New rule: for each m, compare the empirical uniform
miss rate at k = 1 with the mean of the per-episode predictions 1 − n_i/N (n_i = realised cost of that episode),
using z = (empirical − mean prediction) / SE with the Poisson-binomial SE √(Σ p_i(1 − p_i)) / n. Agreement if
|z| < 2.50 (two-sided 5% with Bonferroni over the 4 values of m). H1d is supported if all four agree.

**A2. H1b scope aligned with its decision rule.** v1 stated "for every k" in §2 but "for a majority of k" in §4.
At k = 20 both methods' miss rates are ≈ 0 (floor), so Δ(20, 20) − Δ(1, 20) is not an informative contrast. H1b is
tested over k ∈ {1, 2, 5, 10} and supported if the fold-clustered 95% CI of Δ(20, k) − Δ(1, k) lies above 0 for at
least 3 of these 4 values; k = 20 is reported descriptively.

**A3. H1c sample size.** The equivalence margin stays ±0.05. From the corrected exploratory data at m = 10, k = 1
(difference −0.027, discordance 0.319), the projected 95% CI half-width is 0.035 at 1,000 episodes, 0.025 at 2,000
and 0.020 at 3,000, so equivalence could only be shown with ≈ 3,000 episodes if the true difference is near the
exploratory estimate. The cells (m = 10, k = 1) and (m = 20, k = 1) therefore get **600 episodes per fold (3,000 per
cell)**; all other cells keep 200 per fold. Implemented as
`--cell_episodes 10:1:600 20:1:600`.

**A4. Seed families.** Main confirmatory run: 20260401. The four pre-specified sensitivity analyses (§6) use
20260402 (max pooling), 20260403 (single-cell pooling), 20260404 (τ sweep) and 20260405 (single-host pool), so no
analysis shares random streams with another (the exploratory τ sweep reused the main run's streams, which made its
deviations from 1 − n/N correlated). Part III uses its own family (20260501+), defined in amendment 2.

**A5. Analysis code frozen with this amendment.** The confirmatory decisions are computed by
`scripts/09_confirmatory_analysis.py`, committed together with this amendment and not modified afterwards. It was
test-run on the exploratory data only to check that it executes. For transparency: on those exploratory data
(1,000 episodes per cell) it returns H1a, H1b, H1d, KT1 supported and **H1c not supported** — the CI is too wide
at that sample size, which is the reason for A3. The margin was not changed, and the confirmatory outcome of H1c
is not guaranteed by the larger sample: if the true difference is below −0.05, equivalence will still fail.

**Confirmatory command (to be run once, after tagging prereg-v1.1):**
```
python scripts/05_killtests.py --out results/confirm_v1 --seed 20260401 --episodes 200 --cell_episodes 10:1:600 20:1:600
python scripts/09_confirmatory_analysis.py --dir results/confirm_v1
```

Amendment frozen: 2026-09-26 (git tag prereg-v1.1)
