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

---
# Amendment 2 (v2) — Part III (RQ3) and RQ4, 2026-09-27, BEFORE their confirmatory data exist
No episode with the seeds below has been generated, and no RV-PBS embedding exists, at the time of writing.
Exploratory Part III results (`certification_v1.md`, `certification_v2.md`) and the KU-Optofil RQ4 analyses
(`rq4_v1.md`) motivated these hypotheses.

## Part III — certified threshold decisions (RQ3)
- **Population**: one specimen = N = 100 fields × m = 5 cells = M = 500 cells. **Estimand**: blast proportion p over
  those 500 cells. **Assertion**: p > θ. **Thresholds, run and reported separately, never merged**: θ = 0.10 (ICC 2022,
  MDS/AML boundary) and θ = 0.20 (WHO 5th ed. 2022, AML in most contexts). The study certifies statements about finite
  populations of cells; it establishes nothing about diagnostic validity.
- **p grid**: θ = 0.10 → {0.05, 0.08, 0.10, 0.12, 0.15, 0.20}; θ = 0.20 → {0.15, 0.18, 0.20, 0.22, 0.25, 0.30}.
  **1,000 episodes per (θ, p)**.
- **Change from exploratory runs**: blasts are drawn **without replacement from all 489 KU-Optofil blasts** and background
  without replacement from all blast-free patients (no within-fold restriction). Valid because every cell's score is
  cross-fitted; removes the blast reuse of the exploratory runs.
- **Designs, all run on the same specimen (paired)**: `alpha_uniform` (ALPHA, uniform without replacement);
  `bet_uniform`, `bet_prop` (classifier-driven, prop-MS), `bet_prop_cv` (+ control variate) — Shekhar et al.
  (arXiv:2305.06884) importance-weighted betting test; `naive` (classifier-prioritised order, analysed as if uniform —
  deliberately invalid); conventional 100- and 200-cell counts. α = 0.05. Labels come only from verification of sampled fields.
- **Seeds**: θ = 0.10 → 20260511; θ = 0.20 → 20260512 (per p: seed + round(1000·p)).
- **Hypotheses and decision rules** (implemented in `scripts/13_part3_analysis.py`, frozen with this amendment):
  - **H3a (validity)**: for each of the four valid designs, the false-certification rate at p = θ is ≤ 0.05 **and** a
    one-sided exact binomial test of H0: rate ≤ 0.05 is not rejected at 0.05. Supported if all four pass.
  - **H3b (uncorrected adaptive selection is invalid)**: at p = θ the `naive` rate's 95% Clopper–Pearson lower bound > 0.05.
  - **H3c (classifier-driven sampling is cheaper)**: at every p > θ, cells to certify under `bet_prop` < `alpha_uniform`,
    paired one-sided Wilcoxon signed-rank, Holm-corrected at 0.05 across the p values above θ. Supported if all reject.
  - **H3d (cost governed by the margin)**: median `alpha_uniform` cells strictly decrease as p increases above θ.
  - Descriptive only: `bet_prop_cv` vs `bet_prop`; rate at which conventional counts call "above" at each p.
- Commands (run once, after tagging prereg-v2):
  `python scripts/12_certify_confirm.py --out results/part3_theta010 --theta 0.10 --p_list 0.05 0.08 0.10 0.12 0.15 0.20 --episodes 1000 --seed 20260511`
  `python scripts/12_certify_confirm.py --out results/part3_theta020 --theta 0.20 --p_list 0.15 0.18 0.20 0.22 0.25 0.30 --episodes 1000 --seed 20260512`
  `python scripts/13_part3_analysis.py --dirs results/part3_theta010 results/part3_theta020`

## RQ4 — where do the assumptions break?
- **KU-Optofil error-structure analyses** (`rq4_v1.md`) were run on the complete data and cannot be re-run on fresh data;
  they are reported as **exploratory/descriptive**, not confirmatory.
- **External source (RV-PBS)**, not yet embedded: 752 annotated cell crops (153 blasts) from 583 images, cropped by
  `scripts/14_crop_rvpbs.py` (square crop, 1.5× the polygon's longer bounding-box side, resized to 368 px — a CHOICE,
  since KU-Optofil's cropping rule is undocumented; differences in cropping are a confound of any shift estimate).
  Embedded with the identical encoder and preprocessing (`notebooks/15_external_embeddings_colab.ipynb`).
  Classifier: logistic probe trained on all KU-Optofil cells, threshold fixed at the KU cross-fitted 95%-specificity
  threshold t = 0.03104 (KU pooled blast recall 0.906, false-positive rate 0.050 at t). CIs by bootstrap over source image.
  - **H4a**: blast recall on RV-PBS is lower than on KU-Optofil: upper 95% CI bound < 0.906.
  - **H4b**: false-positive rate on RV-PBS non-blast cells is higher: lower 95% CI bound > 0.050.
  - **H4c (validity under shifted side information)**: specimens of 500 RV-PBS cells, θ = 0.10, classifier-driven
    betting test; false certification at p = θ (1,000 episodes) satisfies the H3a rule. Cost at p = 0.15, 0.20 descriptive.
  - Seed 20260521. Implemented in `scripts/16_external_analysis.py`, frozen with this amendment.
  - H4a/H4b are directional expectations of acquisition shift; either outcome is reported. PBC-DIB is dropped from the plan.

Amendment 2 frozen: 2026-09-26 (git tag prereg-v2)

# Amendment 3 (v3) — spatial clustering and a defensive mixture (RQ4, RQ3), 2026-09-25, BEFORE their data exist
No episode with the seeds below has been generated at the time of writing. The scripts were smoke-tested on
**synthetic** stand-in data only; they have not been run on KU-Optofil or RV-PBS cells.

**Motivation.**
- **E3.1.** Sensitivity to spatial arrangement is a stated project requirement, and it is untested. RV-PBS cannot
  test it: it has 752 annotated cells over 583 images, about 1.3 per image.
- **E3.2.** The exploratory cost reversal under shift (`confirmatory_v2.md`: classifier-driven sampling needed a
  median of 190 cells vs 110 for uniform sampling at p = 0.20 on RV-PBS) was not a pre-registered hypothesis. E3.2
  tests whether it replicates on fresh episodes, and whether a standard safeguard removes it.

## E3.1 — blasts clustered across fields (RQ4: exchangeability and spatial dependence)
- **Population, estimand, assertion**: as in Amendment 2 (N = 100 fields × m = 5 cells = 500 cells; blast
  proportion p; assertion p > θ). θ = 0.10 (ICC 2022) and θ = 0.20 (WHO 5th ed. 2022) are run and reported
  separately. Blasts come from all 489 KU-Optofil blasts and background from blast-free patients, all sampled
  without replacement. Scores are cross-fitted.
- **Arrangements**:
  - `scattered`: the B = round(pM) blasts occupy B cells chosen uniformly at random (the Part III arrangement).
  - `hot0.6` and `hot1.0`: n_hot = min(N, ⌈B/(m·h)⌉) fields, chosen uniformly at random, receive all B blasts,
    placed uniformly at random among their cells. h ∈ {0.6, 1.0} is the nominal blast density of a hot field.
  - These are **stress models**, not estimates. No available data establish how blasts are distributed across
    fields in real smears.
- **Designs (paired on the same specimen)**:
  - `alpha_field`: ALPHA, uniform field sampling.
  - `bet_field`: Shekhar et al. test, uniform field sampling.
  - `prop_field`: Shekhar et al. test, classifier-driven field sampling (prop-MS).
  - `naive_cell`: the same uniform field order as `alpha_field`, with the cells flattened and analysed by
    cell-level ALPHA (M = 500) as if each cell had been sampled individually.
  - `conv100`: first 100 cells of that order; calls "above" if their blast fraction exceeds θ.
- **p grid**: θ = 0.10 → {0.08, 0.10, 0.12, 0.15, 0.20}; θ = 0.20 → {0.18, 0.20, 0.22, 0.25, 0.30}.
  1,000 episodes at p = θ, 500 otherwise, for each arrangement.
- **Expectation from theory, stated in advance**:
  - The three field-level designs are valid for **any** arrangement, because the sampling unit is the analysis
    unit.
  - `naive_cell` is valid only if the cell sequence is exchangeable. That holds under `scattered` (the flattened
    sequence is then a uniformly random permutation of cells) and not under clustering.
  - E3.1 therefore tests *how large* the failure is under clustering, not whether field-level validity holds in
    theory.
- **Hypotheses and decision rules** (in `scripts/20_amend3_analysis.py`, frozen with this amendment):
  - **H5a (field-level validity under clustering).** For `alpha_field`, `bet_field` and `prop_field`, in every
    arrangement, at p = θ and p = θ − 0.02, for both θ, the Amendment 2 H3a rule holds (rate ≤ 0.05 and the
    one-sided exact binomial test of rate ≤ 0.05 is not rejected at 0.05). Supported if all 36 cells pass.
    Because of the multiplicity, the rule is conservative, and any failing cell will be examined individually.
  - **H5b (mismatched analysis unit breaks validity under clustering).** At p = θ, for both θ: `naive_cell`
    passes the H3a rule under `scattered`, **and** its 95% Clopper–Pearson lower bound exceeds 0.05 under
    `hot1.0`. `hot0.6` is descriptive.
  - **H5c (clustering raises the cost of valid certification).** At p = θ + 0.05, for both θ: the median
    `alpha_field` cells are strictly ordered scattered < hot0.6 < hot1.0, **and** scattered < hot1.0 by a
    one-sided Mann–Whitney test at 0.05 (different specimens, so the test is unpaired).
  - Descriptive: `prop_field` cost by arrangement; `conv100` "above" rates at every p.
- **Seeds**: θ = 0.10 → 20260611; θ = 0.20 → 20260612. Per-episode RNG
  `default_rng([seed, arrangement index, round(1000 p), episode])`.

## E3.2 — defensive mixture under side-information shift (RQ4)
- **Sampler**: q_t(i) = (1 − λ)·S(i)/Σ_remaining S + λ/n_remaining. **λ = 0.50 is primary**; λ = 0.25 is
  descriptive. This is valid under Shekhar et al. Theorem 1 (any predictable q_t that is positive on every
  remaining unit).
- **Sources**:
  - `rvpbs`: RV-PBS crops; probe trained on all KU-Optofil cells, as in `scripts/16`.
  - `ku`: KU-Optofil; cross-fitted scores.
- **Specimen and assertion**: N = 100 × m = 5, cells shuffled (scattered); θ = 0.10 (ICC 2022) only.
- **Designs (paired)**: `uniform`, `prop`, `mix25`, `mix50`, all using the Shekhar et al. betting test.
- **p grid**: rvpbs → {0.10, 0.15, 0.20}; ku → {0.10, 0.12, 0.15, 0.20}. 1,000 episodes at p = 0.10, 500
  otherwise.
- **Hypotheses and decision rules**:
  - **H6a (validity).** `mix25` and `mix50` pass the H3a rule at p = 0.10 on both sources.
  - **H6d (the reversal replicates).** RV-PBS, p = 0.20: `prop` needs more cells than `uniform`, by a paired
    one-sided Wilcoxon signed-rank test at 0.05. This turns the exploratory observation into a confirmatory
    test on fresh episodes.
  - **H6b (the mixture caps the loss under shift).** RV-PBS, at p = 0.15 **and** p = 0.20: the upper 95%
    bootstrap bound (2,000 paired episode resamples, seed 20260631) of median(`mix50`)/median(`uniform`) is
    < 1.20.
  - **H6c (the mixture keeps most of the gain in distribution).** KU-Optofil, at p ∈ {0.12, 0.15, 0.20}: the
    retained saving (med_uniform − med_mix50)/(med_uniform − med_prop) is ≥ 0.5, **and** `mix50` < `uniform` by
    a paired one-sided Wilcoxon test, with Holm correction over the three p at 0.05. All three p must pass.
  - H6b and H6c pull in opposite directions (safety vs efficiency). Either may fail, and both outcomes are
    reported.
- **Seeds**: rvpbs → 20260621; ku → 20260622. Per-episode RNG
  `default_rng([seed, source index, round(1000 p), episode])`.

## Commands (run once, after tagging prereg-v3; the scripts resume if interrupted)
```
python scripts/18_spatial_confirm.py --out results/a3_spatial_theta010 --theta 0.10 --p_list 0.08 0.10 0.12 0.15 0.20 --seed 20260611
python scripts/18_spatial_confirm.py --out results/a3_spatial_theta020 --theta 0.20 --p_list 0.18 0.20 0.22 0.25 0.30 --seed 20260612
python scripts/19_mixture_confirm.py --out results/a3_mixture_rvpbs --source rvpbs --p_list 0.10 0.15 0.20 --seed 20260621
python scripts/19_mixture_confirm.py --out results/a3_mixture_ku --source ku --p_list 0.10 0.12 0.15 0.20 --seed 20260622
python scripts/20_amend3_analysis.py --spatial results/a3_spatial_theta010 results/a3_spatial_theta020 --mixture results/a3_mixture_rvpbs results/a3_mixture_ku --out results/a3_decisions.json
```

Amendment 3 frozen: 2026-09-25 (git tag prereg-v3)

# Amendment 4 (v4) — real-specimen spatial structure (RQ4, bridge to RQ3), 2026-09-30, BEFORE the analyses are run
**Data.** BoneMarrowWSI-PediatricLeukemia (Höfener et al.; NCI Imaging Data Commons collection
`bonemarrowwsi_pediatricleukemia`; source DOI 10.5281/zenodo.14933087; CC BY 4.0). Only the annotation series were
downloaded (`notebooks/21_bmwsi_annotations_colab.ipynb`; `data/bmwsi/bmwsi_ann.zip`, SHA-256 c8446c59…aadc6d) and parsed by
`scripts/22_parse_bmwsi.py`.

The dataset documentation states that in 228 ROIs of 2048 × 2048 px (40×, 0.11 µm/px, about 225 µm), **every
haematological structure was boxed and consensus-labelled** by at least two of three or more experts. The ROIs were
placed by the annotators in the evaluable monolayer.

**Scope change, stated openly.** These are pediatric **bone-marrow aspirate** smears, not peripheral blood. WHO 5th ed.
(2022) and ICC (2022) blast criteria refer to blood or marrow. This amendment measures spatial structure in real smears;
it does not transfer any cost number to peripheral blood.

**What was inspected before writing this amendment.** The following design facts, and nothing about spatial
arrangement:
- 113 patients with labelled ROIs (112 with 2 ROIs, 1 with 1).
- Nucleated cells per ROI: median 68 (IQR 37–109). Total 17,003 nucleated cells, 6,660 blasts.
- Per-patient blast proportion: 46 patients ≤ 0.05, 1 in (0.05, 0.10], 2 in (0.10, 0.20], 8 in (0.20, 0.50],
  56 > 0.50.
- Consensus label frequencies.

Because almost no patient lies near either threshold, real-proportion certification near the boundary cannot be
tested on these data. Hence the design below: measure the clustering on real data, then carry it into a
pre-specified simulation.

## A4.1 — is blast status spatially clustered within real specimens? (`scripts/23_bmwsi_clustering.py`)
**Definitions** (as in the script header):
- **Blast** = {lymphocytic_blast, myelocytic_blast, monocytic_blast, unknown_blast}.
- **Nucleated** = all consensus labels except artefacts, damaged/smudge cells, unfit, platelets and aggregates,
  erythrocytes, spicules, NETs, no-consensus, megakaryocytic cells, mitoses, phagocytosis, pseudo-Gaucher cells and
  macrophages.
- **Eligible ROI**: 1500–3000 px sides, ≥ 20 nucleated cells, ≥ 1 blast and ≥ 1 non-blast.

**Scales:**
- **S1 (primary)**: 2 × 2 quadrants within an ROI (about 113 µm).
- **S2**: 4 × 4 within an ROI (about 56 µm).
- **S3**: between the ROIs of a patient.

**Method:**
- **Estimator**: pooled ANOVA intraclass correlation for binary data, with subfields nested in their parent unit.
- **Uncertainty**: patient bootstrap, 2,000 resamples. Seed 20260701.
- **H7a (clustering exists at the field scale)**: the S1 ICC exceeds its permutation null. Blast labels are permuted
  among the cells of each ROI with positions fixed; 2,000 permutations; one-sided p < 0.05.
- **Estimation, no binary hypothesis**: S1/S2/S3 ICC with 95% CIs; k = 5 nearest-neighbour blast enrichment against
  its permutation null.
- **Sensitivity analyses**: promonocytes counted as blasts; damaged and smudge cells kept as non-blast nucleated cells.

## A4.2 — what does the measured clustering imply for validity? (`scripts/24_bmwsi_bridge_sim.py`; model-based)
- **Specimens**: N = 30 fields × m = 17 cells. Field propensities follow a Beta distribution with mean θ and ICC ρ,
  with the blast total fixed so that p = θ exactly.
- **ρ values**: {0, S1 lower CI bound, S1 estimate, S1 upper CI bound}, floored at 0, taken from A4.1.
- **Designs**: `alpha_field` and `naive_cell`, as in E3.1.
- **Runs**: θ = 0.10 (ICC 2022) and θ = 0.20 (WHO 5th ed.); 1,000 episodes per (θ, ρ); seed 20260711.
- **E3.1 placement**: the ICC of the E3.1 arrangements, computed with the same estimator. Seed 20260712.
- **Pre-specified reading rule**:
  - If `naive_cell` false certification at the S1 **upper** CI bound exceeds 0.05 at either θ, we report that the
    H5b failure is plausible at the clustering level observed in real smears.
  - If it is ≤ 0.05 at both θ, we report that the clustering observed in these monolayer ROIs is too weak to break
    the cell-level analysis at this field size.
  - In either case the statement is model-based. It is not an observed false-certification rate.

## Known limits, stated before the results
- ROIs were chosen by annotators in the evenly dispersed monolayer. This selection may bias measured clustering
  downward relative to the whole smear.
- Two labelled ROIs per patient; spatial scale ≤ 225 µm within an ROI.
- Pediatric bone marrow only.
- Smudge and damaged cells may hide blasts (sensitivity analysis above).
- The simulation's Beta-binomial field model is an assumption.

## Commands (run once, after tagging prereg-v4)
```
python scripts/23_bmwsi_clustering.py --out results/a4_bmwsi/clustering.json
python scripts/24_bmwsi_bridge_sim.py --clustering results/a4_bmwsi/clustering.json --out results/a4_bmwsi/bridge.json
```

Amendment 4 frozen: 2026-09-30 (git tag prereg-v4)
