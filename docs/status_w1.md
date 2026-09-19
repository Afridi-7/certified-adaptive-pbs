# Thesis status — Week 1 (2026-09-19)

Scope this session: W1–2 of the timeline (source verification, derivation checks, pre-registration inputs) plus a first synthetic pilot of Kill-test 2 (RQ1). Code: `sim/pilot.py`, `verify/part2_and_budget.py`, `verify/mcnemar_power.py`.

## 1. Source verification

| Item | Status | Notes |
|---|---|---|
| InfoDPP-PAC, arXiv:2608.23574 | Exists. Mittal, Srivastava, Chauhan, submitted 7 Jul 2026. No code link on the arXiv page. | **Eqs. 1–8 NOT yet transcribed** (fetch was rate-limited all session). The abstract describes: teacher-seeded GP relevance, log-det DPP diversity, submodular greedy, concentration-based adaptive stopping; evaluated on **gastrointestinal WSIs**, not blood smears; "83.7% fewer patches … retaining 97.9% of full-budget composite selection quality". The quotes used in the proposal (certificate "does not bound the generalisation error of a classifier"; stopping uses only the IG term; coverage 0.169 vs 0.266) still have to be checked against the full text. |
| KU-Optofil, Zenodo 10.5281/zenodo.17333317 | Exists. Yarıkan et al., CC BY 4.0, v2 15 Oct 2025. Files: `dataset.zip` (564.7 MB), `metadata.csv` (1.7 MB), `metadata_with_patient_level_splits.csv` (2.3 MB). | Descriptor: Yarıkan et al., *Sci. Data* (2026), doi 10.1038/s41597-026-06761-y. The cell count (31,489), the 13 classes, the crop size, and whether field coordinates exist have **not yet been confirmed** from the files. Zenodo is blocked from the cloud sandbox, so download it on the local machine. A patient-level split file exists: **use it**, because it avoids patient leakage between classifier training and the composite-field evaluation. |
| RV-PBS, Shekhar et al., Spertus et al., Csillag et al., Rümke | Not yet checked | Next session. |

**Design consequence (important):** KU-Optofil appears to be a *single-cell* dataset. "Fields" are therefore **composites** that we build. Composite fields contain no real spatial dependence, stain-drift structure or slide-level clustering, so RQ4 (assumption violations) **cannot** be answered on KU-Optofil. RQ4 has to rely on RV-PBS or other real whole-field data, and that makes the "secondary, underpowered" external set more important than the proposal currently suggests. This is a methodological vulnerability that should be addressed in the pre-registration.

## 2. Derivation checks (Part II, RQ2)

- **k = 1.** Under uniform sampling without replacement, P(miss) = 1 − n/N exactly, so n ≥ ⌈(1−α)N⌉. Checked by Monte Carlo (N ∈ {100, 500, 1000}, α ∈ {0.05, 0.01}); all values agree to 3 decimals.
- **k > 1.** P(all k missed) = C(N−k, n)/C(N, n); this is correct. Minimal n at α = 0.05, N = 1000: k = 1 → 950, 2 → 777, 3 → 631, 5 → 450, 10 → 258, 20 → 138, 50 → 57. A good closed form is n ≈ N(1 − α^{1/k}), which comes from the with-replacement approximation and is slightly conservative.
- **Addition for the thesis (a stronger statement than the proposal makes).** If field features carry no information about location, *every* strategy that examines n fields, adaptive or not, has miss probability 1 − n/N averaged over a uniformly random critical position. Any deterministic order has worst-case miss probability 1. Uniform random sampling is therefore minimax-optimal (a Yao-type argument), which makes "near-exhaustive" a lower bound on all methods, not a property of uniform sampling. **This needs a written proof.**
- **Caveat.** The hypergeometric bound needs k, or a lower bound on k. When k is unknown, certifying *absence* falls back to the k = 1 regime. This is exactly why Part III reformulates the estimand.
- **The α = 0.01 claim.** It is reproduced **if** the calibration uses a conformal-risk-control-type correction (n/(n+1))·R̂ + B/(n+1) ≤ α with B = 1. The usable budget α − (1−α)/n is negative for n < 1/α − 1 = 99. At α = 0.05 the threshold is n ≥ 19. **This is an inference about which correction the proposal means. The procedure should be named explicitly in the proposal/pre-registration.**

## 3. Power analysis

The paired-McNemar range of **226–463 could not be reproduced**, because the proposal does not state the discordance assumptions (p₁₀, p₀₁). Connor (1987) at two-sided α = 0.05 and 80% power gives, for example, (0.10, 0.03) → 206, (0.10, 0.04) → 303, (0.08, 0.03) → 343. The *exact* McNemar test at the Connor n reaches only about 0.77 power in simulation, so add about 10%. **Action:** state p₁₀ and p₀₁ explicitly and justify them from the pilot.

## 4. Pilot: Kill-test 2 on abstract synthetic fields (RQ1). Not InfoDPP-PAC.

Setup: 500 fields (20 × 25). d = 4 features. The critical fields are shifted by *s* along one axis (s = 0 means invisible to the encoder). The infosat rule is a generic GP greedy max-variance selector that stops after 5 consecutive information gains below ε. ε was calibrated on negative specimens to about 16% of fields, matching the "83.7% fewer" operating point. Comparisons are paired at matched cost against uniform and raster, with 200 episodes per cell, Wilson 95% CIs, and exact McNemar tests.

Results (premature-stop = miss rate):

- **The stop time is content-blind.** The median is 83 fields (16.6%) for every k and every s. The halting statistic reacts to the geometry of the feature space, not to the presence of critical evidence. This is the Kill-test-1 mechanism in its simplest form.
- **Critical evidence invisible to the encoder (s = 0).** infosat ≈ uniform at matched cost. At k = 1 the miss rate is 0.86 [0.81, 0.90], compared with the exact uniform value of 0.834. The failure is large in absolute terms but is **not** worse than uniform.
- **Critical fields are outliers in feature space (s = 1.5, 3).** infosat is **much better** than matched-cost uniform, because log-det diversity preferentially picks outliers. At k = 1, s = 3 the miss rate is 0.23 vs 0.87 (McNemar b/c = 7/134). This result survives 2–20% distractor outliers in other directions.
- **Rarity gradient.** For all methods, miss rates rise steeply as k falls, and the absolute gap to exhaustive widens with rarity. The *relative* gap to uniform, however, does **not** widen as the proposal hypothesises. It depends on s, not on k.

**Implication.** The proposal's causal claim ("a rare cell is not representative, so a novelty rule stops before it") holds only if critical cells are **not** outliers in the embedding used by the selector. That is an empirical property of the encoder on real blast morphology, not a logical necessity. The primary hypothesis should be split:
- H1a: stop time is independent of critical-evidence content. Supported in the pilot.
- H1b: miss rate relative to matched-cost uniform increases with rarity. Not supported in the pilot; it depends on feature separability.
- The real-embedding experiment has to measure *s* directly, i.e. the separability of blasts from the other 12 classes in the selector's embedding.

## 5. Next steps
1. Transcribe InfoDPP-PAC Eqs. 1–8 and replace the generic infosat rule with a faithful reimplementation, documenting every discrepancy.
2. Download KU-Optofil locally, confirm the class list/counts and the patient splits, compute embeddings, and estimate blast separability *s* and the rate of distractor classes.
3. Choose the real-field dataset(s) for RQ4 and confirm the RV-PBS details.
4. Write the pre-registration: H1a/H1b, rarity levels k/N, matched-cost definition, McNemar with stated p₁₀/p₀₁, α = 0.05, the named risk-control procedure.
5. Write the proof of the minimax k = 1 lower bound.
