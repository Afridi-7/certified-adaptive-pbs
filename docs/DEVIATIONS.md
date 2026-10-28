# Deviations from published methods

| Date | Method | Paper spec | Our implementation | Reason | Possible consequence |
|---|---|---|---|---|---|
| 2026-09-19 | InfoDPP-PAC halting | Eqs. 1–8, arXiv:2608.23574 | `sim/pilot.py:infosat_order` is a GENERIC GP max-variance + IG-saturation rule | Eqs. not yet transcribed | Pilot results are not evidence about InfoDPP-PAC itself |
| 2026-09-19 | KU-Optofil split | provided patient-level `split` (test: 3 blast patients, val: 0) | patient-grouped 5-fold cross-fitting (`scripts/02_make_folds.py`) | too few blast patients in test for critical-evidence evaluation | the classifier/encoder differs per fold; results pooled over folds; provided split kept for sensitivity |
| 2026-09-19 | KU-Optofil cells | 31,484 labelled | 28,280 used (3,204 `unknown`-patient cells excluded) | cannot be grouped by patient, so leakage risk; source unverified | fewer blasts (489 vs 814); fewer training cells |
| 2026-09-19 | InfoDPP-PAC hyperparameters | τ, δ, ρ, kernel, ℓ, β, m0 unspecified | calibrate τ on negative training-fold specimens to about 16% of budget, then sweep | the paper does not report them | the stopping point reflects our calibration; must be shown with a sensitivity sweep |
| 2026-09-19 | InfoDPP-PAC embedding | Phikon (H&E histopathology, 768-d) | DINOv2 ViT-S/14 frozen (`vit_small_patch14_dinov2.lvd142m`, 384-d, 224 px) | no established open Wright–Giemsa encoder; Phikon is trained on H&E tissue | geometry of the embedding (hence stopping behaviour and blast separability) is encoder-specific; single encoder, no comparisons (scope) |
| 2026-09-21 | InfoDPP-PAC Eq. 7 | ρ is the Assumption-1 bounded-difference constant | ρ = 1 by fiat; τ calibrated to the paper's reported operating point (~16.3% of budget) on blast-free specimens from non-evaluation folds | neither ρ nor τ is reported in the paper | the implemented rule reproduces the paper's *cost*, not its nominal PAC guarantee; the certificate's δ must not be claimed as the paper's |
| 2026-09-21 | Field representation | not applicable (the paper's units are image patches) | a field is the mean of its cells' embeddings | composites have to define a field-level representation | mean pooling causes the dilution effect; max/single pooling are implemented for sensitivity |
| 2026-09-21 | Raster-scan baseline | spatial raster over a slide | fixed index order over composite fields | composites carry no geometry | reported as `fixed_order`, not as a spatial baseline |
| 2026-09-25 | Part III weighted test | Shekhar et al. (arXiv:2305.06884) Sec. 2 | v1 used an invalid HT-into-ALPHA construction; **withdrawn** and replaced by a transcription of Shekhar et al.'s importance-weighted betting martingale (`weighted_betting.py`) | error found on reading the paper | v1 finding 3 ("adaptivity did not help") is withdrawn; v2 shows a 35–58% cost reduction with valid size |

## Erratum (2026-09-25): date in PREREG.md Amendment 2 header
The Amendment 2 header reads "2026-09-27"; the correct date is 2026-09-25. This is a typing error in the frozen text and is
left unedited so that the text still matches tag `prereg-v2`. The authoritative freeze time is the timestamp of the `prereg-v2`
commit, which precedes the creation of `results/part3_theta010`, `results/part3_theta020` and `results/external_v1/decisions.json`.

## Additions (2026-09-25)
| Date | Component | Proposal / source | What we did | Why | Consequence |
|---|---|---|---|---|---|
| 2026-09-25 | Part III side information | Proposal lists prediction-powered e-values (Csillag, Struchiner & Goedert; arXiv identifier to be verified) | **Not implemented.** Side information enters only through Shekhar et al.'s prop-MS sampling and control variate | one valid weighted-sampling framework covers the adaptive design; prediction-powered e-values were not needed to test H3a–H3d | the study reports no comparison with prediction-powered e-values; this must be stated as a limitation, not implied as covered |
| 2026-09-25 | Expert verification | a haematologist verifies the sampled fields | verification is simulated by the dataset's expert labels (verifier error β = 0) | no haematologist annotation available | certificates hold for a perfect verifier; `part2_proof.md` Remark 4 shows that verifier error sets a floor |

**Erratum addendum.** The freeze line of Amendment 2 reads "2026-09-26". The tag was created on 2026-09-25, and the
`prereg-v2` tag timestamp is authoritative.

## Note on document dates (2026-10-03)
Several documents carry dates of 2026-09-26 or 2026-09-27: the PREREG v1/v1.1/v2 headers and freeze lines, confirmatory_v1, rq4_v1, and the killtests_v1 and sensitivity_v1 notes. These were typed by hand and may not match when the files were actually written. The **git commit and tag timestamps are authoritative**. The pre-registration claim rests on the ordering of the tags relative to the result commits, not on the typed dates.

## Audit trail and episode counts (2026-10-19)
`docs/PREREG_AUDIT.md` tabulates, for every tag, the commit time, the GitHub push time and the first confirmatory output.
The output always follows the push, by 43 s to 18 min. It also lists the typed dates that are wrong, and content hashes
of the outputs.

Episode counts reported in `confirmatory_v3.md` and in the paper draft of 2026-10-03 were wrong. The files contain the
counts that Amendment 3 specifies:
- **E3.1: 18,000 episodes.** That is 3 arrangements × 2 θ × (1,000 at p = θ + 4 × 500), against the 15,000 written.
- **E3.2: 4,500 episodes.** RV-PBS has 1,000 + 2 × 500 and KU-Optofil has 1,000 + 3 × 500, against the 5,000 written.

No analysis changes; only the reported totals change.

## Pre-registered sensitivity analyses run late (2026-10-19)
PREREG §6 and Amendment 1 (A4) fix four sensitivity analyses: max pooling (seed 20260402), single-cell pooling (20260403), a
τ sweep (20260404) and a single-host background (20260405). They had not been run with these seeds; the paper cited the
exploratory `sensitivity_v1` runs instead (500 episodes per cell, reusing the main run's random streams). They were run on
2026-10-19 with the frozen `05_killtests.py`, unchanged, and 200 episodes per (fold, m, k).
- **τ sweep scope.** The sweep was restricted to k = 1, the H1c cell, because §6 does not state its grid.
- **Patient-level bootstrap.** §6.4 could not be done from the frozen outputs, because they do not record which blast
  patient supplied each critical field. Fold-level intervals are reported instead. The Kill-test 1 re-analysis of
  2026-10-19 records patients.

## Correction (2026-10-03): RV-PBS citation
Earlier documents attributed RV-PBS to "Sen et al." and "Vergara-Lluri et al."; both attributions were wrong. The correct citation is J. B. Pal, A. Bhattacharyea, D. Banerjee, T. Maharaj, *Expert Systems with Applications* 249:123660 (2024), doi:10.1016/j.eswa.2024.123660. This was checked against the dataset's GitHub README. `data_audit_rv_pbs.md` and `runbook.md` are corrected.

## InfoDPP-PAC review corrections (2026-10-19, second pass)
- **Operating point misdescribed.** Earlier drafts called 16.3% the share of the slide that InfoDPP-PAC examines. In the
  paper it is 16.3% of a 300-patch budget (mean k⋆ = 48.9), about 2.4% of a slide's 2,075 patches. The confirmatory
  calibration (16.3% of all fields) is unchanged and favours InfoDPP-PAC about 7-fold. An exploratory run at 2.4% is in
  `results/sens_op_paperpool` (seed 20261022).
- **Notation.** The paper defines K_S through the posterior kernel (Eq. 3) but κ_ii through the prior kernel (Sec. 3). Our
  implementation uses the posterior kernel throughout; this is our reading, logged as an interpretation.
- **Pure information-gain rule** (`scripts/31_stopping_class.py`): λ2 = 0 gives 0·(−∞) = nan in the blended gain (chosen
  fields are masked with −∞), so the run uses λ2 = 1e-9. A first run with λ2 = 0 was discarded (it re-selected chosen
  fields).
- **ALPHA prior guess.** Every ALPHA run used the default η0 = (μ0 + u)/2 of `src/cadpbs/alpha.py`. This was never
  pre-registered as a choice. `scripts/39_alpha_prior.py` shows it is the best of three fixed guesses for field-level
  draws; for individually sampled cells (`scripts/32`) a guess of θ ± 0.05 is used and the default is reported alongside.

## Real-tile Part I (2026-10-20, exploratory)
- Specimens of 16–32 tiles are too small for the seed rule m0 ≈ 0.01n with our floor of 5; `scripts/36` uses m0_min = 2.
- k-center coverage stopping is not cost-matched; it is compared with uniform sampling at its own cost, h(N, k, n_kc).
- Raster = row by row within a region (tile index i*4 + j, i = column), regions in random order.
- Notebook 35 was revised before its first successful run (2026-10-20) only to rebuild missing input tables; the tables it
  produced match the repo tables (fingerprints 23d6a9d242d1e4a1 / e73ccda91a271f29).

## Background includes non-leukocytes (found 2026-10-20)
All Part I and Part III composites draw background cells from every non-blast KU-Optofil image of blast-free patients,
including giant platelets (2,725), platelet clusters (57) and erythroblasts (17): 10.3% of the pool. A leukocyte
differential would not count them, so the composites' p is the blast share of imaged objects, not of leukocytes. This was
never stated. `scripts/41_leukocyte_background.py` (seed 20261032) re-runs the Part III designs with a leukocyte-only
background: false certification ≤ 0.030; median cells 395/160 (ALPHA, θ = 0.10, p = 0.12/0.15), 225/75 (classifier-driven),
445/235 and 335/140 at θ = 0.20 — within 15 cells of the share-0 condition of `scripts/33` (largest gap: uniform betting at p = 0.15, 267.5 vs 252.5). The paper now states the composition.
- MLL23 (`scripts/38`): the 2026-10-19 exclusion matched "atypical_promyelocyte", but the release's class is
  "promyelocyte_atypical", so atypical promyelocytes would have entered the background. Fixed before any MLL23 run;
  smudge cells and normoblasts are now also excluded, and notebook 37 does not download those three classes.

## MLL23 (2026-10-20, exploratory)
- No patient identifiers in the release: recall/FP intervals bootstrap images, as for RV-PBS, and understate uncertainty.
- Cost comparisons use the Shekhar et al. betting test with uniform sampling as the reference (as in H6d), not ALPHA.
- Specimens: 100 fields × 5 cells, θ = 0.10 only (as RV-PBS H4c); the WHO 20% threshold was not run on MLL23.

## Verifier error (2026-10-20, exploratory)
- Reader = reading round (session 0 = first reader of each cell), not an identified person; the release does not record
  who read which cell.
- A box counts for a reader if that reader's label is not in EXCL; reader counts can therefore differ from consensus counts
  in the denominator as well as the numerator.
- θ_s = consensus p_s (exact null for the truth); the plug-in β and the margins are estimated on the same 83 specimens.
- The IDC v25 consensus series has 28,459 boxes (114 patients) vs 28,395 (113) in `data/bmwsi/cells.csv` from the
  2026-09-30 download; the analysis joins to the repo table, so the extra patient is not used.

## Reference cross-check (2026-10-20)
`scripts/44_crosscheck_shangrla.py` compared our code with SHANGRLA's `NonnegMean` (shangrla/core/NonnegMean.py, sha256
1fb7cdce…; AGPL-3.0, not copied into this repository) on 2,000 random cases (seed 20261035):
- ALPHA supermartingale given the same predictable eta: max |Δ log| 4.3e-14 over 233,765 draws.
- Betting, equal weights, uniform sampling (bets reparametrised λ' = λ(N − t + 1)/N): max |Δ log| 2.3e-13 over 364,373 draws.
- Full ALPHA procedure: same stopping draw in 1,860/2,000 cases; in all 140 others ours stops LATER, every one explained
  by a rule our frozen `alpha.py` lacks: reject as soon as the verified sample total exceeds the null total (the null is
  then false with certainty). Conservative, never a validity issue, cannot occur at p = θ.
- Impact (`scripts/45_alpha_deterministic_rule.py`, seed 20261036): the rule binds in at most 7% of specimens, only at
  M = 500 within 0.02 of θ (θ = 0.10, p = 0.12: 7.0%; θ = 0.20, p = 0.22: 6.5%; M = 1,000: 0.25%), and changes no median
  cost in any cell. The frozen module is left unchanged; reported ALPHA costs stand.
- `src/cadpbs/weighted_betting_w.py` (not frozen) gained an optional `return_path` output for this check; with the
  default it behaves exactly as before (checked on 200 cases).
- The paper's statement "the frozen scripts and src/ are unchanged since their tags" was imprecise (modules were added to
  src/ later). It now says that the frozen scripts and the modules they import are unchanged (checked with `git diff`
  against prereg-v4 for alpha.py, weighted_betting.py, infodpp.py and composites.py).

## Supervisor review of 2026-10-22 (corrections; one exploratory re-analysis)
- **Proposition 1 lower bound.** The paper stated (1 − n/(N−k+1))^k ≤ h(N, k, n) without restriction. It fails for even
  k and n > N − k + 1, where h = 0 (e.g. N = 10, k = 2, n = 10: 1/81 > 0). Corrected to the positive-part form
  [max{0, 1 − n/(N−k+1)}]^k, which holds for all 0 ≤ n ≤ N. `scripts/17_part2_verify.py` (not frozen) previously checked
  only n ≤ N − k; check 2 now covers the full range and check 2b records the counterexample; `tests/test_theory.py`
  gained `test_prop1_bounds_full_range`. No other result used the bound (Theorem 1 uses h exactly; the sufficient
  sample size uses the upper bound). `results/part2_v1/verify.json` gained one key; `nstar_table.csv` is unchanged.
- **H1a reclassified** from V (implementation/theorem check) to E (empirical claim): the halting statistic has no
  relevance term but depends on the embeddings and the selected set, which depend on the cells.
- **InfoDPP-PAC seed accounting** (`scripts/47_infodpp_seed_accounting.py`, exploratory, seed 20261022 as the 2.4% run).
  The m0 = 5 random seed fields are scored by the teacher but were never counted as examined; cost and misses count only
  greedy selections. At 2.4% the greedy rule made a median of 4–5 selections (range 2–10) after the seed and never
  reselected a seed field. The script reproduces all 12,000 stored episodes of `results/sens_op_paperpool` exactly and
  adds the alternative accounting (seed counted, uniform sampling at n + m0): k = 1 misses 0.78/0.92/0.93/0.95 vs
  0.95/0.95/0.96/0.94 at m = 1/5/10/20, the same pattern as before. 2.4% is described as a ratio of reported averages.
- **Wording**: content-blindness stated as an empirical finding; Proposition 3 described as an average-case result under
  an independence model, not a proof that only relevance ranking helps and not a worst-case guarantee; certification
  cost described as governed mainly by the margin, with dependence on M diminishing only at margins ≥ 0.05; the
  Discussion's near-threshold numbers corrected (classifier-driven sampling needed 660 and 798 cells at θ = 0.10,
  p = 0.12, M = 2,000 and 4,000, not more than 1,000); savings stated as expert-verified cells only; verifier margins
  described as exploratory and direction-specific; field counts in the verifier analysis stated (first reader's).
- **Figures**: Fig. 2(b) of main.tex now plots the field-level teacher (exploratory); Fig. 4(b) foregrounds the
  bias-corrected ICC, with the raw estimate hatched as a comparison (`scripts/25_figures.py`, not frozen).

## Further checks (2026-10-22; exploratory)
- **InfoDPP-PAC settings not reported in the paper.** In addition to τ (calibrated) and ρ = 1 (fixed by us), the
  reimplementation fixes the RBF length-scale (median heuristic), the observation-noise variance (0.01) and the seed
  floor (m0 = 5). `scripts/48_infodpp_settings.py` varies each one at a time on 2,000 new specimens (seed 20261040),
  recalibrating τ every time; the frozen module `src/cadpbs/infodpp.py` is unchanged (the length-scale is passed through
  `InfoDPPConfig.length_scale`). Miss rate and difference from uniform sampling changed by at most 0.03, with no
  significant paired change. Consequence: the Part I conclusions do not depend on these four choices; they may depend
  on choices not varied (encoder, teacher).
- **Number audit of the CMPB article** (`paper/audit_cmpb_numbers.py`). It found four statements that did not match
  the stored results (Wald's approximation applied to all designs instead of uniform ALPHA; near-threshold cost said to
  grow "towards" values it already exceeded; "up to a third" for 37%; "≤ 0.014 under shift" for 0.011). All are
  corrected in the text; no analysis changed.


## Full audit of 2026-10-22 (corrections; two exploratory re-analyses; tests)
Applied to `paper/cmpb.tex`, `paper/main.tex` (and the supplement built from it), `paper/jbhi_compact.tex`, the cover
letter, title page and highlights. Nothing frozen was modified.
- **Cost reference (C1).** Earlier text reported classifier-driven savings against uniform *betting* with
  approximate-Kelly bets, which itself needs 1.03–1.83 times the cells of uniform ALPHA. The reference is now uniform
  ALPHA, the best valid design without a classifier in our runs. `scripts/49_baseline_reference.py` (exploratory, seed
  family 20261101) reruns all designs on the same specimens for KU-Optofil, MLL23 and RV-PBS: classifier-driven/ALPHA
  0.45–0.75 in distribution, 1.30/1.46 on MLL23 and 1.67/3.00 on RV-PBS (p = 0.15/0.20); within 0.02 of the
  threshold (p = 0.12) every design approached a census of the 500 cells and the ratio was 0.78 (MLL23) and 0.97 (RV-PBS),
  so "1.3–3.0 under shift" is stated for margins of 0.05–0.10. Confirmatory comparisons
  against uniform betting (H3c, H6b–d) are still reported as pre-registered. Fig. 4(a) now shows this rerun.
- **Hedged bets (new, exploratory).** `src/cadpbs/hedged_betting.py`: the average of fixed-fraction betting processes
  (κ ∈ {0.05, 0.1, 0.2, …, 0.9}, fixed in advance), a test supermartingale because averages of test supermartingales
  are. Uniform hedged betting matched ALPHA (0.98–1.17); the hedged 50% mixture needed 0.61–0.81 times the cells of
  ALPHA in distribution and at most 1.15 under shift. Not pre-registered; reported as exploratory everywhere.
- **RV-PBS probe refit.** Script 49 refits the KU-trained probe; its RV-PBS recall is 0.327, not the confirmatory 0.340
  (H4a keeps 0.340). Reported in main.tex.
- **ICC null bias (C4).** `scripts/51_icc_diagnostics.py` (seed 20261102) shows why the permutation null of the pooled
  ANOVA ICC is positive (between mean squares weighted by k_r − 1, within by n_r − k_r, and label variance falls with
  region size; Spearman −0.55) and gives a null-centred estimator; the bias-corrected values reported earlier stand. In
  the 32 regions with ≤ 30% blasts the 4 × 4 ICC is 0.059 [0.007, 0.108] (exploratory).
- **Verifier result reframed.** No new analysis: at exact nulls examined almost completely, false certification is
  decided by the sign of the reader's error (23 over-counting specimens → "above", 52 under-counting → "below", 8
  without error → neither). Remark 1 now states that β0, β1 are per-specimen bounds on the same cells.
- **Statements corrected.** Positivity of q_t is needed for power and "below" certificates, not for validity of
  "above" (f ≥ 0); Assumption 1 reworded (the class of specimens regarded as possible contains every placement);
  Remark 1 needs the same denominator and per-specimen bounds (counterexample: excluding 30 non-blasts gives
  p̃ = 0.1031 > 0.10 with β0 = 0); a corollary of Theorem 1 for halting rules that ignore verified labels; Wald's
  approximation is within 20% only at M = 4,000 (it fails at M = 1,000 for θ = 0.20, p = 0.25: 303 vs 406); KT2
  (the proposal's rarity hypothesis) listed as not supported on exploratory data and not counted among the 20
  pre-registered hypotheses; episode-level McNemar tests described as Monte Carlo precision only (11 blast patients),
  with fold-level t intervals for miss-rate differences; 200-cell counts added to the decision table (from
  `results/decision_level_v1`, unchanged); survey-sampling framing (Kish design effect, Cochran mean-of-ratios,
  Hansen–Hurwitz/Horvitz–Thompson PPS) added with references.
- **Tests (C5).** `tests/test_validity.py` gained 24 checks: the control-variate branch, the 25% and 50% uniform
  mixtures, hedged bets with uniform, oracle, anti-oracle, mixture and adaptive samplers, a sampler with zero mass on
  half of the unsampled fields ("above" stays valid), an independent recomputation of the hedged wealth, and a power
  check. 50 checks pass; the largest false-certification rate among valid procedures is 0.034 (500 repetitions).
- **Output convention.** Script 49 writes run settings to `config_<source>.json` and numbers only to
  `summary_<source>.json`, so `python reproduce.py verify baseline_reference` compares numbers; the stored summaries were
  split accordingly (no number changed; a full MLL23 rerun matches exactly).
- **Compact (JBHI) version.** To stay within 8 pages after the corrections, the assumptions table and the M-scaling
  table were replaced by text, the proof of Theorem 1 moved from the appendix into the text, and the malaria-protocol
  aside omitted; all three remain in main.tex and the CMPB article.

## JBHI manuscript, 10 pages (2026-10-28)
- The target journal is now IEEE JBHI. Two extra pages were agreed with the supervisor, so `paper/jbhi_compact.tex`
  (10 pages) restores what the 8-page version left out: the M-scaling table (Table I) and the assumptions table
  (Table V); the fixed-order and confidence-based baselines with their confirmatory numbers; recall at the stop point;
  the bound for halting rules that ignore verified labels, stated with the observed mean stop (1 − 32.4/200 ≈ 0.84);
  the field-teacher swap; raster, ranked scan and k-center stopping on real tiles; savings by specimen size, the effect
  of ALPHA's prior guess and the cost at α = 0.01; error patterns by cell class and image quality; the caveats of the
  real arrangements; the unequal-field bias in real regions; and the sensitivity intervals. The abstract follows the
  JBHI format (at most 250 words, no abbreviations, a closing sentence on significance); spelling is American.
- `paper/audit_jbhi_numbers.py` checks 169 statements and table cells of the JBHI text against `results/`; all pass
  (`paper/audit_jbhi_numbers.md`).
- Errors found by this audit, corrected in `main.tex` and the supplement as well:
  - the single-host background intervals were rounded wrongly: the stored fold-level intervals [−0.381, −0.165] and
    [−0.195, −0.099] give [−0.38, −0.16] and [−0.19, −0.10], not [−0.38, −0.17] and [−0.20, −0.10];
  - the sentence on confidence-based stopping in `main.tex` ("median 8–41 fields", "0.78–0.88 at m ≥ 10") matched no
    stored result. It now reports the confirmatory run (k = 1: median 10–33 fields, missed 0.74–0.90 against
    0.77–0.92 for uniform sampling at the same cost) and the fixed order (0.83–0.85);
  - "the saving vanished for M ≥ 2,000" (p = 0.22 against θ = 0.20) now reads "at most 3%" (2.8% at M = 2,000,
    −0.5% at M = 4,000).
- Wording in the JBHI text: the weighted observation has conditional mean at most the unsampled mass (equal when q_t
  is positive on every unsampled field); errors by class and by image quality now give the stored values and include
  brightness and sharpness.
- `refs.bib`: pages of Stark (2023) added (641–679); the InfoDPP-PAC entry no longer carries an ambiguous version
  date; Conformal Risk Control (Angelopoulos et al., ICLR 2024) added for the remark on α = 0.01. The CMPB PDF was
  rebuilt for the reference list only.
- Both audit scripts now count the tests with `--noconftest`, because the conftest hook rewrote
  `results/validity_tests_v1/report.json` during collection. The full suite was rerun: 50 passed.

