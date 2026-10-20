# P0 fixes from the 2026-10-19 review: what was done, and what it showed

All new analyses are **exploratory**: they were added after the pre-registered analyses, on new seed families, and none
changes a pre-registered verdict. Paper edits are in `paper/main.tex`. The previous draft is kept as
`paper/main_2026-10-03_backup.tex`.

## P0-1 Pre-registration record (RQ1–RQ4)
- **Finding.** The record is sound. For every tag, the first confirmatory output was written after the tag was pushed to
  GitHub, by 43 s to 18 min. Every frozen script is unchanged since its tag. Only the dates typed in `PREREG.md` are wrong.
  Evidence: `docs/PREREG_AUDIT.md` and paper Appendix C.
- **Paper.**
  - "Five git tags", not four.
  - A new paragraph says what "confirmatory" means: fresh draws over the same cell pool. Only H4a–c and H7a were tested
    on data not analysed before.
  - Table 2 gains a Type column, E (empirical claim) or V (implementation check). V marks H1a, H1d, H3a, H4c, H5a and H6a.
  - Corrected total of confirmatory specimens: 67,500.
- **Still to do (you).**
  - Archive the repository and the result folders on Zenodo at submission.
  - Make the GitHub repository public, or give reviewers access.
  - Do not re-create the lightweight tags `prereg-v3` and `prereg-v4`.

## P0-2 Kill-test 1 teacher (RQ1)
- **Bug confirmed.** The KT1 relevance AUC scores `bank.p[cells].max(1)`, the largest cell-level score in a field. That
  score is a cell-detector oracle, and it was also InfoDPP-PAC's teacher.
- **Fix.** `scripts/26_kt1_field_teacher.py` adds a field-level teacher. It is a logistic regression on the pooled field
  embeddings, trained on the other four folds only. Run: 10,000 episodes, seed 20261019, recording the blast patient of
  each critical field.
- **Results.**

  | m | information-gain AUC | cell-level teacher AUC | field-level teacher AUC [fold t 95%] | InfoDPP-PAC miss, field teacher (k = 1) | uniform miss | ranked by field teacher |
  |---|---|---|---|---|---|---|
  | 1 | 0.775 | 0.978 | 0.978 [0.961, 0.995] | 0.316 | 0.866 | 0.030 |
  | 5 | 0.589 | 0.938 | 0.925 [0.883, 0.966] | 0.662 | 0.854 | 0.152 |
  | 10 | 0.546 | 0.920 | 0.882 [0.830, 0.933] | 0.764 | 0.830 | 0.242 |
  | 20 | 0.529 | 0.898 | 0.827 [0.779, 0.874] | 0.770 | 0.802 | 0.356 |

- **Reading.** Kill-test 1 holds with a fair teacher. The original numbers overstated relevance at m = 20 by about 0.07.
  InfoDPP-PAC's misses barely move when its teacher changes. A relevance-ranked scan finds far more, but carries no
  guarantee.

## P0-3 Unequal field sizes (RQ3, RQ4)
- **New method.** Proposition 2 in the paper has two parts:
  - (a) Counts known: the weighted Shekhar et al. test with w_i = m_i / M (`src/cadpbs/weighted_betting_w.py`). With equal
    weights it reproduces the frozen code exactly in 400 of 400 test cases. In a check at the exact null with an
    adversarial sampler, the rejection rate was 0.002.
  - (b) Only a bound known: the bounded assorter a_i = (B_i − θ m_i + θ m_max) / m_max.
- **Real data** (BoneMarrowWSI, 758 quadrants, 99 patients):
  - Quadrant size: median 20 cells (IQR 12–29, maximum 91, coefficient of variation 0.62).
  - Mean of quadrant proportions minus the region proportion: −0.004 [−0.008, 0.001]. It exceeded 0.03 in absolute value in
    10% of regions.
  - Within-region Spearman correlation of size with proportion: 0.04.
- **Simulation** (`scripts/27_unequal_fields.py`; 60 fields with real sizes, about 1,300 cells):
  - **Validity.** Averaging field proportions certified falsely in 0.38 (θ = 0.10) and 0.32 (θ = 0.20) of boundary
    specimens with size-independent blasts, and in 1.00 with small-field enrichment. Constructions (a) and (b) stayed at or
    below 0.005 everywhere.
  - **Cost** at θ = 0.10, p = 0.15:
    - counts known, classifier-driven: 297 cells;
    - counts known, uniform: 646;
    - assorter with betting: 631;
    - assorter with ALPHA: 1,031.

## P0-4 Population size (RQ3)
`scripts/28_m_scaling.py`: M ∈ {500, 1,000, 2,000, 4,000}; 17,600 specimens; seed 20261021; α = 0.05 and 0.01.
- **Validity.** Every design at every M: false certification ≤ 0.040 at α = 0.05 and ≤ 0.007 at α = 0.01.
- **Cost depends on the margin, not on M.**
  - At margins ≥ 0.05, M = 4,000 is within 20% of Wald's log(1/α)/KL(p‖θ).
  - Within 0.02 of the threshold, cost grows towards about 1,400 cells (θ = 0.10) and 2,500 cells (θ = 0.20).
  - "Near-census" is a small-specimen effect. Paper Table 1.
- **Classifier-driven sampling** saved 36–65% at every M, except θ = 0.20, p = 0.22, where the saving was 19–23% at M ≤ 1,000
  and zero at M ≥ 2,000.
- **α = 0.01** costs 1.0–1.7 times as much as α = 0.05.

## P0-5 Statistical reporting
Computed by `scripts/29_review_stats.py`; output in `results/review_stats_v1.json`.
- **Fold-level t intervals.** H1b is positive in 5 of 5 folds for every k. For H1c, m = 10 gives [−0.153, −0.021], still
  outside the margin. m = 20 gives [−0.050, 0.003], at the edge of the margin, so equivalence there is borderline.
- **Effect sizes, Part III.**
  - Sampler only (classifier-driven vs uniform betting): 0.28–0.73 times the cells.
  - Test only (uniform betting vs uniform ALPHA): 1.03–1.83 times the cells.
  - The draft's "25–55% saving" mixed the two.
- **Episode counts.** E3.1 has 18,000 episodes, not 15,000. E3.2 has 4,500, not 5,000.
- **Pre-registered sensitivity analyses**, run 2026-10-19 with the frozen script (InfoDPP-PAC minus uniform, k = 1):
  - max pooling: −0.115 at m = 10, −0.056 at m = 20;
  - single-cell: −0.036 at m = 10;
  - τ sweep: never more than 0.13 in absolute value, but −0.08 and −0.10 at m = 20 for the 30% and 50% operating points;
    the 16.3% replicate gives −0.050 at m = 20;
  - **single-host background: −0.27 at m = 10 and −0.15 at m = 20.** The advantage is larger with a one-patient
    background, but that construction duplicates cells heavily. This weakens "the advantage disappears by m = 20". It does
    not weaken "the miss rate stays near the worst case" (InfoDPP-PAC still misses 57–71%).
- **Critical-evidence recall at the stop.** InfoDPP-PAC 0.34–0.69 at m = 1 and 0.19–0.21 at m = 20; uniform 0.15–0.18
  throughout.

## P0-6 Citations and records
- **Bibliography.** InfoDPP-PAC is now marked as a preprint (arXiv v1 of 7 Jul 2026, and Preprints.org), not peer reviewed.
  Added: pages for Shekhar et al. and for Csillag et al.; Hanley & Lippman-Hand (1983); the BoneMarrowWSI descriptor
  (arXiv:2509.15895); Zrnic & Candès (ICML 2024); and Wald (1945).
- **RV-PBS.** The authors were already corrected in the repo (Pal et al.). The project copy of `data_audit_rv_pbs.md` is out
  of date. **RV-PBS states no licence**: the data statement now says so. Consider asking the authors for permission.
- **KU-Optofil.** Added the acquisition details (Sysmex DI-60, May–Grünwald–Giemsa stain, κ > 0.85).
- **Rümke.** A secondary source dates the chapter 1977 (Aspen conference, College of American Pathologists), so the bib entry
  is left as 1977. Its key `rumke1978` is cosmetic.

---

# P1 items (2026-10-19, second pass). All exploratory, with new seed families.

## InfoDPP-PAC status and review (prompted by the question "why unreviewed?")
See `docs/infodpp_pac_review.md`. Short version:
- "Preprint, not peer reviewed" is factually correct (arXiv v1 7 Jul 2026, no journal ref; Preprints.org says "not
  peer-reviewed"). The paper now says neutrally "a recent preprint". **Re-check before submission.**
- The review found that our draft **misdescribed the operating point**: 16.3% is of a 300-patch budget, about 2.4% of a
  slide. Our main runs therefore gave InfoDPP-PAC about 7 times its own examination share, which favours it.
- At 2.4% (seed 20261022, `results/sens_op_paperpool`): InfoDPP-PAC miss 0.81/0.96/0.96/0.97 at m = 1/5/10/20 (k = 1),
  uniform 0.98/0.98/0.98/0.97.
- Theorem 6 rests on Assumption 1, whose target quantity and bounded-difference constant are assumed, and the McDiarmid +
  union-bound argument does not cover the adaptive selection path. One restrained sentence in Related work; our results
  do not depend on the certificate.

## P1-1 Saturation as a class of rules (RQ1) — done
`scripts/31_stopping_class.py`, seed 20261024, 4,000 specimens (m ∈ {1,5,10,20}, k ∈ {1,5}, 500 each), field-level teacher,
each rule calibrated to stop at 16.3% of fields on blast-free specimens.

| k = 1, miss rate | m = 1 | m = 5 | m = 10 | m = 20 |
|---|---|---|---|---|
| InfoDPP-PAC | 0.314 | 0.688 | 0.758 | 0.834 |
| pure information gain (λ2 = 1e-9) | 0.336 | 0.670 | 0.742 | 0.830 |
| k-center coverage stop | 0.338 | 0.708 | 0.766 | 0.850 |
| uniform, same cost | 0.842 | 0.846 | 0.832 | 0.834 |
| ranked by field-level teacher | 0.052 | 0.144 | 0.214 | 0.374 |
| Proposition 3 prediction | 0.051 | 0.146 | 0.221 | 0.394 |

- The three saturation rules are within 0.04 of each other everywhere: the conclusion is about saturation stopping, not
  about one preprint.
- **New Proposition 3 (ranked search):** P(miss) = E_{S_c~F1}[P(Bin(N−1, 1−F0(S_c)) ≥ n)] → 1 − TPR at FPR n/(N−1); it
  reduces to 1 − n/N when F1 = F0 (Theorem 1). Matches observed ranked misses within 0.02.

## P1-2 Part I on real tiles — done 2026-10-20 (exploratory)
Notebook 35 was run in Colab by the author on 2026-10-20 (T4; torch 2.11.0, timm 1.0.30; 190 ROIs, 3,040 tiles). The first
version failed with `missing .../rois.csv`; the fixed version rebuilds the input tables from `bmwsi_ann.zip` and checks
them against fingerprints of the repo tables. `check_roi0.jpg` shows a marrow monolayer, and the consensus boxes overlay
the nuclei, so the annotation coordinates are level-0 pixels. Outputs copied to `data/bmwsi_tiles/` (sha256 of the
embeddings file starts f1b9c7bd14087193).
`scripts/36_partI_real_tiles.py` (seed 20261028; re-run byte-identical) + `scripts/40_partI_real_tiles_analysis.py`
(patient bootstrap, seed 20261031) → `results/partI_real_tiles_v1/`.
- 99 patients, N = 32 tiles (91) or 16 (8); median 5 nucleated cells per tile; 42% of tiles critical; teacher AUC 0.966.
- InfoDPP-PAC stops after 5.4 tiles on average (18% of N). 72 specimens have 0 < k < N.
- **Low-blast specimens (k/N ≤ 0.25, n = 24; 14 with k = 1): InfoDPP-PAC miss 0.80 vs uniform 0.70 (h = 0.70);
  difference +0.105 [0.002, 0.203]; worse in 18/24.** Raster 0.67, ranked 0.54.
- InfoDPP-PAC missed in ≥ 95% of runs in 15/24 low-blast specimens; uniform ranged 0.25–0.90 (≈ h everywhere).
- Teacher AUC within the low-blast specimens 0.69 (single-blast tiles vs negatives 0.67).
- k-center (own cost 8.5 tiles), k = 1: 0.89 vs uniform at the same cost 0.70, +0.19 [0.03, 0.34].
- k/N > 0.5 (42 specimens): every rule ≥ 98.7% detection.
- Paper: new Part I paragraph "Real image tiles", abstract sentence, contribution, data and limitation updates, deviation
  entry (m0 floor 2; k-center not cost-matched).

## P1-3 Observed validity on real arrangements (RQ3, RQ4) — done
`scripts/30_real_arrangement_validity.py`, seed 20261023. 75 patients as specimens, 4×4 subregions as fields (median 5
cells, max 33), θ_s = p_s (exact null), 400 sampling orders per specimen.
- Field-level known counts: 0.001 pooled, max 0.01 in any specimen. Bounded assorter: 0. Cell SRS: 0.027, max 0.055.
- Naive cell-level analysis of field samples: **pooled 0.051, but > 0.05 in 27/75 specimens (21 significant), max 0.335.**
- Mean of field proportions: 0.062 pooled, two specimens at 1.0.
- Caveat: specimens are small (median 105 cells) and mostly high-blast (46/75 have p_s > 0.5; only 2 in 0.05–0.30).

## P1-4 Verifier error from reader labels — done 2026-10-20 (exploratory)
Data: per-session label series of BoneMarrowWSI (IDC v25): Session 0 and 1 for 114 patients; later sessions only where
readers disagreed. The conversion documentation (Zenodo 14933087) says each cell was read in sequence until a majority was
reached and that which person read which cell is not recorded, so sessions are reading rounds, not fixed annotators.
Notebook 42 run in Colab by the author on 2026-10-20: 92,222 labelled boxes (sha256 bcb248bc…). Session boxes match the
consensus boxes (99.98% / 99.95%). Copied to `data/bmwsi_sessions/`.
`scripts/43_verifier_error.py` (seed 20261034) → `results/verifier_error_v1/`.
- Reader 0 vs reader 1 (both counted): second calls blast where first called non-blast 1.1% [0.7, 1.7]; non-blast where
  first called blast 1.2% [0.8, 1.6]. Counted vs not counted (damaged, smudge, …) disagreement: 2.7% of boxes.
- 83 real-arrangement specimens, reader 0 as verifier, truth = consensus, θ_s = p_s:
  - |p_reader − p_s|: median 0.005, 90th percentile 0.033, max 0.075; mean bias −0.007 (reader 0 under-calls).
  - False certification "p > θ_s" 0.28 (23 specimens > α); "p < θ_s" 0.63 (52 specimens).
  - Remark-2 correction with plug-in β (inter-reader rates): 0.16 / 0.42 — not enough.
  - Fixed margin: "above" valid from +0.03 (0.0003 pooled, max 0.013); "below" still 3 specimens > α at −0.05.
- Paper: Remark 2 (imperfect verification of a threshold assertion), RQ4 paragraph, practical implication, limitation,
  contribution 4, abstract clause. Margins chosen in-sample on small, mostly high-blast specimens.

## P1-5 Three-way decisions at fixed budgets (RQ3) — done
`scripts/32_decision_level.py`, seed 20261025, M = 2,000, n ∈ {100, 200, 500}, 2,000 specimens per (θ, p). Paper Table 2.
- Within 0.02 of θ at n = 100: conventional count wrong 18–37%.
- Valid procedures: wrong ≤ 0.7% off the threshold, ≤ 4.7% on it; "undetermined" 61–98% within 0.02 at n ≤ 500.
- Fixed-n exact test is decisive more often than the sequential test (the price of anytime validity).
- ALPHA's default prior guess is poor for cell-level draws (31% vs 86% correct at p = 0.15, θ = 0.10, n = 500).

## P1-6 α = 0.01 — done (text)
The proposal's "α = 0.01 unattainable" concerns a CRC-type calibration correction (budget α − (1−α)/n < 0 for n < 99).
The design-based tests use no calibration set, so it does not apply; α = 0.01 costs 1.0–1.7× α = 0.05 (P0-4). The
proposal's claim is unchanged; the paper now says which procedure it applies to.

## P1-7 Confusable backgrounds — done
`scripts/33_confusable_background.py`, seed 20261026: 0/10/20% of non-blasts from metamyelocytes, myelocytes, reactive
lymphocytes. All designs valid (≤ 0.038). Classifier-driven cost at margin 0.05: 75 → 85 → 105 cells (θ = 0.10) and
140 → 150 → 175 (θ = 0.20); still cheaper than uniform ALPHA (160–240).

## P1-8 Second shift source (MLL23) — done 2026-10-20 (exploratory)
Notebook 37 run in Colab by the author (T4; 15 of 18 class ZIPs, md5-checked; atypical promyelocytes, smudge cells and
normoblasts skipped). 36,529 images, 8,606 myeloblasts; licence CC BY 4.0 (Zenodo 10.5281/zenodo.14277609); descriptor
Shetab Boushehri et al., Sci Data 12:1773 (2025), doi:10.1038/s41597-025-06223-x (verified on PMC). **No patient
identifiers** (file names are class_index), so intervals are image bootstraps and too narrow.
`scripts/38_mll23_shift.py` (seed 20261029; re-run byte-identical after adding effect sizes) → `results/mll23_shift_v1/`.
- KU-trained probe, fixed threshold: AUC 0.915; recall 0.63 [0.62, 0.64] (KU 0.91, RV-PBS 0.34); FP 0.060 [0.057, 0.062].
  FP by class: neoplastic lymphocytes 0.47, promyelocytes 0.22, reactive lymphocytes 0.18, LGL 0.15, lymphocytes 0.11.
- Validity at p = θ = 0.10 (1,000 specimens): uniform 0.011, classifier-driven 0.001, 50% mixture 0.004.
- Cost (medians, 400 specimens): p = 0.15 uniform 282.5, prop 215, mix 220; p = 0.20 uniform 105, prop 95, mix 85.
  Ratios prop/uniform 0.76 [0.70, 0.81] and 0.90 [0.82, 1.05]; mix/uniform 0.78 and 0.81 [0.73, 0.93].
- Reading: the classifier-driven saving shrinks with recall across three sources (in distribution 0.28–0.73×; MLL23
  0.76–0.90×; RV-PBS 1.77×). Three cases, not a dose–response estimate.
- Paper: data bullet, RQ4 paragraph, abstract (iv), contribution 4, limitation, data availability; bib entry added.
- `config.json` of the run (with Zenodo checksums) is in Drive `certified-adaptive-pbs/mll23_v1/`; copy it to
  `data/mll23/` for the record.

## P1-9 Bias-corrected clustering — done
`scripts/34_icc_excess.py`, seed 20261027. Excess ICC over the permutation-null mean: quadrant 0.013 [0.002, 0.028];
4×4 0.043 [0.018, 0.067]. Clustering is real but weak. (Not registered as an amendment: it was computed after the data were
seen, so it stays exploratory.)

## Extra: ALPHA prior guess check
`scripts/39_alpha_prior.py`, seed 20261030, 4,400 specimens. For field-level draws (Part III setting) the default prior
guess is the *best* of three fixed in advance (165 vs 260/238 cells at p = 0.15, θ = 0.10), so the reported ALPHA costs
are not inflated by it. For cell-level draws the reverse holds (P1-5).

## Paper changes in this pass
Abstract shortened (495 → 331 words; a journal version must be ≤ 250 for JBHI-style limits — check the venue's rule);
Part I class-of-rules paragraph + Proposition 3; 2.4% operating-point bullet; Part III α = 0.01 reconciliation, ALPHA
prior, decision table (Table 2), confusable backgrounds; RQ4 excess ICC and observed real-arrangement validity;
discussion and practical implications; limitations; four new deviation entries. 22 pages, compiles cleanly.

## Pre-submission checks (2026-10-20)
- **InfoDPP-PAC status:** Preprints.org 202607.0611 is still version 1 and "not peer-reviewed"; no published version
  linked (checked 2026-10-20). Wording "a recent preprint" stays. Check again at submission.
- **Abstract:** IEEE JBHI allows at most 250 words, without abbreviations, references or equations
  (embs.org/jbhi/prepare-and-submit-your-manuscript). Rewritten to 250 words, no equations (the previous 331-word version
  is kept in the git history).
- **Length:** JBHI regular papers are limited to 14 pages in the IEEE two-column template, including supplementary
  material; pages 9–10 cost $250 each and 11+ $350 each. The current 23-page single-column preprint must be converted and
  probably shortened (appendices and some exploratory analyses to a supplement that still counts, or cut).
- **Repository URL:** github.com/Afridi-7/certified-adaptive-pbs; its visibility could not be checked from here. Make it
  public (or provide reviewer access) and archive a release on Zenodo for a DOI before submission.

## Independent checks of the guarantees (2026-10-20; follow-up to an external review)
1. **Reference implementation.** Our ALPHA and betting statistics agree with SHANGRLA's NonnegMean to within 3e-13 in
   log scale (2,000 random cases). The full procedures disagree only through a conservative rule missing from our frozen
   ALPHA, which changed no median cost (details in DEVIATIONS.md). `results/crosscheck_shangrla_v1`, `results/alpha_detrule_v1`.
2. **Validity test suite** (`tests/`, run `python -m pytest -q tests`; 25 tests, about 30 s, fixed seeds). Exact-null
   false certification for ALPHA (default, aggressive, timid priors) and for weighted betting with uniform, random-score,
   oracle, anti-oracle and adaptive follow-the-blasts samplers, equal and unequal fields with size-dependent blasts; the
   bounded assorter; the verifier-margin remark; identities for h(N, k, n) and the ranked-search formula; four negative
   controls (uncorrected ranking, mean of proportions, cells as draws under clustering, naive threshold under verifier
   error) that must be flagged; a power check. All pass; largest valid-procedure rate 0.029 (1,000 reps). Report:
   `results/validity_tests_v1/report.json`.
3. Still needed: a reading of the proofs and of these checks by a statistician who has not worked on the project.
