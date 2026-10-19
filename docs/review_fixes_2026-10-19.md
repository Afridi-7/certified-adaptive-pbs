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
