# InfoDPP-PAC: status check and critical review (2026-10-19)

Paper: Mittal, Srivastava & Chauhan, *InfoDPP-PAC: Principled Patch Selection for Whole Slide Image Analysis*,
arXiv:2608.23574 (local copy: `papers/InfoDPP-PAC.pdf`). The authors' claims are separated from my assessment below.

## 1. Publication status
| Source | What it says |
|---|---|
| arXiv record | v1 only, submitted 7 July 2026; no journal reference |
| Preprints.org 202607.0611 | v1 (7–8 July 2026): "This version is not peer-reviewed"; no published version linked |
| The PDF | No venue, no "under review" note, no code link. The layout and a row labelled "Gastrointestinal (TMLR test)" in Table 4 suggest it was prepared for TMLR; whether it is under review there is **not confirmed** |

**Conclusion.** Calling it a preprint is accurate. In the paper, however, the plain wording "a recent preprint" is the
standard form; "has not been peer reviewed" reads as a judgement and was removed. **Re-check the status before
submitting.** If it is accepted somewhere, update the citation.

## 2. What the authors claim
- Theorem 1: log det(I + βK_S) is the GP mutual information.
- Theorems 2–4: the objective is monotone submodular, so greedy selection carries (1 − 1/e) guarantees.
- Theorem 6: a PAC-style certificate on residual information gain, under Assumption 1 ("seed-stable" estimates).
- Results: 83.7% fewer patches than a 300-patch budget, keeping 97.9% of full-budget composite selection quality, on 22
  test slides.
- Their own scope statements:
  - the certificate does not bound "the generalisation error of a classifier";
  - it certifies the surrogate's residual information, "not downstream diagnostic performance";
  - the experiments do not establish end-to-end diagnostic improvement.

## 3. Review: findings that matter for our paper
1. **What Theorem 6 concentrates is not stated precisely.**
   - The halting statistic Γ̂ (Eqs. 5–6) is a function of the GP *posterior covariance* (Eq. 2). That covariance depends
     on the seed *locations* and the selected set S, not on the teacher's pseudo-labels.
   - Assumption 1 bounds the effect of "replacing any one seed observation". If the observation is the label, that effect
     is zero and the ε term has nothing to control. If it is the seed location, then neither the "target marginal" that Γ̂
     is said to estimate without bias nor the bounded-difference constant ρ/m₀ is derived; both are assumed.
2. **Adaptivity is not handled.**
   - S_{k−1} is chosen using the same seed set (through μ, σ and q).
   - McDiarmid's inequality and a union bound over n candidates and all k treat each Δ̂_info,k(i) as a fixed function
     of independent seeds. The data-dependent selection path is not covered.
   - So the validity of the certificate is **not established**. This is a gap in the argument, not a demonstration that
     the certificate fails.
3. **What is certified.** Whatever the status of Theorem 6, its object is residual GP information in a teacher-seeded
   surrogate, a property of the embedding geometry. The authors say so themselves. It carries no statement about where
   decision-relevant evidence is, which is the gap Part I measures.
4. **Operating point. Our draft overstated how much InfoDPP-PAC examines in its own experiments.**
   - "83.7% fewer patches" is relative to a 300-patch budget. Slides averaged 2,075 patches after filtering.
   - So the paper's mean k⋆ = 48.9 is about 2.4% of an average slide's patches.
   - Our main runs calibrated τ to 16.3% of all fields, which gives InfoDPP-PAC about 7 times more examination than in its
     own experiments.
   - Re-run at about 2.4% of the fields (median 5 of 200; seed 20261022, 12,000 episodes): InfoDPP-PAC missed the single
     critical field in 0.81, 0.96, 0.96 and 0.97 of specimens at m = 1, 5, 10, 20 (uniform sampling: 0.98, 0.98, 0.98,
     0.97). Our main operating point is therefore favourable to the method.
5. **Teacher.** The teacher is a cosine similarity of the patch embedding to category prototypes, i.e. a function of the
   unit's own embedding. The field-level teacher added on 2026-10-19 matches this better than the pre-registered
   cell-level maximum, so it is the better default for describing the method.
6. **Circular quality metric.**
   - "Relevance" in their evaluation is the mean of q_i, the same GP-derived score that Eq. 3 maximises.
   - "Highest mean teacher-derived relevance among fourteen baselines" is therefore partly true by construction.
   - This does not affect our results, but it limits what their evaluation shows.
7. **Unreported settings.**
   - τ, δ, ρ, the kernel, the length-scale, σ_n and m₀ (beyond "≈ 0.01n") are not given.
   - A "corrected stopping rule" is mentioned in Section 6.2 without stating what was corrected.
   - This confirms why our reimplementation can only reproduce cost, not the configuration.
8. **Notation.**
   - K_S is defined as the posterior kernel (Eq. 3), but κ_ii as the prior kernel (Section 3).
   - Our code uses the posterior kernel throughout, the only consistent reading. Logged as an interpretation.

## 4. What changes in our paper
- Wording: "a recent preprint" (neutral). No code was released, so we reimplemented it from its equations.
- Methods: the operating point is described correctly (16.3% of a 300-patch budget, about 2.4% of a slide's patches). We
  state that our 16.3%-of-specimen choice is favourable to InfoDPP-PAC, and report the 2.4% run.
- One restrained sentence notes that Theorem 6 rests on Assumption 1, whose bounded-difference constant and target
  quantity are not derived. Our results do not depend on whether the certificate holds, and we never claim it.
- Dependence on one preprint method is a vulnerability. It is addressed by P1-1, which tests saturation-type stopping as
  a class: pure information-gain and k-center coverage stops at the same cost (`scripts/31_stopping_class.py`).
  **Outcome:** at k = 1 the three rules missed 0.31–0.34, 0.67–0.71, 0.74–0.77 and 0.83–0.85 at m = 1, 5, 10, 20, within
  0.04 of each other in every cell; uniform at the same cost 0.84, 0.85, 0.83, 0.83. The Part I conclusion therefore does
  not rest on InfoDPP-PAC's particular construction or on whether its certificate holds.
