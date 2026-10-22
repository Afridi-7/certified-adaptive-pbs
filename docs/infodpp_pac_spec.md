# InfoDPP-PAC: transcribed specification and reimplementation notes

Source: Mittal, Srivastava & Chauhan, *InfoDPP-PAC: Principled Patch Selection for Whole Slide Image Analysis*, arXiv:2608.23574v1 (7 Jul 2026). Local copy: `papers/InfoDPP-PAC.pdf`. No code has been released.

## Method (paper equations)
- **Pool:** n patches with embeddings x_i ∈ R^d. **Seed:** D0 is drawn uniformly, |D0| = m0 ≈ 0.01n. Teacher pseudo-labels are y_i = T(x_i) + ε_i.
- **GP posterior (Eqs. 1–2):** μ_D, k_D. Quality is q_i = α1·μ̄_i + α2·σ̄_i, with μ̄ and σ̄ slide-wise affine-normalised to [0, 1].
- **Objective (Eq. 3):** F(S) = λ1·log det(I + βK_S) + λ2·Σ_{i∈S} q_i, where β = σn⁻² and K_S is the **GP posterior kernel** restricted to S.
- **Selection (Eq. 4):** i*_k = argmax of the blended gain Δ_k(i) = λ1·Δ_info,k(i) + λ2·q_i.
- **Information marginal (Eq. 5):** Δ_info,k(i) = log(1 + β·σ²_{i,β}(S_{k−1})), where σ²_{i,β}(S) = κ_ii − β·k_{i,S}ᵀ(I + βK_S)⁻¹k_{i,S}.
- **Stopping statistic (Eq. 6):** Γ̂_k = max over unselected i of Δ_info,k(i). This uses information gain only, never quality.
- **Halting (Eq. 7):** stop before adding a patch when Γ̂_k ≤ τ + ε_{k,n}, where ε_{k,n} = ρ·sqrt(log(π²·n·k²/(3δ)) / (2·m0)).
- **Certificate (Thm 6, McDiarmid + union bound under Assumption 1):** with probability ≥ 1 − δ over the seed set, every unselected patch has *target* information marginal ≤ τ + 2ε.
- **Residual vector (Eq. 8):** r = the q-weighted centroid of the unselected patches.
- **Tuned weights (Sec. 6.2):** λ1 = 1.75, λ2 = 0.25, α1 = 0.5, α2 = 1.5.
- **Setting:** Phikon 768-d embeddings; 256×256 tiles; 300-patch cap in the adaptive runs; k* ranged 6–179 (mean 48.9).

## Claims in our proposal, checked against the paper
| Proposal statement | Status |
|---|---|
| The certificate "does not bound the generalisation error of a classifier" | **Verified** (Sec. 2.4). Sec. 5.4 also says the certificate is "not downstream diagnostic performance". |
| The halting statistic uses only the IG term, not relevance | **Verified** (Eq. 6 and Sec. 7 "The role of PAC stopping"). The authors justify this choice explicitly. |
| Spatial coverage 0.169 vs 0.266 for uniform | **Verified with a qualification:** Table 2 compares at **fixed k = 50**, with the tuned InfoDPP-PAC against the uniform grid and random tissue (both 0.266). It is not a comparison under adaptive stopping. The proposal should say "at matched k = 50". |
| Evaluated on WSIs, not blood smears | Verified: 202 HISTAI GI slides, selection-quality metrics only, no diagnostic endpoint |

## Not specified in the paper (our choices, to be logged in DEVIATIONS.md)
The paper gives no values for **τ, δ, ρ, the kernel type and length-scale, σn (β), m0 beyond "≈0.01n", or the teacher for our domain**.
- ρ is a constant in Assumption 1 that cannot be estimated from the paper.
- **Consequence:** the stopping point of any reimplementation is set by our hyperparameter choices.
- **Mitigation:** calibrate τ on negative (blast-free) specimens from the training folds to reproduce the paper's operating point (~16% of the budget, i.e. 83.7% fewer), then sweep τ as a sensitivity analysis.

## Analytical observations (for Kill-test 1, RQ1)
1. **Given the selected set, the stopping statistic contains no relevance term.** σ²_{i,β}(S) depends on the embeddings, the kernel, the seed *locations* (through the posterior kernel) and S. It does not depend on the teacher labels y or on q. Labels affect stopping only **indirectly**, through which S the blended gain selects; with λ2 = 0.25 that indirect channel is weak. So a critical cell can delay stopping **only by being far from the selected set in embedding space**. Whether blasts are embedding-space outliers is therefore the empirical crux. The synthetic pilot showed that the answer (the separability s) determines the sign of the gap.
2. **Assumption 1 is strong.** The target marginal Δ_info,k(i) is not defined constructively, and Γ̂ is a max of nonlinear functions of the seed. The "unbiased estimate" requirement is stated but not verified in the paper. Our study should treat the PAC certificate as *conditional on Assumption 1*.
3. **The certificate bounds residual information about the teacher-relevance GP, not the probability that a critical cell remains unexamined.** This is the precise gap that Part II addresses. With k = 1 critical field and no informative features, no stopping rule examining n fields can push P(miss) below 1 − n/N.
