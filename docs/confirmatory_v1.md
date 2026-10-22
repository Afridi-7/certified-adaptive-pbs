# Confirmatory results — Part I / Part II (2026-09-26)

Protocol: `docs/PREREG.md` v1 + Amendment 1, frozen as git tags `prereg-v1` and `prereg-v1.1` **before** these data
were generated. Command run exactly as pre-registered (seed 20260401):
`python scripts/05_killtests.py --out results/confirm_v1 --seed 20260401 --episodes 200 --cell_episodes 10:1:600 20:1:600`
Decisions computed by the frozen `scripts/09_confirmatory_analysis.py`, unmodified. Output: `results/confirm_v1/decisions.json`.
**24,000 episodes**, 5 folds (1,000 per cell; 3,000 for the two H1c cells).

## Pre-registered verdicts
| Hypothesis | Verdict | Evidence |
|---|---|---|
| **H1a** content-blind stop point | **SUPPORTED** | median stop = 33 of 200 fields in every cell; range across k ≤ 1 field in each m (rule: ≤ 2); Kruskal–Wallis p = 0.000 (m = 1), 0.29, 0.70, 0.50 |
| **H1b** advantage shrinks with dilution | **SUPPORTED** (4 of 4 tested k) | Δ(20, k) − Δ(1, k), fold-clustered 95% CI: k = 1 [0.392, 0.568]; k = 2 [0.514, 0.590]; k = 5 [0.261, 0.327]; k = 10 [0.066, 0.166]. k = 20 descriptive (floor): [−0.005, 0.019] |
| **H1c** equivalence within ±0.05 at m ≥ 10, k = 1 | **NOT SUPPORTED** | m = 10: difference **−0.087** [−0.107, −0.067] (outside the margin); m = 20: −0.023 [−0.043, −0.004] (inside). Both cells must pass |
| **H1d** uniform miss = 1 − n/N | **SUPPORTED** | z = −1.28, −0.74, +1.07, −1.46 (rule: \|z\| < 2.50) |
| **KT1** signals separable | **SUPPORTED** | AUC(information gain) 0.778 / 0.593 / 0.552 / 0.525 vs AUC(relevance) 0.977 / 0.943 / 0.920 / 0.896 at m = 1 / 5 / 10 / 20; m = 20 value inside [0.45, 0.55] |

## Miss rate at k = 1 (one critical field), matched cost ≈ 16.5% of the specimen
| cells per field | InfoDPP-PAC | uniform | difference |
|---|---|---|---|
| 1 | 0.315 | 0.822 | −0.507 |
| 5 | 0.681 | 0.830 | −0.149 |
| 10 | 0.759 | 0.846 | **−0.087** |
| 20 | 0.805 | 0.828 | −0.023 |

## What H1c's failure means (interpretation — not part of the pre-registered decision)
The pre-registered claim was that at ten or more cells per field, InfoDPP-PAC gives **no** protection beyond random
sampling of the same size. That is false at ten cells per field: it retains an advantage of about **9 points** (miss
0.76 vs 0.85). It is true at twenty (2 points, inside the margin). The correct statement is therefore:
**the advantage decays with dilution and becomes negligible between ten and twenty cells per field.**
The substantive point of Part I is unaffected in absolute terms: at ten cells per field the method still misses
three of four specimens that contain a single critical field, having examined one sixth of the specimen.

**Secondary analysis (not pre-registered), reported to characterise the failure:** the m = 10 difference varies by fold
from −0.033 to −0.160 (SD 0.053), i.e. it depends on which blast patients are in the fold. The fold-clustered 95% CI
is [−0.127, −0.047], also outside the margin, so the verdict does not depend on the CI method. At m = 20 the
clustered CI is [−0.037, −0.004], inside the margin.

**Methodological lesson (to state in the paper):** the pre-registered H1c interval treated episodes as independent,
although episodes within a fold share blast patients; it is narrower than the fold-clustered interval. Here both lead
to the same verdict, but for future equivalence claims the clustered interval is the appropriate one. The
exploratory estimate at m = 10 (−0.027) lay within the fold-to-fold spread now observed; with only 5 folds, that
spread was underestimated when the sample size was planned.

## Status of the study's claims after confirmation
- Content-blind stopping (H1a), separability of the signals (KT1), dilution as the governing axis (H1b) and the Part II
  finite-population formula (H1d) are **confirmed**.
- "No protection under dilution" is **confirmed only at 20 cells per field**; at 10 there is a residual, patient-dependent
  advantage of about 9 points.
