# Blast separability in the selector embedding (RQ1 crux) — 2026-09-21

Run: `results/embeddings_v1/` (DINOv2 ViT-S/14 frozen, 28,280 cells, T4, `config.json`), analysis `scripts/04_separability.py` → `results/separability_v1/`.
All classifier scores are cross-fitted over the patient-grouped folds. CIs are bootstrap over **source patients**, not cells.

## Why this was measured
The synthetic pilot (`docs/status_w1.md`) showed that the sign of the InfoDPP-PAC-vs-uniform gap is set by how distinct
critical evidence is in the embedding used for selection, not by rarity alone. Two different quantities matter:
1. **Discriminability** — can a *classifier* separate blasts? (drives confidence-based stopping, and Part III's predictions)
2. **Outlierness** — do blasts sit *far from the bulk* in embedding space? (this is what a log-det/DPP diversity term responds to;
   the InfoDPP-PAC stopping statistic Γ̂ is a function of posterior variance, i.e. of geometry, never of relevance labels)

## Results
| Quantity | Value (95% CI, patient-clustered) |
|---|---|
| Cell-level discriminability, blast vs rest (cross-fitted logistic probe on frozen embeddings) | AUC **0.979** [0.969, 0.988] |
| Cell-level outlierness (AUC of 10-NN cosine distance to a background pool) | **0.901** [0.862, 0.939] |
| Median 10-NN distance | blast 0.0449 vs other 0.0272 |

**Most confusable classes** (blast vs that class, cross-fitted AUC): metamyelocyte 0.881, myelocyte 0.925, reactive lymphocyte 0.928, basophil 0.942, lymphocyte 0.942. The error concentrates on the immature-granulocyte continuum, i.e. exactly the cells that move a blast-percentage decision.

**Per-source-patient blast recall** at a fixed threshold (95% specificity on non-blasts): median 0.947, range **0.385–1.000** (11 patients). One patient (1604394691, 13 blasts) is at 0.385 while five are at ≥0.95.

**Field-level dilution** (a field = mean embedding of m cells, one of which is a blast in the positive case):

| cells per field m | outlierness AUC at field level |
|---|---|
| 1 | 0.891 |
| 5 | 0.673 |
| 10 | 0.592 |
| 20 | 0.555 |

**Image quality vs. cell-level error** (Pearson r): contrast +0.159, saturation +0.148, blur (Laplacian variance) +0.066, brightness −0.028. Blasts have higher mean contrast (0.212 vs 0.176) and saturation (0.273 vs 0.249) than other cells.

## Interpretation (what this does and does not establish)
1. **A single blast is an outlier; a field containing one blast is not.** Outlierness falls from 0.89 at m = 1 to 0.55 at m = 20. A diversity-driven selector operating on *field* embeddings therefore has little geometric signal to lock onto at realistic cell counts per field, while the same selector operating on *cell* crops would find blasts easily. **This is a mechanism, and it is sharper than the proposal's current claim.** The failure of information-saturation stopping should be attributed to *dilution of critical evidence within the observation unit*, not to rarity by itself; rarity and dilution are separate axes and the experiment can vary them independently.
2. **Refined hypotheses for pre-registration.** H1a (stop time is content-blind) is already supported analytically: the InfoDPP-PAC stopping statistic contains no relevance term (`docs/infodpp_pac_spec.md`). H1b should be restated as: the miss-rate gap against matched-cost uniform grows with **dilution m**, and with rarity only through the number of critical fields. The pilot's null for H1b is consistent with this, because the pilot had m = 1 by construction.
3. **Methodological vulnerability (field representation).** These numbers assume a field is represented by the **mean** of its cell embeddings. Mean-pooling is what causes the dilution. Max-pooling, attention pooling or a detector-style representation would preserve more of the signal. The pooling choice is therefore a hidden assumption that partly determines the headline result, so the thesis must (a) pre-register mean-pooling as the primary specification with a stated justification, and (b) report a sensitivity analysis over pooling.
4. **RQ4 signal already visible on composites.** Blast recall varies from 0.385 to 1.000 across the 11 source patients at a fixed threshold. Classifier error is therefore **patient-dependent**, which is direct evidence against the exchangeability assumed by the audit/inference machinery. Prediction-powered and risk-limiting methods that treat classifier error as identically distributed across sampled units inherit this. The quality correlations are weak (|r| ≤ 0.16) and are confounded with class (blasts differ in contrast and saturation), so they do **not** yet establish a stain-quality effect; the test has to be done **within** class.
5. **Limits.** The 10-NN distance is a proxy for what the InfoDPP-PAC statistic responds to, not the statistic itself (which uses GP posterior variance under an RBF kernel on the seeded posterior). All numbers come from 11 blast patients, so patient-clustered CIs are the honest ones and are wide. These are single-cell crops with no real field geometry.

## Next
Implement the InfoDPP-PAC rule (Eqs. 3–7) on composite fields built from these embeddings, then run Kill-test 1 (is the stopping signal separable from relevance?) and Kill-test 2 across a **rarity × dilution** grid.
