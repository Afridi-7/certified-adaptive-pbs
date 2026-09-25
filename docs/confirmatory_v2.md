# Confirmatory results — Part III (RQ3) and RQ4 external source (2026-09-25)

Protocol: `docs/PREREG.md` Amendment 2, frozen as git tag `prereg-v2` before these data existed. Commands run exactly as
pre-registered; decisions from the frozen `scripts/13_part3_analysis.py` and `scripts/16_external_analysis.py`, unmodified.
Outputs: `results/part3_decisions.json`, `results/external_v1/decisions.json`.
RV-PBS embeddings: Colab T4, `external_crops.zip` SHA-256 42d195b1…9100, identical to the local file.

## Part III — certified threshold decisions (12,000 episodes; M = 500 cells; α = 0.05)

| Hypothesis | θ = 0.10 (ICC 2022) | θ = 0.20 (WHO 5th ed. 2022) |
|---|---|---|
| **H3a** validity of the four designs | **SUPPORTED** | **SUPPORTED** |
| **H3b** uncorrected adaptive selection is invalid | **SUPPORTED** | **SUPPORTED** |
| **H3c** classifier-driven sampling is cheaper | **SUPPORTED** | **SUPPORTED** |
| **H3d** cost governed by the margin | **SUPPORTED** | **SUPPORTED** |

**False certification at p = θ (1,000 episodes each; 95% Clopper–Pearson):**

| design | θ = 0.10 | θ = 0.20 |
|---|---|---|
| ALPHA, uniform | 0.022 [0.014, 0.033] | 0.012 [0.006, 0.021] |
| betting, uniform | 0.006 [0.002, 0.013] | 0.011 [0.006, 0.020] |
| betting, classifier-driven | 0.005 [0.002, 0.012] | 0.001 [0.000, 0.006] |
| + control variate | 0.006 [0.002, 0.013] | 0.005 [0.002, 0.012] |
| **uncorrected classifier-driven order** | **0.993 [0.986, 0.997]** | **0.889 [0.868, 0.908]** |

**Median cells examined to certify (paired, same specimens):**

| θ | true p | ALPHA uniform | classifier-driven | ratio | Wilcoxon p |
|---|---|---|---|---|---|
| 0.10 | 0.12 | 390 | 220 | 0.56 | 1.5e-114 |
| 0.10 | 0.15 | 165 | 75 | 0.45 | 2.9e-85 |
| 0.10 | 0.20 | 60 | 40 | 0.67 | 4.8e-60 |
| 0.20 | 0.22 | 445 | 335 | 0.75 | 1.7e-85 |
| 0.20 | 0.25 | 250 | 145 | 0.58 | 7.0e-65 |
| 0.20 | 0.30 | 105 | 60 | 0.57 | 3.3e-99 |

**Conventional 100-cell count, rate of calling "above":** θ = 0.10: 0.41 at p = θ, 0.139 at p = 0.08, 0.694 at p = 0.12.
θ = 0.20: 0.427 at p = θ, 0.243 at p = 0.18, 0.660 at p = 0.22.

## RQ4 — external source: RV-PBS (752 cells, 153 blasts, 583 images; KU-trained probe, KU-fixed threshold)

| Hypothesis | Verdict | Evidence (95% CI, bootstrap over source image) |
|---|---|---|
| **H4a** blast recall drops | **SUPPORTED** | RV-PBS 0.340 [0.250, 0.421] vs KU-Optofil 0.906 |
| **H4b** false-positive rate rises | **NOT SUPPORTED** | RV-PBS 0.058 [0.037, 0.081] vs KU-Optofil 0.050 |
| **H4c** certificate stays valid with shifted side information | **SUPPORTED** | false certification at p = θ: 3/1,000 = 0.003 |

Descriptive: AUC on RV-PBS 0.786 (KU cross-fitted 0.979). False positives concentrate in promyelocytes (0.18),
lymphocytes (0.15) and basophils (0.14). **Cost under shift (not pre-registered as a hypothesis):** median cells to certify at
p = 0.15: uniform 265, classifier-driven 265; at p = 0.20: uniform **110**, classifier-driven **190**.

## Interpretation
1. **Validity is design-based and held in every condition tested**: two thresholds, four designs, and a classifier that
   lost two thirds of its blast recall under a change of laboratory, stain and camera. False certification stayed between
   0.001 and 0.022 throughout; the only design that failed is the one that ignores how the sample was drawn (0.89–0.99).
2. **Efficiency is model-based, and it can reverse.** With in-distribution side information, classifier-driven sampling
   cut the cells needed by 25–55%. Under acquisition shift it gave **no** saving at p = 0.15 and cost **73% more** cells than
   uniform sampling at p = 0.20 (190 vs 110): a classifier that now misses most blasts steers the sample away from them,
   and the importance weights pay for it in variance. Adaptivity is a wager on the side information; the guarantee is not.
3. **The shift is in sensitivity, not in false alarms.** Recall collapsed (0.91 → 0.34) while the false-positive rate stayed
   within uncertainty of its KU value, so H4b was not supported: the directional expectation was wrong.
4. **Threshold decisions near the boundary** remain close to a census for any valid design (220–445 of 500 cells), while the
   conventional 100-cell count calls "above" 41–43% of the time when p sits exactly on the boundary.

## Threats and limits
- **Cropping confound**: RV-PBS crops use our rule (1.5× the polygon's bounding box, resized to 368 px); KU-Optofil's rule is
  undocumented. Part of the recall drop may come from scale and framing rather than from stain or camera; this study cannot
  separate them.
- **Clustering**: RV-PBS gives no patient identifiers; the bootstrap is over images, which may share patients, so its
  intervals may be too narrow.
- **H4c specimens** reuse the 599 non-blast RV-PBS crops heavily across episodes (each episode draws 450 of them); episodes are
  therefore dependent, which matters for the cost estimates more than for the validity check.
- The cost reversal under shift was **not** a pre-registered hypothesis and should be reported as an exploratory finding. A
  defensive mixture of the classifier-driven and uniform distributions is the natural safeguard; not tested here.
