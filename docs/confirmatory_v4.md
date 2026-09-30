# Confirmatory results — Amendment 4: spatial structure of blasts in real smears (2026-09-30)

Protocol: `docs/PREREG.md` Amendment 4, frozen as git tag `prereg-v4` before these analyses were run.

**Where it ran.** The two frozen scripts were run unmodified in the cloud workspace, because the analysis exceeds the
local shell's 3-minute limit. Script SHA-256:
- `23_bmwsi_clustering.py`: 2223ef5e…d91c
- `24_bmwsi_bridge_sim.py`: 970b3ebc…20e6

Inputs `data/bmwsi/cells.csv` and `rois.csv` were copied from the laptop.

**Outputs.** `results/a4_bmwsi/clustering.json` (SHA-256 68da9986…4908) and `bridge.json` (1690a593…3025).

**Data.** BoneMarrowWSI-PediatricLeukemia, fully labelled monolayer ROIs. **Eligible**: 75 patients, 116 ROIs, 9,375
nucleated cells, 5,841 blasts. An ROI is eligible with ≥ 20 nucleated cells and at least one blast and one non-blast;
ROIs with no blasts are uninformative about clustering and were excluded by rule.

## A4.1 — is blast status spatially clustered? H7a **SUPPORTED**
| Scale | ICC (95% patient-bootstrap CI) | Permutation null mean (95th pct) | One-sided p |
|---|---|---|---|
| **S1** quadrants within ROI, ~113 µm (primary) | **0.032** [0.016, 0.050] | 0.018 (0.030) | **0.024** |
| S2 4 × 4 within ROI, ~56 µm | 0.098 [0.061, 0.138] | 0.056 (0.074) | 0.0005 |
| S3 between the two ROIs of a patient | 0.007 [−0.001, 0.016] | 0.004 (0.012) | 0.20 |

- **k = 5 nearest-neighbour enrichment.** The fraction of a blast's neighbours that are blasts is 0.914, against
  0.907 under permutation (p = 0.0005).
- **Sensitivity.**
  - Promonocytes counted as blasts: S1 0.032, p = 0.030. Unchanged.
  - Damaged and smudge cells kept as non-blasts: S1 0.035, p = 0.0005. S3 becomes significant (0.037, p = 0.0005),
    which suggests these cells are distributed unevenly between ROIs.

**Reading.** Blasts are clustered at the scale of tens of micrometres, within a field-sized region. The effect is
statistically clear but small. There is no detectable difference between the two ROIs of a patient; ROIs lie millimetres
apart within the monolayer.

**Estimator bias (noticed after the run; exploratory).** The permutation null of the ICC is not centred at 0: its
mean is 0.018 at S1 and 0.056 at S2. The pooled ANOVA estimator is biased upward for binary data with small, unequal
subfields. The permutation test (H7a) is unaffected. The **excess** over the null, 0.014 at S1 and 0.042 at S2, is a
better estimate of true clustering than the raw ICC. The raw ICC and its CI overstate it.

## A4.2 — what the measured clustering implies (model-based; pre-specified reading rule)
Naive cell-level analysis of field-sampled cells, false certification at p = θ. Specimens are 30 fields × 17 cells;
1,000 episodes per cell.

| ρ used | θ = 0.10: field-level / **naive** | θ = 0.20: field-level / **naive** |
|---|---|---|
| 0 (no clustering) | 0.008 / 0.033 | 0.001 / 0.035 |
| 0.016 (S1 lower CI) | 0.011 / 0.043 | 0.003 / 0.051 |
| 0.032 (S1 estimate) | 0.012 / **0.061** | 0.004 / **0.095** |
| 0.050 (S1 upper CI) | 0.014 / **0.070** | 0.006 / **0.103** |

- **Pre-specified rule.** Naive false certification at the upper CI bound exceeds 0.05 at both θ (0.070 and 0.103).
  **Reading: the H5b failure is plausible at the clustering level observed in real smears**, in this model-based
  sense.
- **Qualification (exploratory, due to the estimator bias above).** At a bias-corrected ρ of about 0.014, close to
  the S1 lower bound, the naive rate is 0.043–0.051, i.e. at the nominal level. The honest statement is therefore:
  at realistic monolayer clustering, the cell-level mis-analysis inflates false certification **from 0.05 to between
  about 0.05 and 0.10**.
- The field-level test stayed at or below 0.014 throughout.
- **E3.1 placement on the same ICC scale.** Scattered ≈ 0.00, hot0.6 ≈ 0.49–0.54, hot1.0 = 1.0. The E3.1 stress
  arrangements are **one to two orders of magnitude more clustered** than the real monolayer ROIs. The 33–48% false
  certification of E3.1 H5b describes an extreme, not typical, specimen.

## Interpretation
1. Real smears are **not exchangeable at the cell level**. Blasts cluster measurably within field-sized regions, which
   answers the RQ4 spatial-dependence question with real data.
2. At the clustering measured in evenly dispersed monolayer regions, analysing cells as individual draws is **mildly**
   anti-conservative: about 1–2× the nominal error, according to the model. It becomes severely anti-conservative only
   under strong clustering (E3.1). Sampling and analysing by field removes the problem at every level tested.
3. For the paper, the safe recommendation, analyse at the unit that was sampled, holds regardless. The size of the
   danger is modest in the monolayer and larger wherever blasts are more clustered.

## Limits
- Pediatric bone-marrow aspirates, not peripheral blood.
- ROIs were selected in the monolayer, so clustering elsewhere in the smear may be larger.
- Two ROIs per patient.
- The ICC estimator bias was noticed after the run (above).
- A4.2 is a Beta-binomial model, not an observed false-certification rate.
