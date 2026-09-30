# Confirmatory results — Amendment 3: spatial clustering and a defensive mixture (2026-09-25)

Protocol: `docs/PREREG.md` Amendment 3, frozen as git tag `prereg-v3` before these data existed. The commands were
run exactly as pre-registered. Decisions come from the frozen `scripts/20_amend3_analysis.py`, unmodified.
Output: `results/a3_decisions.json`.
Episodes: 15,000 for E3.1 (3 arrangements × 2 θ × 5 p) and 5,000 for E3.2.

## Verdicts
| Hypothesis | Verdict | Evidence |
|---|---|---|
| **H5a** field-level certificates stay valid under clustering | **SUPPORTED** | all 36 cells pass; maximum false-certification rate 0.042 (ALPHA, hot0.6, θ = 0.20) |
| **H5b** analysing field-sampled cells as individual draws breaks validity under clustering | **SUPPORTED** | false certification at p = θ: scattered 0.029 (θ = 0.10) / 0.041 (θ = 0.20); **hot1.0 0.327 [0.298, 0.357] / 0.478 [0.447, 0.509]**; hot0.6 0.236 / 0.278 (descriptive) |
| **H5c** clustering raises the cost of valid certification | **SUPPORTED** | median cells (ALPHA, uniform over fields): θ = 0.10, p = 0.15: 160 → 320 → 345; θ = 0.20, p = 0.25: 235 → 370 → 410 |
| **H6a** the mixture is valid | **SUPPORTED** | false certification at p = θ: 0.001–0.014 across sources and λ |
| **H6d** the cost reversal under shift replicates | **SUPPORTED** | RV-PBS, p = 0.20: classifier-driven 195 vs uniform 110 cells (paired Wilcoxon p = 1.5e-54) |
| **H6b** the mixture caps the loss under shift below 20% | **NOT SUPPORTED** | median(mix50)/median(uniform): p = 0.15: 1.10 [1.04, 1.16] (passes); **p = 0.20: 1.23 [1.11, 1.30]** (fails). Both p must pass |
| **H6c** the mixture keeps most of the gain in distribution | **SUPPORTED** | share of the classifier-driven saving retained: 0.58 / 0.62 / 0.73 at p = 0.12 / 0.15 / 0.20 |

## Descriptive results (not hypotheses)
- **Classifier-driven sampling matters more when blasts cluster.** Median cells at θ = 0.10, p = 0.15 (uniform →
  classifier-driven):
  - scattered: 160 → 80;
  - hot0.6: 320 → 60;
  - hot1.0: 345 → 55.

  The classifier finds the hot fields, whereas uniform sampling pays the design effect.
- **The conventional 100-cell count degrades under clustering.** At p = θ − 0.02 (truly below threshold), it calls
  "above" in:
  - θ = 0.10: 0.158 (scattered), 0.254 (hot0.6), 0.216 (hot1.0);
  - θ = 0.20: 0.212, 0.320, 0.266.

  At p = θ + 0.05 it calls "above" in only 0.61–0.77 of the clustered specimens, vs 0.92 at θ = 0.10 and 0.87 at
  θ = 0.20 when scattered.
- **Mixture, λ = 0.25.** It keeps more of the in-distribution gain (e.g. 115 vs 145 cells at p = 0.15), but it
  gives almost no protection under shift: 1.73× uniform at p = 0.20, the same as pure classifier-driven sampling.

## Interpretation
1. **The dependence between cells is the analysis unit's problem, not the sampling design's.** When the unit that
   is sampled (the field) is the unit that is analysed, validity holds under every arrangement tested. When cells
   from sampled fields are treated as if each had been drawn individually, the false-certification rate reaches
   0.33–0.48 under clustering, 7–10× the nominal 0.05.

   This is the RQ4 assumption violation with the largest measured effect in this study. It is also the implicit
   model behind reading a differential count off consecutive fields.
2. **Clustering moves cost in opposite directions for the two samplers.** It roughly doubles the cost for uniform
   field sampling (the design effect), and it slightly reduces the cost for classifier-driven sampling, which
   exploits it. The value of side information therefore depends on spatial structure that the available data
   cannot measure.
3. **The mixture is a partial safeguard, not a guarantee of efficiency.** At λ = 0.5 it keeps 58–73% of the
   in-distribution saving and cuts the loss under shift from +77% to +23% at p = 0.20. It does not meet the
   pre-registered 20% bound. Validity is never at stake; only cost is.

## Threats and limits
- The hot-field model is a stress model with two density levels. No data here estimate how blasts are actually
  distributed across fields in real smears, so the size of the H5b and H5c effects in practice is unknown. Only
  their direction and mechanism are established.
- Fields hold 5 cells. The design effect grows with cells per field, so larger fields would enlarge both the H5b
  failure and the H5c cost increase.
- The E3.2 RV-PBS episodes reuse the same 599 non-blast crops heavily (as in H4c).
- The expert is simulated with perfect labels (DEVIATIONS.md).
