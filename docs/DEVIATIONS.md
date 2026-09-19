# Deviations from published methods

| Date | Method | Paper spec | Our implementation | Reason | Possible consequence |
|---|---|---|---|---|---|
| 2026-09-19 | InfoDPP-PAC halting | Eqs. 1–8, arXiv:2608.23574 | `sim/pilot.py:infosat_order` is a GENERIC GP max-variance + IG-saturation rule | Eqs. not yet transcribed | Pilot results are not evidence about InfoDPP-PAC itself |
| 2026-09-19 | KU-Optofil split | provided patient-level `split` (test: 3 blast patients, val: 0) | patient-grouped 5-fold cross-fitting (`scripts/02_make_folds.py`) | too few blast patients in test for critical-evidence evaluation | the classifier/encoder differs per fold; results pooled over folds; provided split kept for sensitivity |
| 2026-09-19 | KU-Optofil cells | 31,484 labelled | 28,280 used (3,204 `unknown`-patient cells excluded) | cannot be grouped by patient, so leakage risk; source unverified | fewer blasts (489 vs 814); fewer training cells |
| 2026-09-19 | InfoDPP-PAC hyperparameters | τ, δ, ρ, kernel, ℓ, β, m0 unspecified | calibrate τ on negative training-fold specimens to about 16% of budget, then sweep | the paper does not report them | the stopping point reflects our calibration; must be shown with a sensitivity sweep |
