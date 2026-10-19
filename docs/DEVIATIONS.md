# Deviations from published methods

| Date | Method | Paper spec | Our implementation | Reason | Possible consequence |
|---|---|---|---|---|---|
| 2026-09-19 | InfoDPP-PAC halting | Eqs. 1–8, arXiv:2608.23574 | `sim/pilot.py:infosat_order` is a GENERIC GP max-variance + IG-saturation rule | Eqs. not yet transcribed | Pilot results are not evidence about InfoDPP-PAC itself |
| 2026-09-19 | KU-Optofil split | provided patient-level `split` (test: 3 blast patients, val: 0) | patient-grouped 5-fold cross-fitting (`scripts/02_make_folds.py`) | too few blast patients in test for critical-evidence evaluation | the classifier/encoder differs per fold; results pooled over folds; provided split kept for sensitivity |
| 2026-09-19 | KU-Optofil cells | 31,484 labelled | 28,280 used (3,204 `unknown`-patient cells excluded) | cannot be grouped by patient, so leakage risk; source unverified | fewer blasts (489 vs 814); fewer training cells |
| 2026-09-19 | InfoDPP-PAC hyperparameters | τ, δ, ρ, kernel, ℓ, β, m0 unspecified | calibrate τ on negative training-fold specimens to about 16% of budget, then sweep | the paper does not report them | the stopping point reflects our calibration; must be shown with a sensitivity sweep |
| 2026-09-19 | InfoDPP-PAC embedding | Phikon (H&E histopathology, 768-d) | DINOv2 ViT-S/14 frozen (`vit_small_patch14_dinov2.lvd142m`, 384-d, 224 px) | no established open Wright–Giemsa encoder; Phikon is trained on H&E tissue | geometry of the embedding (hence stopping behaviour and blast separability) is encoder-specific; single encoder, no comparisons (scope) |
| 2026-09-21 | InfoDPP-PAC Eq. 7 | ρ is the Assumption-1 bounded-difference constant | ρ = 1 by fiat; τ calibrated to the paper's reported operating point (~16.3% of budget) on blast-free specimens from non-evaluation folds | neither ρ nor τ is reported in the paper | the implemented rule reproduces the paper's *cost*, not its nominal PAC guarantee; the certificate's δ must not be claimed as the paper's |
| 2026-09-21 | Field representation | not applicable (the paper's units are image patches) | a field is the mean of its cells' embeddings | composites have to define a field-level representation | mean pooling causes the dilution effect; max/single pooling are implemented for sensitivity |
| 2026-09-21 | Raster-scan baseline | spatial raster over a slide | fixed index order over composite fields | composites carry no geometry | reported as `fixed_order`, not as a spatial baseline |
| 2026-09-25 | Part III weighted test | Shekhar et al. (arXiv:2305.06884) Sec. 2 | v1 used an invalid HT-into-ALPHA construction; **withdrawn** and replaced by a transcription of Shekhar et al.'s importance-weighted betting martingale (`weighted_betting.py`) | error found on reading the paper | v1 finding 3 ("adaptivity did not help") is withdrawn; v2 shows a 35–58% cost reduction with valid size |

## Erratum (2026-09-25): date in PREREG.md Amendment 2 header
The Amendment 2 header reads "2026-09-27"; the correct date is 2026-09-25. This is a typing error in the frozen text and is
left unedited so that the text still matches tag `prereg-v2`. The authoritative freeze time is the timestamp of the `prereg-v2`
commit, which precedes the creation of `results/part3_theta010`, `results/part3_theta020` and `results/external_v1/decisions.json`.

## Additions (2026-09-25)
| Date | Component | Proposal / source | What we did | Why | Consequence |
|---|---|---|---|---|---|
| 2026-09-25 | Part III side information | Proposal lists prediction-powered e-values (Csillag, Struchiner & Goedert; arXiv identifier to be verified) | **Not implemented.** Side information enters only through Shekhar et al.'s prop-MS sampling and control variate | one valid weighted-sampling framework covers the adaptive design; prediction-powered e-values were not needed to test H3a–H3d | the thesis reports no comparison with prediction-powered e-values; this must be stated as a limitation, not implied as covered |
| 2026-09-25 | Expert verification | a haematologist verifies the sampled fields | verification is simulated by the dataset's expert labels (verifier error β = 0) | no haematologist annotation available | certificates hold for a perfect verifier; `part2_proof.md` Remark 4 shows that verifier error sets a floor |

**Erratum addendum.** The freeze line of Amendment 2 reads "2026-09-26". The tag was created on 2026-09-25, and the
`prereg-v2` tag timestamp is authoritative.

## Note on document dates (2026-10-03)
Several documents carry dates of 2026-09-26 or 2026-09-27: the PREREG v1/v1.1/v2 headers and freeze lines, confirmatory_v1, rq4_v1, and the killtests_v1 and sensitivity_v1 notes. These were typed by hand and may not match when the files were actually written. The **git commit and tag timestamps are authoritative**. The pre-registration claim rests on the ordering of the tags relative to the result commits, not on the typed dates.

## Audit trail and episode counts (2026-10-19)
`docs/PREREG_AUDIT.md` tabulates, for every tag, the commit time, the GitHub push time and the first confirmatory output.
The output always follows the push, by 43 s to 18 min. It also lists the typed dates that are wrong, and content hashes
of the outputs.

Episode counts reported in `confirmatory_v3.md` and in the paper draft of 2026-10-03 were wrong. The files contain the
counts that Amendment 3 specifies:
- **E3.1: 18,000 episodes.** That is 3 arrangements × 2 θ × (1,000 at p = θ + 4 × 500), against the 15,000 written.
- **E3.2: 4,500 episodes.** RV-PBS has 1,000 + 2 × 500 and KU-Optofil has 1,000 + 3 × 500, against the 5,000 written.

No analysis changes; only the reported totals change.

## Pre-registered sensitivity analyses run late (2026-10-19)
PREREG §6 and Amendment 1 (A4) fix four sensitivity analyses: max pooling (seed 20260402), single-cell pooling (20260403), a
τ sweep (20260404) and a single-host background (20260405). They had not been run with these seeds; the paper cited the
exploratory `sensitivity_v1` runs instead (500 episodes per cell, reusing the main run's random streams). They were run on
2026-10-19 with the frozen `05_killtests.py`, unchanged, and 200 episodes per (fold, m, k).
- **τ sweep scope.** The sweep was restricted to k = 1, the H1c cell, because §6 does not state its grid.
- **Patient-level bootstrap.** §6.4 could not be done from the frozen outputs, because they do not record which blast
  patient supplied each critical field. Fold-level intervals are reported instead. The Kill-test 1 re-analysis of
  2026-10-19 records patients.

## Correction (2026-10-03): RV-PBS citation
Earlier documents attributed RV-PBS to "Sen et al." and "Vergara-Lluri et al."; both attributions were wrong. The correct citation is J. B. Pal, A. Bhattacharyea, D. Banerjee, T. Maharaj, *Expert Systems with Applications* 249:123660 (2024), doi:10.1016/j.eswa.2024.123660. This was checked against the dataset's GitHub README. `data_audit_rv_pbs.md` and `runbook.md` are corrected.
