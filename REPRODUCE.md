# Reproducing the paper

This file explains how to regenerate every number, table and figure of `paper/main.tex` (full preprint) and
`paper/jbhi_compact.tex` (journal version). `reproduce.py` wraps the same commands.

## 1. Three levels

| Level | What it shows | Needs | Time |
|---|---|---|---|
| **A. Checks without data** | The validity claims hold (exact nulls, adversarial and adaptive samplers, unequal fields, verifier margin), the Part II theorems are verified exactly, and our ALPHA and betting code agree with the reference implementation | Python only (+ one SHANGRLA file for the cross-check) | about 5 min |
| **B. Re-run analyses from derived data** | Every Part I, Part III and RQ4 result, from the embeddings and annotation tables | the derived files of Section 3 | minutes to hours per step (Section 4) |
| **C. Full pipeline from raw data** | Also the embeddings themselves | the public datasets and a GPU (Google Colab is enough) | about a day |

```
python reproduce.py check --shangrla <folder with NonnegMean.py>    # level A
python reproduce.py list                                             # all steps
python reproduce.py verify mll23 verifier real_arrangements          # level B: re-run into a temporary folder and compare
python reproduce.py run partIII_confirm                              # level B: re-run in place (then: git diff results/)
```
`verify` compares the regenerated summary files with the stored ones (numbers within 1e-9). It was run on 2026-10-20
for `part2`, `crosscheck`, `detrule`, `decisions` and `alpha_prior`: every compared value was identical (and the
`part2` files byte-identical).

The reference implementation is SHANGRLA (github.com/pbstark/SHANGRLA, file `shangrla/core/NonnegMean.py`,
AGPL-3.0). It is not included here; download that one file into a folder and pass the folder with `--shangrla`. The
version used on 2026-10-20 had SHA-256 `1fb7cdce04324b28f9093c9f8d79f12d2190190103667617ecddb41cfe9efde6`.

## 2. Environment

```
python -m venv .venv
# Windows: .venv\Scripts\Activate.ps1      Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt           # ranges; requirements-lock.txt has the exact versions used most recently
```
GPU notebooks run in Google Colab (T4) and install their own packages; their package versions are recorded in each
output's `config.json`.

The frozen confirmatory runs (2026-09-25 to 2026-09-30) did not record package versions. Their outputs are identified
by the content hashes in `docs/PREREG_AUDIT.md`. Re-running one Part I cell later in a separate environment reproduced
the stored episodes (identical miss indicators; 2 of 1,000 stop steps differed by one field because the calibrated
threshold differed in the sixth significant digit).

## 3. Data

Nothing in `data/` is tracked by git. See `data/README.md` for sources, licences and download steps, and
`data/CHECKSUMS.txt` for the SHA-256 of every derived file used in the paper. Derived files that a level-B
reproduction needs:

| File | Made by | Used by |
|---|---|---|
| `results/embeddings_v1/embeddings_dinov2_vits14.npy`, `cells_features.csv` | `notebooks/03_embeddings_colab.ipynb` | Part I, Part III, RQ4 (KU-Optofil) |
| `results/external_v1/rvpbs_embeddings_dinov2_vits14.npy`, `rvpbs_features.csv` | `scripts/14_crop_rvpbs.py`, `notebooks/15_external_embeddings_colab.ipynb` | RV-PBS shift |
| `data/bmwsi/cells.csv`, `rois.csv`, `series_index.csv` | `notebooks/21_bmwsi_annotations_colab.ipynb`, `scripts/22_parse_bmwsi.py` | clustering, real arrangements, unequal fields |
| `data/bmwsi_tiles/tiles.csv`, `tiles_embeddings_dinov2_vits14.npy` | `notebooks/35_bmwsi_tiles_colab.ipynb` | Part I on real tiles |
| `data/mll23/mll23_cells.csv`, `mll23_embeddings_dinov2_vits14.npy` | `notebooks/37_mll23_embeddings_colab.ipynb` | MLL23 shift |
| `data/bmwsi_sessions/annotations_sessions.csv.gz` | `notebooks/42_bmwsi_sessions_colab.ipynb` | verifier error |

All source datasets are public (KU-Optofil, MLL23 and BoneMarrowWSI under CC BY 4.0; RV-PBS without a stated
licence). The derived embeddings and tables can be shared under the same terms and should be included in the
archived release (Zenodo), so that level B needs no GPU.

## 4. Steps

Status: P pre-registered (confirmatory), X exploratory, C check, D data preparation. Cost: S < 5 min, M 5-60 min,
L > 1 h on a laptop CPU (approximate); GPU = Colab notebook. Seeds are fixed in the commands; exploratory seed families
never overlap with the confirmatory ones (`docs/PREREG.md`). Paper items refer to `paper/main.tex`; the compact
version numbers some propositions and tables differently.

| ID | Status | Cost | Paper item | How |
|---|---|---|---|---|
| `tests` | C | S | Validity claims; Section 'Implementation checks' | `python -m pytest -q tests` |
| `part2` | C | S | Prop. 1, Theorem 1, Corollary (exact and LP checks) | `python scripts/17_part2_verify.py --out results/part2_v1` |
| `crosscheck` | C | S | Implementation checks (reference ALPHA, SHANGRLA) | `python -I scripts/44_crosscheck_shangrla.py --shangrla <SHANGRLA folder> --R 2000 --seed 20261035 --out results/crosscheck_shangrla_v1` |
| `ku_inspect` | D | S | Data (KU-Optofil counts) | `python scripts/01_inspect_kuoptofil.py --data data/ku_optofil` |
| `ku_folds` | D | S | Patient-grouped folds | `python scripts/02_make_folds.py --data data/ku_optofil --out results/folds_v1 --seed 20260919` |
| `ku_embeddings` | D | GPU | Embeddings and cross-fitted probe (AUC 0.979) | Colab: `notebooks/03_embeddings_colab.ipynb` |
| `separability` | X | S | Blast separability by pooling (background for Kill-test 1) | `python scripts/04_separability.py --out results/separability_v1` |
| `partI_confirm` | P | L | Fig. 1; H1a-H1d; Kill-test 1 | `python scripts/05_killtests.py --out results/confirm_v1 --seed 20260401 --episodes 200 --cell_episodes 10:1:600 20:1:600`<br>`python scripts/09_confirmatory_analysis.py --dir results/confirm_v1` |
| `partI_sens_max` | P | L | Sensitivity: max pooling | `python scripts/05_killtests.py --out results/sensconf_pool_max --seed 20260402 --episodes 200 --pooling max` |
| `partI_sens_single` | P | L | Sensitivity: single-cell representation | `python scripts/05_killtests.py --out results/sensconf_pool_single --seed 20260403 --episodes 200 --pooling single` |
| `partI_sens_tau` | P | L | Sensitivity: operating point 5-50% | `python scripts/05_killtests.py --out results/sensconf_tau_0.05 --seed 20260404 --episodes 200 --k_list 1 --target_frac 0.05`<br>`python scripts/05_killtests.py --out results/sensconf_tau_0.10 --seed 20260404 --episodes 200 --k_list 1 --target_frac 0.10`<br>`python scripts/05_killtests.py --out results/sensconf_tau_0.163 --seed 20260404 --episodes 200 --k_list 1 --target_frac 0.163`<br>`python scripts/05_killtests.py --out results/sensconf_tau_0.30 --seed 20260404 --episodes 200 --k_list 1 --target_frac 0.30`<br>`python scripts/05_killtests.py --out results/sensconf_tau_0.50 --seed 20260404 --episodes 200 --k_list 1 --target_frac 0.50` |
| `partI_sens_host` | P | L | Sensitivity: single-host background | `python scripts/05_killtests.py --out results/sensconf_pool_host --seed 20260405 --episodes 200 --pool host` |
| `partI_op24` | X | M | InfoDPP-PAC at its own operating point (2.4%) | `python scripts/05_killtests.py --out results/sens_op_paperpool --seed 20261022 --episodes 200 --k_list 1 2 5 --target_frac 0.024` |
| `partI_field_teacher` | X | M | Kill-test 1 with the field-level teacher | `python scripts/26_kt1_field_teacher.py --out results/kt1_field_teacher --seed 20261019` |
| `partI_rules` | X | M | Other saturation rules; ranked search (Prop. 3) | `python scripts/31_stopping_class.py --out results/stopping_class_v1 --seed 20261024` |
| `review_stats` | X | S | Fold t intervals, effect sizes, sensitivity table | `python scripts/29_review_stats.py` |
| `tiles_embeddings` | D | GPU | Real BoneMarrowWSI tiles | Colab: `notebooks/35_bmwsi_tiles_colab.ipynb` |
| `partI_tiles` | X | S | Part I on real tiles | `python scripts/36_partI_real_tiles.py --tiles data/bmwsi_tiles --out results/partI_real_tiles_v1 --seed 20261028`<br>`python scripts/40_partI_real_tiles_analysis.py` |
| `partIII_confirm` | P | M | Fig. 3; H3a-H3d | `python scripts/12_certify_confirm.py --out results/part3_theta010 --theta 0.10 --p_list 0.05 0.08 0.10 0.12 0.15 0.20 --episodes 1000 --seed 20260511`<br>`python scripts/12_certify_confirm.py --out results/part3_theta020 --theta 0.20 --p_list 0.15 0.18 0.20 0.22 0.25 0.30 --episodes 1000 --seed 20260512`<br>`python scripts/13_part3_analysis.py --dirs results/part3_theta010 results/part3_theta020` |
| `m_scaling` | X | L | Population size table | `python scripts/28_m_scaling.py --out results/m_scaling_v1 --seed 20261021` |
| `alpha_prior` | X | M | ALPHA prior guess | `python scripts/39_alpha_prior.py --out results/alpha_prior_v1 --seed 20261030` |
| `decisions` | X | M | Decisions at fixed budgets (table) | `python scripts/32_decision_level.py --out results/decision_level_v1 --seed 20261025` |
| `confusable` | X | M | Confusable backgrounds | `python scripts/33_confusable_background.py --out results/confusable_bg_v1 --seed 20261026` |
| `leukocyte` | X | M | Leukocyte-only background | `python scripts/41_leukocyte_background.py --out results/leukocyte_bg_v1 --seed 20261032` |
| `detrule` | X | M | Impact of the rule missing from the frozen ALPHA | `python scripts/45_alpha_deterministic_rule.py --out results/alpha_detrule_v1 --seed 20261036` |
| `unequal` | X | M | Unequal field sizes; Prop. 2 | `python scripts/27_unequal_fields.py --out results/unequal_fields_v1 --seed 20261020` |
| `ku_error` | X | S | Classifier error is not exchangeable | `python scripts/10_rq4_error_structure.py`<br>`python scripts/11_rq4_certify_by_error_group.py` |
| `rvpbs_crops` | D | S | RV-PBS crops | `python scripts/14_crop_rvpbs.py` |
| `rvpbs_embeddings` | D | GPU | RV-PBS embeddings | Colab: `notebooks/15_external_embeddings_colab.ipynb` |
| `rvpbs` | P | M | H4a-H4c (acquisition shift) | `python scripts/16_external_analysis.py` |
| `spatial_mixture` | P | L | H5a-H5c, H6a-H6d; Fig. 4 | `python scripts/18_spatial_confirm.py --out results/a3_spatial_theta010 --theta 0.10 --p_list 0.08 0.10 0.12 0.15 0.20 --seed 20260611`<br>`python scripts/18_spatial_confirm.py --out results/a3_spatial_theta020 --theta 0.20 --p_list 0.18 0.20 0.22 0.25 0.30 --seed 20260612`<br>`python scripts/19_mixture_confirm.py --out results/a3_mixture_rvpbs --source rvpbs --p_list 0.10 0.15 0.20 --seed 20260621`<br>`python scripts/19_mixture_confirm.py --out results/a3_mixture_ku --source ku --p_list 0.10 0.12 0.15 0.20 --seed 20260622`<br>`python scripts/20_amend3_analysis.py --spatial results/a3_spatial_theta010 results/a3_spatial_theta020 --mixture results/a3_mixture_rvpbs results/a3_mixture_ku --out results/a3_decisions.json` |
| `bmwsi_annotations` | D | GPU | BoneMarrowWSI annotations (no GPU needed) | Colab: `notebooks/21_bmwsi_annotations_colab.ipynb` |
| `bmwsi_parse` | D | S | cells.csv, rois.csv | `python scripts/22_parse_bmwsi.py` |
| `bmwsi_clustering` | P | M | H7a; Beta-binomial bridge | `python scripts/23_bmwsi_clustering.py --out results/a4_bmwsi/clustering.json`<br>`python scripts/24_bmwsi_bridge_sim.py --clustering results/a4_bmwsi/clustering.json --out results/a4_bmwsi/bridge.json` |
| `icc_excess` | X | M | Excess ICC | `python scripts/34_icc_excess.py --out results/a4_bmwsi/icc_excess.json --seed 20261027` |
| `real_arrangements` | X | M | Observed validity on 75 real specimens | `python scripts/30_real_arrangement_validity.py --out results/real_arrangements_v1 --seed 20261023` |
| `mll23_embeddings` | D | GPU | MLL23 embeddings | Colab: `notebooks/37_mll23_embeddings_colab.ipynb` |
| `mll23` | X | S | Second external source (MLL23) | `python scripts/38_mll23_shift.py --mll data/mll23 --out results/mll23_shift_v1 --seed 20261029` |
| `sessions_download` | D | GPU | Reader labels (no GPU needed) | Colab: `notebooks/42_bmwsi_sessions_colab.ipynb` |
| `verifier` | X | M | Imperfect verification; Remark 2 | `python scripts/43_verifier_error.py --sessions data/bmwsi_sessions/annotations_sessions.csv.gz --out results/verifier_error_v1 --seed 20261034` |
| `benchmark` | X | S | Benchmark table across conditions (main Table 3, compact Table IV) | `python scripts/46_benchmark_table.py --out results/benchmark_v1` |
| `figures` | C | S | Figs. 1-4 (copy figures/*.pdf to paper/figures/) | `python scripts/25_figures.py` |


Figures: `python scripts/25_figures.py` writes `figures/fig1_partI.pdf` ... `fig4_rq4.pdf`; copy them to
`paper/figures/`. Papers: `cd paper`, then `latexmk -pdf main.tex` and `latexmk -pdf jbhi_compact.tex`
(the latter needs the IEEEtran class, included in TeX Live as `texlive-publishers`).

## 5. What is frozen

The pre-registration lives in `docs/PREREG.md`, frozen by the git tags `prereg-v1` ... `prereg-v4`
(`git tag`; push times in `docs/PREREG_AUDIT.md`). The confirmatory scripts (05, 09, 12, 13, 16, 18, 19, 20, 23, 24)
and the modules they import (`src/cadpbs/alpha.py`, `weighted_betting.py`, `infodpp.py`, `composites.py`) are
unchanged since the tag at which each was frozen (checked on 2026-10-20 with `git diff <tag> -- <file>`: 05, 09, infodpp.py and composites.py against prereg-v1.1; 12, 13, 16, alpha.py and weighted_betting.py against prereg-v2; 18, 19, 20 against prereg-v3; 23, 24 against prereg-v4). Everything added later is exploratory and
listed, with every deviation, in `docs/DEVIATIONS.md`.
