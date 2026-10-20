# Data (not tracked by git)

All four sources are public. Derived files and their SHA-256 are listed in `CHECKSUMS.txt`; how each result uses
them is in `../REPRODUCE.md`.

## KU-Optofil (primary platform; CC BY 4.0)
- Zenodo 10.5281/zenodo.17333317; descriptor: Yarıkan et al., *Scientific Data* (2026), doi:10.1038/s41597-026-06761-y.
- Put `dataset.zip`, `metadata.csv`, `metadata_with_patient_level_splits.csv` in `data/ku_optofil/` and extract
  `dataset.zip` to `data/ku_optofil/dataset/`.
- `scripts/02_make_folds.py` builds the patient-grouped folds; `notebooks/03_embeddings_colab.ipynb` (Colab, GPU) writes
  `results/embeddings_v1/` (embeddings, cross-fitted probe scores, image-quality features).

## RV-PBS (acquisition-shift probe; no licence stated)
- Pal, Bhattacharyea, Banerjee & Maharaj, *Expert Systems with Applications* 249:123660 (2024),
  doi:10.1016/j.eswa.2024.123660; data on GitHub (`Jimut123/RV-PBS`). Put the class folders with their
  `annotations.xml` in `data/rv_pbs/`.
- `scripts/14_crop_rvpbs.py` cuts the annotated cells into `data/rv_pbs_crops/`; zip them as `data/external_crops.zip`,
  upload to Drive, and run `notebooks/15_external_embeddings_colab.ipynb` (writes `results/external_v1/`).

## MLL23 (second shift source; CC BY 4.0)
- Shetab Boushehri et al., *Scientific Data* 12:1773 (2025), doi:10.1038/s41597-025-06223-x; Zenodo
  10.5281/zenodo.14277609. `notebooks/37_mll23_embeddings_colab.ipynb` downloads 15 of the 18 class archives
  (atypical promyelocytes, smudge cells and normoblasts are not used), checks their MD5 and writes `data/mll23/`
  (copy from Drive `certified-adaptive-pbs/mll23_v1/`). The release has no patient identifiers.

## BoneMarrowWSI-PediatricLeukemia (real spatial structure, tiles, reader labels; CC BY 4.0)
- NCI Imaging Data Commons collection `bonemarrowwsi_pediatricleukemia`; Zenodo 10.5281/zenodo.14933087 (conversion
  documentation, IDC manifests).
- Annotations: `notebooks/21_bmwsi_annotations_colab.ipynb` → `data/bmwsi/bmwsi_ann.zip`; extract to `data/bmwsi/ann/`
  and run `scripts/22_parse_bmwsi.py` (→ `cells.csv`, `rois.csv`; the series index is in the zip).
- Image tiles of the annotated regions: `notebooks/35_bmwsi_tiles_colab.ipynb` → `data/bmwsi_tiles/`.
- Per-round reader labels: `notebooks/42_bmwsi_sessions_colab.ipynb` → `data/bmwsi_sessions/`.
