# Real-data candidates for the paper (2026-09-25)

**Goal.** Obtain real specimens in which every cell of a field or region is labelled. This addresses three reviewer
objections: specimens are synthetic, only 11 blast patients, and no measured spatial structure.

**Serves.** RQ4 (spatial dependence, exchangeability), and RQ3 on real specimens instead of composites.

## Search result
A 2025 scoping review of leukemia datasets reports that no public set combines whole field-of-view images,
cell-level labels and complete annotation (MethodsX, S2215016125005667). This confirms the gap. Two candidates come
close.

| Dataset | What it is | Fit | Unknowns / caveats |
|---|---|---|---|
| **BoneMarrowWSI-PediatricLeukemia** (Höfener et al.; IDC collection `bonemarrowwsi_pediatricleukemia`; source DOI 10.5281/zenodo.14933087; **CC BY 4.0**, verified in the IDC index) | 245 pediatric ALL/AML/CML cases, bone-marrow aspirate WSIs (40×, 0.11 µm/px). **Verified from the IDC index:** annotation series are monolayer ROIs (245 cases), consensus cell labels (113 cases), unlabeled cell boxes (97 cases), six annotator sessions. Annotations total 17.3 MB; slides 2.3 TB. Per-case lab blast values are listed. | Real WSIs, so cell coordinates are real, allowing **spatial clustering of blasts to be measured directly**. Real per-case blast proportions. About 47k labelled cells, roughly 400 per labelled case (inferred, not yet verified). | **Bone marrow, not peripheral blood.** The WHO/ICC blast criteria cover blood *or* marrow, but it is a different specimen type and the title would need adjusting. Pediatric only. **Whether every cell inside each ROI is labelled is not yet verified**; the "unlabeled boxes" series suggests it can be checked. The sandbox cannot reach S3, so the download runs on Colab. |
| **LeukemiaAttri / Sparse-LeukemiaAttri** (Rehman et al., MICCAI 2024, arXiv:2405.10803; arXiv:2504.02602) | Peripheral-blood leukemia fields: 2,416 fields × 3 magnifications × 2 microscopes × 2 cameras. 14 WBC classes incl. myeloblast, lymphoblast, monoblast. The 2025 paper states **48 patients** (18 ALL, 22 AML, 2 CLL, 4 CML, 2 APML) with patient-level splits. Fully annotated 100× source fields; the 40× test set is fully annotated, the 40× training set covers about 20% of each field. | Peripheral blood, whole fields, real per-field blast fractions, 48 patients. **The same fields were imaged on cheap and expensive microscopes**, giving a controlled acquisition-shift test for RQ4 that is stronger than RV-PBS. | **Licence not stated. Patient IDs in the release not confirmed.** Access via Google Drive link or the authors (im.itu.edu.pk). Field positions on the slide are probably not recorded, so between-field variation can be measured but spatial adjacency cannot. |
| AML-Cytomorphology_MLL_Helmholtz (TCIA, CC BY 4.0) | 189 patients, 99–500 WBCs scanned per patient | A real per-patient cell population | Per-cell labels for all cells **not confirmed**; no positions. Fallback only. |
| ALL-IDB1 | 108 fields, lymphoblast centroids only | — | Other WBCs unlabelled, so no blast proportion. Not usable. |

## Plan
1. **Now.** Download the BoneMarrowWSI annotations only (`notebooks/21_bmwsi_annotations_colab.ipynb`, about
   6 MB). Then verify locally:
   - that the annotation is exhaustive within the ROIs;
   - the number of labelled cells and blasts per case;
   - the ROI geometry.
2. **If exhaustive.** Pre-register Amendment 4 (RQ4/RQ3, real specimens), using labels and coordinates only; no
   slide pixels are needed:
   - (a) measure blast clustering within and between ROIs: intraclass correlation, design effect, Moran's I;
   - (b) run the Part III certification designs on real per-case cell populations, uniform vs field-unit sampling,
     with and without the H5b cell-level mis-analysis.

   Classifier-driven sampling on this set would need image crops; that is deferred.
3. **In parallel.** Request the LeukemiaAttri data and licence from the authors, and confirm patient IDs. It is the
   better peripheral-blood set and gives a paired microscope shift.

## Sources
- MethodsX scoping review: https://www.sciencedirect.com/science/article/pii/S2215016125005667
- IDC collection: https://portal.imaging.datacommons.cancer.gov/collections/bonemarrowwsi_pediatricleukemia
- Zenodo record: https://zenodo.org/records/17486627 (dataset description, v4)
- LeukemiaAttri: https://arxiv.org/abs/2405.10803 ; https://github.com/intelligentMachines-ITU/Blood-Cancer-Dataset-Lukemia-Attri-MICCAI-2024
- Sparse-LeukemiaAttri: https://arxiv.org/html/2504.02602v2
- AML-Cytomorphology_MLL_Helmholtz: https://www.cancerimagingarchive.net/collection/aml-cytomorphology_mll_helmholtz/
