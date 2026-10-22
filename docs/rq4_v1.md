# RQ4 — where do the assumptions break? Exploratory results v1 (2026-09-26)

Scripts: `scripts/10_rq4_error_structure.py`, `scripts/11_rq4_certify_by_error_group.py`.
Results: `results/rq4_error_v1/`, `results/rq4_certify_v1/`. Data: KU-Optofil, cross-fitted classifier scores
(each cell scored by a model that never saw its patient). Operating point: 95% specificity on non-blast cells.
Status: exploratory; not covered by `PREREG.md` (Parts I–II). To be pre-registered in amendment 2 with Part III.

## 1. Classifier error is not exchangeable
| Question | Result |
|---|---|
| Blast recall across the 11 source patients | 0.385–1.000; homogeneity χ² = 61.1, df = 10, Monte-Carlo p = 0.0002 |
| False-positive rate across 237 blast-free patients (≥ 50 cells) | pooled 4.9%, range 0–35%, 90th percentile 12.2%; variance **6.6×** that of a common binomial rate |
| False-positive rate by true cell class | metamyelocyte 33%, myelocyte 26%, reactive lymphocyte 20%, basophil 19%, lymphocyte 17%, monocyte 14% |
| Error vs image quality **within class** (log-odds per SD, patient-bootstrap 95% CI) | contrast +0.83 [0.46, 1.31]; brightness +0.78 [0.31, 1.03]; blur (Laplacian var.) +0.21 [0.09, 0.30]; blue channel −0.98 [−1.35, −0.40] |

**Reading.** Errors cluster by patient (both misses and false alarms), by morphology (the immature granulocyte
continuum and lymphoid cells — exactly the cells that sit near a blast-percentage decision) and by stain/imaging
characteristics even after adjusting for class. A model of classifier error as identically distributed across
cells is contradicted on every axis tested.

**Limits of this reading.** The colour channels and brightness are collinear, so individual channel coefficients are
unstable; the defensible statement is that *stain colour composition and contrast* are associated with error, not
that any single channel is. Stain characteristics are largely patient-level (one smear per patient), so quality and
patient effects cannot be fully separated with these data; the patient bootstrap accounts for the clustering but does
not identify the cause. These are crops from one laboratory; the cross-source shift (RV-PBS, PBC-DIB) is not yet run.

## 2. Does that break the certificate? No — it changes the cost
Specimens of 500 cells, θ = 0.10 (ICC 2022), blasts drawn only from the **low-recall** patients (4 patients, 166 blasts,
recall < 0.90) or only from the **high-recall** ones (7 patients, 323 blasts); 300 episodes per cell.

| Blast source | false certification at p = θ (uniform / classifier-driven) | median cells, p = 0.12 | p = 0.15 | p = 0.20 |
|---|---|---|---|---|
| high-recall patients | 0.003 / **0.000** | 425 → **205** (−52%) | 260 → **60** (−77%) | 105 → **35** (−67%) |
| low-recall patients | 0.007 / **0.000** | 420 → **235** (−44%) | 258 → **140** (−46%) | 110 → **48** (−57%) |

**Validity holds in both groups** (false certification ≤ 0.007), as the theory says it must: the importance-weighted
test does not rely on the classifier being right, only on the sampling probabilities being known and positive, and on
labels coming from verification. **Efficiency does depend on the classifier**: on the patients where it is worst,
certification at p = 0.15 needs 140 cells instead of 60, i.e. 2.3× more.

## 3. What this answers for RQ4
The assumptions separate into two kinds, and the data show that they fail differently:

| Assumption | Needed for | On these data |
|---|---|---|
| Known, positive selection probabilities; labels from verification | **validity** | holds by design (enforced by the sampler); the uncorrected design shows what happens when it does not: 99.7% false certification |
| Bounded, non-negative assorter | **validity** | holds by construction (a = b/2θ ∈ [0, 1/2θ]) |
| Exchangeable / identically distributed classifier error | **efficiency** of classifier-driven designs; **validity** only of methods that impute unverified units from predictions | **violated** (patient, morphology, stain) |
| Spatial independence of fields | validity of uniform-over-the-smear designs if the sampling frame were not random; efficiency otherwise | **not testable** with available data (RV-PBS annotation is class-partial) |

The guarantee therefore survives exactly the violation the data exhibit, because it rests on the sampling design rather
than on the model. Methods that would use predictions *in place of* verification for unexamined cells would not have this
protection; that is the boundary this study can draw with evidence.

## Next (RQ4)
- Cross-source shift: embed RV-PBS crops (752 annotated instances, 153 blasts) and PBC-DIB crops on Colab, apply the
  KU-Optofil-trained probe, and measure how recall and false-positive rates move under a change of laboratory, stain
  and camera.
- Pre-register these as part of amendment 2.
