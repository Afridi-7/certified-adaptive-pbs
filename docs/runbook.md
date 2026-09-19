# Thesis runbook (living plan)

Each step lists **RQ**, **Where**, **How**, **What to send back**, and the **Decision gate** it feeds.
Work location: your Windows PC. Keep one folder, e.g. `C:\Users\<you>\Documents\thesis_pbs`, with this layout:

```
thesis_pbs/
  data/ku_optofil/        # dataset.zip, metadata.csv, metadata_with_patient_level_splits.csv
  data/rv_pbs/            # later
  papers/                 # PDFs of all cited papers
  scripts/                # numbered scripts I give you (01_, 02_, ...)
  sim/                    # simulation code (pilot.py)
  results/<run_name>/     # every run writes config.json + csv here; never overwrite
  PREREG.md               # pre-registration (frozen before kill-tests)
```
Rule: every run gets a new `results/<run_name>` folder containing its config. Send me the console output and the `summary.csv`.

---

## Phase 0: Setup (today, ~1 h)
**Where:** PowerShell on your PC.
```powershell
mkdir $HOME\Documents\thesis_pbs; cd $HOME\Documents\thesis_pbs
mkdir data\ku_optofil, papers, scripts, sim, results
python -m venv .venv; .\.venv\Scripts\Activate.ps1
pip install numpy scipy pandas pillow scikit-learn matplotlib torch torchvision timm tqdm
python -c "import torch;print(torch.__version__, torch.cuda.is_available())"
```
Download from https://zenodo.org/records/17333317 into `data\ku_optofil\`: `dataset.zip`, `metadata.csv`, `metadata_with_patient_level_splits.csv`.
Extract it with `Expand-Archive data\ku_optofil\dataset.zip data\ku_optofil\dataset`.
**Send back:** the output of the torch line, plus your GPU model (or "none").

## Phase 1: Verify the data and sources (W1–2), RQ1–4 foundations
1. `python scripts\01_inspect_kuoptofil.py --data data\ku_optofil` → paste the full output.
   We will confirm: the cell count (proposal says 31,489), the 13 class names, whether a blast class exists and how many blasts there are, the number of patients, how the split is defined, crop size, and whether any field or slide IDs or coordinates exist.
2. Download these PDFs into `papers\` and upload them to the project/chat:
   - InfoDPP-PAC (arXiv:2608.23574): **priority**. I need Eqs. 1–8 verbatim to reimplement it.
   - Shekhar et al. 2023 (arXiv:2305.06884)
   - Spertus, Sridhar & Stark (arXiv:2409.06680)
   - Csillag et al. (arXiv:2502.04294)
   - Stark ALPHA (AOAS 2023)
   - the RV-PBS paper (Vergara-Lluri et al. 2024) plus the RV-PBS dataset link
3. **Decision gate A:** Is the blast count per patient large enough for composites at the rarity levels k ∈ {1, 2, 5, 10, 25, 50} per 500 fields, *using test-split patients only*? If not, we reduce N or the rarity grid, and we document the change.

## Phase 2: Pre-registration (W1–2)
I draft `PREREG.md` from the proposal and the pilot (H1a/H1b split, rarity grid, matched-cost definition, McNemar with stated p₁₀/p₀₁, α = 0.05, the named risk-control procedure, stopping rules, exclusion rules). You review it with your supervisor and **freeze it (commit/date it) before running the kill-tests.**

## Phase 3: Patch environment and baselines (W3–5), RQ1
1. Train a cell classifier (blast vs. rest, plus 13-class) **on train-split patients only**. Report per-class recall with CIs on the test patients.
2. Extract embeddings for all cells with (a) the trained classifier and (b) a frozen generic encoder (the "teacher" in InfoDPP-PAC terms; this follows the paper once transcribed).
3. **Measure separability s:** blast vs. each of the other classes in embedding space (AUC, Mahalanobis distance). The pilot showed that the whole RQ1 outcome depends on this number.
4. Composite generator: a field = m cells sampled from test-split patients; N fields per synthetic specimen; k critical fields; the arrangement is scattered or clustered. Exact ground truth is stored.
5. Baselines: exhaustive, raster, uniform, confidence stopping. Every run logs order, selection probabilities, stop step, and cost.

## Phase 4: Kill-test 1 (W6–8), RQ1
The faithful InfoDPP-PAC reimplementation replaces `infosat` in `sim/pilot.py`. Every deviation from the paper goes into `DEVIATIONS.md`.
Measure the three signals per step (information gain, relevance, true critical-evidence indicator). Test whether the stop time depends on k (H1a) and whether IG is correlated with critical content.

## Phase 5: Kill-test 2, the rarity gap (W9–12), RQ1
Run the full grid on real embeddings, paired against matched-cost uniform/raster, with ≥ the pre-registered number of positive episodes. Report miss-rate curves vs. k/N with CIs, McNemar results, and effect sizes. Report the result whichever way it comes out; if H1b fails, map the regime (s × k) instead.

## Phase 6: Part II write-up (parallel, any time), RQ2
Proof of: (i) 1 − n/N and the hypergeometric formula for k; (ii) the minimax lower bound for any strategy that has no location information; (iii) when the bound does *not* apply (informative features, known k ≥ 2). I draft; you check.

## Phase 7: Certified threshold counting (W13–16), RQ3
Estimand: the blast proportion over all cells in the finite specimen vs. a threshold (WHO-HAEM5 20% or ICC 10%, stated per experiment). Adaptive selection with recorded probabilities, simulated expert verification (ground truth) on a subset, and an anytime-valid test (ALPHA / Shekhar weighted without replacement / prediction-powered e-values). Outcome: fields examined until certification at α = 0.05, compared with fixed conventions (100/200/500 cells) at matched risk.

## Phase 8: Real-data assumptions (W17–20), RQ4
On RV-PBS real fields, test: spatial autocorrelation of blasts (Moran's I), classifier error vs. morphology and stain quality, bounded assorters, exchangeability. For each violation: characterise it, test whether validity breaks (does the empirical coverage of the certificate drop below 1 − α?), and propose a weakened guarantee.

## Phase 9: Analysis and writing (W21–26)
