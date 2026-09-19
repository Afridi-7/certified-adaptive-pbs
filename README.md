# Certified Adaptive Examination of Peripheral Blood Smears

MSc thesis code: *When Information-Gain Stopping Fails, and What Can Be Guaranteed Instead.*
Author: Qudrat Ullah · Supervisor: Balázs Harangi

This is a methodological/empirical study. It is **not** a clinical tool.

## Research questions
- **RQ1** Does information-saturation stopping fail systematically as decision-changing evidence becomes rare? (Part I, kill-tests 1–2)
- **RQ2** Why does certifying the absence of rare evidence require near-exhaustive examination in the worst case? (Part II)
- **RQ3** How many fields must be examined to certify a threshold decision under adaptive sampling? (Part III)
- **RQ4** Where do the assumptions of anytime-valid audit methods break on real smear images? (Part III, real data)

## Layout
```
data/        datasets (not tracked; see data/README.md)
papers/      reference PDFs (not tracked)
scripts/     numbered, runnable steps (01_, 02_, ...)
sim/         simulation experiments (pilot.py = synthetic Kill-test 2 pilot)
src/cadpbs/  reusable package code (samplers, stopping rules, certifiers)
verify/      numerical checks of derivations (Part II, power analysis)
results/     one folder per run: config.json + summary.csv are tracked
docs/        runbook, status notes, pre-registration, deviations log
```

## Setup (Windows / PowerShell)
```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121   # or /cpu
pip install -r requirements.txt
```

## Reproducibility rules
- Every run writes to a new `results/<run_name>/`, with `config.json` (including seed) and `summary.csv`.
- `docs/PREREG.md` is frozen with a git tag (`prereg-v1`) **before** the kill-tests are run.
- Any deviation from InfoDPP-PAC (arXiv:2608.23574) is logged in `docs/DEVIATIONS.md`.
