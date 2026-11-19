# Certified Adaptive Examination of Blood Smears

Research code and paper: *Why Information-Gain Stopping Cannot Certify Rare Evidence, and What Can Be Guaranteed Instead.*
Author: Qudrat Ullah · Supervisor: Balázs Harangi · University of Debrecen

This is a methodological and empirical study. It is **not** a clinical tool, and nothing here establishes
diagnostic validity.

## Research questions
- **RQ1** Does information-saturation stopping fail systematically as decision-changing evidence becomes rare? (Part I, kill-tests 1–2)
- **RQ2** Why does certifying the absence of rare evidence require near-exhaustive examination in the worst case? (Part II)
- **RQ3** How many fields must be examined to certify a threshold decision under adaptive sampling? (Part III)
- **RQ4** Where do the assumptions of anytime-valid audit methods break on real smear images? (Part III, real data)

## Reproducing the paper
See **`REPRODUCE.md`**. In short:
```
pip install -r requirements.txt
python reproduce.py check --shangrla <folder with SHANGRLA's NonnegMean.py>   # tests, Part II proofs, reference cross-check
python reproduce.py list                                                    # every step, which paper item it supports, its cost
python reproduce.py verify <step id> ...                                    # re-run a step and compare with the stored results
```
`REPRODUCE.md` also lists the data sources (with licences), the derived files and their checksums
(`data/CHECKSUMS.txt`), the Colab notebooks that need a GPU, and which scripts are frozen by the pre-registration.

## Tests
`pytest -q tests` (56 checks, about 60 s) checks by simulation that each test used in the paper keeps its error rate at
the exact null (uniform, adversarial and adaptive samplers; the control variate, the uniform mixture and the exploratory
hedged bets; a sampler that is not positive on every unsampled field; equal and unequal fields; imperfect verifier with a
margin), that the four invalid analyses the paper warns against do not, the Part II results, and Proposition 4 of the
JBHI article (`tests/test_limit.py`; proof in `docs/part3_limit_proof.md`). The last report is in
`results/validity_tests_v1/report.json`. `python scripts/54_cv_validity_checks.py` (about 5 min) runs the same kind of
checks for the exploratory control-variate test and the stratified sampler (`results/cv_validity_v1/report.json`).

## Layout
```
paper/        jbhi_compact.tex (10-page IEEE JBHI version) and jbhi_supplement.tex (its supplementary material), with
              audit_jbhi_numbers.py (checks the numbers of both against results/); main.tex (longer version, October
              2026, without the November analyses); refs.bib
src/cadpbs/   package code: samplers, stopping rules (InfoDPP-PAC re-implementation), ALPHA, weighted betting,
              hedged betting (hedged_betting.py, exploratory), control variate and score strata
              (cv_hedged_betting.py, exploratory)
scripts/      numbered analysis steps (01_ ... 59_); each writes results/<run>/ with config.json
notebooks/    Google Colab notebooks for the GPU and download steps (03, 15, 21, 35, 37, 42, 57)
tests/        validity test suite (pytest)
verify/       numerical checks of the power analysis and the calibration budget
sim/          early synthetic pilot
results/      one folder per run: config.json and small summaries are tracked
figures/      paper figures (scripts/25_figures.py)
data/         datasets, not tracked (see data/README.md and data/CHECKSUMS.txt)
docs/         pre-registration, audit, deviations log, reviews, notes
reproduce.py  one entry point for checks, runs and comparisons
```

## Reproducibility rules
- Every run writes to its own `results/<run_name>/`, with `config.json` (including the seed).
- `docs/PREREG.md` and its amendments were frozen with git tags (`prereg-v1` to `prereg-v4`) before the confirmatory
  runs; `docs/PREREG_AUDIT.md` checks the runs against them. Frozen scripts and the modules they import are not edited.
- Every departure from the pre-registration or from InfoDPP-PAC (arXiv:2608.23574) is logged in `docs/DEVIATIONS.md`.

## Licence and citation
The code is released under the MIT licence (`LICENSE`); `CITATION.cff` gives the citation. The manuscript sources in
`paper/` are not covered by the licence: copyright in the article rests with the authors and, on acceptance, the
publisher. The datasets are not redistributed here; their sources and licences are listed in `REPRODUCE.md`.
