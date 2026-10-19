# Pre-registration audit trail (compiled 2026-10-19)

**Purpose.** The paper says each confirmatory analysis was frozen before its data existed. This file is the evidence for
that claim, so that a reviewer does not have to rely on dates typed into `PREREG.md`. Those typed dates are known to be
wrong in several places (see "Typed dates" below). They are left unedited, so that the frozen text still matches its tags.

## Evidence used, from strongest to weakest
1. **GitHub push times** (recorded in `.git/logs/refs/remotes/origin/main`). GitHub holds its own server-side record of
   each push, so these can be checked independently of this machine.
2. **Git commit and tag times.** These are self-reported by the committing machine.
3. **File modification times of result files.** Results are git-ignored (`.gitignore`: `results/**`), so these are the only
   record of when outputs were written. They are weak evidence, because copying or touching a file can change them.

## Freeze versus first output (all times UTC)
| Tag | Commit | Tag type | Commit time | Pushed to GitHub | First confirmatory output (file, mtime) | Output after push by |
|---|---|---|---|---|---|---|
| `prereg-v1` | b3f114d | annotated | 2026-09-25 10:14:49 | 10:14:53 | none (v1.1 amended v1 before any run) | n/a |
| `prereg-v1.1` | db4463d | annotated | 2026-09-25 10:51:15 | 10:51:18 | `confirm_v1/episodes_fold0_m1.csv`, 10:52:01 | 43 s |
| `prereg-v2` | fc4f1e8 | annotated | 2026-09-25 11:36:40 | 11:36:43 | `part3_theta010/episodes_p0.05.csv`, 11:55:08; `external_v1/config.json`, 11:59:27 | 18 min |
| `prereg-v3` | f014352 | **lightweight** | 2026-09-25 12:39:56 | 12:39:59 | `a3_spatial_theta010/ep_scattered_p0.08_e0000.csv`, 12:41:27 | 1.5 min |
| `prereg-v4` | 40c838d | **lightweight** | 2026-09-30 10:22:14 | 10:22:17 | `a4_bmwsi/clustering.json`, 10:26:58 | 4.7 min |

**Frozen code unchanged.** `git diff <tag> HEAD` is empty for every frozen script and for `src/`. Checked for
`09_confirmatory_analysis.py` and `05_killtests.py` (v1.1); `12_certify_confirm.py`, `13_part3_analysis.py` and
`16_external_analysis.py` (v2); `18_spatial_confirm.py`, `19_mixture_confirm.py` and `20_amend3_analysis.py` (v3); and
`23_bmwsi_clustering.py` and `24_bmwsi_bridge_sim.py` (v4).

**Reproducibility check (2026-10-19).** Re-running `05_killtests.py --seed 20260401` for fold 0, m = 1 in a separate
environment (numpy 2.5) reproduced the stored 1,000 episodes. Miss indicators and relevance AUCs were identical. The
calibrated τ differed in the sixth significant digit, which moved the stop step by one field in 2 of the 1,000 episodes and
moved the information-gain AUC by at most 0.0005.

## Content hashes of the confirmatory outputs (SHA-256, first 16 hex digits; computed 2026-10-19)
These fix the content from today onward. They do not prove when the files were created.

| File | SHA-256 (prefix) |
|---|---|
| results/confirm_v1/episodes.csv | 2d3b4466241afd0b |
| results/confirm_v1/decisions.json | 33784ab05dcfae80 |
| results/part3_theta010/episodes.csv | e179712ff582e263 |
| results/part3_theta020/episodes.csv | 7acd719349813174 |
| results/part3_decisions.json | 412d65de93fcbe12 |
| results/external_v1/decisions.json | 1cc6960e47fec202 |
| results/external_v1/h4c_episodes.csv | 544b11f022705d52 |
| results/a3_decisions.json | a5bf9f633470a7e9 |
| results/a4_bmwsi/clustering.json | 68da9986521a8afd (matches confirmatory_v4.md) |
| results/a4_bmwsi/bridge.json | 1690a59356d0d7c3 (matches confirmatory_v4.md) |

## Typed dates in PREREG.md that are wrong
| Location in PREREG.md | Typed | Correct (from git, UTC+2 local) |
|---|---|---|
| Status line (line 3) | "DRAFT — to be frozen" | frozen at `prereg-v1`, 2026-09-25 12:14 |
| v1 freeze line | 2026-09-26 | 2026-09-25 12:14 |
| Amendment 1 header and freeze line | 2026-09-26 | 2026-09-25 12:51 |
| Amendment 2 header | 2026-09-27 | 2026-09-25 13:36 |
| Amendment 2 freeze line | 2026-09-26 | 2026-09-25 13:36 |
| Amendment 3, Amendment 4 | correct | — |

## What "confirmatory" means here (stated in the paper)
- **Fresh draws, same data.** For Parts I and III and for E3.1, "confirmatory" means fresh Monte Carlo draws, on a disjoint
  seed family, over the same fixed pool of 28,280 cells and cross-fitted scores that the exploratory runs used. This
  protects against analytic flexibility (choice of test, margin, decision rule). It does not protect against effects
  specific to this dataset.
- **Data not analysed before.** Only two hypothesis groups were confirmed on such data:
  - H4a–H4c: RV-PBS was first embedded after `prereg-v2`.
  - H7a: the BoneMarrowWSI spatial analysis came after `prereg-v4`, which inspected only counts.
- **Implementation checks.** Several "hypotheses" check the implementation rather than test an empirical claim:
  - H1a follows from the absence of a relevance term in Eq. 6 plus a τ calibrated on blast-free specimens.
  - H1d tests a combinatorial identity through the random-number generator.
  - H3a, H4c, H5a and H6a check validity that a theorem already guarantees.

  The paper marks these as type V in its hypothesis table (Table 2).

## Remaining weaknesses, and what to do
1. `prereg-v3` and `prereg-v4` are **lightweight tags**, which carry no tagger date and can be moved silently. The commit and
   push times above are the evidence. Do not re-create the tags; doing so would look like tampering.
2. **Results are not under version control.** Archive the repository and the confirmatory output folders on Zenodo at
   submission, and cite the DOI. Include this file and the hashes above.
3. The GitHub repository is private. Make it public, or give reviewers access, at submission, so that the push history can
   be checked.
