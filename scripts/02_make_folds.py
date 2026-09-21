"""
Step 2: build the clean cell table and the patient-grouped cross-fitting folds (RQ1, RQ3 foundations).

Why folds instead of the provided split: the provided patient-level test split contains only 3 blast-positive
patients (86 blasts) and the validation split has none. Critical-evidence experiments need every blast patient
to be evaluable by a model that never saw that patient. So we use K patient-grouped folds: each cell is embedded
and scored by the model trained on the other K-1 folds.

Exclusions (logged):
  - patient_id == 'unknown' (3,204 cells, filename prefixes EO_/VLY_/BL_/...; no patient grouping possible,
    so they cannot be placed in folds without leakage risk; source unverified)
  - files on disk not listed in metadata (4 '... - Copy.jpg' exact duplicates)
The folder train/val/test split is NOT patient-level (265/276 patients span folders) and is never used.

Output: data/ku_optofil/cells.csv  (image_name, rel_path, patient_id, cell_type, is_blast, fold, provided_split)
        rel_path is relative to the DATASET ROOT (the directory that directly contains train/ val/ test/),
        e.g. "train/Blast/1604468629_006.jpg". The root itself is detected and logged, so the table stays
        valid whether the zip was extracted to <data>/dataset or <data>/dataset/dataset.
        results/folds_v1/config.json, summary.csv
"""
import os, json, argparse
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--data", default="data/ku_optofil")
ap.add_argument("--k", type=int, default=5)
ap.add_argument("--seed", type=int, default=20260919)
ap.add_argument("--out", default="results/folds_v1")
a = ap.parse_args()
rng = np.random.default_rng(a.seed)

m = pd.read_csv(os.path.join(a.data, "metadata_with_patient_level_splits.csv"), dtype={"patient_id": str})
# locate the dataset root: the directory that directly contains the class folders
def find_root(base):
    cands = [d for d, subs, _ in os.walk(base) if {"train", "test", "val"} <= set(subs)]
    if not cands: raise SystemExit(f"could not find a directory containing train/ val/ test/ under {base}")
    return min(cands, key=len)
ROOT = find_root(os.path.join(a.data, "dataset"))
print("dataset root:", ROOT)

# map image_name -> path relative to ROOT
disk = {}
for r, _, fs in os.walk(ROOT):
    for f in fs:
        if f.endswith(".jpg"):
            disk[f] = os.path.relpath(os.path.join(r, f), ROOT).replace("\\", "/")
m["rel_path"] = m.image_name.map(disk)
assert m.rel_path.notna().all(), "metadata file missing on disk"
n0 = len(m)
excl_unknown = int((m.patient_id == "unknown").sum())
m = m[m.patient_id != "unknown"].copy()
m["is_blast"] = (m.cell_type == "blast").astype(int)

# patient-level fold assignment: blast patients dealt round-robin (largest first), others shuffled round-robin
pt = m.groupby("patient_id").agg(n=("is_blast", "size"), nb=("is_blast", "sum")).reset_index()
pos = pt[pt.nb > 0].sort_values("nb", ascending=False).patient_id.tolist()
neg = pt[pt.nb == 0].patient_id.tolist(); rng.shuffle(neg)
fold = {}
for i, p in enumerate(pos): fold[p] = i % a.k
off = len(pos) % a.k
for i, p in enumerate(neg): fold[p] = (i + off) % a.k
m["fold"] = m.patient_id.map(fold)
m = m.rename(columns={"split": "provided_split"})
cols = ["image_name", "rel_path", "patient_id", "cell_type", "is_blast", "fold", "provided_split"]
m[cols].to_csv(os.path.join(a.data, "cells.csv"), index=False)

os.makedirs(a.out, exist_ok=True)
s = m.groupby("fold").agg(patients=("patient_id", "nunique"), cells=("is_blast", "size"),
                          blasts=("is_blast", "sum"),
                          blast_patients=("patient_id", lambda x: int(m.loc[x.index].groupby("patient_id").is_blast.sum().gt(0).sum())))
s.to_csv(os.path.join(a.out, "summary.csv"))
bp = m[m.is_blast == 1].groupby(["fold", "patient_id"]).size().rename("blasts").reset_index()
bp.to_csv(os.path.join(a.out, "blast_patients.csv"), index=False)
cfg = dict(vars(a), dataset_root=os.path.relpath(ROOT, a.data).replace("\\", "/"), n_metadata=n0, excluded_unknown=excl_unknown, n_cells=len(m), n_patients=int(m.patient_id.nunique()),
           n_blast_patients=len(pos), note="unknown-patient cells excluded; 4 '- Copy' duplicate files ignored (not in metadata)")
json.dump(cfg, open(os.path.join(a.out, "config.json"), "w"), indent=1)
print(json.dumps(cfg, indent=1)); print(s.to_string()); print(bp.to_string(index=False))
