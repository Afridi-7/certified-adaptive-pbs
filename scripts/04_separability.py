"""
Step 4 (RQ1 crux): how distinct are blasts in the selector's embedding space?

The synthetic pilot showed the sign of the InfoDPP-PAC-vs-uniform gap is governed by the separability of
critical evidence in the embedding used for selection, not by rarity alone. This script measures it, at
the cell level and at the field level (a field = mean embedding of m cells), and asks the question that
actually drives a diversity-based selector: ARE BLASTS OUTLIERS? (log-det/DPP selection prefers points far
from what is already selected, so it helps only if blasts are far from the bulk.)

All classifier scores are cross-fitted over the patient-grouped folds, so no patient is scored by a model
that saw it. Reported with bootstrap CIs clustered on SOURCE PATIENT (11 blast patients), because episodes
that reuse the same patient's blasts are not independent.

Outputs: results/separability_v1/{summary.json,by_class.csv,by_patient.csv,field_level.csv,quality.csv}
"""
import json, os, argparse
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score

ap = argparse.ArgumentParser()
ap.add_argument("--emb", default="results/embeddings_v1")
ap.add_argument("--out", default="results/separability_v1")
ap.add_argument("--seed", type=int, default=20260919)
ap.add_argument("--m_list", type=int, nargs="+", default=[1, 5, 10, 20])  # cells per field
ap.add_argument("--boot", type=int, default=400)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
rng = np.random.default_rng(a.seed)

E = np.load(f"{a.emb}/embeddings_dinov2_vits14.npy").astype(np.float32)
df = pd.read_csv(f"{a.emb}/cells_features.csv", dtype={"patient_id": str})
assert len(df) == len(E)
y = df.is_blast.values; folds = df.fold.values; pid = df.patient_id.values
En = E / np.linalg.norm(E, axis=1, keepdims=True)          # cosine geometry
summary = {"n_cells": int(len(df)), "n_blasts": int(y.sum()), "dim": int(E.shape[1])}

_GIDX = None
def cluster_boot_ci(stat_fn, groups, B=a.boot, alpha=0.05):
    """Bootstrap resampling whole source patients (clusters), not individual cells."""
    global _GIDX
    if _GIDX is None:
        _GIDX = {g: np.flatnonzero(groups == g) for g in np.unique(groups)}
    gs = np.array(list(_GIDX)); out = []
    for _ in range(B):
        take = rng.choice(gs, len(gs), replace=True)
        idx = np.concatenate([_GIDX[g] for g in take])
        try: out.append(stat_fn(idx))
        except Exception: pass
    return float(np.quantile(out, alpha / 2)), float(np.quantile(out, 1 - alpha / 2))

# ---------- 1. cross-fitted discriminability (supervised upper bound on "can a classifier see it") -------
p = df.probe_p_oof.values
auc = roc_auc_score(y, p)
lo, hi = cluster_boot_ci(lambda i: roc_auc_score(y[i], p[i]) if y[i].sum() and (~y[i].astype(bool)).sum() else np.nan, pid)
summary["cellwise_auc_blast_vs_rest"] = dict(point=float(auc), ci95=[lo, hi], note="cross-fitted logistic probe, patient-clustered bootstrap")

# per class: blast vs each other class separately (the confusable classes matter for the selector)
rows = []
for c in sorted(df.cell_type.unique()):
    if c == "blast": continue
    m = (y == 1) | (df.cell_type.values == c)
    rows.append(dict(other_class=c, n_other=int((df.cell_type.values == c).sum()),
                     auc_blast_vs_class=float(roc_auc_score(y[m], p[m]))))
by_class = pd.DataFrame(rows).sort_values("auc_blast_vs_class")
by_class.to_csv(f"{a.out}/by_class.csv", index=False)

# per source patient: recall at a threshold with 95% specificity on non-blasts
thr = float(np.quantile(p[y == 0], 0.95))
bp = df[y == 1].assign(p=p[y == 1]).groupby("patient_id").agg(blasts=("p", "size"), recall=("p", lambda s: float((s >= thr).mean())))
bp.to_csv(f"{a.out}/by_patient.csv")
summary["threshold_at_95pct_specificity"] = thr
summary["per_patient_recall"] = dict(median=float(bp.recall.median()), min=float(bp.recall.min()), max=float(bp.recall.max()))

# ---------- 2. UNSUPERVISED outlierness: what a diversity-based selector actually responds to -----------
# distance from each cell to its k-th nearest neighbour among a random background pool (non-blast cells)
def knn_dist(Q, ref, k=10, chunk=2000):
    d = np.empty(len(Q))
    for s in range(0, len(Q), chunk):
        sim = Q[s:s+chunk] @ ref.T
        d[s:s+chunk] = 1 - np.partition(sim, -k, axis=1)[:, -k]
    return d
bg = rng.choice(np.flatnonzero(y == 0), 4000, replace=False)
dk = knn_dist(En, En[bg], k=10)
auc_out = roc_auc_score(y, dk)   # 0.5 = blasts no more outlying than normal cells
lo, hi = cluster_boot_ci(lambda i: roc_auc_score(y[i], dk[i]) if y[i].sum() and (~y[i].astype(bool)).sum() else np.nan, pid)
summary["outlierness_auc_cell"] = dict(point=float(auc_out), ci95=[lo, hi],
    note="AUC of 10-NN distance to a background pool for blast vs rest; >0.5 means blasts sit further from the bulk, which HELPS a log-det/DPP selector find them")
summary["knn_dist_medians"] = dict(blast=float(np.median(dk[y == 1])), other=float(np.median(dk[y == 0])))

# ---------- 3. field level: a field is the mean embedding of m cells; does one blast still stand out? ---
rows = []
pool_neg = np.flatnonzero(y == 0); pool_pos = np.flatnonzero(y == 1)
for m in a.m_list:
    nf = 3000
    neg_idx = rng.choice(pool_neg, (nf, m))
    F_neg = En[neg_idx].mean(1)
    pos_cell = rng.choice(pool_pos, nf)
    F_pos = np.concatenate([En[pos_cell][:, None, :], En[rng.choice(pool_neg, (nf, m - 1))]], 1).mean(1) if m > 1 else En[pos_cell]
    F_neg /= np.linalg.norm(F_neg, axis=1, keepdims=True); F_pos /= np.linalg.norm(F_pos, axis=1, keepdims=True)
    ref = F_neg[:2000]
    d_neg, d_pos = knn_dist(F_neg[2000:], ref, k=10), knn_dist(F_pos, ref, k=10)
    lab = np.r_[np.zeros(len(d_neg)), np.ones(len(d_pos))]
    rows.append(dict(m_cells_per_field=m, outlierness_auc_field=float(roc_auc_score(lab, np.r_[d_neg, d_pos])),
                     median_knn_dist_neg=float(np.median(d_neg)), median_knn_dist_pos=float(np.median(d_pos))))
pd.DataFrame(rows).to_csv(f"{a.out}/field_level.csv", index=False)
summary["field_level"] = rows

# ---------- 4. RQ4 preliminary: does probe error track image quality / morphology? ---------------------
qcols = [c for c in df.columns if c.startswith("q_")]
err = (p >= thr).astype(int) ^ y                      # cell-level error at the fixed operating point
qc = pd.DataFrame({q: [float(np.corrcoef(df[q].values, err)[0, 1])] for q in qcols}, index=["pearson_r_with_error"]).T
qc["blast_mean"] = [float(df.loc[y == 1, q].mean()) for q in qcols]
qc["other_mean"] = [float(df.loc[y == 0, q].mean()) for q in qcols]
qc.to_csv(f"{a.out}/quality.csv")
# error rate by source patient (heterogeneity of classifier error across patients = non-exchangeability signal)
summary["blast_error_rate_by_patient"] = {k: float(v) for k, v in (1 - bp.recall).items()}

json.dump(summary, open(f"{a.out}/summary.json", "w"), indent=1)
print(json.dumps(summary, indent=1))
print("\nblast vs each class (lowest = most confusable):\n", by_class.head(6).to_string(index=False))
print("\nper source patient:\n", bp.to_string())
print("\nfield level:\n", pd.DataFrame(rows).to_string(index=False))
print("\nquality:\n", qc.round(4).to_string())
