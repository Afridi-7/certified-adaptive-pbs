"""
Kill-test 1 re-analysis with a FIELD-LEVEL teacher (P0-2 of the 2026-10-19 review; exploratory, not pre-registered).

Problem fixed. In 05_killtests.py the selector sees a field only through the MEAN of its cell embeddings, but the
"teacher relevance" is bank.p[cells].max(1): the maximum of the cross-fitted CELL-level blast scores. That relevance
signal presupposes a detector that finds and classifies every cell, which the selector's representation does not
assume (PREREG section 1). It is also a deviation from InfoDPP-PAC, whose teacher T(x_i) acts on the patch embedding.

What this script does, for every evaluation fold f and dilution m:
  1. Trains a field-level teacher T_f: logistic regression on POOLED field embeddings (same pooling as the selector),
     using composite fields built ONLY from the other four folds (patient-grouped, so no patient is seen in training
     and evaluation). Positive training field = 1 blast + (m - 1) blast-free cells; negative = m blast-free cells.
  2. Builds evaluation specimens with the unchanged generator (make_specimen, N = 200 fields).
  3. Records, per episode, the AUC for critical vs non-critical fields of
       - the information gain before selection (dinfo0; independent of the teacher),
       - the original max-of-cell-scores relevance (cell-level oracle),
       - the field-level teacher (same information as the selector).
  4. Runs InfoDPP-PAC with the FIELD-level teacher (tau recalibrated with it on blast-free specimens of other folds)
     and records the miss at k critical fields for InfoDPP-PAC, uniform at the same cost, and (exploratory) a
     relevance-ranked scan by the field-level teacher at the same cost.
  5. Stores the source patient of every critical blast, so that patient-level bootstrap intervals are possible.

Usage: python scripts/26_kt1_field_teacher.py --out results/kt1_field_teacher --seed 20261019
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.composites import CellBank, make_specimen
from cadpbs.infodpp import InfoDPPConfig, run_infodpp, calibrate_tau
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score

ap = argparse.ArgumentParser()
ap.add_argument("--emb", default="results/embeddings_v1")
ap.add_argument("--out", default="results/kt1_field_teacher")
ap.add_argument("--N", type=int, default=200)
ap.add_argument("--k_list", type=int, nargs="+", default=[1, 2, 5, 10, 20])
ap.add_argument("--m_list", type=int, nargs="+", default=[1, 5, 10, 20])
ap.add_argument("--folds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
ap.add_argument("--episodes", type=int, default=100, help="per (fold, m, k)")
ap.add_argument("--n_train", type=int, default=3000, help="positive and negative training fields each, per (fold, m)")
ap.add_argument("--target_frac", type=float, default=0.163)
ap.add_argument("--n_calib", type=int, default=8)
ap.add_argument("--seed", type=int, default=20261019)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

E = np.load(f"{a.emb}/embeddings_dinov2_vits14.npy").astype(np.float32)
df = pd.read_csv(f"{a.emb}/cells_features.csv", dtype={"patient_id": str})
bank = CellBank(E, df)
base_cfg = InfoDPPConfig(tau=None, seed=a.seed)
folds_all = np.unique(bank.fold)


def pool_fields(cells):                       # identical to make_specimen(pooling="mean")
    F = bank.X[cells].mean(1)
    return F / np.linalg.norm(F, axis=1, keepdims=True)


def train_field_teacher(f, m, rng):
    other = [g for g in folds_all if g != f]
    bl = np.concatenate([bank.blast_idx[g] for g in other]); bg = np.concatenate([bank.bg_idx[g] for g in other])
    neg = rng.choice(bg, (a.n_train, m), replace=True)
    pos = rng.choice(bg, (a.n_train, m), replace=True); pos[:, 0] = rng.choice(bl, a.n_train, replace=True)
    Xtr = np.vstack([pool_fields(pos), pool_fields(neg)]); ytr = np.r_[np.ones(a.n_train), np.zeros(a.n_train)]
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=3000))
    clf.fit(Xtr, ytr)
    return clf


def calibrate(f, m, clf, rng):
    cal = []
    for g in [x for x in folds_all if x != f][:3]:
        for _ in range(max(1, a.n_calib // 3)):
            s = make_specimen(bank, g, a.N, 0, m, rng)
            cal.append((s["F"], clf.predict_proba(s["F"])[:, 1]))
    return calibrate_tau(cal, base_cfg, a.target_frac, rng=rng)


rows_all = []
for f in a.folds:
    for m in a.m_list:
        path = f"{a.out}/episodes_fold{f}_m{m}.csv"
        if os.path.exists(path):
            rows_all.append(pd.read_csv(path)); continue
        rng = np.random.default_rng([a.seed, f, m]); t0 = time.time()
        clf = train_field_teacher(f, m, rng)
        tau = calibrate(f, m, clf, rng); cfg = InfoDPPConfig(**{**base_cfg.__dict__, "tau": tau})
        rows = []
        for k in a.k_list:
            for ep in range(a.episodes):
                s = make_specimen(bank, f, a.N, k, m, rng)
                crit = s["crit"]; ci = np.flatnonzero(crit)
                T_field = clf.predict_proba(s["F"])[:, 1]
                res = run_infodpp(s["F"], T_field, cfg, np.random.default_rng(rng.integers(1 << 30)))
                n = res["k_star"]
                uni = rng.permutation(a.N)
                ranked = np.argsort(-T_field, kind="stable")
                def miss(order):
                    return int(not np.isin(ci, order[:n]).any())
                blast_cells = s["cells"][crit][bank.is_blast[s["cells"][crit]]]
                rows.append(dict(
                    fold=f, m=m, k=k, ep=ep, tau=tau, infodpp_n=n,
                    auc_info=float(roc_auc_score(crit, res["dinfo0"])),
                    auc_rel_cellmax=float(roc_auc_score(crit, s["teacher"])),
                    auc_rel_field=float(roc_auc_score(crit, T_field)),
                    infodpp_fieldteacher_miss=miss(res["order"]),
                    uniform_miss=miss(uni),
                    ranked_field_miss=miss(ranked),
                    blast_patients=";".join(sorted(set(bank.pid[blast_cells])))))
        d = pd.DataFrame(rows); d.to_csv(path, index=False); rows_all.append(d)
        print(f"fold {f} m {m}: tau={tau:.3f} {len(d)} episodes in {time.time()-t0:.0f}s", flush=True)

all_df = pd.concat(rows_all, ignore_index=True); all_df.to_csv(f"{a.out}/episodes.csv", index=False)
json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print("total episodes:", len(all_df))
