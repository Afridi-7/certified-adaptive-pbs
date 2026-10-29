"""
Clustering as the excess of the pooled ANOVA ICC over its permutation-null mean (RQ4; exploratory).
The pooled ICC is biased upward for binary data in small unequal subfields, so the null mean is recomputed
inside every patient-bootstrap sample. Definitions are imported from 23_bmwsi_clustering.py.

Usage: python scripts/34_icc_excess.py --out results/a4_bmwsi/icc_excess.json --seed 20261027
"""
import argparse, json, os, runpy, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/a4_bmwsi/icc_excess.json"); ap.add_argument("--seed", type=int, default=20261027)
ap.add_argument("--n_perm", type=int, default=500); ap.add_argument("--n_boot", type=int, default=500)
ap.add_argument("--n_perm_boot", type=int, default=50)
a = ap.parse_args()
# load the frozen helper functions without running the frozen analysis (its module body runs analyses on import)
src = open(os.path.join(os.path.dirname(__file__), "23_bmwsi_clustering.py")).read()
head = src.split("\nC = prepare()")[0]
sys.argv = ["23_bmwsi_clustering.py", "--n_perm", "1", "--n_boot", "1"]
ns = {}; exec(compile(head, "23_bmwsi_clustering.py", "exec"), ns)
prepare, icc, scale_ids = ns["prepare"], ns["icc"], ns["scale_ids"]
rng = np.random.default_rng(a.seed)


def perm_within(groups, y):
    out = y.copy()
    for g in groups: out[g] = rng.permutation(y[g])
    return out


C = prepare(); res = {}
for sc in ["S1", "S2"]:
    par, sub = scale_ids(C, sc); y = C.y.values.astype(float); rid = C.rid.values
    groups = [np.flatnonzero(rid == r) for r in np.unique(rid)]
    obs = icc(par, sub, y)[0]; null = np.array([icc(par, sub, perm_within(groups, y))[0] for _ in range(a.n_perm)])
    pats, pidx = np.unique(C.patient.values, return_inverse=True); idx = [np.flatnonzero(pidx == j) for j in range(len(pats))]
    _, pc = np.unique(par, return_inverse=True); _, scd = np.unique(sub, return_inverse=True); big = scd.max() + 1
    boots = []
    for _ in range(a.n_boot):
        pick = rng.integers(0, len(pats), len(pats)); ii = np.concatenate([idx[p] for p in pick])
        tag = np.concatenate([np.full(len(idx[p]), j) for j, p in enumerate(pick)])
        P_, S_, Y_ = pc[ii] + tag * big, scd[ii] + tag * big, y[ii]
        G_ = [np.flatnonzero(P_ == q) for q in np.unique(P_)]
        nb = np.mean([icc(P_, S_, perm_within(G_, Y_))[0] for _ in range(a.n_perm_boot)])
        boots.append(icc(P_, S_, Y_)[0] - nb)
    res[sc] = dict(icc_raw=float(obs), null_mean=float(null.mean()), excess=float(obs - null.mean()),
                   excess_ci95=[float(np.quantile(boots, .025)), float(np.quantile(boots, .975))],
                   perm_p=float((1 + (null >= obs).sum()) / (1 + a.n_perm)))
os.makedirs(os.path.dirname(a.out), exist_ok=True); json.dump(res, open(a.out, "w"), indent=1); print(json.dumps(res, indent=1))
