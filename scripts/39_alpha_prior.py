"""
Sensitivity of ALPHA's cost to the prior guess eta0 (RQ3; exploratory).

Compares the default eta0 = (0.5 + u) / 2 used in the reported runs with guesses p = theta + 0.05 and
theta + 0.10, and with uniform betting, on the same sampling orders. Any eta0 fixed in advance is valid;
only power changes.

Usage: python scripts/39_alpha_prior.py --out results/alpha_prior_v1 --seed 20261030
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import binomtest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/alpha_prior_v1"); ap.add_argument("--seed", type=int, default=20261030)
ap.add_argument("--n_null", type=int, default=1000); ap.add_argument("--n_alt", type=int, default=400)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
df = pd.read_csv("results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})
y = df.is_blast.values.astype(bool); score = df.probe_p_oof.values
BL = np.flatnonzero(y); BG = np.flatnonzero(~y & ~df.patient_id.isin(set(df.loc[y, "patient_id"])).values)
M, m = 500, 5; N = M // m; pi = np.full(N, 1 / N); rows = []
for th in (0.10, 0.20):
    u = 1 / (2 * th)
    priors = dict(alpha_default=None, alpha_m05=(th + 0.05) / (2 * th), alpha_m10=(th + 0.10) / (2 * th))
    for p_true in (th, th + 0.02, th + 0.05, th + 0.10):
        rng = np.random.default_rng([a.seed, int(th * 100), int(round(p_true * 1000))])
        nb = int(round(p_true * M))
        for ep in range(a.n_null if np.isclose(p_true, th) else a.n_alt):
            cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(BG, M - nb, replace=False)]
            rng.shuffle(cells); cells = cells.reshape(N, m); frac = y[cells].mean(1); S = score[cells].mean(1)
            order = rng.permutation(N); r = dict(theta=th, p=round(p_true, 2), ep=ep)
            for nm, e0 in priors.items():
                _, st = alpha_martingale(frac[order] / (2 * th), u, 0.5, N, alpha=0.05, eta0=e0)
                r[f"{nm}_cert"] = int(st > 0); r[f"{nm}_cells"] = (st if st > 0 else N) * m
            o = betting_test_w(frac, pi, S, sampler_uniform_w(N), th, rng=rng)
            r["bet_uniform_cert"] = int(o["stop"] > 0); r["bet_uniform_cells"] = (o["stop"] if o["stop"] > 0 else N) * m
            rows.append(r)
D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)
des = ["alpha_default", "alpha_m05", "alpha_m10", "bet_uniform"]; out = []
for (th, p), g in D.groupby(["theta", "p"]):
    r = dict(theta=th, p=p, n=len(g))
    for d in des:
        if np.isclose(p, th):
            k = int(g[f"{d}_cert"].sum()); r[f"{d}_false_cert"] = k / len(g)
            r[f"{d}_p_excess"] = float(binomtest(k, len(g), 0.05, alternative="greater").pvalue)
        else:
            r[f"{d}_median_cells"] = float(g[f"{d}_cells"].median())
    out.append(r)
S_ = pd.DataFrame(out); S_.to_csv(f"{a.out}/summary.csv", index=False); json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print(S_.to_string())
