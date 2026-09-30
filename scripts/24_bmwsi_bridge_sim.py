"""
PREREG Amendment 4, A4.2 (RQ4 -> RQ3): what does the clustering MEASURED in A4.1 imply for certification validity?
MODEL-BASED extrapolation, reported as such. Also places the E3.1 stress arrangements on the same ICC scale.

Specimen model: N = 30 fields x m = 17 cells (m = median nucleated cells per S1 quadrant is ~17; M = 510).
Field blast propensities w_f ~ Beta(a, b) with mean theta and ICC rho = 1/(a+b+1); counts k_f ~ Binomial(m, w_f), then
blasts are added/removed at random (proportional to w_f / 1 - w_f) until the total equals round(theta*M) exactly, so every
specimen sits exactly on the null boundary p = theta. rho in {0 (binomial), A4.1 S1 lower CI, estimate, upper CI} (floored at 0).
Designs per specimen: alpha_field (ALPHA on field means, uniform field order) and naive_cell (same order, cells flattened,
cell-level ALPHA). theta in {0.10 (ICC 2022), 0.20 (WHO 5th ed.)}; 1,000 episodes per (theta, rho); seed 20260711.
The realised ICC of the simulated specimens is re-estimated with the A4.1 estimator and reported.
E3.1 placement: ICC of E3.1 specimens (scattered / hot0.6 / hot1.0; N = 100 x m = 5) at p = theta, 200 each, seed 20260712.
"""
import argparse, json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.alpha import assorter_bound, alpha_martingale

ap = argparse.ArgumentParser()
ap.add_argument("--clustering", default="results/a4_bmwsi/clustering.json")
ap.add_argument("--out", default="results/a4_bmwsi/bridge.json"); ap.add_argument("--episodes", type=int, default=1000)
a = ap.parse_args()
src = open(os.path.join(os.path.dirname(__file__), "23_bmwsi_clustering.py")).read()
ns = {}; exec(src[src.index("def icc"):src.index("def scale_ids")], {"np": np}, ns); icc = ns["icc"]

cl = json.load(open(a.clustering))["primary"]["S1"]
RHOS = sorted(set(round(max(0.0, r), 4) for r in [0.0, cl["ci95"][0], cl["icc"], cl["ci95"][1]]))
N, m = 30, 17; M = N * m


def specimen(rng, theta, rho):
    B = int(round(theta * M))
    if rho <= 0: w = np.full(N, theta)
    else:
        s = 1 / rho - 1; w = rng.beta(theta * s, (1 - theta) * s, N)
    k = rng.binomial(m, w)
    while k.sum() != B:
        if k.sum() < B:
            cand = np.flatnonzero(k < m); pr = w[cand] + 1e-9; f = rng.choice(cand, p=pr / pr.sum()); k[f] += 1
        else:
            cand = np.flatnonzero(k > 0); pr = 1 - w[cand] + 1e-9; f = rng.choice(cand, p=pr / pr.sum()); k[f] -= 1
    cells = np.zeros((N, m), int)
    for f in range(N): cells[f, rng.choice(m, k[f], replace=False)] = 1
    return cells


res = dict(rhos=RHOS, sims={})
for theta, sd in [(0.10, 0), (0.20, 1)]:
    u = assorter_bound(theta)
    for rho in RHOS:
        rng = np.random.default_rng([20260711, sd, int(round(rho * 10000))])
        fa = nc = 0; realised = []
        for ep in range(a.episodes):
            cells = specimen(rng, theta, rho)
            order = rng.permutation(N)
            _, st = alpha_martingale(cells.mean(1)[order] / (2 * theta), u, 0.5, N, alpha=0.05); fa += st > 0
            _, sc = alpha_martingale(cells[order].reshape(-1) / (2 * theta), u, 0.5, M, alpha=0.05); nc += sc > 0
            if ep < 200:
                par = np.zeros(M, int); sub = np.repeat(np.arange(N), m); realised.append(icc(par, sub, cells.reshape(-1))[0])
        res["sims"][f"theta{theta}_rho{rho}"] = dict(theta=theta, rho=rho, realised_icc=float(np.mean(realised)),
            alpha_field_false_cert=fa / a.episodes, naive_cell_false_cert=nc / a.episodes, n=a.episodes)
        print(res["sims"][f"theta{theta}_rho{rho}"], flush=True)

# E3.1 placement on the ICC scale (same generator as scripts/18, labels only)
rng = np.random.default_rng(20260712); place = {}
for theta in (0.10, 0.20):
    for arr, h in [("scattered", None), ("hot0.6", 0.6), ("hot1.0", 1.0)]:
        v = []
        for _ in range(200):
            Nn, mm = 100, 5; MM = Nn * mm; nb = int(round(theta * MM)); y = np.zeros(MM, int)
            if h is None: y[rng.choice(MM, nb, replace=False)] = 1
            else:
                nhot = min(Nn, int(np.ceil(nb / (mm * h)))); hot = rng.choice(Nn, nhot, replace=False)
                y[rng.choice((hot[:, None] * mm + np.arange(mm)).ravel(), nb, replace=False)] = 1
            v.append(icc(np.zeros(MM, int), np.repeat(np.arange(Nn), mm), y)[0])
        place[f"theta{theta}_{arr}"] = float(np.mean(v))
res["e31_icc_placement"] = place
os.makedirs(os.path.dirname(a.out), exist_ok=True); json.dump(res, open(a.out, "w"), indent=1)
print(json.dumps(place, indent=1))
