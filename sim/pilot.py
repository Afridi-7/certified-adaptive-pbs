"""
Pilot simulation for Kill-test 2 (RQ1) on abstract synthetic fields (no images).

The stopper here is a generic GP information-gain saturation rule, not InfoDPP-PAC (see src/cadpbs/infodpp.py).

Specimen: N fields on a grid, k critical. Features x ~ N(0, I), critical x ~ N(s u, I); scattered or clustered.
Policies: exhaustive, raster, uniform, confidence-based, info-saturation.
Metrics: miss (no critical field examined), recall, fields examined; uniform and raster at matched cost (paired).
"""
import argparse, json, time
import numpy as np
from math import comb
from scipy.stats import binomtest


def make_specimen(rng, R, C, k, s, d, arrangement, q=0.0):
    """q: fraction of NON-critical fields that are 'distractors' (outliers of equal magnitude s in random
    directions orthogonal to u, e.g. artefacts / rare benign morphologies competing for the diversity budget)."""
    N = R * C
    coords = np.array([(r, c) for r in range(R) for c in range(C)], float)
    if k == 0:
        crit = np.zeros(N, bool)
    elif arrangement == "scattered":
        crit = np.zeros(N, bool); crit[rng.choice(N, k, replace=False)] = True
    elif arrangement == "clustered":
        centre = coords[rng.integers(N)]
        dist = np.linalg.norm(coords - centre, axis=1) + 1e-6 * rng.random(N)
        crit = np.zeros(N, bool); crit[np.argsort(dist)[:k]] = True
    else:
        raise ValueError(arrangement)
    u = np.zeros(d); u[0] = 1.0
    X = rng.standard_normal((N, d)) + s * crit[:, None] * u
    if q > 0:
        dist = (~crit) & (rng.random(N) < q)
        dirs = rng.standard_normal((dist.sum(), d)); dirs[:, 0] = 0
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
        X[dist] += s * dirs
    return X, crit


def classifier_probs(rng, X, a=2.0, b=None, s=None):
    # logistic score on the discriminative direction; threshold at s/2 (Bayes boundary for equal variance)
    z = a * (X[:, 0] - (s / 2 if s else 1.0))
    return 1 / (1 + np.exp(-z))


def infosat_order(X, ell, sn2, eps, w, max_steps):
    """Greedy max-posterior-variance selection (RBF GP), stop on w consecutive gains < eps."""
    N = X.shape[0]
    sq = (X ** 2).sum(1)
    K = np.exp(-(sq[:, None] + sq[None, :] - 2 * X @ X.T) / (2 * ell ** 2))
    var = np.diag(K).copy()
    L = np.zeros((N, max_steps))
    order, gains, low = [], [], 0
    chosen = np.zeros(N, bool)
    for t in range(max_steps):
        v = np.where(chosen, -np.inf, var)
        j = int(np.argmax(v))
        g = 0.5 * np.log1p(var[j] / sn2)
        order.append(j); gains.append(g); chosen[j] = True
        low = low + 1 if g < eps else 0
        if low >= w:
            break
        c = (K[:, j] - L[:, :t] @ L[j, :t]) / np.sqrt(var[j] + sn2)
        L[:, t] = c
        var = np.maximum(var - c ** 2, 0.0)
    return np.array(order), np.array(gains)


def median_ell(X):
    D = np.sqrt(np.maximum(((X[:, None] - X[None]) ** 2).sum(-1), 0)); return float(np.median(D[D > 0]))


def calibrate_eps(cfg, target_frac, n_spec=20, seed=12345):
    """Choose eps on NEGATIVE specimens (k=0) so median stop ~= target_frac*N (the paper's reported
    operating point: 83.7% fewer patches -> ~16% of budget). Calibration specimens are disjoint from eval."""
    rng = np.random.default_rng(seed); N = cfg["R"] * cfg["C"]
    Xs = [make_specimen(rng, cfg["R"], cfg["C"], 0, 0.0, cfg["d"], "scattered")[0] for _ in range(n_spec)]
    lo, hi = 1e-4, 2.0
    for _ in range(25):
        mid = np.sqrt(lo * hi)
        med = np.median([len(infosat_order(X, median_ell(X), cfg["sn2"], mid, cfg["w"], N)[0]) for X in Xs])
        if med > target_frac * N: lo = mid
        else: hi = mid
    return float(np.sqrt(lo * hi))


def confidence_stop(order, p, tau, m):
    run = 0
    for t, i in enumerate(order):
        if p[i] >= tau:
            return t + 1
        run = run + 1 if p[i] <= 1 - tau else 0
        if run >= m:
            return t + 1
    return len(order)


def run(cfg):
    rng = np.random.default_rng(cfg["seed"])
    R, C = cfg["R"], cfg["C"]; N = R * C
    rows = []
    for k in cfg["k_list"]:
        for arr in cfg["arrangements"]:
            for s in cfg["s_list"]:
              for q in cfg.get("q_list", [0.0]):
                for ep in range(cfg["episodes"]):
                    X, crit = make_specimen(rng, R, C, k, s, cfg["d"], arr, q)
                    p = classifier_probs(rng, X, s=s)
                    ci = np.flatnonzero(crit)
                    uni = rng.permutation(N)  # uniform order, shared across uniform-based rules (paired)
                    raster = np.arange(N)
                    o_info, _ = infosat_order(X, median_ell(X), cfg["sn2"], cfg["eps"], cfg["w"], N)
                    n_info = len(o_info)
                    n_conf = confidence_stop(uni, p, cfg["tau"], cfg["m"])
                    def ev(order, n):
                        seen = np.isin(ci, order[:n])
                        return int(k > 0 and not seen.any()), (seen.mean() if k else np.nan)
                    rec = dict(k=k, arrangement=arr, s=s, q=q, ep=ep, N=N)
                    for name, (order, n) in {
                        "infosat": (o_info, n_info),
                        "uniform@info": (uni, n_info),
                        "raster@info": (raster, n_info),
                        "confidence": (uni, n_conf),
                        "uniform@conf": (rng.permutation(N), n_conf),
                    }.items():
                        miss, recall = ev(order, n)
                        rec[f"{name}_n"] = int(n); rec[f"{name}_miss"] = miss; rec[f"{name}_recall"] = recall
                    rec["exact_uniform_miss@info"] = comb(N - k, n_info) / comb(N, n_info) if k else np.nan
                    rows.append(rec)
    return rows


def wilson(x, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    ph = x / n; den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den; h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return (c - h, c + h)


def summarise(rows):
    import pandas as pd
    df = pd.DataFrame(rows)
    out = []
    for (k, arr, s, q), g in df[df.k > 0].groupby(["k", "arrangement", "s", "q"]):
        r = dict(k=k, arrangement=arr, s=s, q=q, episodes=len(g), n_info_median=g["infosat_n"].median(),
                 n_conf_median=g["confidence_n"].median())
        for name in ["infosat", "uniform@info", "raster@info", "confidence", "uniform@conf"]:
            x = g[f"{name}_miss"].sum(); lo, hi = wilson(x, len(g))
            r[f"{name}_miss"] = f"{x/len(g):.3f} [{lo:.3f},{hi:.3f}]"
        # paired McNemar (exact) infosat vs uniform at matched cost
        b = int(((g["infosat_miss"] == 1) & (g["uniform@info_miss"] == 0)).sum())
        c = int(((g["infosat_miss"] == 0) & (g["uniform@info_miss"] == 1)).sum())
        r["mcnemar_b_c"] = f"{b}/{c}"; r["mcnemar_p"] = binomtest(b, b + c, 0.5).pvalue if b + c else 1.0
        r["exact_uniform_miss@info_mean"] = g["exact_uniform_miss@info"].mean()
        out.append(r)
    return pd.DataFrame(out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--out", default="results/pilot")
    a = ap.parse_args()
    cfg = dict(seed=20260919, R=20, C=25, d=4, k_list=[1, 2, 5, 10, 25, 50],
               arrangements=["scattered", "clustered"], s_list=[0.0, 1.5, 3.0], episodes=200,
               sn2=0.01, eps=None, target_frac=0.163, w=5, tau=0.9, m=30)
    if a.config: cfg.update(json.load(open(a.config)))
    import os; os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    if cfg.get("eps") is None:
        cfg["eps"] = calibrate_eps(cfg, cfg["target_frac"]); print("calibrated eps =", cfg["eps"])
    rows = run(cfg)
    json.dump(cfg, open(f"{a.out}/config.json", "w"), indent=1)
    import pandas as pd
    pd.DataFrame(rows).to_csv(f"{a.out}/episodes.csv", index=False)
    s = summarise(rows); s.to_csv(f"{a.out}/summary.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(s.to_string(index=False)); print(f"{time.time()-t0:.0f}s")
