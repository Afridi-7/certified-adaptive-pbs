"""
Three-way decisions (above / below / undetermined) at fixed budgets (RQ3; exploratory).

On the same simple random sample of cells (M = 2,000):
  sequential    two one-sided ALPHA tests at alpha/2, eta0 = theta +/- 0.05
  seq_default   the same with alpha.py's default eta0
  fixed_exact   one-sided hypergeometric tests at the fixed n
  conventional  observed fraction vs theta, no error control

Usage: python scripts/32_decision_level.py --out results/decision_level_v1 --seed 20261025
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import hypergeom
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.alpha import alpha_martingale

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/decision_level_v1"); ap.add_argument("--seed", type=int, default=20261025)
ap.add_argument("--M", type=int, default=2000); ap.add_argument("--episodes", type=int, default=2000)
ap.add_argument("--budgets", type=int, nargs="+", default=[100, 200, 500]); ap.add_argument("--alpha", type=float, default=0.05)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
M, h = a.M, a.alpha / 2
GRID = {0.10: [0.05, 0.08, 0.10, 0.12, 0.15], 0.20: [0.15, 0.18, 0.20, 0.22, 0.25]}
nmax = max(a.budgets); rows = []
for th, plist in GRID.items():
    Bnull_hi = int(np.floor(th * M)); Bnull_lo = int(np.ceil(th * M))
    # critical values of the exact tests at each budget (least favourable nulls)
    crit_above = {n: next(x for x in range(n + 2) if hypergeom.sf(x - 1, M, Bnull_hi, n) <= h) for n in a.budgets}
    crit_below = {n: max([x for x in range(-1, n + 1) if x < 0 or hypergeom.cdf(x, M, Bnull_lo, n) <= h]) for n in a.budgets}
    for p in plist:
        B = int(round(p * M)); truth = "above" if B > th * M else ("below" if B < th * M else "none")
        rng = np.random.default_rng([a.seed, int(th * 100), int(round(p * 1000))])
        pop = np.r_[np.ones(B), np.zeros(M - B)]
        for ep in range(a.episodes):
            y = pop[rng.permutation(M)][:nmax]
            # sequential: run both tests over the first nmax draws (null means use the full population size M)
            T = {}
            for tag, e_a, e_b in (("tuned", (th + 0.05) / (2 * th), (1 - th + 0.05) / (2 * (1 - th))), ("default", None, None)):
                la, _ = alpha_martingale(y / (2 * th), 1 / (2 * th), 0.5, M, alpha=h, eta0=e_a)
                lb, _ = alpha_martingale((1 - y) / (2 * (1 - th)), 1 / (2 * (1 - th)), 0.5, M, alpha=h, eta0=e_b)
                ta = np.flatnonzero(la >= np.log(1 / h)); tb = np.flatnonzero(lb >= np.log(1 / h))
                T[tag] = (ta[0] + 1 if len(ta) else 10 ** 9, tb[0] + 1 if len(tb) else 10 ** 9)
            cs = np.cumsum(y)
            for n in a.budgets:
                dec = lambda ta, tb: "above" if ta <= n and ta <= tb else ("below" if tb <= n else "undetermined")
                seq = dec(*T["tuned"]); seqd = dec(*T["default"])
                x = int(cs[n - 1])
                fx = "above" if x >= crit_above[n] else ("below" if x <= crit_below[n] else "undetermined")
                cv = "above" if x / n > th else "below"
                rows.append(dict(theta=th, p=p, truth=truth, n=n, ep=ep, sequential=seq, seq_default=seqd, fixed_exact=fx, conventional=cv))
D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)
summ = []
for (th, p, n), g in D.groupby(["theta", "p", "n"]):
    truth = g.truth.iloc[0]
    for proc in ["sequential", "seq_default", "fixed_exact", "conventional"]:
        dec = g[proc]
        correct = float((dec == truth).mean()) if truth != "none" else 0.0
        wrong = float(((dec != truth) & (dec != "undetermined")).mean()) if truth != "none" else float((dec != "undetermined").mean())
        summ.append(dict(theta=th, p=p, n=n, procedure=proc, correct=correct, wrong=wrong, undetermined=float((dec == "undetermined").mean())))
S = pd.DataFrame(summ); S.to_csv(f"{a.out}/summary.csv", index=False)
json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
pd.set_option("display.width", 200); print(S.round(3).to_string(index=False))
