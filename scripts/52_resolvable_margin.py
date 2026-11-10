"""
Resolvable margin of a fixed-size count (review of 2026-11-10; exact calculation, no simulation, no data). RQ3.

Question. If n cells are drawn uniformly without replacement from the M cells of a specimen and classified without
error, for which true blast proportions p can a valid test decide "p > theta" or "p < theta" reliably?

Procedure (the `fixed_exact` procedure of scripts/32_decision_level.py): two exact one-sided hypergeometric tests at
alpha/2 each, evaluated once at n.
  above: reject H0: B <= floor(theta M) when X >= c_hi, c_hi = min{x : P(X >= x | B = floor(theta M)) <= alpha/2}
  below: reject H0: B >= ceil(theta M)  when X <= c_lo, c_lo = max{x : P(X <= x | B = ceil(theta M))  <= alpha/2}
X ~ Hypergeometric(M, B, n) is the number of blasts among the n cells; B is the number of blasts in the specimen.

Output. For each (theta, M, n) and power level: p_up = the smallest p = B/M > theta at which "above" is declared with
probability >= power, and p_lo = the largest p < theta at which "below" is declared with probability >= power. For
p_lo < p < p_up the test returns "undetermined" (or the correct decision with probability below `power`); this interval
is the UNDETERMINED ZONE of an n-cell count. M = inf uses the binomial distribution (p on a grid of 0.0005).

Assumptions: uniform sampling of individual cells, perfect classification of the sampled cells, the sampling frame is
the M cells. With field sampling the zone is wider by the design effect; with verifier error it shifts.
The numbers are exact given these assumptions (no Monte Carlo), so no seed is used.

Usage: python scripts/52_resolvable_margin.py --out results/resolvable_margin_v1
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy.stats import hypergeom, binom

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/resolvable_margin_v1"); ap.add_argument("--alpha", type=float, default=0.05)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); h = a.alpha / 2

THETAS = [0.10, 0.20, 0.01, 0.02]          # ICC 2022 (10%), WHO 5th ed. (20%); 1% and 2%: DIPSS / MIPSS70 cut-offs
SIZES = [500, 1000, 2000, 4000, None]      # None: binomial limit
BUDGETS = [100, 200, 500, 1000, 2000]
POWERS = [0.8, 0.9]


def crit(theta, M, n):
    """Critical values (c_hi, c_lo) of the two exact tests; c_lo = -1 if "below" can never be declared."""
    if M is None:
        sf = lambda x: binom.sf(x - 1, n, theta); cdf = lambda x: binom.cdf(x, n, theta)
    else:
        sf = lambda x: hypergeom.sf(x - 1, M, int(np.floor(theta * M + 1e-9)), n)
        cdf = lambda x: hypergeom.cdf(x, M, int(np.ceil(theta * M - 1e-9)), n)
    c_hi = next(x for x in range(n + 2) if sf(x) <= h)
    lo = [x for x in range(0, n + 1) if cdf(x) <= h]
    return c_hi, (max(lo) if lo else -1)


rows = []
for th in THETAS:
    for M in SIZES:
        for n in BUDGETS:
            if M is not None and n > M:
                continue
            c_hi, c_lo = crit(th, M, n)
            if M is None:
                ps = np.round(np.arange(0.0, 1.0 + 1e-9, 0.0005), 4)
                p_above = binom.sf(c_hi - 1, n, ps); p_below = binom.cdf(c_lo, n, ps) if c_lo >= 0 else np.zeros_like(ps)
            else:
                B = np.arange(0, M + 1); ps = B / M
                p_above = hypergeom.sf(c_hi - 1, M, B, n)
                p_below = hypergeom.cdf(c_lo, M, B, n) if c_lo >= 0 else np.zeros(M + 1)
            for pw in POWERS:
                up = ps[(ps > th + 1e-12) & (p_above >= pw)]; lo = ps[(ps < th - 1e-12) & (p_below >= pw)]
                rows.append(dict(theta=th, M="inf" if M is None else M, n=n, power=pw, c_hi=int(c_hi), c_lo=int(c_lo),
                                 p_lo=None if len(lo) == 0 else float(lo.max()),
                                 p_up=None if len(up) == 0 else float(up.min())))
D = pd.DataFrame(rows); D.to_csv(f"{a.out}/summary.csv", index=False)
json.dump(dict(alpha=a.alpha, thetas=THETAS, sizes=["inf" if M is None else M for M in SIZES], budgets=BUDGETS, powers=POWERS),
          open(f"{a.out}/config.json", "w"), indent=1)
pd.set_option("display.width", 200, "display.max_rows", 500)
print(D[D.power == 0.8].to_string(index=False))
