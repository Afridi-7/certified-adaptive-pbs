"""
Part II (RQ2): numerical verification of the rare-event certification results in docs/part2_proof.md.
Every check is exact (rational arithmetic or an LP), except check 6, which is a seeded Monte Carlo.
Usage: python scripts/17_part2_verify.py --out results/part2_v1
"""
import argparse, itertools, json, os
from fractions import Fraction
from math import comb, ceil
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/part2_v1")
ap.add_argument("--seed", type=int, default=20260701)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
log = {}

def miss(N, k, n):            # P(all k critical fields unexamined) under SRSWOR of size n, exact
    return Fraction(comb(N - k, n), comb(N, n))

def n_star(N, k, alpha):      # smallest n with miss <= alpha
    alpha = Fraction(alpha).limit_denominator(10**6)
    for n in range(N + 1):
        if miss(N, k, n) <= alpha: return n

# 1. k = 1: miss = 1 - n/N exactly; n* = ceil((1-alpha) N)
ok1 = all(miss(N, 1, n) == 1 - Fraction(n, N) for N in range(1, 301) for n in range(N + 1))
ok1b = all(n_star(N, 1, al) == ceil((1 - Fraction(al).limit_denominator(10**6)) * N)
           for N in range(1, 301) for al in (0.05, 0.01, 0.10))
print(f"[1] miss(N,1,n) == 1 - n/N for all N<=300: {ok1};  n* == ceil((1-alpha)N): {ok1b}")
log["check1_k1_exact"] = ok1 and ok1b

# 2. product form and sandwich bounds  [max(0, 1 - n/(N-k+1))]^k <= miss <= (1-n/N)^k, over the FULL range 0 <= n <= N.
#    (Corrected 2026-10-22: the earlier check covered only n <= N-k, and the paper stated the lower bound without the
#    positive part, which fails for even k and n > N-k+1, e.g. N = 10, k = 2, n = 10: 1/81 > 0 = miss.)
ok2 = True
for N in range(1, 121):
    for k in range(1, N + 1):
        for n in range(0, N + 1):
            m = miss(N, k, n)
            prod = np.prod([Fraction(max(N - n - i, 0), N - i) for i in range(k)])
            lo = max(Fraction(0), 1 - Fraction(n, N - k + 1)) ** k; hi = (1 - Fraction(n, N)) ** k
            ok2 &= (m == prod) and (lo <= m <= hi)
print(f"[2] product form and sandwich bounds hold for all N<=120 and all 0<=n<=N: {ok2}")
log["check2_bounds"] = bool(ok2)
bad = (1 - Fraction(10, 9)) ** 2 > miss(10, 2, 10)                  # the bound without the positive part fails here
print(f"[2b] without the positive part the lower bound fails at N=10, k=2, n=10: {bad}")
log["check2b_unclipped_bound_fails_at_endpoint"] = bool(bad)

# 3. monotone in k (worst case of unknown k >= 1 is k = 1)
ok3 = all(miss(N, k + 1, n) <= miss(N, k, n) for N in range(2, 121) for k in range(1, N) for n in range(N + 1))
print(f"[3] miss non-increasing in k for all N<=120: {ok3}")
log["check3_monotone_k"] = ok3

# 4. table of n* and the sufficient closed form ceil(N(1 - alpha^(1/k)))
rows = []
for N in (200, 500, 1000):
    for al in (0.05, 0.01):
        for k in (1, 2, 3, 5, 10, 20, 50):
            ns = n_star(N, k, al); approx = ceil(N * (1 - al ** (1 / k)))
            rows.append(dict(N=N, alpha=al, k=k, n_star=ns, frac=round(ns / N, 3), closed_form=approx,
                             closed_form_minus_exact=approx - ns))
import csv
with open(f"{a.out}/nstar_table.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print("[4] n* (exact) and closed form ceil(N(1-alpha^(1/k))):")
for r in rows: print("   ", r)
log["check4_closed_form_never_below_exact"] = all(r["closed_form_minus_exact"] >= 0 for r in rows)
print("    closed form is always >= exact (sufficient):", log["check4_closed_form_never_below_exact"])

# 5. minimax value of the search game (no location information) by linear programming, small N
#    Examiner mixes over n-subsets; adversary picks a k-subset C; payoff = 1{C disjoint from subset}.
from scipy.optimize import linprog
ok5 = True; lp_rows = []
for N in range(2, 8):
    for k in range(1, N):
        for n in range(0, N - k + 1):
            S = list(itertools.combinations(range(N), n)); Cs = list(itertools.combinations(range(N), k))
            A = np.array([[1.0 if not set(c) & set(s) else 0.0 for s in S] for c in Cs])  # rows: C, cols: strategy
            # min v  s.t.  A x <= v, sum x = 1, x >= 0
            nS = len(S); c = np.r_[np.zeros(nS), 1.0]
            A_ub = np.c_[A, -np.ones(len(Cs))]; b_ub = np.zeros(len(Cs))
            A_eq = np.r_[np.ones(nS), 0.0][None, :]
            res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[1.0], bounds=[(0, None)] * nS + [(None, None)])
            v = res.fun; exact = float(miss(N, k, n))
            det_worst = A.max(axis=0).min()          # best pure (deterministic) strategy's worst case
            ok5 &= abs(v - exact) < 1e-9 and (det_worst == 1.0 or n > N - k)
            lp_rows.append((N, k, n, round(v, 6), round(exact, 6), det_worst))
print(f"[5] LP minimax value == hypergeometric miss, and every deterministic strategy has worst case 1 (N<=7): {ok5}")
log["check5_minimax_lp"] = bool(ok5)

# 6. Monte Carlo of uniform SRSWOR (sanity, seeded)
rng = np.random.default_rng(a.seed); mc = []
for N, k, n in [(200, 1, 33), (200, 1, 190), (500, 5, 225), (1000, 10, 258)]:
    R = 200_000; hits = 0
    crit = np.arange(k)
    for _ in range(R // 1000):
        samp = np.argsort(rng.random((1000, N)), axis=1)[:, :n]
        hits += (~np.isin(samp, crit).any(axis=1)).sum()
    est = hits / R; ex = float(miss(N, k, n)); se = np.sqrt(ex * (1 - ex) / R)
    mc.append(dict(N=N, k=k, n=n, mc=est, exact=ex, z=(est - ex) / se if se > 0 else 0.0))
    print(f"[6] N={N} k={k} n={n}: MC {est:.4f} vs exact {ex:.4f}  z={(est-ex)/se:+.2f}")
log["check6_mc"] = mc

# 7. CRC-type calibration budget (the proposal's alpha = 0.01 claim):  (n/(n+1)) R + B/(n+1) <= alpha, B = 1
#    => usable budget on empirical risk R: alpha - (1 - alpha)/n, negative iff n < (1 - alpha)/alpha
for al in (0.05, 0.01):
    nmin = min(n for n in range(1, 10**5) if al - (1 - al) / n >= 0)
    print(f"[7] alpha={al}: budget alpha-(1-alpha)/n >= 0 first at n = {nmin}  (negative for n < {nmin})")
    log[f"check7_crc_nmin_alpha{al}"] = nmin

json.dump(dict(log, seed=a.seed), open(f"{a.out}/verify.json", "w"), indent=1, default=float)
print("\nALL EXACT CHECKS PASS:", all(log[k] for k in ["check1_k1_exact", "check2_bounds", "check3_monotone_k",
                                                    "check4_closed_form_never_below_exact", "check5_minimax_lp"]))
