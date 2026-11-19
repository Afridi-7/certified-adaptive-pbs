"""
Proposition 4 (limit for certifying a threshold decision; RQ3). Review of 2026-11-19. Exact computation, no data.

Statement (proof in docs/part3_limit_proof.md). If what the examiner sees before verification carries no information
about which cells are blasts (labels exchangeable given everything it observes), then for every procedure that
certifies "p > theta" with probability at most alpha for every population with p <= theta -- any adaptive choice of
cells, any stopping rule -- the probability of certifying within n verified cells, for a population with B > theta*M
blasts, is at most the power of the randomised one-sided hypergeometric test of floor(theta*M) against B based on n
uniform draws without replacement. The bound is attained by that fixed-n test. With fields as units, the same holds
when fields (not cells) are exchangeable, for every null field population with p <= theta.

Outputs (results/threshold_limit_v1/summary.json):
  zones       : 80%-power range (alpha/2 = 0.025 per side) that no valid procedure resolves within n cells, M = 2,000,
                n = 100, 200, 500, 1,000, theta = 0.10 and 0.20; next to the non-randomised exact test of scripts/52.
  table2      : the bound at every Table II setting off the threshold ("above" and "below" rows) and the share of it
                reached by the exact fixed-n test and by sequential ALPHA (values of results/decision_level_v1/summary.csv).
  fields      : fields of m = 5 cells holding either no blast or c blasts (N = 400, M = 2,000, theta = 0.10): verified
                cells needed for 80% power with fields against individually sampled cells, and the design effect
                (c - p m) / (1 - p) of that arrangement.
  checks      : Monte Carlo check of the reduction (adaptive policy with uninformative features) and of the bound
                against a valid sequential test, seed 20261119.

Usage: python scripts/59_threshold_limit.py --out results/threshold_limit_v1
"""
import argparse, json, os, platform
from math import floor, ceil
import numpy as np, pandas as pd
from scipy.stats import hypergeom, chisquare

ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--seed", type=int, default=20261119)
ap.add_argument("--reps", type=int, default=20000)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")


def mp_power(M, B1, B0, n, alpha):
    """Power at B1 of the randomised most powerful level-alpha test of B0 (reject for large X ~ Hyp(M, B, n))."""
    xs = np.arange(n + 1); p0 = hypergeom.pmf(xs, M, B0, n); p1 = hypergeom.pmf(xs, M, B1, n)
    tail0 = np.cumsum(p0[::-1])[::-1]
    c = next(x for x in range(n + 2) if x == n + 1 or tail0[x] <= alpha)
    size_strict = tail0[c] if c <= n else 0.0
    gamma = 0.0 if c == 0 or p0[c - 1] == 0 else min(max((alpha - size_strict) / p0[c - 1], 0.0), 1.0)
    return float((p1[c:].sum() if c <= n else 0.0) + (gamma * p1[c - 1] if c >= 1 else 0.0))


def exact_power(M, B1, B0, n, alpha):
    """Non-randomised exact one-sided test (scripts/52): reject iff P0(X >= x) <= alpha."""
    xs = np.arange(n + 1); tail0 = np.cumsum(hypergeom.pmf(xs, M, B0, n)[::-1])[::-1]
    c = next(x for x in range(n + 2) if x == n + 1 or tail0[x] <= alpha)
    return float(hypergeom.sf(c - 1, M, B1, n)) if c <= n else 0.0


def zone(M, th, n, power=0.8, a_side=0.025, f=mp_power):
    B0a, B0b = floor(th * M), ceil(th * M)
    up = next(B for B in range(B0a + 1, M + 1) if f(M, B, B0a, n, a_side) >= power)
    lo = [B for B in range(B0b) if f(M, M - B, M - B0b, n, a_side) >= power]      # "below": mirror on non-blasts
    return (max(lo) / M if lo else None), up / M


out = dict(zones={}, table2={}, fields={}, checks={})
M = 2000
for th in (0.10, 0.20):
    for n in (100, 200, 500, 1000):
        out["zones"][f"theta{th:.2f}_n{n}"] = dict(bound=zone(M, th, n), exact_test=zone(M, th, n, f=exact_power))

dl = pd.read_csv(os.path.join(ROOT, "results/decision_level_v1/summary.csv"))
for th, p, n in [(th_, p_, n_) for th_, ps in ((0.10, (0.08, 0.12, 0.15)), (0.20, (0.18, 0.22, 0.25))) for p_ in ps for n_ in (200, 500)]:
    g = dl[np.isclose(dl.theta, th) & np.isclose(dl.p, p) & (dl.n == n)].set_index("procedure").correct
    if p > th:     # certify "above": null floor(theta M) blasts
        b = mp_power(M, round(p * M), floor(th * M), n, 0.025); e = exact_power(M, round(p * M), floor(th * M), n, 0.025)
    else:          # certify "below": mirror on non-blasts, null ceil(theta M) blasts
        b = mp_power(M, M - round(p * M), M - ceil(th * M), n, 0.025); e = exact_power(M, M - round(p * M), M - ceil(th * M), n, 0.025)
    out["table2"][f"theta{th:.2f}_p{p:.2f}_n{n}"] = dict(bound=b, exact_test_exact=e, exact_test_simulated=float(g["fixed_exact"]),
                                                         sequential_simulated=float(g["sequential"]),
                                                         sequential_share_of_bound=float(g["sequential"]) / b)

m, N = 5, 400; th = 0.10
for c in (5, 3, 2, 1):
    for p in (0.15, 0.20):
        K1 = round(p * M / c)
        if K1 > N or K1 * c != round(p * M): continue
        K0 = floor(th * M / c)
        nf = next(k for k in range(1, N + 1) if mp_power(N, K1, K0, k, 0.025) >= 0.8)
        nc = next(k for k in range(1, M + 1) if mp_power(M, round(p * M), floor(th * M), k, 0.025) >= 0.8)
        out["fields"][f"c{c}_p{p:.2f}"] = dict(fields_needed=nf, cells_with_fields=nf * m, cells_individually=nc,
                                               ratio=nf * m / nc, design_effect=(c - p * m) / (1 - p), null_rich_fields=K0)

# Monte Carlo checks -----------------------------------------------------------------------------------------------------
rng = np.random.default_rng(a.seed)
Mr, Br, nr = 30, 9, 8; cnt = np.zeros(nr + 1)
for _ in range(a.reps):
    lab = np.zeros(Mr, int); lab[rng.choice(Mr, Br, replace=False)] = 1; feat = rng.normal(size=Mr)
    avail = list(range(Mr)); seen = []; obs = []
    for t in range(nr):
        j = (min(avail, key=lambda i: abs(feat[i] - feat[seen[-1]])) if obs and obs[-1] == 1 else max(avail, key=lambda i: feat[i]))
        avail.remove(j); seen.append(j); obs.append(lab[j])
    cnt[sum(obs)] += 1
ex = hypergeom.pmf(np.arange(nr + 1), Mr, Br, nr) * a.reps; k = ex > 5
out["checks"]["reduction_chisq_p"] = float(chisquare(cnt[k], ex[k] * cnt[k].sum() / ex[k].sum()).pvalue)


def seq_test(seq, Mp, B0, alpha):   # valid anytime betting test of H0: B <= B0 under uniform sampling without replacement
    W, x = 1.0, 0
    for t, s in enumerate(seq):
        mu = (B0 - x) / (Mp - t)
        if mu < 0: return t
        W *= 1 + (0.5 / mu if mu > 0 else 0.0) * (s - mu); x += s
        if W >= 1 / alpha: return t + 1
    return 0


gaps = []
for (Mp, B0, B1, n) in [(40, 8, 16, 12), (60, 12, 24, 20), (100, 10, 25, 30), (200, 20, 40, 50)]:
    pop1 = np.r_[np.ones(B1), np.zeros(Mp - B1)]; pop0 = np.r_[np.ones(B0), np.zeros(Mp - B0)]
    cert = np.mean([0 < seq_test(rng.permutation(pop1)[:n], Mp, B0, 0.05) <= n for _ in range(a.reps)])
    size = np.mean([0 < seq_test(rng.permutation(pop0), Mp, B0, 0.05) for _ in range(a.reps)])
    bnd = mp_power(Mp, B1, B0, n, 0.05)
    out["checks"][f"seq_M{Mp}_B0{B0}_B1{B1}_n{n}"] = dict(size=float(size), certified_by_n=float(cert), bound=bnd)
    gaps.append(cert - bnd)
out["checks"]["max_rate_minus_bound"] = float(max(gaps))
import scipy
json.dump(dict(config=vars(a), env=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)),
          open(f"{a.out}/config.json", "w"), indent=1)
json.dump(out, open(f"{a.out}/summary.json", "w"), indent=1)
print(json.dumps(out, indent=1))
