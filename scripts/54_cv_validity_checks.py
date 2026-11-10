"""
Validity checks for src/cadpbs/cv_hedged_betting.py (2026-11-10; exploratory; RQ3, RQ4). Deterministic (fixed seeds).

Kept outside tests/ because it takes several minutes; `pytest -q tests` remains the 50-check suite cited in the papers.
Every check builds specimens at an EXACT null (p = theta) or inside the null and counts false certifications of
"p > theta" by the hedged betting test with a prediction control variate, for uniform sampling, the 50% mixture and
score-stratified sampling. Side information is chosen to hurt: an oracle, an anti-oracle, a "liar" that scores blast
fields near zero and blast-free fields high (large negative control-variate payoffs), noise, a constant, and a noisy
oracle. Also checked: with c_grid = (0,) the test reproduces hedged_betting_w draw for draw, and power with good scores.
Per-run RNG: default_rng([20261110, cell seed, repetition]).

Usage: python scripts/54_cv_validity_checks.py --out results/cv_validity_v1
"""
import argparse, json, math, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from cadpbs.hedged_betting import hedged_betting_w
from cadpbs.weighted_betting_w import sampler_uniform_w, sampler_prop_w
from cadpbs.cv_hedged_betting import hedged_cv_betting_w, sampler_strata_w, sampler_mix_w

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/cv_validity_v1"); ap.add_argument("--reps", type=int, default=1000)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); T0 = time.time()
ALPHA, REPS = 0.05, a.reps
LIMIT = ALPHA + 3 * math.sqrt(ALPHA * (1 - ALPHA) / REPS)          # about 3 standard errors above alpha


def specimen(N, m, p, rng, packed=False, sizes=None):
    """Field values f = B_i / m_i and weights m_i / M for a specimen with exactly round(p M) blasts."""
    m_i = np.full(N, m) if sizes is None else np.asarray(sizes)
    M = int(m_i.sum()); nb = int(round(p * M)); assert abs(nb - p * M) < 1e-9
    B = np.zeros(N)
    if packed:                                   # blasts fill whole fields (strongest clustering)
        left = nb
        for i in rng.permutation(N):
            B[i] = min(left, m_i[i]); left -= B[i]
    else:
        fld = np.repeat(np.arange(N), m_i)
        B = np.bincount(fld[rng.choice(M, nb, replace=False)], minlength=N).astype(float)
    return B / m_i, m_i / M


def side(kind, f, rng):
    if kind == "oracle": return np.clip(f, 0, 1)
    if kind == "anti": return 1 - f
    if kind == "liar": return np.where(f > 0, 1e-4, 0.9 + 0.1 * rng.random(len(f)))
    if kind == "noise": return rng.random(len(f))
    if kind == "const": return np.full(len(f), 0.3)
    if kind == "noisy_oracle": return np.clip(f + 0.3 * rng.standard_normal(len(f)), 0, 1)
    raise ValueError(kind)


def rate(theta, p, kind, design, packed=False, sizes=None, seed=0, reps=REPS, N=100, m=5):
    k = 0
    for r in range(reps):
        rng = np.random.default_rng([20261110, seed, r])
        f, w = specimen(N, m, p, rng, packed, sizes(rng) if callable(sizes) else sizes); S = side(kind, f, rng)
        smp = {"uniform": lambda: sampler_uniform_w(N), "mix50": lambda: sampler_mix_w(sampler_prop_w(w, S), 0.5),
               "strata": lambda: sampler_strata_w(w, S)}[design]()
        k += hedged_cv_betting_w(f, w, S, smp, theta, rng=rng)["stop"] > 0
    return k / reps


ROWS = []
def record(group, rate_, limit, **kw):
    ROWS.append(dict(group=group, false_certification_rate=float(rate_), allowed_max=float(limit), ok=bool(rate_ <= limit), **kw))


# 1. c_grid = (0,) reproduces hedged_betting_w exactly
same = 0
for r in range(50):
    rng = np.random.default_rng([7, r]); f, w = specimen(100, 5, 0.15, rng); S = side("noisy_oracle", f, rng)
    x = hedged_betting_w(f, w, sampler_prop_w(w, S), 0.10, rng=np.random.default_rng(r))
    y = hedged_cv_betting_w(f, w, S, sampler_prop_w(w, S), 0.10, c_grid=(0.0,), rng=np.random.default_rng(r))
    same += int(x["stop"] == y["stop"] and np.array_equal(x["order"], y["order"]))
print(f"reduces to hedged_betting_w with c = 0: {same}/50", flush=True)

# 2. exact null, scattered arrangement: 6 kinds of side information x 3 designs x 2 thresholds
for si, kind in enumerate(["oracle", "anti", "liar", "noise", "const", "noisy_oracle"]):
    for di, design in enumerate(["uniform", "mix50", "strata"]):
        for theta in (0.10, 0.20):
            x = rate(theta, theta, kind, design, seed=100 * si + 10 * di + int(theta * 10))
            record("exact_null_scattered", x, LIMIT, side=kind, design=design, theta=theta, p=theta, reps=REPS)
    print(f"  scattered, {kind}: done ({time.time() - T0:.0f}s)", flush=True)
# 3. exact null, blasts packed into whole fields
for si, kind in enumerate(["oracle", "liar", "noisy_oracle"]):
    for di, design in enumerate(["uniform", "mix50", "strata"]):
        x = rate(0.10, 0.10, kind, design, packed=True, seed=5000 + 10 * si + di)
        record("exact_null_packed", x, LIMIT, side=kind, design=design, theta=0.10, p=0.10, reps=REPS)
# 4. inside the null
for di, design in enumerate(["uniform", "mix50", "strata"]):
    x = rate(0.10, 0.08, "oracle", design, seed=7000 + di, reps=REPS // 2)
    record("inside_null", x, 0.02, side="oracle", design=design, theta=0.10, p=0.08, reps=REPS // 2)
# 5. unequal field sizes (2-8 cells), exact null
def sizes(rng):
    s = rng.integers(2, 9, 100); s[-1] += (-s.sum()) % 10; return s      # M divisible by 10
for di, design in enumerate(["uniform", "strata"]):
    x = rate(0.10, 0.10, "noisy_oracle", design, sizes=sizes, seed=8000 + di)
    record("exact_null_unequal_fields", x, LIMIT, side="noisy_oracle", design=design, theta=0.10, p=0.10, reps=REPS)
# 6. power with good scores (p = 0.20 against theta = 0.10, uniform sampling)
k = 0
for r in range(200):
    rng = np.random.default_rng([99, r]); f, w = specimen(100, 5, 0.20, rng); S = side("noisy_oracle", f, rng)
    k += hedged_cv_betting_w(f, w, S, sampler_uniform_w(100), 0.10, rng=rng)["stop"] > 0

null_rows = [r for r in ROWS if r["group"].startswith("exact_null")]
out = dict(alpha=ALPHA, reps=REPS, allowed_max=LIMIT, reduces_to_hedged=dict(identical=same, of=50),
           n_exact_null_cells=len(null_rows), max_exact_null_rate=max(r["false_certification_rate"] for r in null_rows),
           max_inside_null_rate=max(r["false_certification_rate"] for r in ROWS if r["group"] == "inside_null"),
           power_good_scores=k / 200, all_ok=bool(all(r["ok"] for r in ROWS) and same == 50 and k / 200 >= 0.95), cells=ROWS)
json.dump(out, open(os.path.join(a.out, "report.json"), "w"), indent=1)
print({k_: v for k_, v in out.items() if k_ != "cells"}); print(f"DONE in {time.time() - T0:.0f}s", flush=True)
sys.exit(0 if out["all_ok"] else 1)
