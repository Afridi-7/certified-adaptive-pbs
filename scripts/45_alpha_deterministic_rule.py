"""
Impact of a rule missing from our ALPHA implementation (found by scripts/44_crosscheck_shangrla.py, 2026-10-20).

The reference implementation (SHANGRLA) also rejects as soon as the verified sample already contains more blasts than
the null allows in the whole specimen (sample total > theta * M): the null is then false with certainty. The frozen
src/cadpbs/alpha.py rejects in that situation only when a further positive field is drawn, so it can stop later. This
never affects validity (it is conservative), and it cannot occur at p = theta. Here we measure how much it inflated the
reported ALPHA costs: the Part III / population-size setting (m = 5 cells per field, N = M / 5 fields, uniform field
sampling), 400 specimens per cell, same orders for both rules. The frozen module is not changed.

Usage: python scripts/45_alpha_deterministic_rule.py --out results/alpha_detrule_v1 --seed 20261036
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.alpha import alpha_martingale

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/alpha_detrule_v1"); ap.add_argument("--seed", type=int, default=20261036)
ap.add_argument("--reps", type=int, default=400)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
rows = []
for th in (0.10, 0.20):
    u = 1 / (2 * th)
    for M in (500, 1000, 2000, 4000):
        N = M // 5
        for p in (th + 0.02, th + 0.05, th + 0.10):
            rng = np.random.default_rng([a.seed, int(th * 100), M, int(round(p * 1000))]); nb = int(round(p * M))
            for ep in range(a.reps):
                cells = np.zeros(M); cells[:nb] = 1; rng.shuffle(cells); frac = cells.reshape(N, 5).mean(1)
                X = frac[rng.permutation(N)] / (2 * th)
                _, st = alpha_martingale(X, u, 0.5, N, alpha=0.05)
                k = np.flatnonzero(np.cumsum(X) > 0.5 * N + 1e-12)
                det = int(k[0]) + 1 if len(k) else -1                   # reject at the draw that pushes the total over
                both = min(s for s in (st, det) if s > 0) if (st > 0 or det > 0) else -1
                rows.append(dict(theta=th, M=M, p=round(p, 2), ep=ep, alpha_cells=(st if st > 0 else N) * 5,
                                 with_rule_cells=(both if both > 0 else N) * 5, rule_binds=int(det > 0 and (st < 0 or det < st))))
D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)
S = D.groupby(["theta", "M", "p"]).agg(median_frozen=("alpha_cells", "median"), median_with_rule=("with_rule_cells", "median"),
                                       share_rule_binds=("rule_binds", "mean")).reset_index()
S["median_change"] = S.median_with_rule - S.median_frozen
S.to_csv(f"{a.out}/summary.csv", index=False); json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print(S.to_string(index=False))
