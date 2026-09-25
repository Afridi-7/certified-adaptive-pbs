"""
PREREG Amendment 3, E3.1 (RQ4, RQ3): certification when blasts are CLUSTERED across fields, and what happens when
field-sampled cells are analysed as if each cell had been sampled individually.

Specimen: N = 100 fields x m = 5 cells = M = 500 cells, as in Part III. B = round(p*M) blasts are drawn without
replacement from all 489 KU-Optofil blasts; background cells without replacement from blast-free patients.
Scores are cross-fitted (out-of-patient), as in Amendment 2.

Arrangements (fixed in Amendment 3):
  scattered : the B blasts occupy B cells chosen uniformly at random among the M cells (the Part III arrangement)
  hot0.6    : n_hot = min(N, ceil(B / (m*h))) fields, chosen uniformly at random, receive all B blasts, which are
  hot1.0      placed uniformly at random among their n_hot*m cells. h is the nominal blast density of a hot field.
Because every design samples fields uniformly or by classifier score, never by position, only which fields are hot
matters, not whether they are adjacent.

Designs, all run on the SAME specimen (paired):
  alpha_field : uniform field sampling w/o replacement; ALPHA on field means           (valid for any arrangement)
  bet_field   : uniform field sampling; Shekhar et al. betting test on field means      (valid for any arrangement)
  prop_field  : classifier-driven field sampling (prop-MS); Shekhar et al. betting test (valid for any arrangement)
  naive_cell  : the SAME field order as alpha_field, with the cells flattened and analysed by a cell-level ALPHA as
                if cells had been sampled one at a time w/o replacement. This is valid only if the cell sequence is
                exchangeable, which holds for 'scattered' and fails under clustering.
  conv100     : first 100 cells in that order; calls 'above' if their blast fraction > theta (no guarantee)

The per-episode RNG is default_rng([seed, arrangement index, round(1000 p), episode]), so results do not depend on
how the run is chunked.
Usage: python scripts/18_spatial_confirm.py --out results/a3_spatial_theta010 --theta 0.10 --p_list ... --seed 20260611
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.alpha import assorter_bound, alpha_martingale
from cadpbs.weighted_betting import betting_test, sampler_uniform, sampler_prop_side

ARR = [("scattered", None), ("hot0.6", 0.6), ("hot1.0", 1.0)]

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True); ap.add_argument("--theta", type=float, required=True)
ap.add_argument("--p_list", type=float, nargs="+", required=True)
ap.add_argument("--episodes", type=int, default=500); ap.add_argument("--null_episodes", type=int, default=1000)
ap.add_argument("--seed", type=int, required=True)
ap.add_argument("--arr", nargs="+", default=[a for a, _ in ARR])
ap.add_argument("--chunk", type=int, default=250, help="episodes per part file")
ap.add_argument("--max_seconds", type=float, default=1e9, help="stop after this long; rerun to resume")
ap.add_argument("--features", default="results/embeddings_v1/cells_features.csv")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

df = pd.read_csv(a.features, dtype={"patient_id": str})
y = df.is_blast.values.astype(bool); score = df.probe_p_oof.values
blast_pat = set(df.loc[y, "patient_id"])
BL = np.flatnonzero(y); BG = np.flatnonzero(~y & ~df.patient_id.isin(blast_pat).values)
N, m = 100, 5; M = N * m; th = a.theta; u = assorter_bound(th); T0 = time.time()


def specimen(rng, nb, h):
    grid = np.empty(M, int)
    if h is None:
        pos = rng.choice(M, nb, replace=False)
    else:
        nhot = min(N, int(np.ceil(nb / (m * h))))
        hot = rng.choice(N, nhot, replace=False)
        pos = rng.choice((hot[:, None] * m + np.arange(m)).ravel(), nb, replace=False)
    mask = np.zeros(M, bool); mask[pos] = True
    grid[mask] = rng.choice(BL, nb, replace=False)
    grid[~mask] = rng.choice(BG, M - nb, replace=False)
    return grid.reshape(N, m)


for ai, (arr, h) in enumerate(ARR):
    if arr not in a.arr: continue
    for p in a.p_list:
        E = a.null_episodes if abs(p - th) < 1e-9 else a.episodes
        nb = int(round(p * M))
        for s0 in range(0, E, a.chunk):
            path = f"{a.out}/ep_{arr}_p{p:.2f}_e{s0:04d}.csv"
            if os.path.exists(path): continue
            if time.time() - T0 > a.max_seconds: print("time budget reached; rerun to resume"); sys.exit(0)
            rows = []
            for ep in range(s0, min(s0 + a.chunk, E)):
                rng = np.random.default_rng([a.seed, ai, int(round(p * 1000)), ep])
                cells = specimen(rng, nb, h)
                f = y[cells].mean(1) / (2 * th); S = score[cells].mean(1) / (2 * th)
                r = dict(theta=th, arr=arr, p=p, ep=ep, n_hot=int((y[cells].sum(1) > 0).sum()))
                order = rng.permutation(N)
                _, st = alpha_martingale(f[order], u, 0.5, N, alpha=0.05)
                r["alpha_field_cert"] = int(st > 0); r["alpha_field_cells"] = (st if st > 0 else N) * m
                for nm, smp in [("bet_field", sampler_uniform(N)), ("prop_field", sampler_prop_side(S))]:
                    o = betting_test(f, S, smp, 0.5, alpha=0.05, rng=rng)
                    r[f"{nm}_cert"] = int(o["stop"] > 0); r[f"{nm}_cells"] = (o["stop"] if o["stop"] > 0 else N) * m
                flat = y[cells[order]].reshape(-1)
                _, sc = alpha_martingale(flat / (2 * th), u, 0.5, M, alpha=0.05)
                r["naive_cell_cert"] = int(sc > 0); r["naive_cell_cells"] = (int(np.ceil(sc / m)) * m) if sc > 0 else M
                r["conv100_call_above"] = int(flat[:100].mean() > th)
                rows.append(r)
            pd.DataFrame(rows).to_csv(path, index=False)
            print(f"{arr} p={p:.2f} episodes {s0}-{min(s0 + a.chunk, E) - 1} done ({time.time() - T0:.0f}s)", flush=True)
print("ALL DONE")
