"""
Certification cost vs population size, M = 500 to 4,000 cells (RQ3; exploratory).

Designs (paired): ALPHA uniform, betting uniform, betting prop-MS, conventional first-n counts.
alpha = 0.05 (primary) and 0.01, read off the same wealth path.
When a specimen needs more than the 489 KU-Optofil blasts they are drawn with replacement (flag blast_reuse).

Usage: python scripts/28_m_scaling.py --out results/m_scaling_v1 --seed 20261021
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w, sampler_prop_w

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/m_scaling_v1"); ap.add_argument("--seed", type=int, default=20261021)
ap.add_argument("--M_list", type=int, nargs="+", default=[500, 1000, 2000, 4000]); ap.add_argument("--m", type=int, default=5)
ap.add_argument("--thetas", type=float, nargs="+", default=[0.10, 0.20])
ap.add_argument("--margins", type=float, nargs="+", default=[0.02, 0.05, 0.10])
ap.add_argument("--n_null", type=int, default=1000); ap.add_argument("--n_alt", type=int, default=400)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
df = pd.read_csv("results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})
y = df.is_blast.values.astype(bool); score = df.probe_p_oof.values
BL = np.flatnonzero(y); BG = np.flatnonzero(~y & ~df.patient_id.isin(set(df.loc[y, "patient_id"])).values)
LEVELS = (0.05, 0.01)

for th in a.thetas:
    for M in a.M_list:
        N = M // a.m; u = 1 / (2 * th); pi = np.full(N, 1 / N)
        for p_true in [th] + [round(th + d, 4) for d in a.margins]:
            path = f"{a.out}/ep_theta{th:.2f}_M{M}_p{p_true:.2f}.csv"
            if os.path.exists(path): continue
            n_ep = a.n_null if np.isclose(p_true, th) else a.n_alt
            rng = np.random.default_rng([a.seed, int(th * 100), M, int(round(p_true * 1000))]); t0 = time.time()
            nb = int(round(p_true * M)); reuse = nb > len(BL); rows = []
            for ep in range(n_ep):
                cells = np.r_[rng.choice(BL, nb, replace=reuse), rng.choice(BG, M - nb, replace=False)]
                rng.shuffle(cells); cells = cells.reshape(N, a.m)
                frac = y[cells].mean(1); S = score[cells].mean(1)
                r = dict(theta=th, M=M, N=N, p=p_true, ep=ep, blast_reuse=int(reuse))
                order = rng.permutation(N)
                lt, _ = alpha_martingale(frac[order] / (2 * th), u, 0.5, N, alpha=0.05)
                for lv in LEVELS:
                    hit = np.flatnonzero(lt >= np.log(1 / lv)); st = int(hit[0]) + 1 if len(hit) else -1
                    r[f"alpha_uniform_cert_a{lv}"] = int(st > 0); r[f"alpha_uniform_cells_a{lv}"] = (st if st > 0 else N) * a.m
                for nm, smp in (("bet_uniform", sampler_uniform_w(N)), ("bet_prop", sampler_prop_w(pi, S))):
                    o = betting_test_w(frac, pi, S, smp, th, rng=rng, alphas=LEVELS)
                    for lv in LEVELS:
                        st = o["stops"][lv]; r[f"{nm}_cert_a{lv}"] = int(st > 0); r[f"{nm}_cells_a{lv}"] = (st if st > 0 else N) * a.m
                flat = y[cells[order].reshape(-1)]
                for nc in (100, 200, 500):
                    r[f"conv{nc}_call_above"] = int(flat[:nc].mean() > th) if nc <= M else np.nan
                rows.append(r)
            pd.DataFrame(rows).to_csv(path, index=False)
            print(f"theta={th} M={M} p={p_true}: {n_ep} episodes {time.time()-t0:.0f}s reuse={reuse}", flush=True)

import glob
D = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{a.out}/ep_*.csv"))], ignore_index=True)
D.to_csv(f"{a.out}/episodes.csv", index=False); json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print("episodes:", len(D))
