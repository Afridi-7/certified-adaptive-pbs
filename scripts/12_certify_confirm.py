"""
Part III confirmatory data generation (PREREG amendment 2). One specimen per episode; EVERY design is run on the
SAME specimen (paired). Specimen = N = 100 fields x m = 5 cells = 500 cells. Blasts drawn WITHOUT replacement from
all 489 blasts; background WITHOUT replacement from all blast-free patients. Allowed because every cell's classifier
score is cross-fitted (out-of-patient), so mixing folds introduces no leakage; it removes the blast reuse of the
exploratory runs. theta and the p grid are fixed in PREREG amendment 2.
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.alpha import assorter_bound, alpha_martingale, sample_uniform, sample_prioritised
from cadpbs.weighted_betting import betting_test, sampler_uniform, sampler_prop_side

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True); ap.add_argument("--theta", type=float, required=True)
ap.add_argument("--p_list", type=float, nargs="+", required=True)
ap.add_argument("--episodes", type=int, default=1000); ap.add_argument("--seed", type=int, required=True)
ap.add_argument("--N", type=int, default=100); ap.add_argument("--m", type=int, default=5); ap.add_argument("--alpha", type=float, default=0.05)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
df = pd.read_csv("results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})
y = df.is_blast.values.astype(bool); score = df.probe_p_oof.values
blast_pat = set(df.loc[y, "patient_id"])
BL = np.flatnonzero(y); BG = np.flatnonzero(~y & ~df.patient_id.isin(blast_pat).values)
N, m, M, th = a.N, a.m, a.N * a.m, a.theta; u = assorter_bound(th)
for p_true in a.p_list:
    path = f"{a.out}/episodes_p{p_true:.2f}.csv"
    if os.path.exists(path): print("skip", path); continue
    rng = np.random.default_rng(a.seed + int(round(p_true * 1000)))
    nb = int(round(p_true * M)); assert nb <= len(BL), "would require blast reuse"
    rows = []; t0 = time.time()
    for ep in range(a.episodes):
        cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(BG, M - nb, replace=False)]
        rng.shuffle(cells); cells = cells.reshape(N, m)
        f = y[cells].mean(1) / (2 * th); S = score[cells].mean(1) / (2 * th)
        r = dict(theta=th, p=p_true, ep=ep, M=M)
        order, _ = sample_uniform(N, rng)
        _, st = alpha_martingale(f[order], u, 0.5, N, alpha=a.alpha)
        r["alpha_uniform_cert"] = int(st > 0); r["alpha_uniform_cells"] = (st if st > 0 else N) * m
        for nm, smp, cv in [("bet_uniform", sampler_uniform(N), False), ("bet_prop", sampler_prop_side(S), False),
                            ("bet_prop_cv", sampler_prop_side(S), True)]:
            o = betting_test(f, S, smp, 0.5, alpha=a.alpha, use_cv=cv, rng=rng)
            r[f"{nm}_cert"] = int(o["stop"] > 0); r[f"{nm}_cells"] = (o["stop"] if o["stop"] > 0 else N) * m
        o2, _ = sample_prioritised(N, score[cells].max(1), rng, 0.5)
        _, st3 = alpha_martingale(f[o2], u, 0.5, N, alpha=a.alpha)
        r["naive_cert"] = int(st3 > 0); r["naive_cells"] = (st3 if st3 > 0 else N) * m
        flat = y[cells[order].reshape(-1)]
        for nc in (100, 200): r[f"conv{nc}_call_above"] = int(flat[:nc].mean() > th)
        rows.append(r)
    pd.DataFrame(rows).to_csv(path, index=False); print(f"theta={th} p={p_true}: {a.episodes} episodes {time.time()-t0:.0f}s", flush=True)
import glob
parts = sorted(glob.glob(f"{a.out}/episodes_p*.csv"))
pd.concat([pd.read_csv(x) for x in parts]).to_csv(f"{a.out}/episodes.csv", index=False)
json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
