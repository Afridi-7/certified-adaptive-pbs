"""
Background restricted to leukocytes (exploratory; found 2026-10-20). RQ3.

Every Part III composite (scripts/12, 18, 19, 27, 28, 33) draws its non-blast cells from ALL non-blast KU-Optofil images of
blast-free patients. That pool contains giant platelets, platelet clusters and erythroblasts, which a leukocyte
differential does not count. The specimen's p is then blasts / all images, not blasts / leukocytes. Validity of the tests
does not depend on what the units are, but (i) the clinical reading of p shifts and (ii) platelets are easy negatives for
the classifier, which may flatter classifier-driven sampling.

Here the background is restricted to leukocyte classes (giant_platelet, platelet_cluster and erythroblast removed) and the
Part III designs are re-run as in scripts/33 (share 0): M = 500 (100 fields of 5); theta = 0.10 and 0.20; p = theta
(validity, 1,000 specimens) and theta + 0.02, + 0.05 (cost, 400 specimens). Seed family 20261032.

Usage: python scripts/41_leukocyte_background.py --out results/leukocyte_bg_v1 --seed 20261032
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import binomtest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w, sampler_prop_w

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/leukocyte_bg_v1"); ap.add_argument("--seed", type=int, default=20261032)
ap.add_argument("--n_null", type=int, default=1000); ap.add_argument("--n_alt", type=int, default=400)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
df = pd.read_csv("results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})
y = df.is_blast.values.astype(bool); score = df.probe_p_oof.values
NONLEUKO = ["giant_platelet", "platelet_cluster", "erythroblast"]
BL = np.flatnonzero(y); BG_all = np.flatnonzero(~y & ~df.patient_id.isin(set(df.loc[y, "patient_id"])).values)
BG = BG_all[~np.isin(df.cell_type.values[BG_all], NONLEUKO)]
info = dict(bg_all=int(len(BG_all)), bg_leuko=int(len(BG)),
            nonleuko_share_of_bg=float(np.isin(df.cell_type.values[BG_all], NONLEUKO).mean()),
            bg_composition_all={str(k): int(v) for k, v in pd.Series(df.cell_type.values[BG_all]).value_counts().items()},
            mean_score_bg_all=float(score[BG_all].mean()), mean_score_bg_leuko=float(score[BG].mean()),
            fp_rate_bg_all_at_0p031=float((score[BG_all] >= 0.0310392350137264).mean()),
            fp_rate_bg_leuko_at_0p031=float((score[BG] >= 0.0310392350137264).mean()))
N, m = 100, 5; M = N * m; rows = []
for th in (0.10, 0.20):
    u = 1 / (2 * th); pi = np.full(N, 1 / N)
    for p in (th, round(th + 0.02, 4), round(th + 0.05, 4)):
        n_ep = a.n_null if np.isclose(p, th) else a.n_alt
        rng = np.random.default_rng([a.seed, int(th * 100), int(round(p * 1000))]); nb = int(round(p * M))
        for ep in range(n_ep):
            cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(BG, M - nb, replace=False)]
            rng.shuffle(cells); cells = cells.reshape(N, m); frac = y[cells].mean(1); S = score[cells].mean(1)
            r = dict(theta=th, p=p, ep=ep)
            lt, _ = alpha_martingale(frac[rng.permutation(N)] / (2 * th), u, 0.5, N, alpha=0.05)
            hit = np.flatnonzero(lt >= np.log(20)); st = int(hit[0]) + 1 if len(hit) else -1
            r["alpha_uniform_cert"] = int(st > 0); r["alpha_uniform_cells"] = (st if st > 0 else N) * m
            for nm, smp in (("bet_uniform", sampler_uniform_w(N)), ("bet_prop", sampler_prop_w(pi, S))):
                o = betting_test_w(frac, pi, S, smp, th, rng=rng); st = o["stop"]
                r[f"{nm}_cert"] = int(st > 0); r[f"{nm}_cells"] = (st if st > 0 else N) * m
            rows.append(r)
D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)
summ = []
for (th, p), g in D.groupby(["theta", "p"]):
    r = dict(theta=th, p=p, n=len(g))
    for d in ["alpha_uniform", "bet_uniform", "bet_prop"]:
        if np.isclose(p, th):
            k = int(g[f"{d}_cert"].sum()); r[f"{d}_false_cert"] = k / len(g)
            r[f"{d}_p_excess"] = float(binomtest(k, len(g), 0.05, alternative="greater").pvalue)
        else:
            r[f"{d}_median_cells"] = float(g[f"{d}_cells"].median())
    summ.append(r)
S_ = pd.DataFrame(summ); S_.to_csv(f"{a.out}/summary.csv", index=False)
json.dump(dict(vars(a), **info), open(f"{a.out}/config.json", "w"), indent=1)
pd.set_option("display.width", 220); print(json.dumps(info, indent=1)); print(S_.round(3).to_string(index=False))
