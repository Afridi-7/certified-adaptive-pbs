"""
Confusable background (RQ3, RQ4; exploratory): a share s in {0, 0.1, 0.2} of non-blast cells drawn from
metamyelocytes, myelocytes and reactive lymphocytes, the classes most often mistaken for blasts.
Designs as in 28_m_scaling.py; M = 500, theta in {0.10, 0.20}.

Usage: python scripts/33_confusable_background.py --out results/confusable_bg_v1 --seed 20261026
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w, sampler_prop_w

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/confusable_bg_v1"); ap.add_argument("--seed", type=int, default=20261026)
ap.add_argument("--n_null", type=int, default=1000); ap.add_argument("--n_alt", type=int, default=400)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
df = pd.read_csv("results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})
y = df.is_blast.values.astype(bool); score = df.probe_p_oof.values
BL = np.flatnonzero(y); BG = np.flatnonzero(~y & ~df.patient_id.isin(set(df.loc[y, "patient_id"])).values)
CONF = np.flatnonzero(~y & df.cell_type.isin(["metamyelocyte", "myelocyte", "reactive_lymphocyte"]).values)
BG = np.setdiff1d(BG, CONF)
N, m = 100, 5; M = N * m; rows = []
for th in (0.10, 0.20):
    u = 1 / (2 * th); pi = np.full(N, 1 / N)
    for share in (0.0, 0.10, 0.20):
        for p in (th, round(th + 0.02, 4), round(th + 0.05, 4)):
            n_ep = a.n_null if np.isclose(p, th) else a.n_alt
            rng = np.random.default_rng([a.seed, int(th * 100), int(share * 100), int(round(p * 1000))])
            nb = int(round(p * M)); nc = int(round(share * (M - nb)))
            for ep in range(n_ep):
                cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(CONF, nc, replace=False), rng.choice(BG, M - nb - nc, replace=False)]
                rng.shuffle(cells); cells = cells.reshape(N, m); frac = y[cells].mean(1); S = score[cells].mean(1)
                r = dict(theta=th, share=share, p=p, ep=ep, mean_score_nonblast=float(score[cells][~y[cells]].mean()))
                lt, _ = alpha_martingale(frac[rng.permutation(N)] / (2 * th), u, 0.5, N, alpha=0.05)
                hit = np.flatnonzero(lt >= np.log(20)); st = int(hit[0]) + 1 if len(hit) else -1
                r["alpha_uniform_cert"] = int(st > 0); r["alpha_uniform_cells"] = (st if st > 0 else N) * m
                for nm, smp in (("bet_uniform", sampler_uniform_w(N)), ("bet_prop", sampler_prop_w(pi, S))):
                    o = betting_test_w(frac, pi, S, smp, th, rng=rng); st = o["stop"]
                    r[f"{nm}_cert"] = int(st > 0); r[f"{nm}_cells"] = (st if st > 0 else N) * m
                rows.append(r)
D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)
summ = []
for (th, share, p), g in D.groupby(["theta", "share", "p"]):
    r = dict(theta=th, share=share, p=p, n=len(g))
    for d in ["alpha_uniform", "bet_uniform", "bet_prop"]:
        r[f"{d}_false_cert" if np.isclose(p, th) else f"{d}_median_cells"] = float(g[f"{d}_cert"].mean() if np.isclose(p, th) else g[f"{d}_cells"].median())
    summ.append(r)
S_ = pd.DataFrame(summ); S_.to_csv(f"{a.out}/summary.csv", index=False); json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
pd.set_option("display.width", 200); print(S_.round(3).to_string(index=False)); print("confusable pool size:", len(CONF))
