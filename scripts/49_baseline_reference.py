"""
Paired comparison of all Part III designs against uniform ALPHA as the reference (RQ3, RQ4; exploratory).

Designs: ALPHA uniform; approximate-Kelly betting (uniform, prop-MS); hedged betting (uniform, prop-MS,
50% mixture). Sources: KU-Optofil, MLL23, RV-PBS. N = 100 fields x 5 cells; theta = 0.10 (all) and 0.20
(KU-Optofil). Costs are expert-verified cells; uncertified specimens count all 500.

Usage: python scripts/49_baseline_reference.py --source ku --out results/baseline_reference_v1
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import binomtest, beta, wilcoxon
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w, sampler_prop_w
from cadpbs.hedged_betting import hedged_betting_w

SOURCES = ["ku", "mll23", "rvpbs"]
DESIGNS = ["alpha_uniform", "bet_uniform", "bet_prop", "hedge_uniform", "hedge_prop", "hedge_mix50"]
T_KU = 0.0310392350137264

ap = argparse.ArgumentParser()
ap.add_argument("--source", choices=SOURCES, required=True); ap.add_argument("--out", required=True)
ap.add_argument("--seed", type=int, default=20261101); ap.add_argument("--n_null", type=int, default=1000)
ap.add_argument("--n_alt", type=int, default=400); ap.add_argument("--n_boot", type=int, default=4000)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); T0 = time.time()
ROOT = os.path.join(HERE, "..")
ku = pd.read_csv(f"{ROOT}/results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})


def ku_probe():
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    Eku = np.load(f"{ROOT}/results/embeddings_v1/embeddings_dinov2_vits14.npy").astype(np.float32)
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, C=0.1, class_weight="balanced")).fit(Eku, ku.is_blast.values)


if a.source == "ku":
    y = ku.is_blast.values.astype(bool); score = ku.probe_p_oof.values
    bp = set(ku.loc[y, "patient_id"]); BL = np.flatnonzero(y); BG = np.flatnonzero(~y & ~ku.patient_id.isin(bp).values)
    grid = [(0.10, [0.10, 0.12, 0.15, 0.20]), (0.20, [0.20, 0.22, 0.25, 0.30])]
elif a.source == "mll23":
    mm = pd.read_csv(f"{ROOT}/data/mll23/mll23_cells.csv"); Em = np.load(f"{ROOT}/data/mll23/mll23_embeddings_dinov2_vits14.npy").astype(np.float32)
    keep = ~mm.cls.str.contains("atypical|smudge|normoblast", case=False, regex=True)
    Em, mm = Em[keep.values], mm[keep].reset_index(drop=True)
    y = mm.cls.str.contains("myeloblast", case=False).values; score = ku_probe().predict_proba(Em)[:, 1]
    BL = np.flatnonzero(y); BG = np.flatnonzero(~y); grid = [(0.10, [0.10, 0.12, 0.15, 0.20])]
else:
    rv = pd.read_csv(f"{ROOT}/results/external_v1/rvpbs_features.csv")
    Erv = np.load(f"{ROOT}/results/external_v1/rvpbs_embeddings_dinov2_vits14.npy").astype(np.float32)
    y = rv.is_blast.values.astype(bool); score = ku_probe().predict_proba(Erv)[:, 1]
    BL = np.flatnonzero(y); BG = np.flatnonzero(~y); grid = [(0.10, [0.10, 0.12, 0.15, 0.20])]
recall = float((score[y] >= T_KU).mean()); fp = float((score[~y] >= T_KU).mean())
print(f"{a.source}: {len(BL)} blasts, {len(BG)} background cells, recall at KU threshold {recall:.3f}, FP {fp:.3f}", flush=True)

N, m = 100, 5; M = N * m; pi = np.full(N, 1.0 / N); si = SOURCES.index(a.source)


def mix(S, lam):
    base = sampler_prop_w(pi, S)
    return lambda rem, t: (1 - lam) * base(rem, t) + lam / len(rem)


rows = []
for ti, (th, plist) in enumerate(grid):
    for p in plist:
        nb = int(round(p * M)); assert nb <= len(BL) and M - nb <= len(BG)
        E = a.n_null if abs(p - th) < 1e-9 else a.n_alt
        for ep in range(E):
            rng = np.random.default_rng([a.seed, si, ti, int(round(p * 1000)), ep])
            cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(BG, M - nb, replace=False)]
            rng.shuffle(cells); cells = cells.reshape(N, m)
            f = y[cells].mean(1).astype(float); S = score[cells].mean(1)
            r = dict(source=a.source, theta=th, p=p, ep=ep)
            _, st = alpha_martingale(f[rng.permutation(N)] / (2 * th), 1 / (2 * th), 0.5, N, alpha=0.05)
            r["alpha_uniform_stop"] = st
            r["bet_uniform_stop"] = betting_test_w(f, pi, S, sampler_uniform_w(N), th, rng=rng)["stop"]
            r["bet_prop_stop"] = betting_test_w(f, pi, S, sampler_prop_w(pi, S), th, rng=rng)["stop"]
            r["hedge_uniform_stop"] = hedged_betting_w(f, pi, sampler_uniform_w(N), th, rng=rng)["stop"]
            r["hedge_prop_stop"] = hedged_betting_w(f, pi, sampler_prop_w(pi, S), th, rng=rng)["stop"]
            r["hedge_mix50_stop"] = hedged_betting_w(f, pi, mix(S, 0.5), th, rng=rng)["stop"]
            for d in DESIGNS:
                s = r.pop(f"{d}_stop"); r[f"{d}_cert"] = int(s > 0); r[f"{d}_cells"] = (s if s > 0 else N) * m
            rows.append(r)
        print(f"  theta={th} p={p}: {E} specimens ({time.time() - T0:.0f}s)", flush=True)

D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes_{a.source}.csv", index=False)
rb = np.random.default_rng([a.seed, si, 999]); out = dict(source=a.source, recall_at_ku_threshold=recall, fp_at_ku_threshold=fp,
                                                       n_blast_pool=int(len(BL)), n_background_pool=int(len(BG)), validity={}, median_cells={}, ratios={})
for (th, p), g in D.groupby(["theta", "p"]):
    key = f"theta{th:.2f}_p{p:.2f}"
    if abs(p - th) < 1e-9:
        for d in DESIGNS:
            k = int(g[f"{d}_cert"].sum()); n = len(g)
            out["validity"][f"{key}_{d}"] = dict(rate=k / n, cp95_upper=float(beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0,
                                                 p_excess=float(binomtest(k, n, 0.05, alternative="greater").pvalue))
        continue
    out["median_cells"][key] = {d: float(g[f"{d}_cells"].median()) for d in DESIGNS}
    for num, den in [("bet_prop", "alpha_uniform"), ("bet_prop", "bet_uniform"), ("bet_uniform", "alpha_uniform"),
                     ("hedge_uniform", "alpha_uniform"), ("hedge_prop", "alpha_uniform"), ("hedge_prop", "hedge_uniform"),
                     ("hedge_mix50", "alpha_uniform")]:
        x, z = g[f"{num}_cells"].values, g[f"{den}_cells"].values; bs = []
        for _ in range(a.n_boot):
            i = rb.integers(0, len(g), len(g)); bs.append(np.median(x[i]) / np.median(z[i]))
        try: wp = float(wilcoxon(x, z).pvalue)
        except ValueError: wp = 1.0
        out["ratios"][f"{key}_{num}_over_{den}"] = dict(ratio=float(np.median(x) / np.median(z)),
                                                        ci95=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))], wilcoxon_p=wp)
import platform, numpy, scipy, sklearn
cfg = dict(config=vars(a), env=dict(python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__,
                                    pandas=pd.__version__, sklearn=sklearn.__version__))
json.dump(cfg, open(f"{a.out}/config_{a.source}.json", "w"), indent=1)          # run settings (not compared by verify)
json.dump(out, open(f"{a.out}/summary_{a.source}.json", "w"), indent=1)         # numbers only
print(f"DONE {a.source} in {time.time() - T0:.0f}s", flush=True)
