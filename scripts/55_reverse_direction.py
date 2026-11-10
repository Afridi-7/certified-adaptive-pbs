"""
Reverse-direction replication (2026-11-10; exploratory; RQ3, RQ4).

Every earlier result trains the blast probe on KU-Optofil and applies it to KU-Optofil (in distribution) and to MLL23
and RV-PBS (acquisition shift). The in-distribution saving of classifier-driven sampling therefore rests on blasts from
11 patients of one laboratory. This script exchanges the roles of the two large sources:
  mll23_in : probe cross-fitted on MLL23 (5 folds), specimens of MLL23 cells scored out of fold   -- in distribution
  ku_shift : probe trained on all MLL23 cells, applied unchanged to KU-Optofil cells              -- acquisition shift
Same probe as before (standardisation + logistic regression, C = 0.1, balanced class weights, frozen DINOv2 ViT-S/14
embeddings), same specimen construction as scripts/49 and scripts/53 (N = 100 fields x m = 5 cells, scattered), same
designs as scripts/53. theta = 0.10 (ICC 2022) for both sources; theta = 0.20 (WHO 5th ed.) in distribution only, as in
the original direction.

LIMITATION (stated wherever the result is used). MLL23 has no patient identifiers, so its folds are random by cell
(stratified by class): cells of one patient fall in training and test folds, and the in-distribution scores are
optimistic compared with the patient-grouped folds used for KU-Optofil. The run shows whether the PATTERN (saving in
distribution, loss of the saving under shift, validity in both) is specific to KU-Optofil; it does not add independent
patients. Labels also differ between the sources (MLL23: myeloblasts; KU-Optofil: blasts), which is part of the shift.

Per-specimen RNG: default_rng([20261111, source index, theta index, round(1000 p), episode]); folds: seed 20261111.
The seed family is new. Uncertified specimens count all 500 cells. Costs are counts of expert-verified cells.

Usage: python scripts/55_reverse_direction.py --out results/reverse_direction_v1
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import binomtest, beta, wilcoxon
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import sampler_uniform_w, sampler_prop_w
from cadpbs.hedged_betting import hedged_betting_w
from cadpbs.cv_hedged_betting import hedged_cv_betting_w, sampler_strata_w, sampler_mix_w

SOURCES = ["mll23_in", "ku_shift"]
DESIGNS = ["alpha_uniform", "hedge_uniform", "hedge_prop", "hedge_mix50", "cv_uniform", "cv_mix50", "hedge_strata", "cv_strata"]
ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True); ap.add_argument("--seed", type=int, default=20261111)
ap.add_argument("--n_null", type=int, default=1000); ap.add_argument("--n_alt", type=int, default=400)
ap.add_argument("--n_boot", type=int, default=4000); ap.add_argument("--folds", type=int, default=5)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); T0 = time.time(); ROOT = os.path.join(HERE, "..")


def probe():
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, C=0.1, class_weight="balanced"))


def auc(score, y):
    r = pd.Series(score).rank().values; npos = int(y.sum()); nneg = len(y) - npos
    return float((r[y].sum() - npos * (npos + 1) / 2) / (npos * nneg))


# ------------------------------------------------------------------ data and scores
mm = pd.read_csv(f"{ROOT}/data/mll23/mll23_cells.csv"); Em = np.load(f"{ROOT}/data/mll23/mll23_embeddings_dinov2_vits14.npy").astype(np.float32)
keep = ~mm.cls.str.contains("atypical|smudge|normoblast", case=False, regex=True)          # as in scripts/38, 49, 53
Em, mm = Em[keep.values], mm[keep].reset_index(drop=True)
ym = mm.cls.str.contains("myeloblast", case=False).values
ku = pd.read_csv(f"{ROOT}/results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})
Eku = np.load(f"{ROOT}/results/embeddings_v1/embeddings_dinov2_vits14.npy").astype(np.float32)
yk = ku.is_blast.values.astype(bool)

oof = np.zeros(len(mm))                                    # out-of-fold scores on MLL23 (folds random by cell: see LIMITATION)
for tr, te in StratifiedKFold(a.folds, shuffle=True, random_state=a.seed).split(Em, mm.cls.values):
    oof[te] = probe().fit(Em[tr], ym[tr]).predict_proba(Em[te])[:, 1]
score_ku = probe().fit(Em, ym).predict_proba(Eku)[:, 1]    # MLL23-trained probe applied unchanged to KU-Optofil
T = float(np.quantile(oof[~ym], 0.95))                     # operating point: 95% specificity on MLL23 (out of fold)
desc = dict(n_mll23=int(len(mm)), n_mll23_blast=int(ym.sum()), n_ku=int(len(ku)), n_ku_blast=int(yk.sum()), threshold=T,
            mll23_in=dict(auc=auc(oof, ym), sensitivity=float((oof[ym] >= T).mean()), fp=float((oof[~ym] >= T).mean())),
            ku_shift=dict(auc=auc(score_ku, yk), sensitivity=float((score_ku[yk] >= T).mean()), fp=float((score_ku[~yk] >= T).mean())))
print(json.dumps(desc, indent=1), flush=True)

bp = set(ku.loc[yk, "patient_id"])
POOL = {"mll23_in": (ym, oof, np.flatnonzero(ym), np.flatnonzero(~ym)),
        "ku_shift": (yk, score_ku, np.flatnonzero(yk), np.flatnonzero(~yk & ~ku.patient_id.isin(bp).values))}   # background: blast-free patients
GRID = {"mll23_in": [(0.10, [0.10, 0.12, 0.15, 0.20]), (0.20, [0.20, 0.22, 0.25, 0.30])], "ku_shift": [(0.10, [0.10, 0.12, 0.15, 0.20])]}
N, m = 100, 5; M = N * m; pi = np.full(N, 1.0 / N)

# ------------------------------------------------------------------ specimens and designs
rows = []
for si, src in enumerate(SOURCES):
    y, score, BL, BG = POOL[src]
    for ti, (th, plist) in enumerate(GRID[src]):
        for p in plist:
            nb = int(round(p * M)); assert nb <= len(BL) and M - nb <= len(BG)
            E = a.n_null if abs(p - th) < 1e-9 else a.n_alt
            for ep in range(E):
                rng = np.random.default_rng([a.seed, si, ti, int(round(p * 1000)), ep])
                cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(BG, M - nb, replace=False)]
                rng.shuffle(cells); cells = cells.reshape(N, m)
                f = y[cells].mean(1).astype(float); S = score[cells].mean(1)
                prop = sampler_prop_w(pi, S); strata = sampler_strata_w(pi, S)
                run = {
                    "alpha_uniform": lambda: alpha_martingale(f[rng.permutation(N)] / (2 * th), 1 / (2 * th), 0.5, N, alpha=0.05)[1],
                    "hedge_uniform": lambda: hedged_betting_w(f, pi, sampler_uniform_w(N), th, rng=rng)["stop"],
                    "hedge_prop": lambda: hedged_betting_w(f, pi, prop, th, rng=rng)["stop"],
                    "hedge_mix50": lambda: hedged_betting_w(f, pi, sampler_mix_w(prop, 0.5), th, rng=rng)["stop"],
                    "cv_uniform": lambda: hedged_cv_betting_w(f, pi, S, sampler_uniform_w(N), th, rng=rng)["stop"],
                    "cv_mix50": lambda: hedged_cv_betting_w(f, pi, S, sampler_mix_w(prop, 0.5), th, rng=rng)["stop"],
                    "hedge_strata": lambda: hedged_betting_w(f, pi, strata, th, rng=rng)["stop"],
                    "cv_strata": lambda: hedged_cv_betting_w(f, pi, S, strata, th, rng=rng)["stop"],
                }
                r = dict(source=src, theta=th, p=p, ep=ep)
                for d in DESIGNS:                                  # fixed order, so the random stream is reproducible
                    s = run[d](); r[f"{d}_cert"] = int(s > 0); r[f"{d}_cells"] = (s if s > 0 else N) * m
                rows.append(r)
            print(f"  {src} theta={th} p={p}: {E} specimens ({time.time() - T0:.0f}s)", flush=True)

# ------------------------------------------------------------------ summaries
D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)
out = dict(descriptive=desc, validity={}, median_cells={}, ratios={})
for si, src in enumerate(SOURCES):
    rb = np.random.default_rng([a.seed, si, 999])
    for (th, p), g in D[D.source == src].groupby(["theta", "p"]):
        key = f"{src}_theta{th:.2f}_p{p:.2f}"
        if abs(p - th) < 1e-9:
            for d in DESIGNS:
                k = int(g[f"{d}_cert"].sum()); n = len(g)
                out["validity"][f"{key}_{d}"] = dict(rate=k / n, cp95_upper=float(beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0,
                                                     p_excess=float(binomtest(k, n, 0.05, alternative="greater").pvalue))
            continue
        out["median_cells"][key] = {d: float(g[f"{d}_cells"].median()) for d in DESIGNS}
        z = g["alpha_uniform_cells"].values
        for num in DESIGNS[1:]:
            x = g[f"{num}_cells"].values; bs = []
            for _ in range(a.n_boot):
                i = rb.integers(0, len(g), len(g)); bs.append(np.median(x[i]) / np.median(z[i]))
            try: wp = float(wilcoxon(x, z).pvalue)
            except ValueError: wp = 1.0
            out["ratios"][f"{key}_{num}_over_alpha_uniform"] = dict(ratio=float(np.median(x) / np.median(z)),
                                                                   ci95=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))], wilcoxon_p=wp)
import platform, scipy, sklearn
json.dump(dict(config=vars(a), env=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                                        pandas=pd.__version__, sklearn=sklearn.__version__)), open(f"{a.out}/config.json", "w"), indent=1)
json.dump(out, open(f"{a.out}/summary.json", "w"), indent=1)
print(pd.DataFrame(out["median_cells"]).T.to_string()); print("validity:", {k: v["rate"] for k, v in out["validity"].items()})
print(f"DONE in {time.time() - T0:.0f}s", flush=True)
