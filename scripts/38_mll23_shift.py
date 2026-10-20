"""
Second acquisition-shift source: MLL23 (P1-8; exploratory). RQ4. Mirrors the pre-registered RV-PBS analysis
(scripts/16_external_analysis.py) and the defensive mixture (scripts/19_mixture_confirm.py) on a larger external source.

Input : data/mll23/mll23_embeddings_dinov2_vits14.npy + mll23_cells.csv (notebooks/37_mll23_embeddings_colab.ipynb).
Probe : identical to scripts/16: logistic regression (standardised, C = 0.1, balanced) trained on ALL KU-Optofil cells;
        threshold fixed at the KU cross-fitted 95%-specificity value t = 0.0310392350137264.
Labels: blast = the myeloblast class. Excluded from specimens: atypical promyelocytes (class promyelocyte_atypical; blast
        equivalents in acute promyelocytic leukaemia under WHO 2022 and ICC 2022, which our labels do not count), smudge
        cells and normoblasts (not counted in a leukocyte differential). The 2026-10-19 version matched the wrong class
        name ("atypical_promyelocyte") and would have kept atypical promyelocytes as background; fixed 2026-10-20.
Uncertainty: bootstrap over images (the release may carry no patient identifiers; then intervals are too narrow).
Certification: specimens of N = 100 fields x 5 cells from MLL23 cells; theta = 0.10 (ICC 2022); designs uniform,
        classifier-driven (prop-MS) and the 50% mixture, Shekhar et al. betting test; p = 0.10 (validity, 1,000 specimens),
        0.15 and 0.20 (cost, 400 specimens).

Usage: python scripts/38_mll23_shift.py --mll data/mll23 --out results/mll23_shift_v1 --seed 20261029
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import binomtest
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w, sampler_prop_w

ap = argparse.ArgumentParser()
ap.add_argument("--mll", default="data/mll23"); ap.add_argument("--out", default="results/mll23_shift_v1")
ap.add_argument("--seed", type=int, default=20261029); ap.add_argument("--n_boot", type=int, default=2000)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
T_KU, KU_RECALL, KU_FP = 0.0310392350137264, 0.9059, 0.05
Eku = np.load("results/embeddings_v1/embeddings_dinov2_vits14.npy").astype(np.float32)
ku = pd.read_csv("results/embeddings_v1/cells_features.csv")
Em = np.load(f"{a.mll}/mll23_embeddings_dinov2_vits14.npy").astype(np.float32)
mm = pd.read_csv(f"{a.mll}/mll23_cells.csv")
keep = ~mm.cls.str.contains("atypical|smudge|normoblast", case=False, regex=True)     # class names: promyelocyte_atypical, smudge_cell, normoblast
Em, mm = Em[keep.values], mm[keep].reset_index(drop=True)
mm["is_blast"] = mm.cls.str.contains("myeloblast", case=False).astype(int)
clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, C=0.1, class_weight="balanced")).fit(Eku, ku.is_blast.values)
mm["score"] = clf.predict_proba(Em)[:, 1]; mm["call"] = mm.score >= T_KU
rng = np.random.default_rng(a.seed)


def boot(stat):
    v = [stat(mm.iloc[rng.integers(0, len(mm), len(mm))]) for _ in range(a.n_boot)]
    return [float(np.quantile(v, .025)), float(np.quantile(v, .975))]


rec = lambda d: d.call[d.is_blast == 1].mean(); fp = lambda d: d.call[d.is_blast == 0].mean()
res = dict(n=int(len(mm)), n_blast=int(mm.is_blast.sum()), auc=float(roc_auc_score(mm.is_blast, mm.score)),
           recall=float(rec(mm)), recall_ci95=boot(rec), fp=float(fp(mm)), fp_ci95=boot(fp), ku_recall=KU_RECALL, ku_fp=KU_FP,
           fp_by_class=mm[mm.is_blast == 0].groupby("cls").call.mean().round(3).to_dict())

BL = np.flatnonzero(mm.is_blast.values == 1); BG = np.flatnonzero(mm.is_blast.values == 0)
yb = mm.is_blast.values.astype(float); sc = mm.score.values; th = 0.10; N, m = 100, 5; pi = np.full(N, 1 / N); rows = []


def mix(S, lam):
    base = sampler_prop_w(pi, S)
    def q(rem, t): return (1 - lam) * base(rem, t) + lam / len(rem)
    return q


for p, reps in ((0.10, 1000), (0.15, 400), (0.20, 400)):
    nb = int(round(p * N * m))
    for ep in range(reps):
        cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(BG, N * m - nb, replace=False)]; rng.shuffle(cells)
        cells = cells.reshape(N, m); f = yb[cells].mean(1); S = sc[cells].mean(1); r = dict(p=p, ep=ep)
        for nm, smp in (("uniform", sampler_uniform_w(N)), ("prop", sampler_prop_w(pi, S)), ("mix50", mix(S, 0.5))):
            o = betting_test_w(f, pi, S, smp, th, rng=rng); r[f"{nm}_cert"] = int(o["stop"] > 0); r[f"{nm}_cells"] = (o["stop"] if o["stop"] > 0 else N) * m
        rows.append(r)
D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)
nul = D[D.p == 0.10]
res["validity_p_eq_theta"] = {d: dict(rate=float(nul[f"{d}_cert"].mean()),
                                      p_excess=float(binomtest(int(nul[f"{d}_cert"].sum()), len(nul), 0.05, alternative="greater").pvalue))
                              for d in ["uniform", "prop", "mix50"]}
res["median_cells"] = {str(p): {d: float(g[f"{d}_cells"].median()) for d in ["uniform", "prop", "mix50"]} for p, g in D[D.p > 0.10].groupby("p")}
# effect sizes: ratio of median cells (paired bootstrap over specimens, separate RNG) and paired Wilcoxon
from scipy.stats import wilcoxon
rb = np.random.default_rng(a.seed + 1); res["cost_ratio"] = {}
for p, g in D[D.p > 0.10].groupby("p"):
    for x_, y_ in (("prop", "uniform"), ("mix50", "uniform"), ("mix50", "prop")):
        x, y = g[f"{x_}_cells"].values, g[f"{y_}_cells"].values; bs = []
        for _ in range(4000):
            i = rb.integers(0, len(g), len(g)); bs.append(np.median(x[i]) / np.median(y[i]))
        res["cost_ratio"][f"p{p:.2f}_{x_}_over_{y_}"] = dict(ratio=float(np.median(x) / np.median(y)),
                                                              ci95=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))],
                                                              wilcoxon_p=float(wilcoxon(x, y).pvalue))
json.dump(res, open(f"{a.out}/summary.json", "w"), indent=1); json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print(json.dumps(res, indent=1))
