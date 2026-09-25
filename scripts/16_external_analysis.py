"""
RQ4 external-source analysis (PREREG amendment 2, H4a-H4c). Written before the RV-PBS embeddings exist.
Classifier: logistic probe (standardise, C = 0.1, class_weight balanced) trained on ALL KU-Optofil cells; threshold
FIXED at the KU cross-fitted 95%-specificity threshold t = 0.0310392350137264. Applied unchanged to RV-PBS crops.
Uncertainty: bootstrap clustered on RV-PBS source IMAGE (2,000 resamples).
"""
import json, os, sys
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
sys.path.insert(0, "src")
from cadpbs.weighted_betting import betting_test, sampler_prop_side, sampler_uniform

T = 0.0310392350137264; KU_RECALL = 0.9059; KU_FP = 0.05; OUT = "results/external_v1"
Eku = np.load("results/embeddings_v1/embeddings_dinov2_vits14.npy").astype(np.float32)
ku = pd.read_csv("results/embeddings_v1/cells_features.csv")
Erv = np.load(f"{OUT}/rvpbs_embeddings_dinov2_vits14.npy").astype(np.float32)
rv = pd.read_csv(f"{OUT}/rvpbs_features.csv")
clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, C=0.1, class_weight="balanced")).fit(Eku, ku.is_blast.values)
rv["score"] = clf.predict_proba(Erv)[:, 1]; rv["call"] = rv.score >= T
rng = np.random.default_rng(20260521); imgs = rv.image.unique(); idx = {g: np.flatnonzero(rv.image.values == g) for g in imgs}
def boot(stat):
    v = []
    for _ in range(2000):
        ii = np.concatenate([idx[g] for g in rng.choice(imgs, len(imgs), replace=True)]); v.append(stat(rv.iloc[ii]))
    v = np.array(v, float); v = v[np.isfinite(v)]; return float(np.quantile(v, .025)), float(np.quantile(v, .975))
rec = lambda d: d.call[d.is_blast == 1].mean() if (d.is_blast == 1).any() else np.nan
fp = lambda d: d.call[d.is_blast == 0].mean() if (d.is_blast == 0).any() else np.nan
res = {}
r, ci = rec(rv), boot(rec); res["H4a_recall_drops"] = dict(rv_recall=float(r), ci95=ci, ku_recall=KU_RECALL, supported=bool(ci[1] < KU_RECALL))
f, ci = fp(rv), boot(fp); res["H4b_fp_rises"] = dict(rv_fp=float(f), ci95=ci, ku_fp=KU_FP, supported=bool(ci[0] > KU_FP))
res["descriptive"] = dict(auc=float(roc_auc_score(rv.is_blast, rv.score)),
                          fp_by_label=rv[rv.is_blast == 0].groupby("label").call.mean().round(3).to_dict(),
                          n=int(len(rv)), n_blast=int(rv.is_blast.sum()), n_images=int(len(imgs)))
# H4c: certificate validity with SHIFTED side information. Specimens from RV-PBS crops, N=100 x m=5, theta = 0.10.
BL = np.flatnonzero(rv.is_blast.values == 1); BG = np.flatnonzero(rv.is_blast.values == 0)
yb = rv.is_blast.values.astype(float); sc = rv.score.values; th = 0.10; rows = []
for p_true, reps in [(0.10, 1000), (0.15, 300), (0.20, 300)]:
    nb = int(round(p_true * 500))
    for ep in range(reps):
        cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(BG, 500 - nb, replace=False)]; rng.shuffle(cells)
        cells = cells.reshape(100, 5); f_ = yb[cells].mean(1) / (2 * th); S = sc[cells].mean(1) / (2 * th)
        row = dict(p=p_true)
        for nm, smp in [("uniform", sampler_uniform(100)), ("prop", sampler_prop_side(S))]:
            o = betting_test(f_, S, smp, 0.5, rng=rng); row[f"{nm}_cert"] = int(o["stop"] > 0); row[f"{nm}_cells"] = (o["stop"] if o["stop"] > 0 else 100) * 5
        rows.append(row)
c = pd.DataFrame(rows); c.to_csv(f"{OUT}/h4c_episodes.csv", index=False)
at = c[c.p == 0.10]; x = int(at.prop_cert.sum()); n = len(at)
from scipy.stats import binomtest
pv = binomtest(x, n, 0.05, alternative="greater").pvalue
res["H4c_validity_under_shift"] = dict(false_cert=x, n=n, rate=x / n, p_excess=float(pv), supported=bool(x / n <= 0.05 and pv >= 0.05))
res["H4c_cost_descriptive"] = {str(pt): dict(uniform=float(g.uniform_cells.median()), prop=float(g.prop_cells.median())) for pt, g in c[c.p > 0.10].groupby("p")}
json.dump(res, open(f"{OUT}/decisions.json", "w"), indent=1, default=float)
for k in ["H4a_recall_drops", "H4b_fp_rises", "H4c_validity_under_shift"]: print(k, "SUPPORTED" if res[k]["supported"] else "NOT SUPPORTED")
print(json.dumps(res, indent=1, default=float))
