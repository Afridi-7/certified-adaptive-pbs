"""
Confirmatory analysis for Part I / Part II, implementing EXACTLY the decision rules of docs/PREREG.md
(v1, as amended by Amendment 1). Written and committed BEFORE the confirmatory data exist.

Input : results/confirm_v1/episodes.csv  (produced by scripts/05_killtests.py --seed 20260401 ...)
Output: results/confirm_v1/decisions.json and a printed verdict per hypothesis.
"""
import argparse, json
import numpy as np, pandas as pd
from scipy.stats import kruskal, binomtest, norm

ap = argparse.ArgumentParser(); ap.add_argument("--dir", default="results/confirm_v1"); ap.add_argument("--N", type=int, default=200)
a = ap.parse_args()
d = pd.read_csv(f"{a.dir}/episodes.csv")
assert sorted(d.fold.unique()) == [0, 1, 2, 3, 4], f"confirmatory file must contain all 5 folds, got {sorted(d.fold.unique())}"
rng = np.random.default_rng(20260401)
out = {}

# ---- H1a: median stop range across k <= 2 fields within each m ------------------------------------------
h1a = {}
for m, g in d.groupby("m"):
    med = g.groupby("k").infodpp_n.median(); rngw = float(med.max() - med.min())
    p = float(kruskal(*[x.infodpp_n.values for _, x in g.groupby("k")])[1])
    h1a[int(m)] = dict(median_by_k=med.to_dict(), range=rngw, kruskal_p=p, supported=rngw <= 2)
out["H1a"] = dict(cells=h1a, supported=all(v["supported"] for v in h1a.values()))

# ---- H1b: Delta(20,k) - Delta(1,k) > 0, fold-clustered bootstrap, k in {1,2,5,10} (A2); k=20 descriptive ---
def delta(g, m, k):
    x = g[(g.m == m) & (g.k == k)]
    return x.infodpp_miss.mean() - x["uniform@info_miss"].mean()
folds = np.sort(d.fold.unique()); h1b = {}
for k in [1, 2, 5, 10, 20]:
    point = delta(d, 20, k) - delta(d, 1, k)
    boots = []
    for _ in range(2000):
        take = rng.choice(folds, len(folds), replace=True)
        bd = pd.concat([d[d.fold == f] for f in take])
        boots.append(delta(bd, 20, k) - delta(bd, 1, k))
    lo, hi = np.quantile(boots, [0.025, 0.975])
    h1b[int(k)] = dict(delta_m1=delta(d, 1, k), delta_m20=delta(d, 20, k), contrast=point, ci95=[float(lo), float(hi)],
                       positive=bool(lo > 0), tested=k != 20)
n_pos = sum(v["positive"] for k, v in h1b.items() if v["tested"])
out["H1b"] = dict(cells=h1b, n_positive_of_4=n_pos, supported=n_pos >= 3)

# ---- H1c: equivalence within +-0.05, paired, k=1, m in {10,20} -------------------------------------------
h1c = {}
for m in [10, 20]:
    g = d[(d.m == m) & (d.k == 1)]; n = len(g)
    b = int(((g.infodpp_miss == 1) & (g["uniform@info_miss"] == 0)).sum())
    c = int(((g.infodpp_miss == 0) & (g["uniform@info_miss"] == 1)).sum())
    diff = (b - c) / n; var = (b + c) / n**2 - (b - c)**2 / n**3
    hw = 1.96 * np.sqrt(max(var, 0)); lo, hi = diff - hw, diff + hw
    h1c[int(m)] = dict(episodes=n, b=b, c=c, diff=diff, ci95=[lo, hi],
                      mcnemar_p=float(binomtest(b, b + c, 0.5).pvalue) if b + c else 1.0,
                      equivalent=bool(lo >= -0.05 and hi <= 0.05))
out["H1c"] = dict(cells=h1c, supported=all(v["equivalent"] for v in h1c.values()))

# ---- H1d: uniform miss at k=1 vs per-episode 1 - n_i/N; |z| < 2.50 (Bonferroni over 4 m) (A1) -------------
h1d = {}
for m, g in d[d.k == 1].groupby("m"):
    p_i = 1 - g.infodpp_n.values / a.N; emp = g["uniform@info_miss"].mean(); pred = p_i.mean()
    se = np.sqrt((p_i * (1 - p_i)).sum()) / len(g); z = (emp - pred) / se
    h1d[int(m)] = dict(empirical=emp, predicted=pred, z=float(z), agree=bool(abs(z) < 2.50))
out["H1d"] = dict(cells=h1d, supported=all(v["agree"] for v in h1d.values()))

# ---- KT1: AUC(info) < AUC(relevance) at every m; AUC(info) at m=20 in [0.45, 0.55] ---------------------------
kt = d.groupby("m")[["kt1_auc_info_step1", "kt1_auc_relevance"]].mean()
out["KT1"] = dict(by_m=kt.round(4).to_dict(),
                  supported=bool((kt.kt1_auc_info_step1 < kt.kt1_auc_relevance).all()
                                 and 0.45 <= kt.loc[20, "kt1_auc_info_step1"] <= 0.55))

json.dump(out, open(f"{a.dir}/decisions.json", "w"), indent=1, default=float)
for h in ["H1a", "H1b", "H1c", "H1d", "KT1"]:
    print(f"{h}: {'SUPPORTED' if out[h]['supported'] else 'NOT SUPPORTED'}")
print(json.dumps(out, indent=1, default=float)[:4000])
