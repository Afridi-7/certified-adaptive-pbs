"""Analysis of results/killtests_v1: KT1 signal separability, H1a (content-blind stopping), KT2 gap vs
rarity k and dilution m. CIs are clustered on FOLD (each fold has its own blast patients), and the paired
InfoDPP-vs-uniform comparison uses the exact McNemar test at matched cost."""
import argparse, json
import numpy as np, pandas as pd
from scipy.stats import binomtest

ap = argparse.ArgumentParser(); ap.add_argument("--dir", default="results/killtests_v1"); a = ap.parse_args()
d = pd.read_csv(f"{a.dir}/episodes.csv")
out = {}

def wilson(x, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = x / n; den = 1 + z*z/n; c = (p + z*z/(2*n))/den
    h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den; return (max(0, c-h), min(1, c+h))

print("=== H1a: does the stop point depend on the specimen's content? ===")
h = d.groupby(["m", "k"]).infodpp_n.agg(["median", "mean", "std"]).round(2)
print(h.to_string())
# regression-free test: Kruskal-Wallis across k within each m
from scipy.stats import kruskal
for m, g in d.groupby("m"):
    groups = [gg.infodpp_n.values for _, gg in g.groupby("k")]
    st, p = kruskal(*groups)
    print(f"m={m:2d}: stop-step vs rarity k -> Kruskal-Wallis H={st:.2f}, p={p:.3f}, "
          f"median range {min(np.median(x) for x in groups):.0f}-{max(np.median(x) for x in groups):.0f}")
out["H1a_stop_vs_k"] = {int(m): float(kruskal(*[gg.infodpp_n.values for _, gg in g.groupby('k')])[1]) for m, g in d.groupby('m')}

print("\n=== KT1: are the three stopping signals separable? (AUC for critical vs non-critical fields) ===")
kt1 = d.groupby("m")[["kt1_auc_info_step1", "kt1_auc_relevance", "kt1_corr_info_relevance"]].agg(["mean", "std"]).round(3)
print(kt1.to_string())
out["KT1"] = d.groupby("m")[["kt1_auc_info_step1","kt1_auc_relevance","kt1_corr_info_relevance"]].mean().round(4).to_dict()

print("\n=== KT2: premature-stop rate at matched observation cost ===")
rows = []
for (m, k), g in d.groupby(["m", "k"]):
    n = len(g); mi = int(g["infodpp_miss"].sum()); un = int(g["uniform@info_miss"].sum())
    b = int(((g["infodpp_miss"] == 1) & (g["uniform@info_miss"] == 0)).sum())
    c = int(((g["infodpp_miss"] == 0) & (g["uniform@info_miss"] == 1)).sum())
    pv = binomtest(b, b + c, 0.5).pvalue if b + c else 1.0
    lo1, hi1 = wilson(mi, n); lo2, hi2 = wilson(un, n)
    # fold-clustered spread of the difference
    per_fold = g.groupby("fold").apply(lambda x: x["infodpp_miss"].mean() - x["uniform@info_miss"].mean())
    rows.append(dict(m=m, k=k, episodes=n, cost_median=g.infodpp_n.median(),
                     infodpp=f"{mi/n:.3f} [{lo1:.3f},{hi1:.3f}]", uniform=f"{un/n:.3f} [{lo2:.3f},{hi2:.3f}]",
                     diff=round(mi/n - un/n, 3), fold_range=f"[{per_fold.min():+.2f},{per_fold.max():+.2f}]",
                     mcnemar_b_c=f"{b}/{c}", p=f"{pv:.1e}",
                     recall_infodpp=round(g["infodpp_recall"].mean(), 3),
                     recall_uniform=round(g["uniform@info_recall"].mean(), 3),
                     conf_cost=g.conf_n.median(), conf_miss=round(g["confidence_miss"].mean(), 3)))
kt2 = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(kt2.to_string(index=False))
kt2.to_csv(f"{a.dir}/summary_kt2.csv", index=False)
json.dump(out, open(f"{a.dir}/summary.json", "w"), indent=1, default=float)
