"""
Patient-level analysis of the Part I real-tile run (P1-2; exploratory). RQ1, RQ2.

Input : results/partI_real_tiles_v1/runs.csv (scripts/36_partI_real_tiles.py, seed 20261028).
Unit of inference: the patient (= specimen). Each patient's miss rate is averaged over its 200 runs; differences between
rules are paired within patient; 95% intervals are percentile bootstraps over patients (4,000 resamples, seed 20261031).
k-center stops at its own cost, so it is compared with uniform sampling at that cost, h(N, k, n_kcenter) (Theorem 1).
Also: how often a rule's miss rate in a single specimen is at least 0.95 (near-deterministic failure), and the teacher's
tile-level AUC within the rare specimens.

Usage: python scripts/40_partI_real_tiles_analysis.py
"""
import json
from math import comb
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

OUT = "results/partI_real_tiles_v1"
D = pd.read_csv(f"{OUT}/runs.csv", dtype={"patient": str})
h = lambda N, k, n: comb(N - k, n) / comb(N, n) if n <= N - k else 0.0
D["h_kc"] = [h(N, k, n) for N, k, n in zip(D.N, D.k, D.kcenter_n)]
P = D.groupby("patient").agg(N=("N", "first"), k=("k", "first"), p_s=("p_s", "first"), fold=("fold", "first"),
                             n=("n", "mean"), kcenter_n=("kcenter_n", "mean"), infodpp=("infodpp", "mean"),
                             uniform=("uniform", "mean"), raster=("raster", "mean"), ranked=("ranked", "mean"),
                             kcenter=("kcenter", "mean"), h=("h_uniform", "mean"), h_kc=("h_kc", "mean"))
P["kN"] = P.k / P.N
P.to_csv(f"{OUT}/patient_summary.csv")
rng = np.random.default_rng(20261031)


def boot(g, a, b, B=4000):
    d = (g[a] - g[b]).values; m = d[rng.integers(0, len(d), (B, len(d)))].mean(1)
    return dict(diff=float(d.mean()), ci95=[float(np.quantile(m, .025)), float(np.quantile(m, .975))],
                n_pos=int((d > 0).sum()), n_neg=int((d < 0).sum()))


res = dict(n_patients_0_lt_k_lt_N=int(len(P)), N_values=P.N.value_counts().to_dict(),
           stop_fraction_mean=float((P.n / P.N).mean()), stop_n_mean=float(P.n.mean()), kcenter_n_mean=float(P.kcenter_n.mean()))
groups = {"all": P, "rare_kN_le_0.25": P[P.kN <= 0.25], "k_eq_1": P[P.k == 1], "mid_0.25_0.5": P[(P.kN > 0.25) & (P.kN <= 0.5)],
          "common_gt_0.5": P[P.kN > 0.5]}
for name, g in groups.items():
    r = dict(n_patients=int(len(g)), mean_miss={c: float(g[c].mean()) for c in ["infodpp", "uniform", "raster", "ranked", "h", "kcenter", "h_kc"]})
    for a, b in [("infodpp", "uniform"), ("raster", "uniform"), ("ranked", "uniform"), ("kcenter", "h_kc"), ("infodpp", "raster")]:
        r[f"{a}_minus_{b}"] = boot(g, a, b)
    r["n_specimens_miss_ge_0.95"] = {c: int((g[c] >= 0.95).sum()) for c in ["infodpp", "uniform", "raster", "ranked"]}
    r["uniform_miss_range"] = [float(g.uniform.min()), float(g.uniform.max())]
    res[name] = r
TS = pd.read_csv(f"{OUT}/tiles_scored.csv", dtype={"patient": str})
for name, g in [("all", P), ("rare_kN_le_0.25", P[P.kN <= 0.25])]:
    t = TS[TS.patient.isin(g.index)]
    res[name]["teacher_auc_tiles"] = float(roc_auc_score(t.crit, t.teacher))
    one = t[(t.crit == 0) | (t.n_blast == 1)]
    res[name]["teacher_auc_single_blast_tiles_vs_negative"] = float(roc_auc_score(one.crit, one.teacher))
res["note"] = "Ranked-scan misses are 0 or 1 per specimen because the ranking is deterministic."
json.dump(res, open(f"{OUT}/analysis.json", "w"), indent=1)
print(json.dumps(res, indent=1))
