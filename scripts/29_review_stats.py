"""
Additional statistics from stored outputs (no new data): fold-level intervals for the Part I contrasts,
critical-evidence recall at the stop point, Part III cost ratios with bootstrap intervals, sensitivity analyses,
and the field-level-teacher kill-test.

Output: results/review_stats_v1.json
"""
import json, os
import numpy as np, pandas as pd
from scipy import stats

T4 = stats.t.ppf(0.975, 4)
rng = np.random.default_rng(20261019)
out = {}


def fold_t(values):
    v = np.asarray(values, float); m = v.mean(); se = v.std(ddof=1) / np.sqrt(len(v))
    return dict(mean=m, ci95=[m - T4 * se, m + T4 * se], fold_estimates=v.round(4).tolist())


def dlt(g):
    return g.infodpp_miss.mean() - g["uniform@info_miss"].mean()


d = pd.read_csv("results/confirm_v1/episodes.csv")
out["H1b_fold_t"] = {int(k): fold_t([dlt(g[(g.m == 20) & (g.k == k)]) - dlt(g[(g.m == 1) & (g.k == k)]) for _, g in d.groupby("fold")])
                     for k in [1, 2, 5, 10]}
out["H1c_fold_t"] = {int(m): fold_t([dlt(g[(g.m == m) & (g.k == 1)]) for _, g in d.groupby("fold")]) for m in [10, 20]}
rec = []
for (m, k), g in d.groupby(["m", "k"]):
    rec.append(dict(m=int(m), k=int(k), n=len(g), median_cost=float(g.infodpp_n.median()),
                    **{f"{c}_recall": float(g[f"{c}_recall"].mean()) for c in ["infodpp", "uniform@info", "confidence"]}))
out["critical_evidence_recall"] = rec

eff = []
for f in ["results/part3_theta010/episodes.csv", "results/part3_theta020/episodes.csv"]:
    e = pd.read_csv(f); th = float(e.theta.iloc[0])
    for p, g in e[e.p > th + 1e-9].groupby("p"):
        x = {c: g[f"{c}_cells"].values for c in ["alpha_uniform", "bet_uniform", "bet_prop"]}; n = len(g)
        def ratio(a, b):
            bs = [np.median(x[a][i]) / np.median(x[b][i]) for i in (rng.integers(0, n, n) for _ in range(2000))]
            return dict(ratio=float(np.median(x[a]) / np.median(x[b])), ci95=np.quantile(bs, [.025, .975]).tolist())
        eff.append(dict(theta=th, p=float(p), medians={k: float(np.median(v)) for k, v in x.items()},
                        prop_vs_alpha=ratio("bet_prop", "alpha_uniform"), sampler_only=ratio("bet_prop", "bet_uniform"),
                        test_only=ratio("bet_uniform", "alpha_uniform")))
out["part3_effect_sizes"] = eff

sens = {}
for run in ["sensconf_pool_max", "sensconf_pool_single", "sensconf_pool_host", "sensconf_tau_0.05", "sensconf_tau_0.10",
            "sensconf_tau_0.163", "sensconf_tau_0.30", "sensconf_tau_0.50"]:
    p = f"results/{run}/episodes.csv"
    if not os.path.exists(p): continue
    s = pd.read_csv(p); sens[run] = {}
    for m, g in s[s.k == 1].groupby("m"):
        sens[run][int(m)] = dict(n=len(g), median_cost=float(g.infodpp_n.median()), infodpp_miss=float(g.infodpp_miss.mean()),
                                 uniform_miss=float(g["uniform@info_miss"].mean()),
                                 diff=fold_t([dlt(x) for _, x in g.groupby("fold")]),
                                 auc_info=float(s[s.m == m].kt1_auc_info_step1.mean()))
out["sensitivity_confirmatory"] = sens

if os.path.exists("results/kt1_field_teacher/episodes.csv"):
    k1 = pd.read_csv("results/kt1_field_teacher/episodes.csv", dtype={"blast_patients": str}); kt = {}
    for m, g in k1.groupby("m"):
        r = {c: fold_t(g.groupby("fold")[c].mean()) for c in ["auc_info", "auc_rel_cellmax", "auc_rel_field"]}
        h = g[g.k == 1]; diff = (h.infodpp_fieldteacher_miss - h.uniform_miss).values; pats = h.blast_patients.values
        up = np.unique(pats); byp = {u: diff[pats == u] for u in up}
        bs = [np.concatenate([byp[u] for u in rng.choice(up, len(up))]).mean() for _ in range(2000)]
        r.update(k1_n=len(h), k1_infodpp_miss=float(h.infodpp_fieldteacher_miss.mean()), k1_uniform_miss=float(h.uniform_miss.mean()),
                 k1_ranked_field_miss=float(h.ranked_field_miss.mean()), k1_diff=float(diff.mean()),
                 k1_diff_patient_boot95=np.quantile(bs, [.025, .975]).tolist(), n_patients=int(len(up)))
        kt[int(m)] = r
    out["kt1_field_teacher"] = kt

json.dump(out, open("results/review_stats_v1.json", "w"), indent=1, default=float)
print("written results/review_stats_v1.json")
