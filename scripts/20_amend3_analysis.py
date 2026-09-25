"""
PREREG Amendment 3 decision rules (H5a-H5c for E3.1; H6a-H6d for E3.2). Frozen with the amendment, before any data.
Usage: python scripts/20_amend3_analysis.py --spatial results/a3_spatial_theta010 results/a3_spatial_theta020 \
                                            --mixture results/a3_mixture_rvpbs results/a3_mixture_ku --out results/a3_decisions.json
"""
import argparse, glob, json
import numpy as np, pandas as pd
from scipy.stats import binomtest, beta, wilcoxon, mannwhitneyu

ap = argparse.ArgumentParser()
ap.add_argument("--spatial", nargs="+", required=True); ap.add_argument("--mixture", nargs="+", required=True)
ap.add_argument("--out", default="results/a3_decisions.json"); ap.add_argument("--boot_seed", type=int, default=20260631)
a = ap.parse_args()
ALPHA = 0.05; res = {}


def load(dirs, expect):
    d = pd.concat([pd.read_csv(f) for x in dirs for f in sorted(glob.glob(f"{x}/ep_*.csv"))], ignore_index=True)
    key = [c for c in ["theta", "arr", "source"] if c in d.columns] + ["p"]
    cnt = d.groupby(key).ep.nunique()
    for k, n in cnt.items():
        need = expect(k)
        assert n == need, f"incomplete: {k} has {n} episodes, expected {need}"
    assert not d.duplicated(key + ["ep"]).any(), "duplicate episodes"
    return d


def cp(x, n):
    lo = beta.ppf(ALPHA / 2, x, n - x + 1) if x > 0 else 0.0
    hi = beta.ppf(1 - ALPHA / 2, x + 1, n - x) if x < n else 1.0
    return float(lo), float(hi)


def valid_rule(x, n):     # the H3a rule of Amendment 2
    pv = binomtest(x, n, ALPHA, alternative="greater").pvalue
    return dict(x=int(x), n=int(n), rate=x / n, ci95=cp(x, n), p_excess=float(pv), pass_=bool(x / n <= ALPHA and pv >= 0.05))


def holm(pvals):
    order = np.argsort(pvals); k = len(pvals); rej = np.zeros(k, bool)
    for r, i in enumerate(order):
        if pvals[i] <= ALPHA / (k - r): rej[i] = True
        else: break
    return rej

# ---------------- E3.1 spatial ----------------
S = load(a.spatial, lambda k: 1000 if abs(k[2] - k[0]) < 1e-9 else 500)
thetas = sorted(S.theta.unique()); arrs = ["scattered", "hot0.6", "hot1.0"]
# H5a: the three field-level designs are valid in every arrangement, at p = theta and at p = theta - 0.02
h5a = {}; ok = True
for th in thetas:
    for arr in arrs:
        for p in (th - 0.02, th):
            g = S[(S.theta == th) & (S.arr == arr) & (np.isclose(S.p, p))]
            for d in ["alpha_field", "bet_field", "prop_field"]:
                v = valid_rule(g[f"{d}_cert"].sum(), len(g)); h5a[f"theta{th}_{arr}_p{p:.2f}_{d}"] = v; ok &= v["pass_"]
res["H5a_field_level_valid_all_arrangements"] = dict(cells=h5a, supported=bool(ok))
# H5b: naive cell-level analysis is valid when scattered, invalid (CP lower bound > 0.05) under hot1.0, at p = theta
h5b = {}; ok = True
for th in thetas:
    g = S[(S.theta == th) & np.isclose(S.p, th)]
    sc = g[g.arr == "scattered"]; v1 = valid_rule(sc.naive_cell_cert.sum(), len(sc))
    hc = g[g.arr == "hot1.0"]; x, n = int(hc.naive_cell_cert.sum()), len(hc); lo, hi = cp(x, n)
    md = g[g.arr == "hot0.6"]
    h5b[f"theta{th}"] = dict(scattered=v1, hot1_0=dict(x=x, n=n, rate=x / n, ci95=(lo, hi), invalid=bool(lo > ALPHA)),
                             hot0_6_descriptive=dict(rate=float(md.naive_cell_cert.mean()), ci95=cp(int(md.naive_cell_cert.sum()), len(md))))
    ok &= v1["pass_"] and lo > ALPHA
res["H5b_naive_cell_analysis_breaks_under_clustering"] = dict(cells=h5b, supported=bool(ok))
# H5c: alpha_field cost at p = theta + 0.05 increases with clustering: medians strictly ordered and
#      scattered < hot1.0 by one-sided Mann-Whitney (different specimens) at 0.05, for both thetas
h5c = {}; ok = True
for th in thetas:
    p = th + 0.05; med = {}
    for arr in arrs:
        med[arr] = float(S[(S.theta == th) & (S.arr == arr) & np.isclose(S.p, p)].alpha_field_cells.median())
    x = S[(S.theta == th) & (S.arr == "scattered") & np.isclose(S.p, p)].alpha_field_cells
    z = S[(S.theta == th) & (S.arr == "hot1.0") & np.isclose(S.p, p)].alpha_field_cells
    pv = float(mannwhitneyu(x, z, alternative="less").pvalue)
    strict = med["scattered"] < med["hot0.6"] < med["hot1.0"]
    h5c[f"theta{th}_p{p:.2f}"] = dict(medians=med, mw_p=pv, strict=bool(strict)); ok &= strict and pv < ALPHA
res["H5c_cost_rises_with_clustering"] = dict(cells=h5c, supported=bool(ok))
# descriptive
desc = []
for (th, arr, p), g in S.groupby(["theta", "arr", "p"]):
    row = dict(theta=th, arr=arr, p=p, n=len(g), mean_fields_with_blast=float(g.n_hot.mean()))
    for d in ["alpha_field", "bet_field", "prop_field", "naive_cell"]:
        row[f"{d}_cert_rate"] = float(g[f"{d}_cert"].mean()); row[f"{d}_median_cells"] = float(g[f"{d}_cells"].median())
    row["conv100_call_above"] = float(g.conv100_call_above.mean()); desc.append(row)
res["E31_descriptive"] = desc

# ---------------- E3.2 mixture ----------------
X = load(a.mixture, lambda k: 1000 if abs(k[1] - 0.10) < 1e-9 else 500)
# H6a: mix25 and mix50 valid at p = theta = 0.10 on both sources (H3a rule)
h6a = {}; ok = True
for src in ["rvpbs", "ku"]:
    g = X[(X.source == src) & np.isclose(X.p, 0.10)]
    for d in ["mix25", "mix50", "uniform", "prop"]:
        v = valid_rule(g[f"{d}_cert"].sum(), len(g)); h6a[f"{src}_{d}"] = v
        if d.startswith("mix"): ok &= v["pass_"]
res["H6a_mixture_valid"] = dict(cells=h6a, supported=bool(ok))
# H6d: the reversal replicates on fresh episodes: RV-PBS p = 0.20, prop cells > uniform cells (paired one-sided Wilcoxon)
g = X[(X.source == "rvpbs") & np.isclose(X.p, 0.20)]
pv = float(wilcoxon(g.prop_cells, g.uniform_cells, alternative="greater", zero_method="wilcox").pvalue)
res["H6d_reversal_replicates"] = dict(median_prop=float(g.prop_cells.median()), median_uniform=float(g.uniform_cells.median()),
                                      wilcoxon_p=pv, supported=bool(pv < ALPHA))
# H6b: under shift the mixture's loss vs uniform is < 20%: upper 95% bootstrap bound of median(mix50)/median(uniform) < 1.20,
#      at p = 0.15 and p = 0.20 (episode bootstrap, 2,000 resamples, paired rows)
rng = np.random.default_rng(a.boot_seed); h6b = {}; ok = True
for p in (0.15, 0.20):
    g = X[(X.source == "rvpbs") & np.isclose(X.p, p)]; mu, mx = g.uniform_cells.values, g.mix50_cells.values
    rat = []
    for _ in range(2000):
        ii = rng.integers(0, len(g), len(g)); rat.append(np.median(mx[ii]) / np.median(mu[ii]))
    hi = float(np.quantile(rat, 0.975)); est = float(np.median(mx) / np.median(mu))
    h6b[f"p{p:.2f}"] = dict(ratio=est, ci95=(float(np.quantile(rat, 0.025)), hi),
                            mix25_ratio_descriptive=float(np.median(g.mix25_cells) / np.median(mu)), pass_=bool(hi < 1.20))
    ok &= hi < 1.20
res["H6b_mixture_caps_loss_under_shift"] = dict(cells=h6b, supported=bool(ok))
# H6c: in distribution (KU) the mixture keeps >= 50% of prop's saving at p in {0.12, 0.15, 0.20}:
#      retained = (med_u - med_mix50)/(med_u - med_prop) >= 0.5 AND mix50 < uniform by paired one-sided Wilcoxon, Holm over 3 p
h6c = {}; pvs = []; keys = []
for p in (0.12, 0.15, 0.20):
    g = X[(X.source == "ku") & np.isclose(X.p, p)]
    mu, mp, mm = g.uniform_cells.median(), g.prop_cells.median(), g.mix50_cells.median()
    ret = float((mu - mm) / (mu - mp)) if mu - mp > 0 else float("nan")
    pv = float(wilcoxon(g.mix50_cells, g.uniform_cells, alternative="less", zero_method="wilcox").pvalue)
    h6c[f"p{p:.2f}"] = dict(med_uniform=float(mu), med_prop=float(mp), med_mix50=float(mm), med_mix25=float(g.mix25_cells.median()),
                            retained=ret, wilcoxon_p=pv); pvs.append(pv); keys.append(f"p{p:.2f}")
rej = holm(np.array(pvs)); ok = True
for k, r in zip(keys, rej):
    h6c[k]["holm_reject"] = bool(r); h6c[k]["pass_"] = bool(r and np.isfinite(h6c[k]["retained"]) and h6c[k]["retained"] >= 0.5)
    ok &= h6c[k]["pass_"]
res["H6c_mixture_keeps_gain_in_distribution"] = dict(cells=h6c, supported=bool(ok))

json.dump(res, open(a.out, "w"), indent=1, default=float)
for k in ["H5a_field_level_valid_all_arrangements", "H5b_naive_cell_analysis_breaks_under_clustering", "H5c_cost_rises_with_clustering",
          "H6a_mixture_valid", "H6d_reversal_replicates", "H6b_mixture_caps_loss_under_shift", "H6c_mixture_keeps_gain_in_distribution"]:
    print(f"{k}: {'SUPPORTED' if res[k]['supported'] else 'NOT SUPPORTED'}")
