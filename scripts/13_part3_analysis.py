"""Part III confirmatory decisions, implementing PREREG amendment 2 exactly. Written before the data exist."""
import argparse, json
import numpy as np, pandas as pd
from scipy.stats import beta, binomtest, wilcoxon

ap = argparse.ArgumentParser(); ap.add_argument("--dirs", nargs="+", required=True); a = ap.parse_args()
VALID = ["alpha_uniform", "bet_uniform", "bet_prop", "bet_prop_cv"]
def cp(x, n):  # Clopper-Pearson 95%
    return (0.0 if x == 0 else beta.ppf(0.025, x, n - x + 1), 1.0 if x == n else beta.ppf(0.975, x + 1, n - x))
out = {}
for d in a.dirs:
    e = pd.read_csv(f"{d}/episodes.csv"); th = float(e.theta.iloc[0]); key = f"theta_{th:.2f}"; res = {}
    at = e[np.isclose(e.p, th)]; n = len(at)
    # H3a validity: each valid design, false certification at p = theta: point <= 0.05 and one-sided exact test of
    # H0: rate <= 0.05 not rejected at 0.05.
    h3a = {}
    for dsg in VALID:
        x = int(at[f"{dsg}_cert"].sum()); pv = binomtest(x, n, 0.05, alternative="greater").pvalue
        h3a[dsg] = dict(n=n, false_cert=x, rate=x / n, cp95=cp(x, n), p_excess=pv, ok=bool(x / n <= 0.05 and pv >= 0.05))
    res["H3a_validity"] = dict(cells=h3a, supported=all(v["ok"] for v in h3a.values()))
    # H3b naive design is invalid: CP lower bound of false-certification rate at p = theta exceeds 0.05
    x = int(at.naive_cert.sum()); lo, hi = cp(x, n)
    res["H3b_naive_invalid"] = dict(rate=x / n, cp95=(lo, hi), supported=bool(lo > 0.05))
    # H3c classifier-driven sampling is cheaper than uniform ALPHA: paired Wilcoxon on cells, per p > theta, Holm
    above = sorted(e.p[e.p > th + 1e-9].unique()); tests = []
    for pv_ in above:
        g = e[np.isclose(e.p, pv_)]
        stat = wilcoxon(g.bet_prop_cells, g.alpha_uniform_cells, alternative="less", zero_method="wilcox")
        tests.append(dict(p=pv_, median_bet_prop=float(g.bet_prop_cells.median()), median_alpha_uniform=float(g.alpha_uniform_cells.median()),
                          median_ratio=float(g.bet_prop_cells.median() / g.alpha_uniform_cells.median()), p_value=float(stat.pvalue)))
    order = np.argsort([t["p_value"] for t in tests]); k = len(tests); rej = [False] * k
    for rank, i in enumerate(order):
        if tests[i]["p_value"] <= 0.05 / (k - rank): rej[i] = True
        else: break
    for i, t in enumerate(tests): t["holm_reject"] = rej[i]
    res["H3c_adaptive_cheaper"] = dict(cells=tests, supported=all(rej))
    # H3d cost grows as p approaches theta: median uniform-ALPHA cells strictly decreasing in p above theta
    meds = [float(e[np.isclose(e.p, pv_)].alpha_uniform_cells.median()) for pv_ in above]
    res["H3d_margin"] = dict(p=above, median_cells=meds, supported=bool(all(np.diff(meds) < 0)))
    # descriptive: conventional 100-cell count, rate of calling "above" per p
    res["conventional_100_call_above"] = e.groupby("p").conv100_call_above.mean().round(4).to_dict()
    out[key] = res
json.dump(out, open("results/part3_decisions.json", "w"), indent=1, default=float)
for k, r in out.items():
    print(k, {h: ("SUPPORTED" if v["supported"] else "NOT SUPPORTED") for h, v in r.items() if isinstance(v, dict) and "supported" in v})
