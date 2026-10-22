"""
Sensitivity of the Part I result to InfoDPP-PAC settings that the original paper does not report (exploratory). RQ1.

Our reimplementation (src/cadpbs/infodpp.py, frozen) fixes several choices the paper leaves open: the constant rho in the
concentration term of the halting rule (rho = 1), the RBF kernel length-scale (median heuristic), the observation-noise
variance (sn2 = 0.01) and the floor of the random seed set (m0 = 5). Here each choice is varied one at a time:
    rho in {0, 0.5, 2};  length-scale x {0.5, 2};  sn2 in {0.001, 0.1};  m0 in {2, 10}
and the halting threshold tau is re-calibrated for every setting so that the median stop on blast-free specimens of
other folds is 16.3% of the fields, exactly as in the confirmatory run. All settings are run on the SAME specimens and
the same uniform-sampling orders (paired), with one critical field (k = 1) and m in {1, 5, 10, 20} cells per field.
The frozen module is not changed: a different length-scale is passed through InfoDPPConfig.length_scale per specimen,
and tau is calibrated with the same bisection as calibrate_tau.

Question: do the conclusions survive these choices? (i) the advantage over uniform sampling of the same size shrinks
with dilution; (ii) it is small at m = 20; (iii) the information-gain term loses its association with the critical
field as m grows.

Each setting is also compared with the main setting on the same specimens: exact McNemar tests of the InfoDPP-PAC miss
indicator for every (setting, m), Holm-corrected over all of them.

Usage: python scripts/48_infodpp_settings.py --out results/infodpp_settings_v1 --seed 20261040 --episodes 100
       (--analyse_only recomputes the summary from an existing episodes.csv)
"""
import argparse, json, os, sys, time
from dataclasses import asdict
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.composites import CellBank, make_specimen
from cadpbs.infodpp import InfoDPPConfig, run_infodpp, median_length_scale
from scipy.stats import binomtest

ap = argparse.ArgumentParser()
ap.add_argument("--emb", default="results/embeddings_v1")
ap.add_argument("--out", default="results/infodpp_settings_v1")
ap.add_argument("--N", type=int, default=200)
ap.add_argument("--m_list", type=int, nargs="+", default=[1, 5, 10, 20])
ap.add_argument("--folds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
ap.add_argument("--episodes", type=int, default=100, help="episodes per (fold, m), k = 1")
ap.add_argument("--target_frac", type=float, default=0.163)
ap.add_argument("--n_calib", type=int, default=8)
ap.add_argument("--seed", type=int, default=20261040)
ap.add_argument("--analyse_only", action="store_true")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)

SETTINGS = {                       # name -> (InfoDPPConfig overrides, length-scale factor)
    "base":      ({}, 1.0),
    "rho0":      ({"rho": 0.0}, 1.0),
    "rho0.5":    ({"rho": 0.5}, 1.0),
    "rho2":      ({"rho": 2.0}, 1.0),
    "ell0.5":    ({}, 0.5),
    "ell2":      ({}, 2.0),
    "sn2_0.001": ({"sn2": 0.001}, 1.0),
    "sn2_0.1":   ({"sn2": 0.1}, 1.0),
    "m0_2":      ({"m0_min": 2}, 1.0),
    "m0_10":     ({"m0_min": 10}, 1.0),
}


def simulate():
    from sklearn.metrics import roc_auc_score
    E = np.load(f"{a.emb}/embeddings_dinov2_vits14.npy").astype(np.float32)
    df = pd.read_csv(f"{a.emb}/cells_features.csv", dtype={"patient_id": str})
    bank = CellBank(E, df)
    BASE = InfoDPPConfig(tau=None, seed=a.seed)

    def cfg_for(over, ell_factor, X, tau):
        c = {**asdict(BASE), **over, "tau": tau}
        if ell_factor != 1.0:
            c["length_scale"] = ell_factor * median_length_scale(X)
        return InfoDPPConfig(**c)

    def calibrate(specs, over, ell_factor, iters=25):
        """Same bisection as cadpbs.infodpp.calibrate_tau, with the length-scale set per specimen."""
        lo, hi = 1e-6, 5.0
        for _ in range(iters):
            mid = np.sqrt(lo * hi)
            ks = [run_infodpp(X, t, cfg_for(over, ell_factor, X, mid), np.random.default_rng(1000 + i))["k_star"]
                  for i, (X, t) in enumerate(specs)]
            if np.median(ks) > a.target_frac * len(specs[0][0]): lo = mid
            else: hi = mid
        return float(np.sqrt(lo * hi))

    rows, taus = [], []
    t0 = time.time()
    for fold in a.folds:
        for m in a.m_list:
            rng = np.random.default_rng([a.seed, fold, m])
            cal = []
            for f in [x for x in np.unique(bank.fold) if x != fold][:3]:
                for _ in range(max(1, a.n_calib // 3)):
                    s = make_specimen(bank, f, a.N, 0, m, rng); cal.append((s["F"], s["teacher"]))
            specs = []
            for ep in range(a.episodes):
                s = make_specimen(bank, fold, a.N, 1, m, rng)
                specs.append((s, rng.permutation(a.N), int(rng.integers(1 << 30))))
            for name, (over, ef) in SETTINGS.items():
                tau = calibrate(cal, over, ef); taus.append(dict(setting=name, fold=fold, m=m, tau=tau))
                for ep, (s, uni, sd) in enumerate(specs):
                    X = s["F"]; res = run_infodpp(X, s["teacher"], cfg_for(over, ef, X, tau), np.random.default_rng(sd))
                    n = res["k_star"]; ci = np.flatnonzero(s["crit"])
                    rows.append(dict(setting=name, fold=fold, m=m, ep=ep, tau=tau, infodpp_n=n,
                                     infodpp_miss=int(not np.isin(ci, res["order"][:n]).any()),
                                     uniform_miss=int(not np.isin(ci, uni[:n]).any()),
                                     auc_info=float(roc_auc_score(s["crit"], res["dinfo0"]))))
            print(f"fold {fold} m {m}: {time.time() - t0:.0f}s", flush=True)

    D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)
    pd.DataFrame(taus).to_csv(f"{a.out}/taus.csv", index=False)


    return D


D = pd.read_csv(f"{a.out}/episodes.csv") if a.analyse_only else simulate()


def fold_t(g):
    per = g.groupby("fold").apply(lambda x: x.infodpp_miss.mean() - x.uniform_miss.mean())
    k = len(per); se = per.std(ddof=1) / np.sqrt(k); tq = {5: 2.776, 4: 3.182, 3: 4.303}.get(k, 2.0)
    return float(per.mean()), [float(per.mean() - tq * se), float(per.mean() + tq * se)]


summ = []
for (name, m), g in D.groupby(["setting", "m"], sort=False):
    d, ci = fold_t(g)
    summ.append(dict(setting=name, m=m, episodes=len(g), median_stop=float(g.infodpp_n.median()),
                     infodpp_miss=float(g.infodpp_miss.mean()), uniform_miss=float(g.uniform_miss.mean()),
                     diff=d, diff_ci95_lo=ci[0], diff_ci95_hi=ci[1], auc_info=float(g.auc_info.mean())))
S = pd.DataFrame(summ); S.to_csv(f"{a.out}/summary.csv", index=False)
checks = {}
for name, g in S.groupby("setting", sort=False):
    g = g.set_index("m")
    checks[name] = dict(advantage_shrinks=bool(g.loc[20, "diff"] > g.loc[1, "diff"]),
                        diff_m1=float(g.loc[1, "diff"]), diff_m20=float(g.loc[20, "diff"]),
                        auc_m1=float(g.loc[1, "auc_info"]), auc_m20=float(g.loc[20, "auc_info"]),
                        miss_m20=float(g.loc[20, "infodpp_miss"]))

# paired comparison with the main setting: exact McNemar on the miss indicator, Holm over all (setting, m)
B = D[D.setting == "base"].set_index(["fold", "m", "ep"]).infodpp_miss
mc = []
for (name, m), g in D[D.setting != "base"].groupby(["setting", "m"], sort=False):
    x = g.set_index(["fold", "m", "ep"]).infodpp_miss; y = B.loc[x.index]
    b01, b10 = int(((x == 1) & (y == 0)).sum()), int(((x == 0) & (y == 1)).sum())
    mc.append(dict(setting=name, m=int(m), more_misses=b01, fewer_misses=b10,
                   p=float(binomtest(b01, b01 + b10).pvalue) if b01 + b10 else 1.0))
order = np.argsort([r["p"] for r in mc]); run = 0.0
for rank, i in enumerate(order):
    run = max(run, min(1.0, (len(mc) - rank) * mc[i]["p"])); mc[i]["p_holm"] = run
others = S[S.setting != "base"].merge(S[S.setting == "base"][["m", "diff", "infodpp_miss"]], on="m", suffixes=("", "_base"))
overall = dict(n_specimens=int(len(D[D.setting == "base"])),
               max_abs_change_miss=float((others.infodpp_miss - others.infodpp_miss_base).abs().max()),
               max_abs_change_diff=float((others["diff"] - others.diff_base).abs().max()),
               min_p_holm=float(min(r["p_holm"] for r in mc)),
               diff_m1_range=[float(S[S.m == 1]["diff"].min()), float(S[S.m == 1]["diff"].max())],
               diff_m20_range=[float(S[S.m == 20]["diff"].min()), float(S[S.m == 20]["diff"].max())],
               auc_m1_range=[float(S[S.m == 1].auc_info.min()), float(S[S.m == 1].auc_info.max())],
               auc_m20_range=[float(S[S.m == 20].auc_info.min()), float(S[S.m == 20].auc_info.max())],
               median_stop_range=[float(S.median_stop.min()), float(S.median_stop.max())],
               all_advantage_shrinks=bool(all(c["advantage_shrinks"] for c in checks.values())))
json.dump(dict(summary=S.to_dict(orient="records"), checks=checks, mcnemar_vs_base=mc, overall=overall),
          open(f"{a.out}/summary.json", "w"), indent=1)
if not a.analyse_only:
    json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
pd.set_option("display.width", 200)
print(S.round(3).to_string(index=False)); print(json.dumps(overall, indent=1))
