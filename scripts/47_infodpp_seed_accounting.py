"""
InfoDPP-PAC seed-set accounting at the 2.4% operating point (RQ1; exploratory).

At this operating point the stop is close to the seed-set size (m0 = 5). Compares the paper's accounting
(seed fields not counted as examined) with counting them, uniform sampling at n + |D0| on the same permutation.
The episode loop reproduces 05_killtests.py's random stream (verify with --check).

Usage: python scripts/47_infodpp_seed_accounting.py --out results/infodpp_seed_v1 --seed 20261022 --episodes 200
       --k_list 1 2 5 --target_frac 0.024 [--check results/sens_op_paperpool/episodes.csv]
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.composites import CellBank, make_specimen
from cadpbs.infodpp import InfoDPPConfig, run_infodpp, calibrate_tau

ap = argparse.ArgumentParser()
ap.add_argument("--emb", default="results/embeddings_v1")
ap.add_argument("--out", default="results/infodpp_seed_v1")
ap.add_argument("--N", type=int, default=200)
ap.add_argument("--k_list", type=int, nargs="+", default=[1, 2, 5])
ap.add_argument("--m_list", type=int, nargs="+", default=[1, 5, 10, 20])
ap.add_argument("--folds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
ap.add_argument("--episodes", type=int, default=200)
ap.add_argument("--pool", default="pool", choices=["pool", "host"])
ap.add_argument("--pooling", default="mean", choices=["mean", "max", "single"])
ap.add_argument("--target_frac", type=float, default=0.024)
ap.add_argument("--n_calib", type=int, default=8)
ap.add_argument("--tau_conf", type=float, default=0.9)
ap.add_argument("--m_conf", type=int, default=30)
ap.add_argument("--seed", type=int, default=20261022)
ap.add_argument("--check", default=None, help="stored episodes.csv of the same run, to confirm identical episodes")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

E = np.load(f"{a.emb}/embeddings_dinov2_vits14.npy").astype(np.float32)
df = pd.read_csv(f"{a.emb}/cells_features.csv", dtype={"patient_id": str})
bank = CellBank(E, df)
base_cfg = InfoDPPConfig(tau=None, seed=a.seed)


def confidence_stop(order, teacher, tau, mrun):                     # as in scripts/05 (kept for the RNG stream)
    run = 0
    for t, i in enumerate(order):
        if teacher[i] >= tau: return t + 1
        run = run + 1 if teacher[i] <= 1 - tau else 0
        if run >= mrun: return t + 1
    return len(order)


def calibrate(fold, m, rng):
    cal = []
    for f in [x for x in np.unique(bank.fold) if x != fold][:3]:
        for _ in range(max(1, a.n_calib // 3)):
            s = make_specimen(bank, f, a.N, 0, m, rng, a.pool, a.pooling)
            cal.append((s["F"], s["teacher"]))
    return calibrate_tau(cal, base_cfg, a.target_frac, rng=rng)


rows = []
for fold in a.folds:
    for m in a.m_list:
        rng = np.random.default_rng(a.seed + 1000 * fold + m)
        t0 = time.time(); tau = calibrate(fold, m, rng)
        cfg = InfoDPPConfig(**{**base_cfg.__dict__, "tau": tau})
        for k in a.k_list:
            for ep in range(a.episodes):
                s = make_specimen(bank, fold, a.N, k, m, rng, a.pool, a.pooling)
                res = run_infodpp(s["F"], s["teacher"], cfg, np.random.default_rng(rng.integers(1 << 30)))
                n = res["k_star"]; ci = np.flatnonzero(s["crit"]); D0 = np.asarray(res["seed_idx"]); m0 = len(D0)
                uni = rng.permutation(a.N)
                n_conf = confidence_stop(uni, s["teacher"], a.tau_conf, a.m_conf); _ = rng.permutation(a.N)
                S = np.asarray(res["order"][:n], int)
                miss = lambda examined: int(not np.isin(ci, examined).any())
                rows.append(dict(fold=fold, m=m, k=k, ep=ep, N=a.N, tau=tau, infodpp_n=n, conf_n=n_conf, m0=m0,
                                 seed_reselected=int(np.isin(S, D0).sum()), seed_has_crit=int(np.isin(ci, D0).any()),
                                 infodpp_miss=miss(S), uniform_at_n_miss=miss(uni[:n]),
                                 infodpp_seedincl_miss=miss(np.union1d(S, D0)),
                                 uniform_at_n_plus_m0_miss=miss(uni[:n + m0]),
                                 n_examined_seedincl=int(len(np.union1d(S, D0)))))
        print(f"fold {fold} m {m}: tau={tau:.3f} in {time.time() - t0:.0f}s", flush=True)

D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)


def wilson(x, n, z=1.96):
    p = x / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return [(c - h) / d, (c + h) / d]


def fold_t(d, col_a, col_b):
    """difference of miss rates a - b with a t interval on the five fold estimates (as in the paper)."""
    from math import sqrt
    per = d.groupby("fold").apply(lambda g: g[col_a].mean() - g[col_b].mean())
    k = len(per); se = per.std(ddof=1) / sqrt(k); tq = {5: 2.776, 4: 3.182, 3: 4.303}.get(k, 2.0)
    return float(per.mean()), [float(per.mean() - tq * se), float(per.mean() + tq * se)]


summ = {}
for (m, k), g in D.groupby(["m", "k"]):
    n_ = len(g); key = f"m{m}_k{k}"
    summ[key] = dict(episodes=n_, selections_mean=float(g.infodpp_n.mean()), selections_median=float(g.infodpp_n.median()),
                     selections_min=int(g.infodpp_n.min()), selections_max=int(g.infodpp_n.max()), m0=int(g.m0.iloc[0]),
                     share_seed_reselected=float((g.seed_reselected > 0).mean()), share_seed_has_crit=float(g.seed_has_crit.mean()),
                     infodpp_miss=float(g.infodpp_miss.mean()), uniform_at_n_miss=float(g.uniform_at_n_miss.mean()),
                     diff_seed_excluded=fold_t(g, "infodpp_miss", "uniform_at_n_miss"),
                     infodpp_seedincl_miss=float(g.infodpp_seedincl_miss.mean()),
                     infodpp_seedincl_miss_wilson=wilson(g.infodpp_seedincl_miss.sum(), n_),
                     uniform_at_n_plus_m0_miss=float(g.uniform_at_n_plus_m0_miss.mean()),
                     uniform_at_n_plus_m0_miss_wilson=wilson(g.uniform_at_n_plus_m0_miss.sum(), n_),
                     diff_seed_included=fold_t(g, "infodpp_seedincl_miss", "uniform_at_n_plus_m0_miss"))
check = None
if a.check:
    old = pd.read_csv(a.check); old = old[old.k.isin(a.k_list) & old.m.isin(a.m_list) & old.fold.isin(a.folds)]
    key = ["fold", "m", "k", "ep"]; mrg = D.merge(old[key + ["infodpp_n", "infodpp_miss", "uniform@info_miss"]], on=key, suffixes=("", "_old"))
    check = dict(stored_episodes=int(len(old)), matched=int(len(mrg)),
                 identical_n=bool((mrg.infodpp_n == mrg.infodpp_n_old).all()),
                 identical_infodpp_miss=bool((mrg.infodpp_miss == mrg.infodpp_miss_old).all()),
                 identical_uniform_miss=bool((mrg.uniform_at_n_miss == mrg["uniform@info_miss"]).all()))
json.dump(dict(summary=summ), open(f"{a.out}/summary.json", "w"), indent=1)
if check is not None: json.dump(check, open(f"{a.out}/check.json", "w"), indent=1)
json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print(json.dumps(dict(check=check), indent=1))
for key, v in summ.items():
    print(key, {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items()
                if kk in ("selections_mean", "selections_median", "selections_min", "selections_max", "infodpp_miss",
                          "uniform_at_n_miss", "infodpp_seedincl_miss", "uniform_at_n_plus_m0_miss", "share_seed_has_crit")})
