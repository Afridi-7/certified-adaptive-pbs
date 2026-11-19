"""
Leave-one-blast-patient-out (exploratory; RQ1, RQ3). Review of 2026-11-19.

All 489 KU-Optofil blasts come from 11 patients, and five of them contribute 75% of the blasts. This script asks how
much the in-distribution results rest on single patients. For each blast-positive patient P (and once with nobody
left out) it
  1. removes every cell of P from the cell bank;
  2. refits the patient-grouped 5-fold probe of notebook 03 (standardised logistic regression, C = 0.1, balanced,
     same fold assignment) on the remaining cells, so P influences neither the specimens nor the classifier;
  3. Part III: reruns the paired design comparison of scripts/49_baseline_reference.py on KU-Optofil (N = 100 fields x
     5 cells; theta = 0.10 and 0.20; same designs, same seed family 20261101, so with nobody left out and the stored
     probe the run equals results/baseline_reference_v1);
  4. Part I: InfoDPP-PAC against uniform sampling at matched cost, one critical field, N = 200, m in {1, 5, 10, 20},
     settings of scripts/05_killtests.py; seed family 20261120.
Validity of the tests does not depend on which patients are present; the check is about costs and miss rates.

Usage: python scripts/58_leave_one_patient_out.py --out results/lopo_v1
"""
import argparse, json, os, sys, time, platform
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w, sampler_prop_w
from cadpbs.hedged_betting import hedged_betting_w
from cadpbs.composites import CellBank, make_specimen
from cadpbs.infodpp import InfoDPPConfig, run_infodpp, calibrate_tau
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

DESIGNS = ["alpha_uniform", "bet_uniform", "bet_prop", "hedge_uniform", "hedge_prop", "hedge_mix50"]
ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--seed3", type=int, default=20261101); ap.add_argument("--seed1", type=int, default=20261120)
ap.add_argument("--n_null", type=int, default=1000); ap.add_argument("--n_alt", type=int, default=400)
ap.add_argument("--episodes", type=int, default=100, help="Part I episodes per fold and m")
ap.add_argument("--m_list", type=int, nargs="+", default=[1, 5, 10, 20])
ap.add_argument("--only", nargs="*", default=None, help="subset of runs (patient ids or 'none')")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); T0 = time.time()

E_all = np.load(f"{ROOT}/results/embeddings_v1/embeddings_dinov2_vits14.npy").astype(np.float32)
ku_all = pd.read_csv(f"{ROOT}/results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})
blast_patients = sorted(ku_all.loc[ku_all.is_blast == 1, "patient_id"].unique())
RUNS = ["none", "none_refit"] + blast_patients            # 'none' = stored probe; 'none_refit' = refit, nobody left out
if a.only: RUNS = [r for r in RUNS if r in a.only]


def refit_oof(E, df):
    p = np.zeros(len(df)); y = df.is_blast.values
    for f in sorted(df.fold.unique()):
        tr, te = df.fold.values != f, df.fold.values == f
        clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=0.1, class_weight="balanced")).fit(E[tr], y[tr])
        p[te] = clf.predict_proba(E[te])[:, 1]
    return p


def part3(df, tag):
    y = df.is_blast.values.astype(bool); score = df.probe_p_oof.values
    bp = set(df.loc[y, "patient_id"]); BL = np.flatnonzero(y); BG = np.flatnonzero(~y & ~df.patient_id.isin(bp).values)
    N, m = 100, 5; M = N * m; pi = np.full(N, 1.0 / N)
    def mix(S, lam):
        base = sampler_prop_w(pi, S)
        return lambda rem, t: (1 - lam) * base(rem, t) + lam / len(rem)
    rows = []
    for ti, (th, plist) in enumerate([(0.10, [0.10, 0.12, 0.15, 0.20]), (0.20, [0.20, 0.22, 0.25, 0.30])]):
        for p in plist:
            nb = int(round(p * M)); assert nb <= len(BL) and M - nb <= len(BG)
            for ep in range(a.n_null if abs(p - th) < 1e-9 else a.n_alt):
                rng = np.random.default_rng([a.seed3, 0, ti, int(round(p * 1000)), ep])     # source index 0 = ku, as script 49
                cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(BG, M - nb, replace=False)]
                rng.shuffle(cells); cells = cells.reshape(N, m)
                f = y[cells].mean(1).astype(float); S = score[cells].mean(1)
                r = dict(run=tag, theta=th, p=p, ep=ep)
                _, st = alpha_martingale(f[rng.permutation(N)] / (2 * th), 1 / (2 * th), 0.5, N, alpha=0.05)
                r["alpha_uniform_stop"] = st
                r["bet_uniform_stop"] = betting_test_w(f, pi, S, sampler_uniform_w(N), th, rng=rng)["stop"]
                r["bet_prop_stop"] = betting_test_w(f, pi, S, sampler_prop_w(pi, S), th, rng=rng)["stop"]
                r["hedge_uniform_stop"] = hedged_betting_w(f, pi, sampler_uniform_w(N), th, rng=rng)["stop"]
                r["hedge_prop_stop"] = hedged_betting_w(f, pi, sampler_prop_w(pi, S), th, rng=rng)["stop"]
                r["hedge_mix50_stop"] = hedged_betting_w(f, pi, mix(S, 0.5), th, rng=rng)["stop"]
                for d in DESIGNS:
                    s = r.pop(f"{d}_stop"); r[f"{d}_cert"] = int(s > 0); r[f"{d}_cells"] = (s if s > 0 else N) * m
                rows.append(r)
    return pd.DataFrame(rows)


def part1(E, df, tag):
    bank = CellBank(E, df); N, k = 200, 1; base_cfg = InfoDPPConfig(tau=None, seed=a.seed1); rows = []
    for fold in sorted(np.unique(bank.fold)):
        if len(bank.blast_idx[fold]) == 0: continue
        for m in a.m_list:
            rng = np.random.default_rng(a.seed1 + 1000 * int(fold) + m); cal = []
            for f in [x for x in np.unique(bank.fold) if x != fold][:3]:
                for _ in range(2):
                    s = make_specimen(bank, f, N, 0, m, rng, "pool", "mean"); cal.append((s["F"], s["teacher"]))
            tau = calibrate_tau(cal, base_cfg, 0.163, rng=rng); cfg = InfoDPPConfig(**{**base_cfg.__dict__, "tau": tau})
            for ep in range(a.episodes):
                s = make_specimen(bank, fold, N, k, m, rng, "pool", "mean")
                res = run_infodpp(s["F"], s["teacher"], cfg, np.random.default_rng(rng.integers(1 << 30)))
                n = res["k_star"]; ci = np.flatnonzero(s["crit"]); uni = rng.permutation(N)
                rows.append(dict(run=tag, fold=int(fold), m=m, ep=ep, infodpp_n=n, infodpp_miss=int(not np.isin(ci, res["order"][:n]).any()),
                                 uniform_miss=int(not np.isin(ci, uni[:n]).any())))
    return pd.DataFrame(rows)


summary = dict(blast_counts=ku_all.loc[ku_all.is_blast == 1].groupby("patient_id").size().astype(int).to_dict(), runs={})
for tag in RUNS:
    t1 = time.time()
    if tag in ("none", "none_refit"):
        keep = np.ones(len(ku_all), bool)
    else:
        keep = (ku_all.patient_id != tag).values
    df = ku_all[keep].reset_index(drop=True).copy(); E = E_all[keep]
    if tag != "none":
        old = df.probe_p_oof.values.copy(); df["probe_p_oof"] = refit_oof(E, df)
        drift = dict(corr=float(np.corrcoef(old, df.probe_p_oof)[0, 1]), max_abs=float(np.abs(old - df.probe_p_oof).max())) if tag == "none_refit" else None
    else:
        drift = None
    D3 = part3(df, tag); D3.to_csv(f"{a.out}/part3_episodes_{tag}.csv", index=False)
    D1 = part1(E, df, tag); D1.to_csv(f"{a.out}/part1_episodes_{tag}.csv", index=False)
    r = dict(blasts_left=int(df.is_blast.sum()), probe_drift_vs_stored=drift, part3={}, part3_validity={}, part1={})
    for (th, p), g in D3.groupby(["theta", "p"]):
        key = f"theta{th:.2f}_p{p:.2f}"
        if abs(p - th) < 1e-9:
            for d in DESIGNS:
                k_ = int(g[f"{d}_cert"].sum()); n_ = len(g)
                r["part3_validity"][f"{key}_{d}"] = dict(rate=k_ / n_, cp95_upper=float(beta.ppf(0.975, k_ + 1, n_ - k_)))
            continue
        med = {d: float(g[f"{d}_cells"].median()) for d in DESIGNS}
        r["part3"][key] = dict(median_cells=med, **{f"{d}_over_alpha": med[d] / med["alpha_uniform"] for d in DESIGNS if d != "alpha_uniform"})
    for m_, g in D1.groupby("m"):
        r["part1"][int(m_)] = dict(n=int(len(g)), infodpp_miss=float(g.infodpp_miss.mean()), uniform_miss=float(g.uniform_miss.mean()),
                                   mean_stop=float(g.infodpp_n.mean()))
    summary["runs"][tag] = r
    json.dump(summary, open(f"{a.out}/summary.json", "w"), indent=1)
    print(f"run {tag}: {r['blasts_left']} blasts left, {time.time() - t1:.0f}s (total {time.time() - T0:.0f}s)", flush=True)

import numpy, scipy, sklearn
json.dump(dict(config=vars(a), runs=RUNS, env=dict(python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__,
               pandas=pd.__version__, sklearn=sklearn.__version__)), open(f"{a.out}/config.json", "w"), indent=1)
print(f"DONE in {time.time() - T0:.0f}s")
