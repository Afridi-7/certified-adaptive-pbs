"""
Kill-test 1 and Kill-test 2 on composite specimens built from KU-Optofil embeddings (RQ1).

KT1 (are the three stopping signals separable?): at every greedy step we record
    - Gamma_k = max residual information gain (the InfoDPP-PAC stopping statistic, Eq. 6)
    - the relevance q of the field actually selected (Eq. 3)
    - whether that field is critical (ground truth)
  and per episode: AUC of the step-1 information gain for critical vs non-critical fields, the same for
  relevance, and the stop step k*. H1a says k* does not depend on rarity k.

KT2 (does a gap against matched-cost baselines appear?): premature-stop (miss) and critical-evidence recall
  for InfoDPP-PAC vs uniform-without-replacement and a fixed-order scan, all evaluated at the SAME number of
  fields that InfoDPP-PAC used in that episode (paired -> exact McNemar). Confidence stopping is reported at
  its own cost. Axes: rarity k and dilution m.

Note on the "raster" baseline: composites have no spatial geometry, so a raster scan is just a fixed
arbitrary order. It is reported as "fixed_order" and must not be read as a spatial baseline.

Resumable: one CSV per (fold, m); re-running skips finished ones.
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.composites import CellBank, make_specimen
from cadpbs.infodpp import InfoDPPConfig, run_infodpp, calibrate_tau
from sklearn.metrics import roc_auc_score

ap = argparse.ArgumentParser()
ap.add_argument("--emb", default="results/embeddings_v1")
ap.add_argument("--out", default="results/killtests_v1")
ap.add_argument("--N", type=int, default=200)
ap.add_argument("--k_list", type=int, nargs="+", default=[1, 2, 5, 10, 20])
ap.add_argument("--m_list", type=int, nargs="+", default=[1, 5, 10, 20])
ap.add_argument("--folds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
ap.add_argument("--episodes", type=int, default=100)
ap.add_argument("--pool", default="pool", choices=["pool", "host"])
ap.add_argument("--pooling", default="mean", choices=["mean", "max", "single"])
ap.add_argument("--target_frac", type=float, default=0.163)   # paper's operating point
ap.add_argument("--n_calib", type=int, default=8)
ap.add_argument("--tau_conf", type=float, default=0.9)        # confidence-stopping threshold
ap.add_argument("--m_conf", type=int, default=30)             # consecutive confident-negative fields
ap.add_argument("--seed", type=int, default=20260919)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

E = np.load(f"{a.emb}/embeddings_dinov2_vits14.npy").astype(np.float32)
df = pd.read_csv(f"{a.emb}/cells_features.csv", dtype={"patient_id": str})
bank = CellBank(E, df)
base_cfg = InfoDPPConfig(tau=None, seed=a.seed)


def confidence_stop(order, teacher, tau, mrun):
    run = 0
    for t, i in enumerate(order):
        if teacher[i] >= tau: return t + 1
        run = run + 1 if teacher[i] <= 1 - tau else 0
        if run >= mrun: return t + 1
    return len(order)


def calibrate(fold, m, rng):
    """tau is calibrated on NEGATIVE specimens from the OTHER folds (never the evaluation fold)."""
    cal = []
    for f in [x for x in np.unique(bank.fold) if x != fold][:3]:
        for _ in range(max(1, a.n_calib // 3)):
            s = make_specimen(bank, f, a.N, 0, m, rng, a.pool, a.pooling)
            cal.append((s["F"], s["teacher"]))
    return calibrate_tau(cal, base_cfg, a.target_frac, rng=rng)


rows_all = []
for fold in a.folds:
    for m in a.m_list:
        path = f"{a.out}/episodes_fold{fold}_m{m}.csv"
        if os.path.exists(path):
            rows_all.append(pd.read_csv(path)); print("skip (done):", path); continue
        rng = np.random.default_rng(a.seed + 1000 * fold + m)
        t0 = time.time(); tau = calibrate(fold, m, rng)
        cfg = InfoDPPConfig(**{**base_cfg.__dict__, "tau": tau})
        rows = []
        for k in a.k_list:
            for ep in range(a.episodes):
                s = make_specimen(bank, fold, a.N, k, m, rng, a.pool, a.pooling)
                res = run_infodpp(s["F"], s["teacher"], cfg, np.random.default_rng(rng.integers(1 << 30)))
                n_info = res["k_star"]; crit = s["crit"]; ci = np.flatnonzero(crit)
                uni = rng.permutation(a.N); fixed = np.arange(a.N)
                n_conf = confidence_stop(uni, s["teacher"], a.tau_conf, a.m_conf)

                def ev(order, n):
                    seen = np.isin(ci, order[:n])
                    return int(not seen.any()), float(seen.mean())
                r = dict(fold=fold, m=m, k=k, ep=ep, N=a.N, tau=tau, blast_fraction=s["blast_fraction"],
                         infodpp_n=n_info, conf_n=n_conf, stop_reason=res["stop_reason"])
                for nm, (o, n) in {"infodpp": (res["order"], n_info), "uniform@info": (uni, n_info),
                                   "fixed_order@info": (fixed, n_info), "confidence": (uni, n_conf),
                                   "uniform@conf": (rng.permutation(a.N), n_conf)}.items():
                    miss, rec = ev(o, n); r[f"{nm}_miss"] = miss; r[f"{nm}_recall"] = rec
                # KT1: are the three signals separable? AUCs for critical vs non-critical fields, using
                # (a) the information gain before any selection, (b) the teacher relevance.
                r["kt1_auc_info_step1"] = float(roc_auc_score(crit, res["dinfo0"]))
                r["kt1_auc_relevance"] = float(roc_auc_score(crit, s["teacher"]))
                r["kt1_corr_info_relevance"] = float(np.corrcoef(res["dinfo0"], s["teacher"])[0, 1])
                r["kt1_rank_of_first_critical"] = int(np.argmax(np.isin(res["order"], ci))) + 1 if np.isin(res["order"], ci).any() else -1
                rows.append(r)
        d = pd.DataFrame(rows); d.to_csv(path, index=False); rows_all.append(d)
        print(f"fold {fold} m {m}: tau={tau:.3f} {len(d)} episodes in {time.time()-t0:.0f}s", flush=True)

if rows_all:
    all_df = pd.concat(rows_all, ignore_index=True)
    all_df.to_csv(f"{a.out}/episodes.csv", index=False)
    json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
    print("total episodes:", len(all_df))
