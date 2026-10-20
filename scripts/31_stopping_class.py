"""
Information-saturation stopping as a CLASS of rules, and a bridge to Theorem 1 (P1-1; exploratory). RQ1, RQ2.

The draft's Part I tests one rule (InfoDPP-PAC). Here every rule below runs on the same composite specimens (the
unchanged generator, N = 200 fields, mean-pooled), with a field-level teacher trained only on the other folds
(as in scripts/26_kt1_field_teacher.py):

  infodpp      InfoDPP-PAC with the field-level teacher (tau calibrated to the target cost on blast-free specimens)
  gp_var       pure information-gain greedy (InfoDPP-PAC with lambda2 -> 0, set to 1e-9): uncertainty sampling, same
               halting rule
  kcenter      farthest-point (k-center greedy) selection, stops when the covering radius falls below r*, a coverage
               (novelty) saturation rule; r* calibrated to the same target cost
  uniform      uniform sampling at InfoDPP-PAC's cost in the same episode
  ranked       fields in decreasing field-teacher score, at InfoDPP-PAC's cost (no guarantee; location information)

Bridge (Proposition 3 in the paper). Under the model that the critical field's teacher score S_c ~ F1 and the N - 1
non-critical scores are i.i.d. ~ F0, independent of S_c, a scan of the top n fields misses the single critical field with
probability
        P(miss) = E_{S_c ~ F1} [ P( Binomial(N - 1, 1 - F0(S_c)) >= n ) ]           (exact under the model)
                ~ 1 - TPR(t*),   where FPR(t*) = n / (N - 1)                        (large-N approximation)
The script stores, per episode, the critical score and the non-critical scores so the prediction can be compared
with the observed miss of the ranked scan.

Usage: python scripts/31_stopping_class.py --out results/stopping_class_v1 --seed 20261024
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import binom
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.composites import CellBank, make_specimen
from cadpbs.infodpp import InfoDPPConfig, run_infodpp, calibrate_tau
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

ap = argparse.ArgumentParser()
ap.add_argument("--emb", default="results/embeddings_v1"); ap.add_argument("--out", default="results/stopping_class_v1")
ap.add_argument("--N", type=int, default=200); ap.add_argument("--k_list", type=int, nargs="+", default=[1, 5])
ap.add_argument("--m_list", type=int, nargs="+", default=[1, 5, 10, 20]); ap.add_argument("--episodes", type=int, default=100)
ap.add_argument("--n_train", type=int, default=3000); ap.add_argument("--target_frac", type=float, default=0.163)
ap.add_argument("--n_calib", type=int, default=8); ap.add_argument("--seed", type=int, default=20261024)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
E = np.load(f"{a.emb}/embeddings_dinov2_vits14.npy").astype(np.float32)
df = pd.read_csv(f"{a.emb}/cells_features.csv", dtype={"patient_id": str})
bank = CellBank(E, df); folds_all = np.unique(bank.fold)
cfg_info = InfoDPPConfig(tau=None, seed=a.seed)
cfg_var = InfoDPPConfig(tau=None, seed=a.seed, lambda2=1e-9)   # exactly 0 gives 0 * (-inf) = nan in Eq. 4 and breaks argmax


def pool_fields(cells):
    F = bank.X[cells].mean(1); return F / np.linalg.norm(F, axis=1, keepdims=True)


def train_teacher(f, m, rng):
    other = [g for g in folds_all if g != f]
    bl = np.concatenate([bank.blast_idx[g] for g in other]); bg = np.concatenate([bank.bg_idx[g] for g in other])
    neg = rng.choice(bg, (a.n_train, m)); pos = rng.choice(bg, (a.n_train, m)); pos[:, 0] = rng.choice(bl, a.n_train)
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=3000))
    clf.fit(np.vstack([pool_fields(pos), pool_fields(neg)]), np.r_[np.ones(a.n_train), np.zeros(a.n_train)]); return clf


def kcenter_stop(F, r_star, rng):
    D = np.sqrt(np.maximum(2 - 2 * F @ F.T, 0)); n = len(F)
    order = [int(rng.integers(n))]; mind = D[order[0]].copy(); mind[order[0]] = -1
    while len(order) < n:
        rad = mind.max()
        if rad <= r_star: break
        j = int(np.argmax(mind)); order.append(j); mind = np.minimum(mind, D[j]); mind[order] = -1
    return np.array(order)


def calib_kcenter(specs, rng, iters=30):
    lo, hi = 0.0, 2.0
    for _ in range(iters):
        mid = (lo + hi) / 2
        ks = [len(kcenter_stop(F, mid, np.random.default_rng(i))) for i, F in enumerate(specs)]
        if np.median(ks) > a.target_frac * a.N: lo = mid
        else: hi = mid
    return (lo + hi) / 2


rows = []; scores = {}
for f in folds_all:
    for m in a.m_list:
        rng = np.random.default_rng([a.seed, int(f), m]); t0 = time.time()
        clf = train_teacher(f, m, rng)
        cal = []
        for g in [x for x in folds_all if x != f][:3]:
            for _ in range(max(1, a.n_calib // 3)):
                s = make_specimen(bank, g, a.N, 0, m, rng); cal.append((s["F"], clf.predict_proba(s["F"])[:, 1]))
        tau_i = calibrate_tau(cal, cfg_info, a.target_frac, rng=rng); tau_v = calibrate_tau(cal, cfg_var, a.target_frac, rng=rng)
        r_star = calib_kcenter([F for F, _ in cal], rng)
        ci_ = InfoDPPConfig(**{**cfg_info.__dict__, "tau": tau_i}); cv_ = InfoDPPConfig(**{**cfg_var.__dict__, "tau": tau_v})
        for k in a.k_list:
            for ep in range(a.episodes):
                s = make_specimen(bank, f, a.N, k, m, rng); crit = s["crit"]; ci = np.flatnonzero(crit)
                T = clf.predict_proba(s["F"])[:, 1]
                oi = run_infodpp(s["F"], T, ci_, np.random.default_rng(rng.integers(1 << 30)))
                ov = run_infodpp(s["F"], T, cv_, np.random.default_rng(rng.integers(1 << 30)))
                ok = kcenter_stop(s["F"], r_star, np.random.default_rng(rng.integers(1 << 30)))
                n = oi["k_star"]; uni = rng.permutation(a.N); ranked = np.argsort(-T, kind="stable")
                miss = lambda order, nn: int(not np.isin(ci, order[:nn]).any())
                rows.append(dict(fold=int(f), m=m, k=k, ep=ep,
                                 infodpp_n=n, infodpp_miss=miss(oi["order"], n),
                                 gp_var_n=ov["k_star"], gp_var_miss=miss(ov["order"], ov["k_star"]),
                                 kcenter_n=len(ok), kcenter_miss=miss(ok, len(ok)),
                                 uniform_miss=miss(uni, n), ranked_miss=miss(ranked, n),
                                 uniform_at_kcenter_miss=miss(rng.permutation(a.N), len(ok))))
                if k == 1:
                    scores.setdefault(m, []).append(dict(sc=float(T[ci[0]]), neg=T[~crit].astype(np.float32), n=n))
        print(f"fold {f} m {m}: tau_i={tau_i:.3f} tau_v={tau_v:.3f} r*={r_star:.4f} {time.time()-t0:.0f}s", flush=True)

D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)
# bridge: predicted vs observed miss of the ranked scan at k = 1
bridge = {}
for m, L in scores.items():
    neg_all = np.sort(np.concatenate([x["neg"] for x in L])); sc = np.array([x["sc"] for x in L]); nn = np.array([x["n"] for x in L])
    fpr_at = lambda s: 1 - np.searchsorted(neg_all, s, side="right") / len(neg_all)    # 1 - F0(s)
    exact = float(np.mean([binom.sf(n - 1, a.N - 1, fpr_at(s)) for s, n in zip(sc, nn)]))   # P(Bin >= n)
    nbar = int(np.median(nn)); t_star = np.quantile(neg_all, 1 - nbar / (a.N - 1))
    approx = float(np.mean(sc <= t_star))                                                  # 1 - TPR(t*)
    within = float(np.mean([(x["neg"] >= x["sc"]).sum() >= x["n"] for x in L]))            # observed, same episodes
    bridge[int(m)] = dict(observed_ranked_miss=within, predicted_exact_iid=exact, predicted_1_minus_TPR=approx,
                          median_n=nbar, fpr_star=nbar / (a.N - 1))
summ = dict(by_m_k=D.groupby(["m", "k"])[[c for c in D.columns if c.endswith("_miss") or c.endswith("_n")]].mean().round(4).reset_index().to_dict("records"),
            bridge=bridge)
json.dump(summ, open(f"{a.out}/summary.json", "w"), indent=1); json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
