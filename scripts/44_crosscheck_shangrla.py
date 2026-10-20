"""
Cross-check of our test statistics against the reference implementation of ALPHA (SHANGRLA, P. B. Stark et al.).

Why: every validity claim in the paper rests on two pieces of code, src/cadpbs/alpha.py (ALPHA, Stark 2023) and
src/cadpbs/weighted_betting_w.py (Shekhar et al. 2023). This script compares them, on identical inputs, with
`NonnegMean` from SHANGRLA (github.com/pbstark/SHANGRLA, shangrla/core/NonnegMean.py; AGPL-3.0, so it is NOT copied
into this repository). Download that one file and pass its folder with --shangrla.

Checks (R random cases each; populations, sizes, thresholds and orders drawn at random; null and alternative cases):
  A  ALPHA, our full procedure vs SHANGRLA alpha_mart with the shrink-trunc estimator set to our parameters
     (d = 20, f = 0, start eta0, c ~ 0): agreement of the first draw at which p <= alpha, and of the p-value path.
  B  ALPHA supermartingale update alone: SHANGRLA alpha_mart fed with the same predictable eta sequence as ours
     (re-derived from the documented formula), so that only the term and null-mean update are compared.
  C  Betting (Shekhar et al.), equal weights and uniform sampling: our log-wealth path vs SHANGRLA betting_mart with
     the bets lambda'_t = lambda_t (N - t + 1) / N (the same martingale in SHANGRLA's parametrisation, because
     Z_t - mu_t = (N - t + 1)/N * (x_t - m_t) when q_t is uniform).
Steps at which either code enters a special case (null mean of the remainder <= 0 or >= u) are excluded from path
comparisons and counted separately.

Usage: python -I scripts/44_crosscheck_shangrla.py --shangrla /path/to/folder_with_NonnegMean.py --R 2000 --seed 20261035
"""
import argparse, hashlib, json, os, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--shangrla", required=True); ap.add_argument("--R", type=int, default=2000)
ap.add_argument("--seed", type=int, default=20261035); ap.add_argument("--out", default="results/crosscheck_shangrla_v1")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))
sys.path.insert(0, os.path.abspath(a.shangrla))
import warnings; warnings.filterwarnings("ignore")
from NonnegMean import NonnegMean                                   # reference implementation
from cadpbs.alpha import alpha_martingale                          # ours (frozen module)
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w
ref_sha = hashlib.sha256(open(os.path.join(a.shangrla, "NonnegMean.py"), "rb").read()).hexdigest()
rng = np.random.default_rng(a.seed); ALPHA_ = 0.05; D = 20.0


def population(rng):
    """Random finite population of assorter values in [0, u] for a random threshold; mean at, below or above 1/2."""
    th = float(rng.uniform(0.05, 0.30)); u = 1 / (2 * th); m = int(rng.choice([1, 5, 10])); N = int(rng.integers(20, 400))
    target = float(rng.choice([th, th, th - 0.03, th + 0.03, th + 0.08]))
    nb = int(round(np.clip(target, 0, 1) * N * m)); cells = np.zeros(N * m); cells[:nb] = 1; rng.shuffle(cells)
    frac = cells.reshape(N, m).mean(1)
    return th, u, N, frac, frac / (2 * th)


def our_eta(X, u, mu0, N, eta0, eps=1e-6):
    """Predictable eta_j of src/cadpbs/alpha.py, re-derived from its documented formula."""
    s = 0.0; out = []
    for j, x in enumerate(X, start=1):
        mu = (N * mu0 - s) / (N - j + 1)
        eta = (D * eta0 + s) / (D + j - 1); eta = min(max(eta, mu + eps), u * (1 - eps)); out.append(eta); s += x
    return np.array(out)


def null_means(X, mu0, N):
    S = np.concatenate([[0.0], np.cumsum(X)[:-1]]); j = np.arange(1, len(X) + 1)
    return (N * mu0 - S) / (N - j + 1)


def first_le(p, lv):
    h = np.flatnonzero(p <= lv); return int(h[0]) + 1 if len(h) else -1


resA = dict(cases=0, same_stop=0, stop_diffs=[]); resB = dict(cases=0, max_abs_log_diff=0.0, n_steps=0, special_steps=0)
resC = dict(cases=0, max_abs_log_diff=0.0, n_steps=0, special_steps=0)
for r in range(a.R):
    th, u, N, frac, X = population(rng)
    order = rng.permutation(N); n = int(rng.integers(5, N)) if rng.random() < 0.8 else N       # partial sample or census
    x = X[order][:n]; mu0 = 0.5
    eta0 = float(rng.uniform(mu0 + 0.05 * (u - mu0), u * 0.99)) if rng.random() < 0.7 else (mu0 + u) / 2
    m = null_means(x, mu0, N); ok = (m > 1e-9) & (m < u - 1e-9)
    # --- A: full procedures
    logT, st = alpha_martingale(x, u, mu0, N, alpha=ALPHA_, eta0=eta0)
    ref = NonnegMean(test=NonnegMean.alpha_mart, estim=NonnegMean.shrink_trunc, u=u, N=N, t=mu0, eta=eta0, d=D, f=0, c=1e-12)
    _, p_ref = ref.test(x); st_ref = first_le(p_ref, ALPHA_)
    resA["cases"] += 1; resA["same_stop"] += int(st == st_ref)
    # classify a disagreement: does the reference stop because the sample total already exceeds the null total?
    # first draw k (1-based) after which the sample total exceeds the null total; SHANGRLA rejects at k if k is the last
    # draw supplied, otherwise at k + 1 (null mean of the remainder < 0)
    tot = np.cumsum(x); det = np.flatnonzero(tot > N * mu0 + 1e-12); det = int(det[0]) + 1 if len(det) else -1
    if st != st_ref:
        eta_ref = ref.shrink_trunc(x); eta_o = our_eta(x, u, mu0, N, eta0)
        resA["stop_diffs"].append(dict(case=r, ours=st, ref=st_ref, n=n, N=N, deterministic_stop=det,
                                       max_eta_diff=float(np.abs(eta_ref - eta_o).max()),
                                       p_ours_at_ref_stop=float(min(1, np.exp(-logT[st_ref - 1]))) if 0 < st_ref <= len(logT) else None,
                                       explained=bool(det > 0 and st_ref in (det, det + 1))))
    # --- B: same predictable eta, compare the martingale itself
    eta = our_eta(x, u, mu0, N, eta0)
    refB = NonnegMean(test=NonnegMean.alpha_mart, estim=lambda self, xx, **k: eta, u=u, N=N, t=mu0)
    _, pB = refB.test(x)
    lT = np.full(len(x), np.inf); lT[:len(logT)] = np.asarray(logT)[:len(x)]                  # ours stops early after a certain rejection
    T_ours = np.exp(lT); p_ours = np.minimum(1, 1 / T_ours)
    # exclude: steps from the first special case on; the final step of a census (SHANGRLA's "total exceeds the null"
    # rule); steps where either statistic is ~0 (SHANGRLA resets those to 1)
    bad = ~ok | (T_ours < 1e-12) | (pB >= 1 - 1e-15) & (T_ours < 1e-12)
    bad |= np.cumsum(x) > N * mu0 + 1e-12          # the sample total already exceeds the null total (deterministic rule)
    if n == N: bad[-1] = True
    upto = np.flatnonzero(bad); upto = upto[0] if len(upto) else len(x)
    if upto > 0:
        d = np.abs(np.log(np.maximum(p_ours[:upto], 1e-300)) - np.log(np.maximum(pB[:upto], 1e-300)))
        resB["max_abs_log_diff"] = max(resB["max_abs_log_diff"], float(d.max())); resB["n_steps"] += int(upto)
    resB["special_steps"] += int(len(x) - upto); resB["cases"] += 1
    # --- C: betting, equal weights, uniform sampling
    f = frac; S = np.zeros(N); w = np.full(N, 1 / N); seed_c = int(rng.integers(1 << 30))
    ob = betting_test_w(f, w, S, sampler_uniform_w(N), th, alpha=1e-300, rng=np.random.default_rng(seed_c), return_path=True)
    oc = ob["order"]; L = len(ob["lam"]); xc = f[oc][:L]; t_idx = np.arange(1, L + 1)
    lam_ref = ob["lam"] * (N - t_idx + 1) / N
    refC = NonnegMean(test=NonnegMean.betting_mart, bet=lambda self, xx, **k: lam_ref, u=1.0, N=N, t=th)
    mart = refC.test(xc, out="mart") if "out" in NonnegMean.betting_mart.__code__.co_varnames else None
    if mart is None:
        _, pC = refC.test(xc); logref = -np.log(np.maximum(pC, 1e-300)); logours = np.maximum(ob["logW"], 0)
    else:
        mart = mart[1] if isinstance(mart, tuple) else mart
        logref = np.log(np.maximum(np.asarray(mart, float), 1e-300)); logours = ob["logW"]
    mC = null_means(xc, th, N); okC = (mC > 1e-9) & (mC < 1 - 1e-9) & (logours > np.log(1e-12))
    if L and f[oc][:L].sum() / N > th + 1e-12: okC[-1] = False      # our deterministic stop = SHANGRLA's final-term rule
    uptoC = np.flatnonzero(~okC); uptoC = uptoC[0] if len(uptoC) else L
    if uptoC > 0:
        dC = np.abs(logours[:uptoC] - logref[:uptoC])
        resC["max_abs_log_diff"] = max(resC["max_abs_log_diff"], float(dC.max())); resC["n_steps"] += int(uptoC)
    resC["special_steps"] += int(L - uptoC); resC["cases"] += 1

res = dict(reference="SHANGRLA shangrla/core/NonnegMean.py", reference_sha256=ref_sha, R=a.R, seed=a.seed,
           A_full_alpha=dict(cases=resA["cases"], identical_stopping_draw=resA["same_stop"],
                             rate=resA["same_stop"] / resA["cases"], disagreements=resA["stop_diffs"][:20],
                             n_disagreements=len(resA["stop_diffs"]),
                             n_explained_by_deterministic_rule=sum(d["explained"] for d in resA["stop_diffs"]),
                             unexplained=[d for d in resA["stop_diffs"] if not d["explained"]],
                             ours_later_in_all_disagreements=all((d["ours"] < 0 or d["ours"] > d["ref"]) for d in resA["stop_diffs"])),
           B_alpha_martingale_same_eta=resB, C_betting_uniform_equal_weights=resC)
json.dump(res, open(f"{a.out}/summary.json", "w"), indent=1, default=float); json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print(json.dumps(res, indent=1, default=float))
