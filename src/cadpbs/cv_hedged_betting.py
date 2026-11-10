"""
Hedged betting with a prediction control variate, and a score-stratified sampler (added 2026-11-10, exploratory;
no frozen module is changed). RQ3, RQ4.

Purpose. The paper uses the classifier only through the SAMPLING design (prop-MS, Shekhar et al., UAI 2023, Sec. 3),
which loses to uniform sampling when the scores are poor. Two classical alternatives use the scores without letting
them drive the design completely:
  (1) a control variate under uniform (or mixed) sampling: the design-based analogue of a difference / prediction-
      powered estimator (Shekhar et al., Sec. 4; cf. Csillag et al., ICML 2025, for the superpopulation setting);
  (2) stratification by score: fields are grouped into H score strata and sampled uniformly within a stratum.

Setting (as in hedged_betting.py). Weights pi(i) > 0 summing to 1, verified field values f(i) in [0, 1], scores S(i)
known for every field before sampling, target m* = sum_i pi(i) f(i), H0: m* <= m0. At step t a field I_t is drawn from
a predictable distribution q_t over the unsampled fields.
    Z_t  = pi(I_t) f(I_t) / q_t(I_t)                    E[Z_t | F_{t-1}] <= sum over unsampled i of pi(i) f(i)
    g_t(i) = pi(i) S(i) / q_t(i)                        the same transform applied to the score
    U_t  = g_t(I_t) - sum_{i: q_t(i) > 0} pi(i) S(i)    E[U_t | F_{t-1}] = 0 for every q_t (no condition on S)
    mu_t = m0 - sum_{s<t} pi(I_s) f(I_s)                null mass that may remain
For a predictable coefficient c_t >= 0 the payoff Y_t = Z_t - c_t U_t - mu_t has E[Y_t | F_{t-1}] <= m* - m0 <= 0
under H0. Because Z_t >= 0 and U_t <= Umax_t := max_i g_t(i) - sum pi S (known before the draw),
    Y_t >= -(mu_t + c_t Umax_t),
so with the predictable bet lambda_t = kappa / (mu_t + c_t Umax_t), kappa in (0, 1),
    W_t^(kappa) = prod_{s<=t} (1 + lambda_s Y_s) >= prod (1 - kappa) > 0
is a non-negative test supermartingale under H0, WHATEVER the scores are (they enter only through predictable
quantities and a zero-mean term). The average over a fixed grid of kappa is again a test supermartingale; rejecting
when it reaches 1/alpha has false-certification probability <= alpha at any stopping time (Ville).

Choice of c (default c_rule = "mix"). The test statistic is the equally weighted average of the processes for every
pair (c, kappa) of two fixed grids; each is a test supermartingale, so the average is one. With c_grid = (0,) the test
equals hedged_betting_w exactly. The alternative c_rule = "ftl" ("follow the leader": at each step use the grid value
whose fixed-c wealth is largest so far, a predictable choice) is also valid but was clearly worse on synthetic
specimens (it switches on noise) and is kept only for the record.
Note the price of a control variate in a betting test: the bet is capped by 1 / (mu_t + c_t Umax_t) instead of
1 / mu_t, because a high-scoring field without blasts gives a large negative payoff. Under prop-MS sampling g_t is
constant, U_t = 0, and the control variate has no effect.
"""
import numpy as np

KAPPAS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
C_GRID = (0.0, 0.25, 0.5)


def _mix(logw):
    mx = logw.max(); return mx + np.log(np.mean(np.exp(logw - mx)))


def hedged_cv_betting_w(f, pi, S, sampler, m0, alpha=0.05, kappas=KAPPAS, c_grid=C_GRID, rng=None, return_c=False, c_rule="mix"):
    """Returns dict(stop = first draw at which H0: m* <= m0 is rejected, or -1; order; c_path if return_c)."""
    rng = rng or np.random.default_rng()
    f = np.asarray(f, float); pi = np.asarray(pi, float); S = np.asarray(S, float); N = len(f)
    assert np.all(f >= -1e-12) and np.all(f <= 1 + 1e-12) and np.all(pi > 0) and abs(pi.sum() - 1) < 1e-9
    assert np.all(S >= 0) and min(c_grid) >= 0 and c_rule in ("mix", "ftl")
    kap = np.asarray(kappas, float); cg = np.asarray(c_grid, float); thr = np.log(1.0 / alpha)
    logW = np.zeros(len(kap))                      # wealth components of the "ftl" rule
    hyp = np.zeros((len(cg), len(kap)))            # wealth components for each fixed c (the "mix" rule averages them)
    remaining = list(range(N)); sampled = 0.0; order = []; c_path = []
    for t in range(1, N + 1):
        rem = np.array(remaining); q = sampler(rem, t)
        supp = q > 0
        g = np.where(supp, pi[rem] * S[rem] / np.where(supp, q, 1.0), 0.0)
        ES = float((pi[rem] * S[rem])[supp].sum()); Umax = max(float(g[supp].max()) - ES, 0.0)
        mu = m0 - sampled
        c_t = float(cg[int(np.argmax([_mix(hyp[a]) for a in range(len(cg))]))]) if t > 1 else float(cg.min())
        j = rng.choice(len(rem), p=q); i = int(rem[j]); qi = float(q[j])
        Z = pi[i] * f[i] / qi; U = float(g[j]) - ES
        if mu > 0:
            logW += np.log1p(kap * (Z - c_t * U - mu) / (mu + c_t * Umax))
            for a_, c in enumerate(cg):
                hyp[a_] += np.log1p(kap * (Z - c * U - mu) / (mu + c * Umax))
        sampled += pi[i] * f[i]; order.append(i); c_path.append(c_t); remaining.pop(j)
        if (_mix(hyp.ravel()) if c_rule == "mix" else _mix(logW)) >= thr or sampled > m0 + 1e-12:
            out = dict(stop=t, order=np.array(order))
            if return_c: out["c_path"] = np.array(c_path)
            return out
    out = dict(stop=-1, order=np.array(order))
    if return_c: out["c_path"] = np.array(c_path)
    return out


def sampler_strata_w(pi, S, H=4, power=0.5, floor=1e-3):
    """Score-stratified sampling. Fields are split into H strata at the quantiles of S (fixed before sampling); a field
    of stratum h is drawn with probability proportional to pi(i) * max(mean score of h, floor) ** power, i.e. uniformly
    within a stratum, with strata of higher mean score over-sampled. power = 0.5 follows the Neyman-type rule
    (probability proportional to the root of the expected squared value); power = 0 is uniform sampling. The largest
    ratio of two selection probabilities is bounded by (max stratum mean / min stratum mean) ** power, whereas prop-MS
    is bounded only by 1 / floor."""
    pi = np.asarray(pi, float); S = np.asarray(S, float)
    edges = np.quantile(S, np.linspace(0, 1, H + 1)[1:-1]); lab = np.searchsorted(edges, S, side="right")
    sbar = np.array([S[lab == h].mean() if np.any(lab == h) else 0.0 for h in range(H)])
    w_all = pi * np.maximum(sbar[lab], floor) ** power
    def q(remaining, t):
        w = w_all[remaining]; return w / w.sum()
    return q


def sampler_mix_w(base, lam):
    """Mixture of a sampler with the uniform distribution on the unsampled fields (uniform share lam)."""
    return lambda rem, t: (1 - lam) * base(rem, t) + lam / len(rem)
