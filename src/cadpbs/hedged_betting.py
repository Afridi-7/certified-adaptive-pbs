"""
Hedged betting test for a finite-population weighted mean under arbitrary predictable sampling without replacement
(added 2026-10-22, exploratory; the frozen modules weighted_betting.py and weighted_betting_w.py are unchanged).

Setting as in Shekhar, Xu, Lipton, Liang & Ramdas (UAI 2023, Thm. 1): weights pi(i) > 0 summing to 1, field values
f(i) in [0, 1], target m* = sum_i pi(i) f(i). At step t a field I_t is drawn from a predictable distribution q_t (positivity on
every unsampled field is needed for power and for lower-tail tests, not for the validity of this test, because f >= 0;
see tests/test_validity.py::test_above_certificate_valid_without_positivity), and
    Z_t  = pi(I_t) f(I_t) / q_t(I_t)                 E[Z_t | F_{t-1}] = sum over unsampled i of pi(i) f(i)
    mu_t = m0 - sum_{s<t} pi(I_s) f(I_s)             (null mass that may remain under H0: m* <= m0)
Under H0, E[Z_t - mu_t | F_{t-1}] = m* - m0 <= 0.

Each component k uses the fixed fraction kappa_k in (0, 1) and the predictable bet lambda_t = kappa_k / mu_t, so
    W_t^(k) = prod_{s<=t} (1 + kappa_k (Z_s / mu_s - 1)) >= prod (1 - kappa_k) > 0
is a non-negative test supermartingale under H0 (Z_s >= 0). The equally weighted average W_t = mean_k W_t^(k) is again a
test supermartingale, so rejecting H0 when W_t >= 1/alpha has false-certification probability <= alpha at any stopping
time (Ville). If the sampled mass already exceeds m0, H0 is false with certainty (logical rule, as in the frozen code).
The grid of fractions is fixed in advance and does not depend on the data.
"""
import numpy as np

KAPPAS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


def hedged_betting_w(f, pi, sampler, m0, alpha=0.05, kappas=KAPPAS, rng=None):
    """Returns dict(stop = first draw at which H0: m* <= m0 is rejected, or -1; order of draws)."""
    rng = rng or np.random.default_rng()
    f = np.asarray(f, float); pi = np.asarray(pi, float); N = len(f)
    assert np.all(f >= -1e-12) and np.all(f <= 1 + 1e-12) and np.all(pi > 0) and abs(pi.sum() - 1) < 1e-9
    kap = np.asarray(kappas, float); logW = np.zeros(len(kap)); thr = np.log(1.0 / alpha)
    remaining = list(range(N)); sampled = 0.0; order = []
    for t in range(1, N + 1):
        rem = np.array(remaining); q = sampler(rem, t)
        j = rng.choice(len(rem), p=q); i = int(rem[j]); qi = float(q[j])
        Z = pi[i] * f[i] / qi; mu = m0 - sampled
        if mu > 0:
            logW += np.log1p(kap * (Z / mu - 1.0))
        sampled += pi[i] * f[i]; order.append(i); remaining.pop(j)
        mx = logW.max(); wealth = mx + np.log(np.mean(np.exp(logW - mx)))
        if wealth >= thr or sampled > m0 + 1e-12:
            return dict(stop=t, order=np.array(order))
    return dict(stop=-1, order=np.array(order))
