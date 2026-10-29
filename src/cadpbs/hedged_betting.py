"""
Hedged betting test for a finite-population weighted mean under predictable sampling without replacement.

Setting of Shekhar, Xu, Lipton, Liang & Ramdas (UAI 2023, Thm. 1): weights pi(i) summing to 1, values f(i) in [0, 1],
H0: sum_i pi(i) f(i) <= m0. With
    Z_t  = pi(I_t) f(I_t) / q_t(I_t),    mu_t = m0 - sum_{s<t} pi(I_s) f(I_s),
each component W^(k) = prod_s (1 + kappa_k (Z_s / mu_s - 1)) is a non-negative test supermartingale under H0,
and so is their average. Reject when the average reaches 1/alpha (Ville). The grid kappa_k is fixed in advance.
Positivity of q_t is needed for power, not for validity of this upper-tail test (f >= 0).
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
