"""
Unequal-weight version of weighted_betting.py (added 2026-10-19 for P0-3; the frozen module is left unchanged).

Shekhar, Xu, Lipton, Liang & Ramdas (UAI 2023, PMLR 216:1932-1941), Sec. 2, Thm. 1, with GENERAL weights pi(i) > 0,
sum_i pi(i) = 1, and values f(i) in [0, 1]:
  m*      = sum_i pi(i) f(i)                       (here: the cell-level blast proportion p when pi(i) = m_i / M and
                                                    f(i) = B_i / m_i, i.e. field i holds m_i cells, B_i of them blasts)
  Z_t     = pi(I_t) f(I_t) / q_t(I_t)              E[Z_t | F_{t-1}] = sum over unsampled i of pi(i) f(i)
  mu_t(m) = m - sum_{s<t} pi(I_s) f(I_s)           null mass that may remain
  W_t     = prod_s (1 + lambda_s (Z_s + beta_s U_s - mu_s))
One-sided H0: m* <= m0. With lambda_s >= 0, predictable, and 1 + lambda_s * (lower bound of the bracket) >= 0, W is a
non-negative test supermartingale under H0 for ANY predictable q_t that is positive on every unsampled unit (Ville).
Equal weights pi = 1/N recover weighted_betting.betting_test exactly.

Two regimes for unequal field sizes (see scripts/27_unequal_fields.py):
  (a) cell counts m_i of all fields are known before verification (e.g. a nucleus detector counts cells, only their
      classes need an expert): use pi(i) = m_i / M, f(i) = B_i / m_i. This module.
  (b) counts of unsampled fields are unknown, only an upper bound m_max is known: use the bounded assorter
      a_i = (B_i - theta m_i + theta m_max) / m_max in [0, 1] with equal weights and null mean theta.
"""
import numpy as np


def sampler_uniform_w(N):
    def q(remaining, t):
        return np.full(len(remaining), 1.0 / len(remaining))
    return q


def sampler_prop_w(pi, S, floor=1e-3):
    """prop-MS: q_t(i) proportional to pi(i) * S(i) over the unsampled units; the floor keeps q_t > 0."""
    w_all = np.asarray(pi, float) * np.maximum(np.asarray(S, float), floor)
    def q(remaining, t):
        w = w_all[remaining]; return w / w.sum()
    return q


def betting_test_w(f, pi, S, sampler, m0, alpha=0.05, c=0.75, use_cv=False, rng=None, alphas=()):
    """f in [0,1]^N (revealed by verification), pi > 0 summing to 1, S side information.
    Returns dict(stop=first draw at which H0 is rejected at `alpha` or -1, stops={a: draw} for extra levels, cells
    order). Extra levels in `alphas` reuse the same wealth path (the bets do not depend on alpha)."""
    rng = rng or np.random.default_rng()
    f = np.asarray(f, float); pi = np.asarray(pi, float); S = np.asarray(S, float); N = len(f)
    assert np.all(f >= -1e-12) and np.all(f <= 1 + 1e-12) and np.all(pi > 0) and abs(pi.sum() - 1) < 1e-9
    remaining = list(range(N)); sampled_mass = 0.0
    A = V = num_b = den_b = 0.0; logW = 0.0; order = []
    levels = sorted(set([alpha, *alphas])); stops = {lv: -1 for lv in levels}
    for t in range(1, N + 1):
        rem = np.array(remaining); q = sampler(rem, t)
        ES = float(q @ S[rem]); mu = m0 - sampled_mass
        beta = float(np.clip(-num_b / den_b, -1, 1)) if (use_cv and den_b > 0) else 0.0
        LB = float((beta * (S[rem] - ES)).min()) - mu
        j = rng.choice(len(rem), p=q); i = int(rem[j]); qi = float(q[j])
        Z = pi[i] * f[i] / qi; U = S[i] - ES; Y = Z + beta * U - mu
        if LB >= 0:
            lam = 0.0 if V == 0 else max(0.0, c * A / V)
        else:
            cap = c / (-LB); lam = cap * 0.5 if V == 0 else float(np.clip(c * A / V, 0.0, cap))
        logW += np.log(max(1.0 + lam * Y, 1e-300))
        A += Y; V += Y * Y
        if use_cv:
            num_b += (Z - mu) * U; den_b += U * U
        sampled_mass += pi[i] * f[i]; order.append(i); remaining.pop(j)
        logical = sampled_mass > m0 + 1e-12
        for lv in levels:
            if stops[lv] < 0 and (logW >= np.log(1 / lv) or logical):
                stops[lv] = t
        if all(v > 0 for v in stops.values()):
            break
    return dict(stop=stops[alpha], stops=stops, order=np.array(order))
