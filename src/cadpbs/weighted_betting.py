"""
One-sided anytime-valid test of a finite-population mean under ARBITRARY (adaptive) sampling without
replacement, following Shekhar, Xu, Lipton, Liang & Ramdas, "Risk-limiting Financial Audits via Weighted
Sampling without Replacement", arXiv:2305.06884v1 (2023), Section 2 (Theorem 1), Section 3 (prop-MS
sampling) and Section 4 (control variates, Theorem 2). Transcription of their construction:

  population values f(i) >= 0, weights pi(i) (here pi = 1/N for equal-size fields), m* = sum_i pi(i) f(i)
  q_t   : any predictable distribution over the not-yet-sampled units (must be > 0 on all of them)
  Z_t   = pi(I_t) f(I_t) / q_t(I_t)                                    (importance-weighted observation)
  mu_t(m) = m - sum_{j<t} pi(I_j) f(I_j)                               (null mass REMAINING, raw values)
  W_t(m) = W_{t-1}(m) * (1 + lambda_t (Z_t + beta_t U_t - mu_t(m)))   (wealth; beta_t U_t optional)
  U_t   = S(I_t) - E_{I~q_t}[S(I)]                                     (zero-mean control variate)
Under m = m*, E[Z_t | F_{t-1}] = mu_t(m*) and E[U_t | F_{t-1}] = 0, so W_t(m*) is a test martingale for any
predictable q_t, lambda_t, beta_t that keep it non-negative (Prop. 1 / Ville). 

Our use (one-sided): H0: m* <= m0. With lambda_t >= 0, E[Z_t - mu_t(m0) | F] = m* - m0 <= 0, so W_t(m0) is a
test SUPERmartingale under every null value; reject when W >= 1/alpha. Non-negativity needs
1 + lambda * LB_t >= 0 where LB_t is the smallest possible value of (Z + beta U - mu), so lambda <= c / (-LB_t).
Because Z >= 0, the size of Z (small q_t) cannot break validity in the one-sided case; q_t only needs to be
positive on every remaining unit.

Betting: ApproxKelly (their Sec. 2.1) lambda_t = c * A_{t-1} / V_{t-1}, truncated to [0, c/(-LB_t)].
Logical bound (their Sec. 2.2): if sum_{sampled} pi f > m0, H0 is false with certainty.
"""
import numpy as np


def sampler_uniform(N):
    def q(remaining, t):
        return np.full(len(remaining), 1.0 / len(remaining))
    return q


def sampler_prop_side(S, floor=1e-3):
    """prop-MS of Shekhar et al.: q_t(i) proportional to pi(i) S(i) over the remaining units (pi uniform here).
    A floor keeps every q_t(i) > 0, which validity requires."""
    S = np.maximum(np.asarray(S, float), floor)
    def q(remaining, t):
        w = S[remaining]; return w / w.sum()
    return q


def betting_test(f, S, sampler, m0, alpha=0.05, c=0.75, use_cv=False, rng=None):
    """Run the sequential test on one population.
    f: true field values (revealed only when a field is sampled, i.e. expert verification), shape (N,)
    S: side information for every field (classifier-based prediction of f), shape (N,)
    Returns dict(stop=first draw index at which H0 is rejected or -1, order, probs, logW trajectory)."""
    rng = rng or np.random.default_rng()
    f = np.asarray(f, float); S = np.asarray(S, float); N = len(f); pi = 1.0 / N
    remaining = list(range(N)); sampled_mass = 0.0
    A = 0.0; V = 0.0; num_b = 0.0; den_b = 0.0
    logW = 0.0; traj = []; order = []; probs = []; stop = -1
    for t in range(1, N + 1):
        rem = np.array(remaining)
        q = sampler(rem, t)
        ES = float(q @ S[rem])
        mu = m0 - sampled_mass
        beta = float(np.clip(-num_b / den_b, -1, 1)) if (use_cv and den_b > 0) else 0.0
        # smallest possible value of (Z + beta*U - mu): Z >= 0, beta*U ranges over beta*(S_rem - ES)
        u_rng = beta * (S[rem] - ES)
        LB = float(u_rng.min()) - mu
        j = rng.choice(len(rem), p=q)
        i = int(rem[j]); qi = float(q[j])
        Z = pi * f[i] / qi
        U = S[i] - ES
        Y = Z + beta * U - mu
        if LB >= 0:                      # the null cannot lose; any bet is safe -> use the cap-free Kelly guess
            lam = 0.0 if V == 0 else max(0.0, c * A / V)
        else:
            cap = c / (-LB)
            lam = cap * 0.5 if V == 0 else float(np.clip(c * A / V, 0.0, cap))
        logW += np.log(max(1.0 + lam * Y, 1e-300))
        # bookkeeping (updated AFTER the bet so everything above is predictable)
        A += Y; V += Y * Y
        if use_cv:
            num_b += (Z - mu) * U; den_b += U * U
        sampled_mass += pi * f[i]
        order.append(i); probs.append(qi); traj.append(logW)
        remaining.pop(j)
        if stop < 0 and (logW >= np.log(1 / alpha) or sampled_mass > m0 + 1e-12):   # Ville, or logical bound
            stop = t
    return dict(stop=stop, order=np.array(order), probs=np.array(probs), logW=np.array(traj))
