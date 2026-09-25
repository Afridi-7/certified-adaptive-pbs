"""
Anytime-valid certification of a finite-population threshold assertion (Part III, RQ3).

ESTIMAND. A specimen is a finite population of M = N*m cells (N fields of m cells). Let b_i in {0,1} mark a
blast. The quantity of interest is the blast proportion p = (1/M) sum_i b_i. The ASSERTION to certify is
p > theta for a stated diagnostic boundary theta. theta is an input, never a default: the WHO 5th edition
(2022) and the ICC (2022) use different boundaries, and which applies depends on the clinical context. This
code certifies a statement about a finite population of cells; it does not establish diagnostic validity.

ASSORTER (Stark's device). a_i = b_i / (2*theta) in [0, u], u = 1/(2*theta). Then
        mean(a) > 1/2   <=>   p > theta,
so the assertion becomes a bounded-mean assertion with null mean mu0 = 1/2.

TEST. ALPHA (Stark, *Ann. Appl. Stat.* 17(1), 2023, DOI 10.1214/22-AOAS1646): a non-negative supermartingale
under the null for sampling WITHOUT replacement from a finite population,
        T_j = T_{j-1} * ( X_j * eta_j/mu_j + (u - X_j) * (u - eta_j)/(u - mu_j) ) / u,
where mu_j = (M*mu0 - sum_{i<j} X_i) / (M - j + 1) is the null mean of the remaining population and eta_j is a
prediction of the alternative mean formed from past draws only (shrink-trunc estimator). By Ville's
inequality, P(sup_j T_j >= 1/alpha) <= alpha under the null, so "reject when T_j >= 1/alpha" is anytime-valid:
the analyst may stop whenever they like.

UNITS. We sample FIELDS, and an expert verifies every cell in a sampled field, so the observation is the field
mean of the cell assorters, which is bounded by the same u. With equal-size fields the population mean over
fields equals the cell-level proportion.

ADAPTIVE SELECTION. Uniform sampling without replacement satisfies ALPHA's conditions directly. For adaptive
(classifier-driven) sampling use cadpbs.weighted_betting, which transcribes Shekhar et al. (arXiv:2305.06884).
WITHDRAWN (2026-09-25): an earlier `ht_assorter` construction fed importance-weighted draws into ALPHA's
UNWEIGHTED null-mean update mu_j = (M*mu0 - sum X)/(M - j + 1). That is not a valid supermartingale
construction; its results (docs/certification_v1.md, finding 3) are withdrawn. `sample_prioritised` is kept
only to build the deliberately invalid "naive" design.
"""
import numpy as np


def assorter_bound(theta):
    return 1.0 / (2.0 * theta)


def alpha_martingale(X, u, mu0, M, d=20.0, eps=1e-6, alpha=0.05, eta0=None):
    """ALPHA supermartingale. X: observed assorter values in draw order (each in [0, u]).
    Returns (log T trajectory, first index (1-based) at which T >= 1/alpha, or -1)."""
    X = np.asarray(X, float)
    logT, s, stop = [], 0.0, -1
    logt = 0.0
    eta0 = eta0 if eta0 is not None else (mu0 + u) / 2
    for j, x in enumerate(X, start=1):
        mu = (M * mu0 - s) / (M - j + 1)                       # null mean of what remains
        if mu <= 0:
            # Under the null the remaining units must all be 0. A positive draw contradicts it outright;
            # a zero draw is fully consistent and must leave the statistic unchanged.
            if x > 0:
                stop = j if stop < 0 else stop; logt = np.inf; logT.append(logt); break
            logT.append(logt); s += x; continue
        if mu >= u:
            # Under the null the remaining units must all equal u; nothing further can be learned.
            logT.append(logt); s += x; continue
        # shrink-trunc prediction of the alternative mean, using past draws only
        eta = (d * eta0 + s) / (d + j - 1)
        eta = min(max(eta, mu + eps), u * (1 - eps))
        term = (x * eta / mu + (u - x) * (u - eta) / (u - mu)) / u
        logt += np.log(max(term, 1e-300))
        logT.append(logt)
        s += x
        if stop < 0 and logt >= np.log(1.0 / alpha):
            stop = j
    return np.array(logT), stop


def sample_uniform(N, rng):
    """Uniform sampling without replacement: pi_t = 1/(N-t+1) at draw t."""
    order = rng.permutation(N)
    probs = 1.0 / (N - np.arange(N))
    return order, probs


def sample_prioritised(N, score, rng, lam=0.5):
    """Classifier-prioritised sampling WITHOUT replacement, probabilities recorded.
    At each draw the remaining fields are sampled with probability proportional to
        w_i = (1 - lam) * score_i / sum(score_remaining) + lam / n_remaining,
    i.e. mixed with the uniform distribution so no probability can vanish (needed for bounded HT weights)."""
    remaining = list(range(N)); order, probs = [], []
    sc = np.asarray(score, float); sc = np.clip(sc, 1e-6, None)
    while remaining:
        w = sc[remaining] / sc[remaining].sum()
        p = (1 - lam) * w + lam / len(remaining)
        j = rng.choice(len(remaining), p=p)
        order.append(remaining[j]); probs.append(float(p[j])); remaining.pop(j)
    return np.array(order), np.array(probs)


# ht_assorter removed: see WITHDRAWN note in the module docstring.
