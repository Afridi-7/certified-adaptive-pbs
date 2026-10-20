"""Helpers for the validity tests. Rates are Monte Carlo estimates with fixed seeds, so every test is deterministic."""
import math
import numpy as np
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w, sampler_prop_w

ALPHA = 0.05


def max_allowed(reps, alpha=ALPHA, z=3.0):
    """Largest Monte Carlo rejection rate compatible with a true rate <= alpha (one-sided, about 3 standard errors)."""
    return alpha + z * math.sqrt(alpha * (1 - alpha) / reps)


def exact_null_fields(N, m, theta, rng, sizes=None, small_enriched=False):
    """Field blast counts B and sizes m_i for a specimen whose cell-level proportion is exactly theta (an exact null)."""
    m_i = np.full(N, m) if sizes is None else np.asarray(sizes)
    M = int(m_i.sum()); nb = int(round(theta * M))
    assert abs(nb - theta * M) < 1e-9, "choose sizes so that theta * M is an integer"
    field_of_cell = np.repeat(np.arange(N), m_i)
    w = 1.0 / m_i[field_of_cell] if small_enriched else np.ones(M)
    blasts = rng.choice(M, nb, replace=False, p=w / w.sum())
    B = np.bincount(field_of_cell[blasts], minlength=N)
    return B, m_i


def alpha_rejects(X, u, N, eta0=None):
    _, st = alpha_martingale(X, u, 0.5, N, alpha=ALPHA, eta0=eta0)
    return st > 0


def betting_rejects(f, w, S, sampler, theta, rng):
    return betting_test_w(f, w, S, sampler, theta, alpha=ALPHA, rng=rng)["stop"] > 0


RESULTS = []


def record(name, rate, reps, limit):
    """Collected by conftest.py into results/validity_tests_v1/report.json."""
    RESULTS.append(dict(test=name, false_certification_rate=float(rate), reps=int(reps),
                        allowed_max=None if limit is None else float(limit)))
