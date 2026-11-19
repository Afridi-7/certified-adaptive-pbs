"""Checks of Proposition 4 (docs/part3_limit_proof.md): limit for certifying a threshold decision."""
import itertools
from fractions import Fraction
from math import comb, floor

import numpy as np
import pytest
from scipy.stats import hypergeom


def mp_power(M, B1, B0, n, alpha):
    xs = np.arange(n + 1); p0 = hypergeom.pmf(xs, M, B0, n); p1 = hypergeom.pmf(xs, M, B1, n)
    tail0 = np.cumsum(p0[::-1])[::-1]
    c = next(x for x in range(n + 2) if x == n + 1 or tail0[x] <= alpha)
    size_strict = tail0[c] if c <= n else 0.0
    gamma = 0.0 if c == 0 or p0[c - 1] == 0 else min(max((alpha - size_strict) / p0[c - 1], 0.0), 1.0)
    return float((p1[c:].sum() if c <= n else 0.0) + (gamma * p1[c - 1] if c >= 1 else 0.0))


def exact_power(M, B1, B0, n, alpha):
    xs = np.arange(n + 1); tail0 = np.cumsum(hypergeom.pmf(xs, M, B0, n)[::-1])[::-1]
    c = next(x for x in range(n + 2) if x == n + 1 or tail0[x] <= alpha)
    return float(hypergeom.sf(c - 1, M, B1, n)) if c <= n else 0.0


@pytest.mark.parametrize("M,th,n", [(50, 0.1, 10), (200, 0.2, 40), (2000, 0.1, 200), (2000, 0.2, 500)])
def test_bound_has_size_alpha_and_dominates_exact_test(M, th, n):
    B0 = floor(th * M)
    assert abs(mp_power(M, B0, B0, n, 0.025) - 0.025) < 1e-12           # randomised test has size exactly alpha
    for B1 in range(B0 + 1, min(M, B0 + M // 5) + 1, max(1, M // 100)):
        assert mp_power(M, B1, B0, n, 0.025) >= exact_power(M, B1, B0, n, 0.025) - 1e-12


def test_reduction_exact_enumeration():
    """An adaptive deterministic policy that uses uninformative features and past labels sees exactly the
    hypergeometric distribution of blasts, by enumeration of all arrangements (rational arithmetic)."""
    M, B, n = 7, 3, 4
    feat = [0.3, 2.0, -1.0, 0.7, 1.5, -0.2, 0.9]
    def policy(seen, obs):
        avail = [i for i in range(M) if i not in seen]
        if obs and obs[-1] == 1:                                     # after a blast: the nearest feature
            return min(avail, key=lambda i: abs(feat[i] - feat[seen[-1]]))
        return max(avail, key=lambda i: feat[i])                     # otherwise the largest feature
    counts = {}
    arrangements = list(itertools.combinations(range(M), B))
    for blasts in arrangements:
        lab = [1 if i in blasts else 0 for i in range(M)]; seen, obs = [], []
        for _ in range(n):
            j = policy(seen, obs); seen.append(j); obs.append(lab[j])
        counts[sum(obs)] = counts.get(sum(obs), 0) + 1
    for x in range(n + 1):
        assert Fraction(counts.get(x, 0), len(arrangements)) == Fraction(comb(B, x) * comb(M - B, n - x), comb(M, n))


def test_sequential_valid_test_never_beats_bound():
    rng = np.random.default_rng(20261119)
    def seq_test(seq, Mp, B0, alpha):
        W, x = 1.0, 0
        for t, s in enumerate(seq):
            mu = (B0 - x) / (Mp - t)
            if mu < 0: return t
            W *= 1 + (0.5 / mu if mu > 0 else 0.0) * (s - mu); x += s
            if W >= 1 / alpha: return t + 1
        return 0
    for Mp, B0, B1, n in [(40, 8, 16, 12), (100, 10, 25, 30)]:
        pop = np.r_[np.ones(B1), np.zeros(Mp - B1)]
        rate = np.mean([0 < seq_test(rng.permutation(pop)[:n], Mp, B0, 0.05) <= n for _ in range(4000)])
        assert rate <= mp_power(Mp, B1, B0, n, 0.05) + 3 * np.sqrt(0.25 / 4000)
