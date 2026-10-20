"""Checks of the combinatorial results (Proposition 1, Theorem 1 attainment, ranked-search proposition)."""
from fractions import Fraction
from math import comb
import numpy as np


def h(N, k, n):
    return Fraction(comb(N - k, n), comb(N, n)) if n <= N - k else Fraction(0)


def test_h_product_form_and_k1():
    for N in range(2, 60):
        for k in range(1, min(N, 6)):
            for n in range(0, N + 1):
                prod = Fraction(1)
                for i in range(k):
                    prod *= Fraction(max(N - n - i, 0), N - i)
                assert h(N, k, n) == prod
        for n in range(N + 1):
            assert h(N, 1, n) == Fraction(N - n, N)


def test_uniform_sampling_attains_h():
    rng = np.random.default_rng(10); N, k, n, R = 40, 3, 12, 40000
    C = np.arange(k); miss = sum(not np.isin(rng.permutation(N)[:n], C).any() for _ in range(R)) / R
    exact = float(h(N, k, n)); se = (exact * (1 - exact) / R) ** 0.5
    assert abs(miss - exact) < 4 * se


def test_ranked_search_formula():
    """P(miss) = E[ P(Bin(N-1, 1 - F0(S_c)) >= n) ] for a top-n scan, Gaussian scores with a shift."""
    from scipy.stats import norm, binom
    rng = np.random.default_rng(11); N, n, shift, R = 50, 8, 1.0, 20000
    s_c = rng.normal(shift, 1, R); others = rng.normal(0, 1, (R, N - 1))
    mc = ((others > s_c[:, None]).sum(1) >= n).mean()
    grid = np.linspace(-6, 8, 4001); dens = norm.pdf(grid, shift, 1)
    formula = np.trapezoid(binom.sf(n - 1, N - 1, 1 - norm.cdf(grid)) * dens, grid)
    assert abs(mc - formula) < 4 * (formula * (1 - formula) / R) ** 0.5
