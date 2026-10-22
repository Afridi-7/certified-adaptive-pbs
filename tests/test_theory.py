"""Checks of the combinatorial results (Proposition 1 including its bounds, Theorem 1 attainment, ranked-search proposition)."""
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


def test_prop1_bounds_full_range():
    """Proposition 1 sandwich bounds over the FULL range 0 <= n <= N (corrected 2026-10-22).
    The lower bound needs the positive part: without it, (1 - n/(N-k+1))^k is positive for even k and n > N-k+1,
    where h = 0 (e.g. N = 10, k = 2, n = 10: 1/81 > 0)."""
    assert h(10, 2, 10) == 0 and (1 - Fraction(10, 9)) ** 2 == Fraction(1, 81)       # the reported counterexample
    for N in range(2, 60):
        for k in range(1, N + 1):
            for n in range(0, N + 1):
                lo = max(Fraction(0), 1 - Fraction(n, N - k + 1)) ** k
                hi = (1 - Fraction(n, N)) ** k
                assert lo <= h(N, k, n) <= hi, (N, k, n)


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
