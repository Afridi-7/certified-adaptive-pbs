"""
Attempts to break the validity claims of the paper (false certification <= alpha at the exact null, for any stopping
time), plus negative controls showing that the suite detects invalid procedures. Each case uses a fixed seed.
"""
import numpy as np
import pytest
from helpers import ALPHA, max_allowed, exact_null_fields, alpha_rejects, betting_rejects, record
from cadpbs.weighted_betting_w import sampler_uniform_w, sampler_prop_w

REPS = 1000


# ---------------------------------------------------------------- ALPHA, uniform sampling of equal fields
@pytest.mark.parametrize("theta,N,m,eta_mode", [(0.10, 100, 5, "default"), (0.20, 100, 5, "default"), (0.10, 60, 1, "default"),
                                                (0.10, 100, 5, "aggressive"), (0.20, 50, 10, "timid")])
def test_alpha_exact_null(theta, N, m, eta_mode):
    rng = np.random.default_rng([1, int(theta * 100), N, m]); u = 1 / (2 * theta)
    eta0 = {"default": None, "aggressive": 0.99 * u, "timid": 0.5 + 1e-3}[eta_mode]
    rej = 0
    for _ in range(REPS):
        B, m_i = exact_null_fields(N, m, theta, rng)
        X = (B / m_i)[rng.permutation(N)] / (2 * theta)
        rej += alpha_rejects(X, u, N, eta0)
    record(f"ALPHA exact null theta={theta} N={N} m={m} eta={eta_mode}", rej / REPS, REPS, max_allowed(REPS))
    assert rej / REPS <= max_allowed(REPS), rej / REPS


# ---------------------------------------------------------------- weighted betting (Shekhar et al.), equal and unequal fields
def _sizes(rng, N):
    """Field sizes like BoneMarrowWSI quadrants (heavy right tail), adjusted so that theta * M is an integer for theta = 0.1, 0.2."""
    s = np.maximum(1, rng.gamma(2.0, 10.0, N).round().astype(int)); s[0] += (-s.sum()) % 10
    return s


@pytest.mark.parametrize("sampler_kind", ["uniform", "prop_random_scores", "prop_oracle", "prop_antioracle"])
@pytest.mark.parametrize("unequal", [False, True])
def test_betting_exact_null(sampler_kind, unequal):
    """Prop. 2(a): weights m_i / M. 'oracle' puts sampling mass on the fields that really have the most blasts (the
    sampler an adversary would choose); 'antioracle' on the fields with the fewest. Validity must hold for any
    predictable sampler that is positive on every unsampled field."""
    theta = 0.10; N = 60; rng = np.random.default_rng([2, int(unequal), ["uniform", "prop_random_scores", "prop_oracle", "prop_antioracle"].index(sampler_kind)])
    rej = 0
    for _ in range(REPS // 2):
        sizes = _sizes(rng, N) if unequal else None
        B, m_i = exact_null_fields(N, 5, theta, rng, sizes=sizes, small_enriched=unequal)
        w = m_i / m_i.sum(); f = B / m_i
        S = {"uniform": np.zeros(N), "prop_random_scores": rng.random(N), "prop_oracle": f + 1e-3,
             "prop_antioracle": 1 - f + 1e-3}[sampler_kind]
        smp = sampler_uniform_w(N) if sampler_kind == "uniform" else sampler_prop_w(w, S)
        rej += betting_rejects(f, w, S, smp, theta, rng)
    record(f"Betting exact null sampler={sampler_kind} unequal={unequal}", rej / (REPS // 2), REPS // 2, max_allowed(REPS // 2))
    assert rej / (REPS // 2) <= max_allowed(REPS // 2), rej / (REPS // 2)


def test_bounded_assorter_exact_null_and_identity():
    """Prop. 2(b): a_i = (B_i - theta m_i + theta m_max) / m_max in [0, 1], mean(a) > theta iff p > theta; ALPHA on a_i
    (null mean theta, i.e. assorter a_i / (2 theta) against 1/2) is valid with unequal sizes and size-dependent blasts."""
    theta = 0.10; N = 60; rng = np.random.default_rng(3); rej = 0
    for _ in range(REPS):
        sizes = _sizes(rng, N); B, m_i = exact_null_fields(N, 5, theta, rng, sizes=sizes, small_enriched=True)
        m_max = int(m_i.max()) + int(rng.integers(0, 10))                       # a valid but possibly loose bound
        a = (B - theta * m_i + theta * m_max) / m_max
        assert a.min() >= -1e-12 and a.max() <= 1 + 1e-12
        assert abs(a.mean() - theta) < 1e-12                                   # exact null in the assorter scale
        X = a[rng.permutation(N)] / (2 * theta); rej += alpha_rejects(X, 1 / (2 * theta), N)
    record("Bounded assorter (Prop. 2b) exact null, size-dependent blasts", rej / REPS, REPS, max_allowed(REPS))
    assert rej / REPS <= max_allowed(REPS), rej / REPS
    for _ in range(2000):                                                      # identity on random populations
        sizes = rng.integers(1, 40, 30); B = rng.binomial(sizes, rng.random()); m_max = int(sizes.max())
        p = B.sum() / sizes.sum(); a = (B - theta * sizes + theta * m_max) / m_max
        assert (a.mean() > theta + 1e-12) == (p > theta + 1e-12)


# ---------------------------------------------------------------- imperfect verification (Remark: threshold inflation)
def test_verifier_margin_remark():
    """A verifier that adds false blasts to at most beta0 of the true non-blasts. At the exact null for the truth, testing
    the verified proportion against theta + beta0 (1 - theta) is valid; testing it against theta is not."""
    theta, beta0, N, m = 0.10, 0.04, 100, 5; rng = np.random.default_rng(4); rej_corr = rej_naive = 0
    for _ in range(REPS):
        B, m_i = exact_null_fields(N, m, theta, rng); M = N * m
        cells = np.zeros(M); start = np.r_[0, np.cumsum(m_i)[:-1]]
        for i in range(N): cells[start[i]:start[i] + B[i]] = 1
        neg = np.flatnonzero(cells == 0); flip = rng.choice(neg, int(beta0 * len(neg)), replace=False)
        ver = cells.copy(); ver[flip] = 1                                      # worst case: all allowed errors are false blasts
        f = ver.reshape(N, m).mean(1); w = np.full(N, 1 / N); smp = sampler_uniform_w(N)
        rej_corr += betting_rejects(f, w, np.zeros(N), smp, theta + beta0 * (1 - theta), rng)
        rej_naive += betting_rejects(f, w, np.zeros(N), smp, theta, rng)
    record("Verifier margin (Remark): corrected threshold", rej_corr / REPS, REPS, max_allowed(REPS))
    record("Verifier margin (Remark): naive threshold [negative control]", rej_naive / REPS, REPS, None)
    assert rej_corr / REPS <= max_allowed(REPS), rej_corr / REPS
    assert rej_naive / REPS > 0.5, rej_naive / REPS                            # negative control


# ---------------------------------------------------------------- negative controls: the suite must catch invalid designs
def test_negative_control_uncorrected_ranking():
    """Ranking fields by a classifier score and analysing them as if sampled uniformly (the paper's 'uncorrected' design)."""
    theta, N, m = 0.10, 100, 5; rng = np.random.default_rng(5); rej = 0
    for _ in range(200):
        B, m_i = exact_null_fields(N, m, theta, rng); f = B / m_i
        S = f + 0.1 * rng.random(N)                                            # an informative classifier
        X = f[np.argsort(-S)] / (2 * theta); rej += alpha_rejects(X, 1 / (2 * theta), N)
    record("Uncorrected ranking [negative control]", rej / 200, 200, None)
    assert rej / 200 > 0.5, rej / 200


def test_negative_control_mean_of_proportions():
    """Averaging field proportions when blasts sit preferentially in small fields (wrong estimand)."""
    theta, N = 0.10, 60; rng = np.random.default_rng(6); rej = 0
    for _ in range(200):
        sizes = _sizes(rng, N); B, m_i = exact_null_fields(N, 5, theta, rng, sizes=sizes, small_enriched=True)
        X = (B / m_i)[rng.permutation(N)] / (2 * theta); rej += alpha_rejects(X, 1 / (2 * theta), N)
    record("Mean of field proportions, small-field enrichment [negative control]", rej / 200, 200, None)
    assert rej / 200 > 0.5, rej / 200


def test_negative_control_cells_as_independent_draws_under_clustering():
    """Cells of uniformly sampled fields analysed as individual draws, with blasts concentrated in a few fields."""
    theta, N, m = 0.10, 100, 5; rng = np.random.default_rng(7); rej = 0
    for _ in range(200):
        cells = np.zeros((N, m)); hot = rng.choice(N, int(theta * N), replace=False); cells[hot] = 1   # all blasts in 10 full fields
        flat = cells[rng.permutation(N)].reshape(-1) / (2 * theta)
        rej += alpha_rejects(flat, 1 / (2 * theta), N * m)
    record("Cells as independent draws, hot fields [negative control]", rej / 200, 200, None)
    assert rej / 200 > 0.2, rej / 200


def test_determinism():
    theta, N = 0.10, 50; f = np.random.default_rng(8).integers(0, 6, N) / 5; w = np.full(N, 1 / N)
    from cadpbs.weighted_betting_w import betting_test_w
    r1 = betting_test_w(f, w, np.zeros(N), sampler_uniform_w(N), theta, rng=np.random.default_rng(9))
    r2 = betting_test_w(f, w, np.zeros(N), sampler_uniform_w(N), theta, rng=np.random.default_rng(9))
    assert r1["stop"] == r2["stop"] and np.array_equal(r1["order"], r2["order"])


def _follow_the_blasts(N, f, boost=20.0):
    """A predictable ADAPTIVE sampler: after a field turns out to contain blasts, its neighbours (fields are on a line)
    get boost-times more sampling mass. It uses only fields already drawn and verified, so it is predictable."""
    def q(remaining, t):
        drawn = np.setdiff1d(np.arange(N), remaining); hot = drawn[f[drawn] > 0]
        near = np.zeros(N, bool)
        for d in (-1, 1): near[np.clip(hot + d, 0, N - 1)] = True
        w = 1.0 + boost * near[remaining]; return w / w.sum()
    return q


@pytest.mark.parametrize("unequal", [False, True])
def test_betting_exact_null_adaptive_sampler_clustered(unequal):
    """Exact null with blasts in contiguous runs of fields (strong spatial clustering) and a sampler that follows them."""
    theta, N = 0.10, 60; rng = np.random.default_rng([12, int(unequal)]); rej = 0; reps = REPS // 2
    for _ in range(reps):
        m_i = _sizes(rng, N) if unequal else np.full(N, 5); M = int(m_i.sum()); nb = int(round(theta * M))
        start = int(rng.integers(0, N)); order = np.r_[np.arange(start, N), np.arange(0, start)]     # fill a contiguous run
        B = np.zeros(N, int); left = nb
        for i in order:
            take = min(left, m_i[i]); B[i] = take; left -= take
            if left == 0: break
        w = m_i / M; f = B / m_i
        rej += betting_rejects(f, w, np.zeros(N), _follow_the_blasts(N, f), theta, rng)
    record(f"Betting exact null, adaptive follow-the-blasts sampler, clustered, unequal={unequal}", rej / reps, reps, max_allowed(reps))
    assert rej / reps <= max_allowed(reps), rej / reps


def test_power_sanity():
    """The valid tests must not pass trivially by never certifying: at p = theta + 0.05 (M = 500) they certify."""
    theta, N, m = 0.10, 100, 5; rng = np.random.default_rng(13); ra = rb = 0; reps = 200
    for _ in range(reps):
        cells = np.zeros(N * m); cells[:int(round(0.15 * N * m))] = 1; rng.shuffle(cells); f = cells.reshape(N, m).mean(1)
        ra += alpha_rejects(f[rng.permutation(N)] / (2 * theta), 1 / (2 * theta), N)
        rb += betting_rejects(f, np.full(N, 1 / N), np.zeros(N), sampler_uniform_w(N), theta, rng)
    record("Power sanity at p = 0.15, theta = 0.10: ALPHA", ra / reps, reps, None)
    record("Power sanity at p = 0.15, theta = 0.10: betting", rb / reps, reps, None)
    assert ra / reps > 0.9 and rb / reps > 0.9, (ra / reps, rb / reps)


# ================================================================ added 2026-10-22 (audit item C5)
# The control-variate branch of the frozen betting test, the defensive mixture sampler and the exploratory hedged bets
# had no direct validity test. Same protocol as above: exact nulls, fixed seeds, Monte Carlo limit max_allowed(reps).
from cadpbs.hedged_betting import hedged_betting_w, KAPPAS
from cadpbs.weighted_betting_w import betting_test_w

SCORES = ("random", "oracle", "antioracle")


def _scores(kind, f, rng):
    N = len(f)
    return {"random": rng.random(N), "oracle": f + 1e-3, "antioracle": 1 - f + 1e-3}[kind]


def _mixture(pi, S, lam):
    """The 50% (or lam) uniform mixture of scripts/19 and scripts/49: q_t = (1 - lam) prop-MS + lam uniform."""
    base = sampler_prop_w(pi, S)
    return lambda rem, t: (1 - lam) * base(rem, t) + lam / len(rem)


def _specimen(rng, N, unequal, theta=0.10):
    sizes = _sizes(rng, N) if unequal else None
    B, m_i = exact_null_fields(N, 5, theta, rng, sizes=sizes, small_enriched=unequal)
    return m_i / m_i.sum(), B / m_i


@pytest.mark.parametrize("score_kind", SCORES)
@pytest.mark.parametrize("unequal", [False, True])
def test_betting_control_variate_exact_null(score_kind, unequal):
    """Frozen betting test with use_cv=True (the 'CV' design of Fig. 3): beta_t is predictable and the bet cap keeps
    1 + lambda Y > 0, so validity must hold for informative, uninformative and adversarial side information."""
    theta, N, reps = 0.10, 60, REPS // 2
    rng = np.random.default_rng([20, int(unequal), SCORES.index(score_kind)]); rej = 0
    for _ in range(reps):
        w, f = _specimen(rng, N, unequal, theta); S = _scores(score_kind, f, rng)
        rej += betting_test_w(f, w, S, sampler_prop_w(w, S), theta, alpha=ALPHA, use_cv=True, rng=rng)["stop"] > 0
    record(f"Betting with control variate, exact null, scores={score_kind}, unequal={unequal}", rej / reps, reps, max_allowed(reps))
    assert rej / reps <= max_allowed(reps), rej / reps


@pytest.mark.parametrize("lam,score_kind,unequal", [(0.5, "oracle", False), (0.5, "antioracle", False),
                                                    (0.25, "oracle", True), (0.25, "antioracle", True)])
def test_mixture_sampler_exact_null(lam, score_kind, unequal):
    """Defensive mixture (H6a): a predictable q_t bounded below by lam / n_remaining, with the frozen betting test."""
    theta, N, reps = 0.10, 60, REPS // 2
    rng = np.random.default_rng([21, int(100 * lam), SCORES.index(score_kind), int(unequal)]); rej = 0
    for _ in range(reps):
        w, f = _specimen(rng, N, unequal, theta); S = _scores(score_kind, f, rng)
        rej += betting_rejects(f, w, S, _mixture(w, S, lam), theta, rng)
    record(f"Betting, {lam:.0%} uniform mixture, exact null, scores={score_kind}, unequal={unequal}", rej / reps, reps, max_allowed(reps))
    assert rej / reps <= max_allowed(reps), rej / reps


HEDGE_SAMPLERS = ("uniform", "prop_oracle", "prop_antioracle", "mix50_antioracle", "follow_the_blasts")


@pytest.mark.parametrize("sampler_kind", HEDGE_SAMPLERS)
@pytest.mark.parametrize("unequal", [False, True])
def test_hedged_betting_exact_null(sampler_kind, unequal):
    """Exploratory hedged bets (average of fixed-fraction betting processes, src/cadpbs/hedged_betting.py)."""
    theta, N, reps = 0.10, 60, REPS // 2
    rng = np.random.default_rng([22, int(unequal), HEDGE_SAMPLERS.index(sampler_kind)]); rej = 0
    for _ in range(reps):
        w, f = _specimen(rng, N, unequal, theta)
        if sampler_kind == "uniform":
            smp = sampler_uniform_w(N)
        elif sampler_kind == "prop_oracle":
            smp = sampler_prop_w(w, f + 1e-3)
        elif sampler_kind == "prop_antioracle":
            smp = sampler_prop_w(w, 1 - f + 1e-3)
        elif sampler_kind == "mix50_antioracle":
            smp = _mixture(w, 1 - f + 1e-3, 0.5)
        else:
            smp = _follow_the_blasts(N, f)
        rej += hedged_betting_w(f, w, smp, theta, alpha=ALPHA, rng=rng)["stop"] > 0
    record(f"Hedged betting, exact null, sampler={sampler_kind}, unequal={unequal}", rej / reps, reps, max_allowed(reps))
    assert rej / reps <= max_allowed(reps), rej / reps


def _zero_mass_oracle(f):
    """A predictable sampler that is NOT positive on all unsampled fields: it puts all mass on the unsampled fields with
    the most blasts (the top half). Z_t then has conditional mean <= the unsampled mass, because f >= 0."""
    def q(rem, t):
        k = max(1, len(rem) // 2); top = np.argsort(-f[rem], kind="stable")[:k]
        out = np.zeros(len(rem)); out[top] = 1.0 / k; return out
    return q


@pytest.mark.parametrize("test_kind", ["kelly", "hedged"])
def test_above_certificate_valid_without_positivity(test_kind):
    """The papers state that positivity of q_t is needed for power and for 'p < theta' certificates but not for the
    validity of 'p > theta'. Check: exact null, unequal fields, sampler with zero mass on half of the unsampled fields."""
    theta, N, reps = 0.10, 60, REPS // 2; rng = np.random.default_rng([23, int(test_kind == "hedged")]); rej = 0
    for _ in range(reps):
        w, f = _specimen(rng, N, True, theta); smp = _zero_mass_oracle(f)
        if test_kind == "kelly":
            rej += betting_rejects(f, w, np.zeros(N), smp, theta, rng)
        else:
            rej += hedged_betting_w(f, w, smp, theta, alpha=ALPHA, rng=rng)["stop"] > 0
    record(f"'Above' certificate without positivity of q_t ({test_kind} bets), exact null", rej / reps, reps, max_allowed(reps))
    assert rej / reps <= max_allowed(reps), rej / reps


def test_hedged_wealth_matches_direct_computation():
    """Implementation check: the hedged stop equals the first time the average of the fixed-fraction products
    prod_s (1 + kappa (Z_s / mu_s - 1)), recomputed independently from the recorded draw order, reaches 1/alpha
    (or the sampled mass exceeds m0)."""
    rng = np.random.default_rng(24)
    for _ in range(200):
        N = int(rng.integers(5, 40)); theta = float(rng.choice([0.1, 0.2]))
        m_i = rng.integers(1, 12, N); w = m_i / m_i.sum(); f = rng.binomial(m_i, rng.uniform(0, 0.4)) / m_i
        S = rng.random(N); seed = int(rng.integers(1 << 30))
        out = hedged_betting_w(f, w, sampler_prop_w(w, S), theta, rng=np.random.default_rng(seed))
        q_fn = sampler_prop_w(w, S); rem = list(range(N)); sampled = 0.0; W = np.ones(len(KAPPAS)); stop = -1
        for t, i in enumerate(out["order"], start=1):
            q = q_fn(np.array(rem), t); qi = q[rem.index(i)]; mu = theta - sampled
            if mu > 0:
                W *= 1 + np.array(KAPPAS) * (w[i] * f[i] / qi / mu - 1)
            sampled += w[i] * f[i]; rem.remove(i)
            if W.mean() >= 1 / ALPHA * (1 - 1e-12) or sampled > theta + 1e-12:
                stop = t; break
        assert stop == out["stop"], (stop, out["stop"])


def test_hedged_power_sanity():
    """Hedged bets must certify a clear alternative: p = 0.15 against theta = 0.10, 100 fields of 5 cells."""
    theta, N, m = 0.10, 100, 5; rng = np.random.default_rng(25); reps = 200; ru = rm = 0
    for _ in range(reps):
        cells = np.zeros(N * m); cells[:int(round(0.15 * N * m))] = 1; rng.shuffle(cells)
        f = cells.reshape(N, m).mean(1); w = np.full(N, 1 / N); S = f + 0.1 * rng.random(N)
        ru += hedged_betting_w(f, w, sampler_uniform_w(N), theta, rng=rng)["stop"] > 0
        rm += hedged_betting_w(f, w, _mixture(w, S, 0.5), theta, rng=rng)["stop"] > 0
    record("Hedged power sanity at p = 0.15, theta = 0.10: uniform", ru / reps, reps, None)
    record("Hedged power sanity at p = 0.15, theta = 0.10: 50% mixture", rm / reps, reps, None)
    assert ru / reps > 0.9 and rm / reps > 0.9, (ru / reps, rm / reps)
