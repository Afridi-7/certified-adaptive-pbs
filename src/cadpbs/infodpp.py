"""
Faithful reimplementation of the InfoDPP-PAC selection + halting rule.

Source: Mittal, Srivastava & Chauhan, arXiv:2608.23574v1 (2026), Sections 4.1-4.4, Algorithm 1.
No code was released; every choice not fixed by the paper is marked CHOICE and logged in docs/DEVIATIONS.md.

Notation (paper -> code)
  G            pool of n units (here: fields of one specimen), embeddings X (n x d)
  D0           seed set, |D0| = m0, drawn uniformly at random               (paper: m0 ~ 0.01 n)
  T            teacher, giving pseudo-labels y_i = T(x_i) + noise           (CHOICE: cross-fitted classifier score)
  k(.,.)       PSD kernel                                                    (CHOICE: RBF, median heuristic length-scale)
  beta = 1/sn2 inverse observation-noise variance
  K            GP POSTERIOR kernel after conditioning on the seed (Eq. 2), used in Eqs. 3 and 5
  mu_i, sig_i  GP posterior mean / sd (Eqs. 1-2); slide-wise affine-normalised to [0,1] -> mubar, sigbar
  q_i          alpha1*mubar_i + alpha2*sigbar_i                              (Eq. 3)
  F(S)         lambda1 * logdet(I + beta K_S) + lambda2 * sum_S q_i          (Eq. 3)
  Delta_k(i)   blended marginal gain, used for SELECTION                     (Eq. 4)
  Dinfo_k(i)   log(1 + beta * sigma^2_{i,beta}(S_{k-1}))                     (Eq. 5)
  Gamma_k      max over unselected i of Dinfo_k(i), used for STOPPING ONLY   (Eq. 6)
  halt when    Gamma_k <= tau + eps_{k,n},  eps = rho*sqrt(log(pi^2 n k^2/(3 delta))/(2 m0))   (Eq. 7)
  r            q-weighted centroid of unselected units                       (Eq. 8)

Key structural fact (used by Kill-test 1): Gamma_k depends on X, the seed LOCATIONS and S only.
It contains no term in y or q, so the stopping time cannot respond to relevance except through which
units the blended gain has already selected.

Tuned weights from the paper (Sec. 6.2): lambda1=1.75, lambda2=0.25, alpha1=0.5, alpha2=1.5.
Not reported anywhere in the paper: tau, delta, rho, the kernel, the length-scale, sn2, and m0 beyond "~0.01n".
"""
from dataclasses import dataclass, asdict, field
import numpy as np


@dataclass
class InfoDPPConfig:
    lambda1: float = 1.75      # paper Sec. 6.2 (tuned)
    lambda2: float = 0.25
    alpha1: float = 0.5
    alpha2: float = 1.5
    sn2: float = 0.01          # CHOICE: observation-noise variance; beta = 1/sn2
    tau: float = None          # CHOICE: calibrated on negative specimens (see calibrate_tau)
    delta: float = 0.05        # CHOICE: PAC failure probability
    rho: float = 1.0           # CHOICE: Assumption 1 bounded-difference constant; not estimable from the paper
    m0_frac: float = 0.01      # paper: "typically m0 ~ 0.01 n"
    m0_min: int = 5            # CHOICE: floor, since 0.01n < 1 for small specimens
    kernel: str = "rbf"        # CHOICE
    length_scale: float = None # CHOICE: None -> median heuristic on the pool
    max_steps: int = None      # cap; None -> n
    seed: int = 0


def _rbf(X, ell):
    sq = (X ** 2).sum(1)
    d2 = np.maximum(sq[:, None] + sq[None, :] - 2 * X @ X.T, 0.0)
    return np.exp(-d2 / (2 * ell ** 2))


def median_length_scale(X, rng=None, cap=2000):
    idx = np.arange(len(X)) if len(X) <= cap else np.random.default_rng(0).choice(len(X), cap, replace=False)
    Z = X[idx]; sq = (Z ** 2).sum(1)
    D = np.sqrt(np.maximum(sq[:, None] + sq[None, :] - 2 * Z @ Z.T, 0.0))
    v = D[np.triu_indices(len(Z), 1)]
    return float(np.median(v[v > 0]))


def gp_posterior(X, teacher, cfg, rng):
    """Seed phase + inference phase (paper Sec. 4.2). Returns posterior kernel K, mean mu, sd sig, seed idx."""
    n = len(X)
    m0 = max(cfg.m0_min, int(round(cfg.m0_frac * n)))
    m0 = min(m0, n)
    D0 = rng.choice(n, m0, replace=False)
    ell = cfg.length_scale or median_length_scale(X)
    Kf = _rbf(X, ell)
    A = Kf[np.ix_(D0, D0)] + cfg.sn2 * np.eye(m0)
    Ainv_kt = np.linalg.solve(A, Kf[D0, :])            # (m0, n)
    K_post = Kf - Kf[:, D0] @ Ainv_kt                  # Eq. 2, posterior covariance
    K_post = (K_post + K_post.T) / 2
    y = np.asarray(teacher, float)[D0]                 # y_i = T(x_i) (+ noise, already in the teacher)
    mu = Kf[:, D0] @ np.linalg.solve(A, y)             # Eq. 1
    sig = np.sqrt(np.maximum(np.diag(K_post), 0.0))
    return K_post, mu, sig, D0, ell


def _affine01(v):
    lo, hi = float(np.min(v)), float(np.max(v))
    return np.zeros_like(v) if hi - lo < 1e-12 else (v - lo) / (hi - lo)


def run_infodpp(X, teacher, cfg: InfoDPPConfig, rng=None):
    """Algorithm 1. Returns a dict with the selected order, per-step signals and the stop step."""
    rng = rng or np.random.default_rng(cfg.seed)
    n = len(X)
    K, mu, sig, D0, ell = gp_posterior(X, teacher, cfg, rng)
    q = cfg.alpha1 * _affine01(mu) + cfg.alpha2 * _affine01(sig)     # Eq. 3
    beta = 1.0 / cfg.sn2
    m0 = len(D0)
    max_steps = cfg.max_steps or n

    var = np.diag(K).copy()                    # sigma^2_{i,beta}(S) with S = {} is K_ii
    L = np.zeros((n, max_steps))               # incremental Cholesky factors
    chosen = np.zeros(n, bool)
    order, gam, dinfo_sel, blended_sel, eps_hist, q_sel = [], [], [], [], [], []
    stop_reason = "budget"; dinfo0 = None
    for k in range(1, max_steps + 1):
        dinfo = np.log1p(beta * np.maximum(var, 0.0))                 # Eq. 5
        dinfo[chosen] = -np.inf
        Gamma = float(dinfo.max())                                    # Eq. 6 (stopping statistic)
        if k == 1: dinfo0 = np.log1p(beta * np.maximum(var, 0.0)).copy()   # per-field IG before any selection
        eps = cfg.rho * np.sqrt(np.log(np.pi ** 2 * n * k ** 2 / (3 * cfg.delta)) / (2 * m0))  # Eq. 7
        gam.append(Gamma); eps_hist.append(float(eps))
        if cfg.tau is not None and Gamma <= cfg.tau + eps:
            stop_reason = "pac"
            break                                                     # halt BEFORE adding (Alg. 1 line 11)
        blended = cfg.lambda1 * dinfo + cfg.lambda2 * np.where(chosen, -np.inf, q)   # Eq. 4
        j = int(np.argmax(blended))
        order.append(j); chosen[j] = True
        dinfo_sel.append(float(dinfo[j])); blended_sel.append(float(blended[j])); q_sel.append(float(q[j]))
        # rank-one update of the noisy posterior variance
        c = (K[:, j] - L[:, :k - 1] @ L[j, :k - 1]) / np.sqrt(var[j] + cfg.sn2)
        L[:, k - 1] = c
        var = np.maximum(var - c ** 2, 0.0)
        if chosen.all():
            break
    S = np.array(order, int)
    unsel = np.flatnonzero(~chosen)
    r = (q[unsel] @ X[unsel] / q[unsel].sum()) if len(unsel) and q[unsel].sum() > 0 else np.zeros(X.shape[1])  # Eq. 8
    return dict(order=S, k_star=len(S), dinfo0=dinfo0, stop_reason=stop_reason, seed_idx=D0, length_scale=ell,
                gamma=np.array(gam), eps=np.array(eps_hist), dinfo_selected=np.array(dinfo_sel),
                blended_selected=np.array(blended_sel), q_selected=np.array(q_sel), q=q, mu=mu, sig=sig,
                residual=r, config=asdict(cfg))


def calibrate_tau(specimens, cfg: InfoDPPConfig, target_frac=0.163, iters=25, rng=None):
    """CHOICE: the paper never reports tau. We calibrate it on NEGATIVE (blast-free) calibration specimens so
    that the median stop lands at the paper's reported operating point (83.7% fewer than the full budget,
    i.e. ~16.3%). Calibration specimens must be disjoint from the evaluation specimens."""
    rng = rng or np.random.default_rng(cfg.seed)
    lo, hi = 1e-6, 5.0
    for _ in range(iters):
        mid = np.sqrt(lo * hi)
        c = InfoDPPConfig(**{**asdict(cfg), "tau": mid})
        ks = [run_infodpp(X, t, c, np.random.default_rng(1000 + i))["k_star"] for i, (X, t) in enumerate(specimens)]
        if np.median(ks) > target_frac * len(specimens[0][0]):
            lo = mid            # stopping too late -> raise tau
        else:
            hi = mid
    return float(np.sqrt(lo * hi))
