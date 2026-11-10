"""
Review of 2026-11-10 (exploratory; RQ3, RQ4): designs that use the classifier without letting it drive the sampling.

The paper compares uniform sampling (no classifier) with prop-MS sampling (classifier decides the design) and a 50%
mixture. A reviewer will ask for the two classical alternatives. On the SAME specimens (paired) this script runs
  alpha_uniform   : ALPHA, uniform field sampling (frozen cadpbs.alpha)                       -- cost reference
  hedge_uniform   : hedged betting, uniform sampling            (cadpbs.hedged_betting)       -- as in scripts/49
  hedge_prop      : hedged betting, prop-MS sampling                                          -- as in scripts/49
  hedge_mix50     : hedged betting, 50% uniform / 50% prop-MS                                 -- as in scripts/49
  cv_uniform      : hedged betting with a prediction control variate, uniform sampling        -- NEW
  cv_mix50        : the same with the 50% mixture                                             -- NEW
  hedge_strata    : hedged betting, score-stratified sampling (4 strata, root allocation)     -- NEW
  cv_strata       : the same with the control variate                                         -- NEW
(cadpbs.cv_hedged_betting; validity checks in scripts/54_cv_validity_checks.py).

Sources.
  ku, mll23, rvpbs : real cells, loaded and assembled exactly as in scripts/49_baseline_reference.py.
  synth_ku, synth_mll23, synth_rvpbs : a SYNTHETIC score model, for checking the code and previewing the comparison
      where the cell scores are not available. Latent score z ~ N(b, 1) for non-blasts and N(d, s^2) for blasts, with
      (b, d, s) chosen to match the reported AUC, the sensitivity at the KU-Optofil 95%-specificity threshold and the
      false-positive rate of each source; probability score = logistic(k (z - t)) with t fixed by the reported
      threshold value (0.031) and one free steepness k (--k; default 3, chosen so that synth_ku reproduces the stored
      KU-Optofil medians of hedged prop-MS and the hedged 50% mixture in results/baseline_reference_v1). No patient
      effects, no morphology: these runs are NOT evidence about the real sources and must not be reported as such.
Specimen: N = 100 fields x m = 5 cells, scattered; theta = 0.10 (ICC 2022) for every source, theta = 0.20 (WHO 5th
ed.) for ku / synth_ku. Exact null p = theta: --n_null specimens; p > theta: --n_alt each.
Per-specimen RNG: default_rng([20261110, source index, theta index, round(1000 p), episode]); seed family 20261110 is
new. Uncertified specimens count all 500 cells. Costs are counts of expert-verified cells.

Usage: python scripts/53_shift_robust_designs.py --source ku --out results/shift_robust_v1
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import binomtest, beta, wilcoxon, norm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import sampler_uniform_w, sampler_prop_w
from cadpbs.hedged_betting import hedged_betting_w
from cadpbs.cv_hedged_betting import hedged_cv_betting_w, sampler_strata_w, sampler_mix_w

SOURCES = ["ku", "mll23", "rvpbs", "synth_ku", "synth_mll23", "synth_rvpbs"]
DESIGNS = ["alpha_uniform", "hedge_uniform", "hedge_prop", "hedge_mix50", "cv_uniform", "cv_mix50", "hedge_strata", "cv_strata"]
T_KU = 0.0310392350137264
# synthetic model: (non-blast mean b, blast mean d, blast sd s), solved from the reported AUC / sensitivity / FP rate:
#   AUC = Phi((d - b) / sqrt(1 + s^2)), sensitivity = Phi((d - 1.645) / s), FP = Phi(b - 1.645)
#   ku: 0.979, 0.906, 0.050 (the wider of two solutions; recall varied 0.39-1.00 across blast patients);
#   mll23: 0.915, 0.631, 0.060; rvpbs: 0.786, 0.327, 0.057
SYNTH = {"synth_ku": (0.0, 3.446, 1.368), "synth_mll23": (0.087, 1.953, 0.922), "synth_rvpbs": (0.063, 1.191, 1.013)}

ap = argparse.ArgumentParser()
ap.add_argument("--source", choices=SOURCES, required=True); ap.add_argument("--out", required=True)
ap.add_argument("--seed", type=int, default=20261110); ap.add_argument("--n_null", type=int, default=1000)
ap.add_argument("--n_alt", type=int, default=400); ap.add_argument("--n_boot", type=int, default=4000)
ap.add_argument("--k", type=float, default=3.0, help="steepness of the synthetic probability score (synthetic sources only)")
ap.add_argument("--designs", nargs="+", default=DESIGNS)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); T0 = time.time()
ROOT = os.path.join(HERE, ".."); synthetic = a.source.startswith("synth")

if not synthetic:
    ku = pd.read_csv(f"{ROOT}/results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})

    def ku_probe():
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler
        from sklearn.pipeline import make_pipeline
        Eku = np.load(f"{ROOT}/results/embeddings_v1/embeddings_dinov2_vits14.npy").astype(np.float32)
        return make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, C=0.1, class_weight="balanced")).fit(Eku, ku.is_blast.values)

    if a.source == "ku":
        y = ku.is_blast.values.astype(bool); score = ku.probe_p_oof.values
        bp = set(ku.loc[y, "patient_id"]); BL = np.flatnonzero(y); BG = np.flatnonzero(~y & ~ku.patient_id.isin(bp).values)
    elif a.source == "mll23":
        mm = pd.read_csv(f"{ROOT}/data/mll23/mll23_cells.csv"); Em = np.load(f"{ROOT}/data/mll23/mll23_embeddings_dinov2_vits14.npy").astype(np.float32)
        keep = ~mm.cls.str.contains("atypical|smudge|normoblast", case=False, regex=True)
        Em, mm = Em[keep.values], mm[keep].reset_index(drop=True)
        y = mm.cls.str.contains("myeloblast", case=False).values; score = ku_probe().predict_proba(Em)[:, 1]
        BL = np.flatnonzero(y); BG = np.flatnonzero(~y)
    else:
        rv = pd.read_csv(f"{ROOT}/results/external_v1/rvpbs_features.csv")
        Erv = np.load(f"{ROOT}/results/external_v1/rvpbs_embeddings_dinov2_vits14.npy").astype(np.float32)
        y = rv.is_blast.values.astype(bool); score = ku_probe().predict_proba(Erv)[:, 1]
        BL = np.flatnonzero(y); BG = np.flatnonzero(~y)
    recall = float((score[y] >= T_KU).mean()); fp = float((score[~y] >= T_KU).mean())
else:
    b0, d0, s0 = SYNTH[a.source]; t0 = 1.645 - np.log(T_KU / (1 - T_KU)) / a.k      # logistic(k (1.645 - t0)) = T_KU
    recall = float(norm.cdf((d0 - 1.645) / s0)); fp = float(norm.cdf(b0 - 1.645))
grid = [(0.10, [0.10, 0.12, 0.15, 0.20])] + ([(0.20, [0.20, 0.22, 0.25, 0.30])] if a.source in ("ku", "synth_ku") else [])
print(f"{a.source}: sensitivity at the KU threshold {recall:.3f}, FP {fp:.3f}", flush=True)

N, m = 100, 5; M = N * m; pi = np.full(N, 1.0 / N); si = SOURCES.index(a.source)


def build(nb, rng):
    """Blast indicator and probability score of the M cells of one specimen, as (N, m) arrays, cells shuffled."""
    if synthetic:
        lab = np.r_[np.ones(nb, bool), np.zeros(M - nb, bool)]
        z = np.r_[d0 + s0 * rng.standard_normal(nb), b0 + rng.standard_normal(M - nb)]
        sc = 1.0 / (1.0 + np.exp(-a.k * (z - t0))); o = rng.permutation(M)
        return lab[o].reshape(N, m), sc[o].reshape(N, m)
    cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(BG, M - nb, replace=False)]
    rng.shuffle(cells); cells = cells.reshape(N, m)
    return y[cells], score[cells]


rows = []
for ti, (th, plist) in enumerate(grid):
    for p in plist:
        nb = int(round(p * M)); E = a.n_null if abs(p - th) < 1e-9 else a.n_alt
        for ep in range(E):
            rng = np.random.default_rng([a.seed, si, ti, int(round(p * 1000)), ep])
            lab, sc = build(nb, rng); f = lab.mean(1).astype(float); S = sc.mean(1)
            prop = sampler_prop_w(pi, S); strata = sampler_strata_w(pi, S)
            run = {
                "alpha_uniform": lambda: alpha_martingale(f[rng.permutation(N)] / (2 * th), 1 / (2 * th), 0.5, N, alpha=0.05)[1],
                "hedge_uniform": lambda: hedged_betting_w(f, pi, sampler_uniform_w(N), th, rng=rng)["stop"],
                "hedge_prop": lambda: hedged_betting_w(f, pi, prop, th, rng=rng)["stop"],
                "hedge_mix50": lambda: hedged_betting_w(f, pi, sampler_mix_w(prop, 0.5), th, rng=rng)["stop"],
                "cv_uniform": lambda: hedged_cv_betting_w(f, pi, S, sampler_uniform_w(N), th, rng=rng)["stop"],
                "cv_mix50": lambda: hedged_cv_betting_w(f, pi, S, sampler_mix_w(prop, 0.5), th, rng=rng)["stop"],
                "hedge_strata": lambda: hedged_betting_w(f, pi, strata, th, rng=rng)["stop"],
                "cv_strata": lambda: hedged_cv_betting_w(f, pi, S, strata, th, rng=rng)["stop"],
            }
            r = dict(source=a.source, theta=th, p=p, ep=ep)
            for d in a.designs:                                   # fixed order, so the random stream is reproducible
                s = run[d](); r[f"{d}_cert"] = int(s > 0); r[f"{d}_cells"] = (s if s > 0 else N) * m
            rows.append(r)
        print(f"  theta={th} p={p}: {E} specimens ({time.time() - T0:.0f}s)", flush=True)

D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes_{a.source}.csv", index=False)
rb = np.random.default_rng([a.seed, si, 999])
out = dict(source=a.source, synthetic=synthetic, recall_at_ku_threshold=recall, fp_at_ku_threshold=fp, validity={}, median_cells={}, ratios={})
for (th, p), g in D.groupby(["theta", "p"]):
    key = f"theta{th:.2f}_p{p:.2f}"
    if abs(p - th) < 1e-9:
        for d in a.designs:
            k = int(g[f"{d}_cert"].sum()); n = len(g)
            out["validity"][f"{key}_{d}"] = dict(rate=k / n, cp95_upper=float(beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0,
                                                 p_excess=float(binomtest(k, n, 0.05, alternative="greater").pvalue))
        continue
    out["median_cells"][key] = {d: float(g[f"{d}_cells"].median()) for d in a.designs}
    if "alpha_uniform" not in a.designs: continue
    z = g["alpha_uniform_cells"].values
    for num in [d for d in a.designs if d != "alpha_uniform"]:
        x = g[f"{num}_cells"].values; bs = []
        for _ in range(a.n_boot):
            i = rb.integers(0, len(g), len(g)); bs.append(np.median(x[i]) / np.median(z[i]))
        try: wp = float(wilcoxon(x, z).pvalue)
        except ValueError: wp = 1.0
        out["ratios"][f"{key}_{num}_over_alpha_uniform"] = dict(ratio=float(np.median(x) / np.median(z)),
                                                               ci95=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))], wilcoxon_p=wp)
import platform, scipy
cfg = dict(config=vars(a), env=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__))
json.dump(cfg, open(f"{a.out}/config_{a.source}.json", "w"), indent=1)          # run settings (not compared by verify)
json.dump(out, open(f"{a.out}/summary_{a.source}.json", "w"), indent=1)         # numbers only
print(pd.DataFrame(out["median_cells"]).T.to_string())
print("validity:", {k: v["rate"] for k, v in out["validity"].items()})
print(f"DONE {a.source} in {time.time() - T0:.0f}s", flush=True)
