"""
PREREG Amendment 3, E3.2 (RQ4, RQ3): a defensive mixture for classifier-driven sampling under side-information shift.

Motivation (exploratory, confirmatory_v2.md): on RV-PBS, classifier-driven (prop-MS) sampling needed MORE cells than
uniform sampling at p = 0.20 (median 190 vs 110). Mixture sampler:
    q_t(i) = (1 - lam) * S(i) / sum_{j remaining} S(j) + lam / n_remaining,   lam in {0.25, 0.50}
This is a predictable distribution that is positive on every remaining field, so the Shekhar et al. test stays valid
(their Theorem 1 allows any such q_t). The mixture bounds how far q_t can move from uniform: q_t(i) >= lam / n_remaining.

Sources:
  rvpbs : RV-PBS crops; scores from the logistic probe trained on ALL KU-Optofil cells (as in scripts/16), i.e. shifted
  ku    : KU-Optofil; cross-fitted scores; blasts from all 489, background from blast-free patients (as scripts/12)
Specimen: N = 100 fields x m = 5 cells, cells shuffled (scattered). theta = 0.10 (ICC 2022). All designs are paired on
the same specimen: uniform, prop, mix25, mix50 (Shekhar betting test, c = 0.75, no control variate).
The per-episode RNG is default_rng([seed, source index, round(1000 p), episode]).
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.weighted_betting import betting_test, sampler_uniform, sampler_prop_side

SOURCES = ["rvpbs", "ku"]


def sampler_mixture(S, lam, floor=1e-3):
    S = np.maximum(np.asarray(S, float), floor)
    def q(remaining, t):
        w = S[remaining]; return (1 - lam) * w / w.sum() + lam / len(remaining)
    return q


ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True); ap.add_argument("--source", choices=SOURCES, required=True)
ap.add_argument("--p_list", type=float, nargs="+", required=True)
ap.add_argument("--episodes", type=int, default=500); ap.add_argument("--null_episodes", type=int, default=1000)
ap.add_argument("--seed", type=int, required=True); ap.add_argument("--chunk", type=int, default=250)
ap.add_argument("--max_seconds", type=float, default=1e9)
ap.add_argument("--ku_features", default="results/embeddings_v1/cells_features.csv")
ap.add_argument("--ku_emb", default="results/embeddings_v1/embeddings_dinov2_vits14.npy")
ap.add_argument("--rv_dir", default="results/external_v1")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
th = 0.10; N, m = 100, 5; M = N * m; T0 = time.time()

ku = pd.read_csv(a.ku_features, dtype={"patient_id": str})
if a.source == "ku":
    y = ku.is_blast.values.astype(bool); score = ku.probe_p_oof.values
    bp = set(ku.loc[y, "patient_id"])
    BL = np.flatnonzero(y); BG = np.flatnonzero(~y & ~ku.patient_id.isin(bp).values)
else:
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    Eku = np.load(a.ku_emb).astype(np.float32)
    Erv = np.load(f"{a.rv_dir}/rvpbs_embeddings_dinov2_vits14.npy").astype(np.float32)
    rv = pd.read_csv(f"{a.rv_dir}/rvpbs_features.csv")
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, C=0.1, class_weight="balanced"))
    clf.fit(Eku, ku.is_blast.values)
    y = rv.is_blast.values.astype(bool); score = clf.predict_proba(Erv)[:, 1]
    BL = np.flatnonzero(y); BG = np.flatnonzero(~y)
si = SOURCES.index(a.source)

for p in a.p_list:
    E = a.null_episodes if abs(p - th) < 1e-9 else a.episodes
    nb = int(round(p * M)); assert nb <= len(BL) and M - nb <= len(BG)
    for s0 in range(0, E, a.chunk):
        path = f"{a.out}/ep_{a.source}_p{p:.2f}_e{s0:04d}.csv"
        if os.path.exists(path): continue
        if time.time() - T0 > a.max_seconds: print("time budget reached; rerun to resume"); sys.exit(0)
        rows = []
        for ep in range(s0, min(s0 + a.chunk, E)):
            rng = np.random.default_rng([a.seed, si, int(round(p * 1000)), ep])
            cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(BG, M - nb, replace=False)]
            rng.shuffle(cells); cells = cells.reshape(N, m)
            f = y[cells].mean(1) / (2 * th); S = score[cells].mean(1) / (2 * th)
            r = dict(source=a.source, p=p, ep=ep)
            for nm, smp in [("uniform", sampler_uniform(N)), ("prop", sampler_prop_side(S)),
                            ("mix25", sampler_mixture(S, 0.25)), ("mix50", sampler_mixture(S, 0.50))]:
                o = betting_test(f, S, smp, 0.5, alpha=0.05, rng=rng)
                r[f"{nm}_cert"] = int(o["stop"] > 0); r[f"{nm}_cells"] = (o["stop"] if o["stop"] > 0 else N) * m
            rows.append(r)
        pd.DataFrame(rows).to_csv(path, index=False)
        print(f"{a.source} p={p:.2f} episodes {s0}-{min(s0 + a.chunk, E) - 1} done ({time.time() - T0:.0f}s)", flush=True)
print("ALL DONE")
