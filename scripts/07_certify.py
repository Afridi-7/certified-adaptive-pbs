"""
Part III (RQ3): how many cells must be examined to CERTIFY a threshold decision?

Assertion: blast proportion p > theta, over the finite population of M = N*m cells of one specimen.
theta is stated per run (WHO-HAEM5 2022 uses 20% for most AML contexts; the ICC 2022 uses 10% for the
MDS/AML boundary). This certifies a statement about a finite population of cells; it does not establish
diagnostic validity, and the two classification systems are never merged.

Sampling designs (fields are the sampling unit; an expert verifies every cell of a sampled field):
  uniform   : SRSWOR, pi_t = 1/(N-t+1)                                   -> ALPHA applies directly
  prio_ht   : classifier-prioritised without replacement, probabilities recorded, HT-corrected assorter
              (our construction, mixed with uniform at rate lam so weights stay bounded)
  prio_naive: the same prioritised order, but the assorter used AS IF the sample were uniform
              (the error we want to measure: adaptive selection without probability correction)

Also reported: the fixed-convention baselines (examine the first 100 / 200 / 500 cells and call
p_hat > theta), at their own error rates, which is what the certified procedure has to be compared against.
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.composites import CellBank, make_proportion_specimen
from cadpbs.alpha import assorter_bound, alpha_martingale, sample_uniform, sample_prioritised, ht_assorter

ap = argparse.ArgumentParser()
ap.add_argument("--emb", default="results/embeddings_v1")
ap.add_argument("--out", default="results/certify_v1")
ap.add_argument("--N", type=int, default=100)
ap.add_argument("--m", type=int, default=5)
ap.add_argument("--theta", type=float, default=0.10)
ap.add_argument("--p_list", type=float, nargs="+", default=[0.05, 0.08, 0.10, 0.12, 0.15, 0.20])
ap.add_argument("--folds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
ap.add_argument("--episodes", type=int, default=100)
ap.add_argument("--alpha", type=float, default=0.05)
ap.add_argument("--lam", type=float, default=0.5)
ap.add_argument("--seed", type=int, default=20260401)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

E = np.load(f"{a.emb}/embeddings_dinov2_vits14.npy").astype(np.float32)
df = pd.read_csv(f"{a.emb}/cells_features.csv", dtype={"patient_id": str})
bank = CellBank(E, df)
u = assorter_bound(a.theta)

rows = []
t0 = time.time()
for fold in a.folds:
    rng = np.random.default_rng(a.seed + fold)
    for p_true in a.p_list:
        for ep in range(a.episodes):
            s = make_proportion_specimen(bank, fold, a.N, a.m, p_true, rng)
            a_field = (s["n_blast_per_field"] / a.m) / (2 * a.theta)      # field-level assorter, in [0, u]
            r = dict(fold=fold, p_target=p_true, p_true=s["blast_fraction"], theta=a.theta, ep=ep,
                     N=a.N, m=a.m, M=s["M"], blast_reuse=round(s["blast_reuse"], 2), replaced=s["replaced"])
            # --- uniform SRSWOR ---
            order, _ = sample_uniform(a.N, rng)
            _, stop = alpha_martingale(a_field[order], u, 0.5, a.N, alpha=a.alpha)
            r["uniform_cert"] = int(stop > 0); r["uniform_cells"] = (stop if stop > 0 else a.N) * a.m
            # --- prioritised, HT-corrected ---
            o2, pr = sample_prioritised(a.N, s["teacher_max"], rng, a.lam)
            xw, uw = ht_assorter(a_field[o2], pr, a.N, u, a.lam)
            _, stop2 = alpha_martingale(xw, uw, 0.5, a.N, alpha=a.alpha)
            r["prio_ht_cert"] = int(stop2 > 0); r["prio_ht_cells"] = (stop2 if stop2 > 0 else a.N) * a.m
            # --- prioritised, NO correction (deliberately invalid) ---
            _, stop3 = alpha_martingale(a_field[o2], u, 0.5, a.N, alpha=a.alpha)
            r["prio_naive_cert"] = int(stop3 > 0); r["prio_naive_cells"] = (stop3 if stop3 > 0 else a.N) * a.m
            # --- fixed conventions: count the first n cells uniformly and compare p_hat with theta ---
            cells_flat = s["cells"][order].reshape(-1)
            b_flat = bank.is_blast[cells_flat]
            for n_cells in (100, 200, 500):
                n_cells = min(n_cells, s["M"])
                r[f"conv{n_cells}_call"] = int(b_flat[:n_cells].mean() > a.theta)
            rows.append(r)
    print(f"fold {fold} done {time.time()-t0:.0f}s", flush=True)

d = pd.DataFrame(rows)
d.to_csv(f"{a.out}/episodes.csv", index=False)
json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print("episodes:", len(d))
