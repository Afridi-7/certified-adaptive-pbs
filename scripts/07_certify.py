"""
Part III (RQ3): cells examined to certify p > theta over the finite population of one specimen.

theta is set per run (WHO 5th ed. 2022: 20%; ICC 2022: 10%); the two systems are not merged.
Fields are the sampling unit; an expert verifies every cell of a sampled field.

Designs:
  uniform     SRSWOR + ALPHA
  bet_*       Shekhar et al. (2023) weighted betting: uniform, prop-MS, prop-MS + control variate
  prio_naive  prioritised order analysed as if uniform (invalid; measures the error)
  conv*       first 100 / 200 / 500 cells, call p_hat > theta (no guarantee)
A former 'prio_ht' design was withdrawn as invalid; see docs/certification_v1.md.
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.composites import CellBank, make_proportion_specimen
from cadpbs.alpha import assorter_bound, alpha_martingale, sample_uniform, sample_prioritised
from cadpbs.weighted_betting import betting_test, sampler_uniform, sampler_prop_side

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
ap.add_argument("--seed", type=int, default=20260501)   # Part III seed family (disjoint from Part I confirmatory 20260401)
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
            # --- Shekhar et al. (2023) betting test: uniform, classifier-driven prop-MS, and prop-MS + control variate ---
            S = s["teacher"] / (2 * a.theta)          # side information: classifier-predicted field assorter
            for nm, smp, cv in [("bet_uniform", sampler_uniform(a.N), False),
                                ("bet_prop", sampler_prop_side(S), False),
                                ("bet_prop_cv", sampler_prop_side(S), True)]:
                res = betting_test(a_field, S, smp, 0.5, alpha=a.alpha, use_cv=cv, rng=rng)
                r[f"{nm}_cert"] = int(res["stop"] > 0); r[f"{nm}_cells"] = (res["stop"] if res["stop"] > 0 else a.N) * a.m
            o2, pr = sample_prioritised(a.N, s["teacher_max"], rng, a.lam)
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
