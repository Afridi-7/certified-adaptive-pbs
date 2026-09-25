"""
RQ4 (exploratory): does patient-dependent classifier error break the certificate, or only raise its cost?
Specimens (N=100 fields x m=5 cells, theta=0.10) whose blasts all come from EITHER the low-recall source patients
OR the high-recall ones (grouped by the per-patient recall in results/rq4_error_v1). Background cells from all
blast-free patients. All scores are cross-fitted, so each cell's score comes from a model that never saw its patient.
Compared: uniform sampling vs classifier-driven sampling (Shekhar et al. prop-MS betting test).
Validity is checked at p = theta; cost above it. Seed family 20260601.
"""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from cadpbs.weighted_betting import betting_test, sampler_uniform, sampler_prop_side

OUT = "results/rq4_certify_v1"; os.makedirs(OUT, exist_ok=True)
df = pd.read_csv("results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})
rec = json.load(open("results/rq4_error_v1/summary.json"))["blast_recall_by_patient"]
low = [k for k, v in rec.items() if v["recall"] < 0.90]; high = [k for k, v in rec.items() if v["recall"] >= 0.90]
y = df.is_blast.values.astype(bool); p = df.probe_p_oof.values
blast_pat = set(df.loc[y, "patient_id"]); bg = np.flatnonzero(~y & ~df.patient_id.isin(blast_pat).values)
pools = {"low_recall": np.flatnonzero(y & df.patient_id.isin(low).values),
         "high_recall": np.flatnonzero(y & df.patient_id.isin(high).values)}
N, m, theta = 100, 5, 0.10; M = N * m
rng = np.random.default_rng(20260601); rows = []
for grp, pool in pools.items():
    for p_true in [0.10, 0.12, 0.15, 0.20]:
        nbl = int(round(p_true * M))
        for ep in range(300):
            cells = np.r_[rng.choice(pool, nbl, replace=nbl > len(pool)), rng.choice(bg, M - nbl, replace=False)]
            rng.shuffle(cells); cells = cells.reshape(N, m)
            f = y[cells].mean(1) / (2 * theta); S = p[cells].mean(1) / (2 * theta)
            r = dict(group=grp, p=p_true, ep=ep, reuse=nbl / len(pool))
            for nm, smp in [("uniform", sampler_uniform(N)), ("prop", sampler_prop_side(S))]:
                out = betting_test(f, S, smp, 0.5, rng=rng)
                r[f"{nm}_cert"] = int(out["stop"] > 0); r[f"{nm}_cells"] = (out["stop"] if out["stop"] > 0 else N) * m
            rows.append(r)
d = pd.DataFrame(rows); d.to_csv(f"{OUT}/episodes.csv", index=False)
print("low-recall patients:", low, "| blasts", len(pools["low_recall"]), "| high-recall blasts", len(pools["high_recall"]))
for (g, pt), x in d.groupby(["group", "p"]):
    side = "AT theta" if pt == theta else "above"
    if side == "AT theta":
        print(f"{g:11s} p={pt:.2f} {side}: false-certification uniform {x.uniform_cert.mean():.3f}  prop {x.prop_cert.mean():.3f}  (n={len(x)})")
    else:
        print(f"{g:11s} p={pt:.2f} {side}: median cells uniform {x.uniform_cells.median():.0f}  prop {x.prop_cells.median():.0f}  "
              f"saving {1 - x.prop_cells.median() / x.uniform_cells.median():.0%}")
