"""
Unequal field sizes (P0-3 of the 2026-10-19 review; exploratory, not pre-registered). RQ3, RQ4.

Problem. The paper's field-level tests average the field blast PROPORTIONS. With unequal field sizes that is the mean of
field proportions, not the cell-level blast proportion p = sum_i B_i / sum_i m_i, unless proportion and size are
unrelated. Real fields are unequal (BoneMarrowWSI quadrants: see Part A).

Two valid constructions, both checked here:
  (a) COUNTS KNOWN for every field before verification (e.g. a nucleus detector counts cells; only classes need the
      expert): Shekhar et al. weighted test with pi_i = m_i / M and f_i = B_i / m_i (src/cadpbs/weighted_betting_w.py).
  (b) COUNTS UNKNOWN, only an a-priori bound m_max >= m_i: the bounded assorter
          a_i = (B_i - theta * m_i + theta * m_max) / m_max  in [0, 1],
          mean_i a_i > theta  <=>  sum_i B_i > theta * sum_i m_i  <=>  p > theta,
      with equal weights, so ALPHA and the equal-weight betting test apply unchanged.
Compared against the NAIVE analysis: ALPHA on the field proportions B_i / m_i with null 1/2 after scaling, i.e. treating
the mean of field proportions as p.

Part A (real data, BoneMarrowWSI): distribution of nucleated-cell counts per 2x2 quadrant of an ROI, and the bias of
the mean of quadrant proportions relative to the ROI proportion (patient bootstrap).
Part B (simulation): N fields with sizes drawn from the Part A quadrant distribution, filled with KU-Optofil cells
(blasts without replacement, cross-fitted scores). Blasts are placed with cell weight m_field^(-gamma):
gamma = 0 (independent of size), gamma = 1 (small fields enriched: about equal blasts per field regardless of size).

Usage: python scripts/27_unequal_fields.py --out results/unequal_fields_v1 --seed 20261020
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import beta as beta_dist, spearmanr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w, sampler_prop_w

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/unequal_fields_v1"); ap.add_argument("--seed", type=int, default=20261020)
ap.add_argument("--N", type=int, default=60); ap.add_argument("--n_null", type=int, default=1000)
ap.add_argument("--n_alt", type=int, default=300); ap.add_argument("--n_boot", type=int, default=2000)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
EXCL = {"artifact", "damaged_cell", "smudge_cell", "technically_unfit", "thrombocyte", "thrombocyte_aggregate", "giant_platelet",
        "erythrocyte", "spicule", "neutrophil_extracellular_trap", "no_consensus_found", "megakaryocyte", "micromegakaryocyte",
        "promegakaryocyte", "mitosis", "phagocytosis", "pseudo_gaucher_cell", "macrophage"}
BLASTL = {"lymphocytic_blast", "myelocytic_blast", "monocytic_blast", "unknown_blast"}

# ---------------- Part A: real field sizes and the naive-estimand bias -------------------------------------------
C = pd.read_csv("data/bmwsi/cells.csv"); R = pd.read_csv("data/bmwsi/rois.csv")
C = C[(C.source == "consensus") & (C.roi >= 0) & (~C.label.isin(EXCL))].copy(); C["y"] = C.label.isin(BLASTL).astype(int)
R = R[((R.x1 - R.x0).between(1500, 3000)) & ((R.y1 - R.y0).between(1500, 3000))]
C = C.merge(R[["patient", "roi", "x0", "y0", "x1", "y1"]], on=["patient", "roi"], suffixes=("", "_r"))
u = (C.cx - C.x0_r) / (C.x1_r - C.x0_r); v = (C.cy - C.y0_r) / (C.y1_r - C.y0_r)
C["q"] = np.minimum((u * 2).astype(int), 1) * 2 + np.minimum((v * 2).astype(int), 1)
g = C.groupby(["patient", "roi"]).y.agg(["size", "sum"]); keep = g[g["size"] >= 20].index
C = C.set_index(["patient", "roi"]).loc[keep].reset_index()
Q = C.groupby(["patient", "roi", "q"]).y.agg(m="size", B="sum").reset_index()
sizes = Q.m.values.astype(int)
m_max = int(sizes.max())
partA = dict(n_rois=int(len(keep)), n_patients=int(C.patient.nunique()), n_quadrants=int(len(Q)),
             quadrant_cells=dict(median=float(np.median(sizes)), q25=float(np.quantile(sizes, .25)),
                                 q75=float(np.quantile(sizes, .75)), mean=float(sizes.mean()), max=m_max,
                                 cv=float(sizes.std() / sizes.mean())))
# bias of mean-of-proportions vs ROI proportion, in ROIs with >= 1 blast and >= 1 non-blast
E = C.groupby(["patient", "roi"]).y.agg(["size", "sum"]); E = E[(E["sum"] >= 1) & (E["sum"] <= E["size"] - 1)].index
QE = Q.set_index(["patient", "roi"]).loc[E].reset_index()
QE["prop"] = QE.B / QE.m
per_roi = QE.groupby(["patient", "roi"]).apply(lambda d: pd.Series(dict(
    mean_prop=d.prop.mean(), p=d.B.sum() / d.m.sum(), n_q=len(d))), include_groups=False).reset_index()
per_roi["bias"] = per_roi.mean_prop - per_roi.p
QE["prop_c"] = QE.prop - QE.groupby(["patient", "roi"]).prop.transform("mean")
QE["m_c"] = QE.m - QE.groupby(["patient", "roi"]).m.transform("mean")
rho = spearmanr(QE.m_c, QE.prop_c).statistic
rng = np.random.default_rng(a.seed); pats = per_roi.patient.unique(); boots = []
by_p = {p: d.bias.values for p, d in per_roi.groupby("patient")}
for _ in range(a.n_boot):
    take = rng.choice(pats, len(pats), replace=True); boots.append(np.concatenate([by_p[p] for p in take]).mean())
partA.update(eligible_rois=int(len(per_roi)), bias_mean=float(per_roi.bias.mean()),
             bias_ci95=[float(x) for x in np.quantile(boots, [.025, .975])],
             abs_bias_median=float(per_roi.bias.abs().median()), abs_bias_q90=float(per_roi.bias.abs().quantile(.9)),
             spearman_size_vs_prop_within_roi=float(rho))
print("Part A:", json.dumps(partA, indent=1))

# ---------------- Part B: simulation --------------------------------------------------------------------------
df = pd.read_csv("results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})
yb = df.is_blast.values.astype(bool); score = df.probe_p_oof.values
BL = np.flatnonzero(yb); BG = np.flatnonzero(~yb & ~df.patient_id.isin(set(df.loc[yb, "patient_id"])).values)
GRID = {0.10: [0.10, 0.12, 0.15, 0.20], 0.20: [0.20, 0.22, 0.25]}
LEVELS = (0.05, 0.01)


def alpha_stops(X, u_, Mfields):
    lt, _ = alpha_martingale(X, u_, 0.5, Mfields, alpha=0.05)
    out = {}
    for lv in LEVELS:
        hit = np.flatnonzero(lt >= np.log(1 / lv)); out[lv] = int(hit[0]) + 1 if len(hit) else -1
    return out


rows = []
for th, plist in GRID.items():
    for gamma in (0.0, 1.0):
        for p_true in plist:
            n_ep = a.n_null if np.isclose(p_true, th) else a.n_alt
            r_ = np.random.default_rng([a.seed, int(th * 100), int(gamma * 10), int(round(p_true * 1000))]); t0 = time.time()
            for ep in range(n_ep):
                m = r_.choice(sizes[sizes > 0], a.N, replace=True); M = int(m.sum())
                B = int(np.floor(p_true * M)); assert B <= len(BL), "would need blast reuse"
                field_of = np.repeat(np.arange(a.N), m)
                w = m[field_of].astype(float) ** (-gamma); w /= w.sum()
                is_bl = np.zeros(M, bool); is_bl[r_.choice(M, B, replace=False, p=w)] = True
                cell_id = np.empty(M, int)
                cell_id[is_bl] = r_.choice(BL, B, replace=False); cell_id[~is_bl] = r_.choice(BG, M - B, replace=False)
                Bf = np.bincount(field_of, is_bl, a.N); Sf = np.bincount(field_of, score[cell_id], a.N) / m
                prop = Bf / m; pi = m / M
                order = r_.permutation(a.N); cum = np.cumsum(m[order])
                cells_at = lambda st, cm=cum: int(cm[st - 1]) if st > 0 else M
                rec = dict(theta=th, gamma=gamma, p=p_true, ep=ep, M=M, p_real=B / M, meanprop=float(prop.mean()))
                # naive: ALPHA on field proportions (estimand = mean of proportions)
                st = alpha_stops(prop / (2 * th), 1 / (2 * th), a.N)
                for lv in LEVELS: rec[f"naive_a{lv}"] = cells_at(st[lv]); rec[f"naive_cert_a{lv}"] = int(st[lv] > 0)
                # (b) bounded assorter, ALPHA, bound m_max (tight) and 2 m_max (loose)
                for nm, mm in (("assort", m_max), ("assort_loose", 2 * m_max)):
                    x = (Bf - th * m + th * mm) / mm
                    st = alpha_stops(x / (2 * th), 1 / (2 * th), a.N)
                    for lv in LEVELS: rec[f"{nm}_a{lv}"] = cells_at(st[lv]); rec[f"{nm}_cert_a{lv}"] = int(st[lv] > 0)
                # (b) bounded assorter, equal-weight betting
                x = (Bf - th * m + th * m_max) / m_max
                o = betting_test_w(x, np.full(a.N, 1 / a.N), Sf, sampler_uniform_w(a.N), th, rng=r_, alphas=LEVELS)
                for lv in LEVELS:
                    s_ = o["stops"][lv]; rec[f"assort_bet_a{lv}"] = int(np.cumsum(m[o["order"]])[s_ - 1]) if s_ > 0 else M
                    rec[f"assort_bet_cert_a{lv}"] = int(s_ > 0)
                # (a) counts known: weighted betting, uniform and classifier-driven sampling
                for nm, smp in (("known_unif", sampler_uniform_w(a.N)), ("known_prop", sampler_prop_w(pi, Sf))):
                    o = betting_test_w(prop, pi, Sf, smp, th, rng=r_, alphas=LEVELS)
                    for lv in LEVELS:
                        s_ = o["stops"][lv]; rec[f"{nm}_a{lv}"] = int(np.cumsum(m[o["order"]])[s_ - 1]) if s_ > 0 else M
                        rec[f"{nm}_cert_a{lv}"] = int(s_ > 0)
                rows.append(rec)
            print(f"theta={th} gamma={gamma} p={p_true}: {n_ep} episodes {time.time()-t0:.0f}s", flush=True)

D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)
DES = ["naive", "assort", "assort_loose", "assort_bet", "known_unif", "known_prop"]
def cp(x, n): return [0.0 if x == 0 else float(beta_dist.ppf(.025, x, n - x + 1)), 1.0 if x == n else float(beta_dist.ppf(.975, x + 1, n - x))]
summ = {"partA": partA, "m_max": m_max, "N_fields": a.N, "validity": {}, "cost": {}}
for (th, gamma), G in D.groupby(["theta", "gamma"]):
    key = f"theta{th:.2f}_gamma{gamma:.0f}"
    nul = G[np.isclose(G.p, th)]
    summ["validity"][key] = {d: {f"a{lv}": dict(rate=float(nul[f"{d}_cert_a{lv}"].mean()), cp95=cp(int(nul[f"{d}_cert_a{lv}"].sum()), len(nul)))
                                 for lv in LEVELS} for d in DES}
    summ["validity"][key]["mean_of_props_minus_p"] = float((nul.meanprop - nul.p_real).mean())
    summ["cost"][key] = {f"p{p:.2f}": {d: {f"a{lv}": float(H[f"{d}_a{lv}"].median()) for lv in LEVELS} | {"median_M": float(H.M.median())}
                                       for d in DES} for p, H in G[~np.isclose(G.p, th)].groupby("p")}
json.dump(summ, open(f"{a.out}/summary.json", "w"), indent=1)
json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print(json.dumps(summ["validity"], indent=1)[:3000])
