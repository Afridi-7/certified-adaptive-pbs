"""
Observed validity on REAL cell arrangements (P1-3; exploratory, not pre-registered). RQ3, RQ4.

Amendment 4 measured blast clustering in real smears and then used a Beta-binomial MODEL to infer what it does to the
cell-level mis-analysis. Here the model is replaced by the real specimens themselves.

Specimen  = one BoneMarrowWSI patient: all nucleated cells, with their real positions and consensus labels, in the
            patient's eligible fully annotated ROIs (definitions as in scripts/23_bmwsi_clustering.py).
Field     = a 4 x 4 subregion of an ROI (about 56 um; median about 4-5 cells). Fields have UNEQUAL sizes; empty
            subregions are not fields.
Threshold = the specimen's own blast proportion, theta_s = p_s. Every specimen is then an exact null of H0: p <= theta,
            so ANY certification is false. Validity of a design-based test does not depend on what theta means
            clinically, so this tests validity on real arrangements. (theta_s is not a diagnostic threshold.)
Repeats   = n_orders independent random sampling orders per specimen: the sampling design is the only randomness.

Designs (each stops at certification or at census):
  field_known    Shekhar et al. weighted betting, weights m_i / M, uniform field sampling (Prop. 2a). Valid.
  field_assort   ALPHA on the bounded assorter with m_max = largest 4x4 count in the data (Prop. 2b). Valid.
  cell_srs       ALPHA on individually sampled cells (simple random sampling of cells). Valid, but not what a
                 microscope does; included as the reference.
  naive_cell     uniform FIELD sampling, cells of sampled fields analysed as if drawn individually (H5b mis-analysis).
  naive_fmean    ALPHA on field proportions (mean of proportions; Section 5 / Prop. 2 warning).
Also: cost (fraction of cells examined) at theta_s = p_s - 0.05 for specimens with p_s >= 0.10.

Usage: python scripts/30_real_arrangement_validity.py --out results/real_arrangements_v1 --seed 20261023
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import beta as beta_dist
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/real_arrangements_v1"); ap.add_argument("--seed", type=int, default=20261023)
ap.add_argument("--n_orders", type=int, default=400); ap.add_argument("--grid", type=int, default=4)
ap.add_argument("--alpha", type=float, default=0.05)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
EXCL = {"artifact", "damaged_cell", "smudge_cell", "technically_unfit", "thrombocyte", "thrombocyte_aggregate", "giant_platelet",
        "erythrocyte", "spicule", "neutrophil_extracellular_trap", "no_consensus_found", "megakaryocyte", "micromegakaryocyte",
        "promegakaryocyte", "mitosis", "phagocytosis", "pseudo_gaucher_cell", "macrophage"}
BLAST = {"lymphocytic_blast", "myelocytic_blast", "monocytic_blast", "unknown_blast"}

C = pd.read_csv("data/bmwsi/cells.csv"); R = pd.read_csv("data/bmwsi/rois.csv")
C = C[(C.source == "consensus") & (C.roi >= 0) & (~C.label.isin(EXCL))].copy(); C["y"] = C.label.isin(BLAST).astype(int)
R = R[((R.x1 - R.x0).between(1500, 3000)) & ((R.y1 - R.y0).between(1500, 3000))]
C = C.merge(R[["patient", "roi", "x0", "y0", "x1", "y1"]], on=["patient", "roi"], suffixes=("", "_r"))
g = C.groupby(["patient", "roi"]).y.agg(["size", "sum"])
ok = g[(g["size"] >= 20) & (g["sum"] >= 1) & (g["sum"] <= g["size"] - 1)].index
C = C.set_index(["patient", "roi"]).loc[ok].reset_index()
G = a.grid
u_ = np.minimum(((C.cx - C.x0_r) / (C.x1_r - C.x0_r) * G).astype(int), G - 1)
v_ = np.minimum(((C.cy - C.y0_r) / (C.y1_r - C.y0_r) * G).astype(int), G - 1)
C["field"] = C.patient + "_" + C.roi.astype(str) + "_" + (u_ * G + v_).astype(str)
m_max = int(C.groupby("field").size().max())


def alpha_stop(X, u, n):
    _, st = alpha_martingale(X, u, 0.5, n, alpha=a.alpha); return st


rng = np.random.default_rng(a.seed); rows = []; t0 = time.time()
for pat, P in C.groupby("patient"):
    F = P.groupby("field").y.agg(m="size", B="sum").reset_index()
    m = F.m.values; B = F.B.values; N = len(F); M = int(m.sum()); p_s = B.sum() / M
    cell_y = P.y.values; field_of_cell = pd.Categorical(P.field, categories=F.field).codes
    for scen, th in (("null", p_s), ("alt", p_s - 0.05)):
        if scen == "alt" and p_s < 0.10: continue
        if not (0 < th < 1): continue
        u = 1 / (2 * th)
        for o in range(a.n_orders):
            order = rng.permutation(N); cum = np.cumsum(m[order])
            r = dict(patient=pat, scen=scen, theta=th, p_s=p_s, N=N, M=M, order=o)
            # valid (a): known counts, weighted betting, uniform field sampling
            ob = betting_test_w(B / m, m / M, np.zeros(N), sampler_uniform_w(N), th, alpha=a.alpha, rng=rng)
            st = ob["stop"]; r["field_known_cert"] = int(st > 0)
            r["field_known_frac"] = (np.cumsum(m[ob["order"]])[st - 1] if st > 0 else M) / M
            # valid (b): bounded assorter with global m_max, ALPHA, uniform field order
            x = (B - th * m + th * m_max) / m_max
            st = alpha_stop(x[order] / (2 * th), u, N); r["field_assort_cert"] = int(st > 0)
            r["field_assort_frac"] = (cum[st - 1] if st > 0 else M) / M
            # reference: simple random sampling of individual cells
            cy = cell_y[rng.permutation(M)]
            st = alpha_stop(cy / (2 * th), u, M); r["cell_srs_cert"] = int(st > 0); r["cell_srs_frac"] = (st if st > 0 else M) / M
            # mis-analysis 1: cells of uniformly sampled fields treated as individual draws
            flat = np.concatenate([cell_y[field_of_cell == i] for i in order])
            st = alpha_stop(flat / (2 * th), u, M); r["naive_cell_cert"] = int(st > 0)
            # mis-analysis 2: mean of field proportions
            st = alpha_stop((B / m)[order] / (2 * th), u, N); r["naive_fmean_cert"] = int(st > 0)
            rows.append(r)
D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)


def cp(x, n):
    return [0.0 if x == 0 else float(beta_dist.ppf(.025, x, n - x + 1)), 1.0 if x == n else float(beta_dist.ppf(.975, x + 1, n - x))]


DES = ["field_known", "field_assort", "cell_srs", "naive_cell", "naive_fmean"]
nul = D[D.scen == "null"]
per_pat = nul.groupby("patient")[[f"{d}_cert" for d in DES]].mean()
summ = dict(n_patients=int(nul.patient.nunique()), n_orders=a.n_orders, grid=G, m_max=m_max,
            field_cells=dict(median=float(C.groupby("field").size().median()), max=m_max),
            cells_per_specimen=dict(median=float(nul.groupby("patient").M.first().median())),
            fields_per_specimen=dict(median=float(nul.groupby("patient").N.first().median())),
            p_s=dict(n_below_0_2=int((nul.groupby("patient").p_s.first() < 0.2).sum()),
                     n_0_2_to_0_5=int(nul.groupby("patient").p_s.first().between(0.2, 0.5).sum()),
                     n_above_0_5=int((nul.groupby("patient").p_s.first() > 0.5).sum())),
            validity={})
for d in DES:
    x = int(nul[f"{d}_cert"].sum()); n = len(nul)
    pr = per_pat[f"{d}_cert"]
    # per-patient one-sided exact test of rate <= alpha
    from scipy.stats import binomtest
    sig = sum(binomtest(int(round(v * a.n_orders)), a.n_orders, a.alpha, alternative="greater").pvalue < 0.05 for v in pr)
    summ["validity"][d] = dict(pooled_rate=x / n, cp95=cp(x, n), max_patient_rate=float(pr.max()),
                               median_patient_rate=float(pr.median()), n_patients_rate_gt_alpha=int((pr > a.alpha).sum()),
                               n_patients_signif_gt_alpha=int(sig))
alt = D[D.scen == "alt"]
if len(alt):
    summ["cost_at_ps_minus_0.05"] = dict(n_patients=int(alt.patient.nunique()),
                                         **{d: dict(median_fraction=float(alt[f"{d}_frac"].median()),
                                                    cert_rate=float(alt[f"{d}_cert"].mean()))
                                            for d in ["field_known", "field_assort", "cell_srs"]})
json.dump(summ, open(f"{a.out}/summary.json", "w"), indent=1); json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print(json.dumps(summ, indent=1)); print(f"{time.time()-t0:.0f}s")
