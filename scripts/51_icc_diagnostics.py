"""
Why the permutation null of the pooled ANOVA ICC is not centred at 0, and clustering in low-blast ROIs (RQ4; exploratory).

1. Under the within-ROI null, E0[MSB] and E0[MSW] weight the ROI variances S2_r by (k_r - 1) and (n_r - k_r);
   if S2_r falls with n_r the null mean is positive.
2. Null-centred estimator: MSW* = sum (k_r - 1) MSW_r / sum (k_r - 1).
3. The same for ROIs with p_r <= 0.30, with the excess over the permutation null.
Definitions as in 23_bmwsi_clustering.py.
"""
import json, os, sys
import numpy as np, pandas as pd
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
EXCL = {"artifact", "damaged_cell", "smudge_cell", "technically_unfit", "thrombocyte", "thrombocyte_aggregate", "giant_platelet",
        "erythrocyte", "spicule", "neutrophil_extracellular_trap", "no_consensus_found", "megakaryocyte", "micromegakaryocyte",
        "promegakaryocyte", "mitosis", "phagocytosis", "pseudo_gaucher_cell", "macrophage"}
BLAST = {"lymphocytic_blast", "myelocytic_blast", "monocytic_blast", "unknown_blast"}
rng = np.random.default_rng(20261102); NPERM = 1000; NBOOT = 2000


def prepare():
    C = pd.read_csv(f"{ROOT}/data/bmwsi/cells.csv"); R = pd.read_csv(f"{ROOT}/data/bmwsi/rois.csv")
    C = C[(C.source == "consensus") & (C.roi >= 0)].copy(); C = C[~C.label.isin(EXCL)].copy()
    C["y"] = C.label.isin(BLAST).astype(int)
    R = R[((R.x1 - R.x0).between(1500, 3000)) & ((R.y1 - R.y0).between(1500, 3000))]
    C = C.merge(R[["patient", "roi", "x0", "y0", "x1", "y1"]], on=["patient", "roi"], suffixes=("", "_r"))
    C["u"] = (C.cx - C.x0_r) / (C.x1_r - C.x0_r); C["v"] = (C.cy - C.y0_r) / (C.y1_r - C.y0_r)
    g = C.groupby(["patient", "roi"]).y.agg(["size", "sum"])
    ok = g[(g["size"] >= 20) & (g["sum"] >= 1) & (g["sum"] <= g["size"] - 1)].index
    C = C.set_index(["patient", "roi"]).loc[ok].reset_index(); C["rid"] = C.patient + "_" + C.roi.astype(str)
    return C


def sub_ids(C, k):
    return (C.rid + "_" + (np.minimum((C.u * k).astype(int), k - 1) * k + np.minimum((C.v * k).astype(int), k - 1)).astype(str)).values


def roi_stats(par, sub, y):
    """Per-ROI SSB, SSW, k (occupied subfields), n, sum of squared subfield sizes, S2."""
    df = pd.DataFrame(dict(par=par, sub=sub, y=y))
    s = df.groupby(["par", "sub"]).y.agg(["size", "sum"]).reset_index()
    s["pm"] = s["sum"] / s["size"]
    r = df.groupby("par").y.agg(["size", "mean"]); s = s.join(r, on="par", rsuffix="_r")
    s["ssb"] = s["size"] * (s.pm - s["mean"]) ** 2; s["ssw"] = s["size"] * s.pm * (1 - s.pm); s["n2"] = s["size"] ** 2
    a = s.groupby("par").agg(ssb=("ssb", "sum"), ssw=("ssw", "sum"), k=("sub", "size"), n=("size", "sum"), n2=("n2", "sum"), p=("mean", "first"))
    a = a[a.k >= 2]; a["S2"] = a.p * (1 - a.p) * a.n / (a.n - 1)
    return a


def estimators(a):
    dfb = (a.k - 1).sum(); dfw = (a.n - a.k).sum(); n0 = (a.n - a.n2 / a.n).sum() / dfb
    msb = a.ssb.sum() / dfb; msw = a.ssw.sum() / dfw
    msw_star = ((a.k - 1) * a.ssw / (a.n - a.k)).sum() / dfb
    pooled = (msb - msw) / (msb + (n0 - 1) * msw); star = (msb - msw_star) / (msb + (n0 - 1) * msw_star)
    e_msb0 = ((a.k - 1) * a.S2).sum() / dfb; e_msw0 = ((a.n - a.k) * a.S2).sum() / dfw
    return dict(pooled=float(pooled), null_centred=float(star), n0=float(n0), analytic_null_E_MSB=float(e_msb0),
                analytic_null_E_MSW=float(e_msw0), analytic_null_icc=float((e_msb0 - e_msw0) / (e_msb0 + (n0 - 1) * e_msw0)))


def permute(C):
    y = C.y.values.copy()
    for idx in C.groupby("rid").indices.values(): y[idx] = rng.permutation(y[idx])
    return y


def analyse(C, k, boot=True):
    sub = sub_ids(C, k); a = roi_stats(C.rid.values, sub, C.y.values); est = estimators(a)
    nul = [estimators(roi_stats(C.rid.values, sub, permute(C))) for _ in range(NPERM)]
    est["perm_null_mean_pooled"] = float(np.mean([z["pooled"] for z in nul]))
    est["perm_null_mean_null_centred"] = float(np.mean([z["null_centred"] for z in nul]))
    est["excess_pooled"] = est["pooled"] - est["perm_null_mean_pooled"]
    est["corr_S2_n_spearman"] = float(pd.Series(a.S2.values).corr(pd.Series(a.n.values), method="spearman"))
    est["n_rois"] = int(len(a)); est["n_patients"] = int(C.patient.nunique()); est["n_cells"] = int(len(C))
    if boot:
        pats = C.patient.unique(); groups = {p: g for p, g in C.groupby("patient")}; bs = []
        for _ in range(NBOOT):
            pick = rng.choice(pats, len(pats))
            B = pd.concat([groups[p].assign(rid=groups[p].rid + f"#{j}") for j, p in enumerate(pick)])
            bs.append(estimators(roi_stats(B.rid.values, sub_ids(B, k), B.y.values))["null_centred"])
        est["null_centred_ci95"] = [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]
    return est


C = prepare(); out = dict(note="exploratory; scripts/51_icc_diagnostics.py; seed 20261102")
pr = C.groupby("rid").y.mean()
for name, CC in (("all_rois", C), ("low_blast_rois_p_le_0.30", C[C.rid.isin(pr[pr <= 0.30].index)])):
    out[name] = {}
    for k, sc in ((2, "S1"), (4, "S2")):
        out[name][sc] = analyse(CC, k); print(name, sc, json.dumps(out[name][sc]), flush=True)
os.makedirs(f"{ROOT}/results/icc_diagnostics_v1", exist_ok=True)
json.dump(out, open(f"{ROOT}/results/icc_diagnostics_v1/summary.json", "w"), indent=1)
