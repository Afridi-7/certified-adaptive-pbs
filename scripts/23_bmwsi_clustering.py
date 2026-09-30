"""
PREREG Amendment 4, A4.1 (RQ4): are blasts spatially clustered in real smears? BoneMarrowWSI-PediatricLeukemia,
classification subset: 228 ROIs (2048 x 2048 px at 40x, 0.11 um/px, i.e. about 225 um) in which EVERY haematological
structure was boxed and consensus-labelled (dataset documentation, "Documentation of conversion process", Zenodo 14933087).

Definitions (fixed in Amendment 4):
  blast      : lymphocytic_blast, myelocytic_blast, monocytic_blast, unknown_blast   (sensitivity: + promonocyte)
  nucleated  : all consensus labels except EXCL below                              (sensitivity: damaged/smudge kept as non-blast)
  eligible ROI: width and height in [1500, 3000] px, >= 20 nucleated cells, >= 1 blast and >= 1 non-blast
  scales     : S1 2x2 quadrants within ROI (primary, ~113 um); S2 4x4 within ROI (~56 um);
               S3 the ROIs of a patient (between-ROI, patients with >= 2 eligible ROIs)
Estimator: pooled one-way ANOVA intraclass correlation for binary data (Fleiss), subfields nested in their parent unit
(ROI for S1/S2, patient for S3), deviations taken from the parent-unit mean so larger-scale differences are removed.
Test (H7a, S1): one-sided permutation test, blast labels permuted among the cells of each ROI (positions fixed), 2,000 permutations.
CI: bootstrap over patients, 2,000 resamples. Descriptive: k = 5 nearest-neighbour blast enrichment within ROI vs permutation (same number of permutations).
"""
import argparse, json
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("--cells", default="data/bmwsi/cells.csv"); ap.add_argument("--rois", default="data/bmwsi/rois.csv")
ap.add_argument("--out", default="results/a4_bmwsi/clustering.json"); ap.add_argument("--seed", type=int, default=20260701)
ap.add_argument("--n_perm", type=int, default=2000); ap.add_argument("--n_boot", type=int, default=2000)
a = ap.parse_args()
import os; os.makedirs(os.path.dirname(a.out), exist_ok=True)
EXCL = {"artifact", "damaged_cell", "smudge_cell", "technically_unfit", "thrombocyte", "thrombocyte_aggregate", "giant_platelet",
        "erythrocyte", "spicule", "neutrophil_extracellular_trap", "no_consensus_found", "megakaryocyte", "micromegakaryocyte",
        "promegakaryocyte", "mitosis", "phagocytosis", "pseudo_gaucher_cell", "macrophage"}
BLAST = {"lymphocytic_blast", "myelocytic_blast", "monocytic_blast", "unknown_blast"}
rng = np.random.default_rng(a.seed)


def prepare(extra_blast=(), keep_damaged=False):
    C = pd.read_csv(a.cells); R = pd.read_csv(a.rois)
    C = C[(C.source == "consensus") & (C.roi >= 0)].copy()
    excl = EXCL - ({"damaged_cell", "smudge_cell"} if keep_damaged else set())
    C = C[~C.label.isin(excl)].copy()
    C["y"] = C.label.isin(BLAST | set(extra_blast)).astype(int)
    R = R[((R.x1 - R.x0).between(1500, 3000)) & ((R.y1 - R.y0).between(1500, 3000))]
    C = C.merge(R[["patient", "roi", "x0", "y0", "x1", "y1"]], on=["patient", "roi"], suffixes=("", "_r"))
    C["u"] = (C.cx - C.x0_r) / (C.x1_r - C.x0_r); C["v"] = (C.cy - C.y0_r) / (C.y1_r - C.y0_r)
    g = C.groupby(["patient", "roi"]).y.agg(["size", "sum"])
    ok = g[(g["size"] >= 20) & (g["sum"] >= 1) & (g["sum"] <= g["size"] - 1)].index
    C = C.set_index(["patient", "roi"]).loc[ok].reset_index()
    C["rid"] = C.patient + "_" + C.roi.astype(str)
    return C


def icc(parent, sub, y):
    """Pooled ANOVA ICC for binary y; sub = subfield id nested in parent. Vectorised."""
    _, P = np.unique(parent, return_inverse=True); _, S = np.unique(sub, return_inverse=True)
    y = np.asarray(y, float)
    ns = np.bincount(S); ys = np.bincount(S, y); sp = np.zeros(S.max() + 1, int); sp[S] = P   # parent of each subfield
    npar = np.bincount(P); ypar = np.bincount(P, y); kpar = np.bincount(sp, minlength=len(npar))
    use_p = kpar >= 2; use_s = use_p[sp]
    ps = ys / ns; pbar = (ypar / npar)[sp]
    ssb = (ns * (ps - pbar) ** 2)[use_s].sum(); dfb = (kpar[use_p] - 1).sum()
    ssw = (ns * ps * (1 - ps))[use_s].sum(); dfw = (npar[use_p] - kpar[use_p]).sum()
    sumsq = np.bincount(sp, ns.astype(float) ** 2, minlength=len(npar))
    n0 = (npar[use_p] - sumsq[use_p] / npar[use_p]).sum() / dfb
    msb, msw = ssb / dfb, ssw / dfw
    return float((msb - msw) / (msb + (n0 - 1) * msw)), dict(n_parents=int(use_p.sum()), n0=float(n0))


def scale_ids(C, scale):
    if scale == "S1": return C.rid.values, (C.rid + "_" + (np.minimum((C.u * 2).astype(int), 1) * 2 + np.minimum((C.v * 2).astype(int), 1)).astype(str)).values
    if scale == "S2": return C.rid.values, (C.rid + "_" + (np.minimum((C.u * 4).astype(int), 3) * 4 + np.minimum((C.v * 4).astype(int), 3)).astype(str)).values
    if scale == "S3": return C.patient.values, C.rid.values


def permute_within(groups, y):
    y = y.copy(); out = np.empty_like(y)
    for idx in groups: out[idx] = rng.permutation(y[idx])
    return out


def analyse(C, perm=True, boot=True):
    res = {}
    for sc in ["S1", "S2", "S3"]:
        par, sub = scale_ids(C, sc); y = C.y.values
        est, info = icc(par, sub, y); r = dict(icc=est, **info)
        if perm:
            groups = [np.flatnonzero(par == p) for p in np.unique(par)]
            null = np.array([icc(par, sub, permute_within(groups, y))[0] for _ in range(a.n_perm)])
            r.update(perm_p=float((1 + (null >= est).sum()) / (1 + a.n_perm)), null_mean=float(null.mean()),
                     null_q95=float(np.quantile(null, 0.95)))
        if boot:
            _, pc = np.unique(par, return_inverse=True); _, scd = np.unique(sub, return_inverse=True)
            pats, pidx = np.unique(C.patient.values, return_inverse=True)
            idx = [np.flatnonzero(pidx == j) for j in range(len(pats))]; bs = []
            big = scd.max() + 1
            for _ in range(a.n_boot):
                pick = rng.integers(0, len(pats), len(pats))
                ii = np.concatenate([idx[p] for p in pick])
                tag = np.concatenate([np.full(len(idx[p]), j) for j, p in enumerate(pick)])
                # duplicated patients count as distinct parents/subfields
                bs.append(icc(pc[ii] + tag * big, scd[ii] + tag * big, y[ii])[0])
            r["ci95"] = [float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))]
        res[sc] = r
    return res


C = prepare()
out = dict(n_patients=int(C.patient.nunique()), n_rois=int(C.rid.nunique()), n_cells=int(len(C)), n_blasts=int(C.y.sum()))
out["primary"] = analyse(C)
out["H7a_clustering_S1"] = dict(icc=out["primary"]["S1"]["icc"], perm_p=out["primary"]["S1"]["perm_p"],
                                supported=bool(out["primary"]["S1"]["perm_p"] < 0.05))
# k-NN enrichment (descriptive): mean fraction of each blast's 5 nearest same-ROI neighbours that are blasts, vs permutation
C = C.reset_index(drop=True); groups = [np.flatnonzero(C.rid.values == r) for r in C.rid.unique()]
NB = []
for g in groups:
    if len(g) < 6: continue
    xy = C.loc[g, ["cx", "cy"]].values; _, nn = cKDTree(xy).query(xy, k=6); NB.append((g, g[nn[:, 1:]]))
def knn_stat(y):
    v = [y[nb[y[g] == 1]].mean(axis=1) for g, nb in NB if (y[g] == 1).any()]
    return float(np.concatenate(v).mean())
obs = knn_stat(C.y.values); null = np.array([knn_stat(permute_within(groups, C.y.values)) for _ in range(a.n_perm)])
out["knn5_enrichment"] = dict(observed=obs, null_mean=float(null.mean()), perm_p=float((1 + (null >= obs).sum()) / (1 + a.n_perm)))
# sensitivity analyses (point estimates + permutation, no bootstrap)
out["sens_promonocyte_as_blast"] = analyse(prepare(extra_blast=("promonocyte",)), boot=False)
out["sens_damaged_smudge_kept"] = analyse(prepare(keep_damaged=True), boot=False)
json.dump(out, open(a.out, "w"), indent=1)
print(json.dumps({k: out[k] for k in ["n_patients", "n_rois", "n_cells", "n_blasts", "H7a_clustering_S1"]}, indent=1))
