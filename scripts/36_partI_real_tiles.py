"""
Part I on REAL image tiles with real geometry (P1-2; exploratory). RQ1, RQ2.

Input  : data/bmwsi_tiles/tiles.csv + tiles_embeddings_dinov2_vits14.npy, produced by notebooks/35_bmwsi_tiles_colab.ipynb
         (4 x 4 tiles of every fully annotated BoneMarrowWSI ROI, embedded with the same frozen DINOv2 ViT-S/14).
         data/bmwsi/cells.csv for the consensus labels.
Unit   : a tile (about 56 um). Critical tile = holds at least one blast (definitions as in scripts/23_bmwsi_clustering.py).
Specimen: one patient = all tiles of their eligible ROIs (N = 16 per ROI). Real arrangement, real dilution (each tile image
         contains its nucleated cells, erythrocytes and background), no compositing.
Teacher: logistic regression on tile embeddings, label = critical, cross-fitted over 5 patient-grouped folds.
Rules at matched cost (InfoDPP-PAC's stop k* in each run): InfoDPP-PAC (tau calibrated to 16.3% of N on blast-free specimens
         of other folds), uniform, raster (row by row from the top-left of a random ROI, then the next ROI), ranked by the
         teacher, and k-center coverage stopping (own cost). The uniform miss is also compared with h(N, k, n) (Theorem 1).

Usage: python scripts/36_partI_real_tiles.py --tiles data/bmwsi_tiles --out results/partI_real_tiles_v1 --seed 20261028
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
from math import comb
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.infodpp import InfoDPPConfig, run_infodpp, calibrate_tau
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score

ap = argparse.ArgumentParser()
ap.add_argument("--tiles", default="data/bmwsi_tiles"); ap.add_argument("--cells", default="data/bmwsi/cells.csv")
ap.add_argument("--out", default="results/partI_real_tiles_v1"); ap.add_argument("--seed", type=int, default=20261028)
ap.add_argument("--runs", type=int, default=200); ap.add_argument("--target_frac", type=float, default=0.163)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
EXCL = {"artifact", "damaged_cell", "smudge_cell", "technically_unfit", "thrombocyte", "thrombocyte_aggregate", "giant_platelet",
        "erythrocyte", "spicule", "neutrophil_extracellular_trap", "no_consensus_found", "megakaryocyte", "micromegakaryocyte",
        "promegakaryocyte", "mitosis", "phagocytosis", "pseudo_gaucher_cell", "macrophage"}
BLAST = {"lymphocytic_blast", "myelocytic_blast", "monocytic_blast", "unknown_blast"}
T = pd.read_csv(f"{a.tiles}/tiles.csv", dtype={"patient": str})
E = np.load(f"{a.tiles}/tiles_embeddings_dinov2_vits14.npy").astype(np.float32); E /= np.linalg.norm(E, axis=1, keepdims=True)
RO = pd.read_csv(os.path.join(os.path.dirname(a.cells), "rois.csv"), dtype={"patient": str})   # guard: tiles must lie in the local ROIs
chk = T.merge(RO, on=["patient", "roi"], suffixes=("", "_r"), how="left")
bad = chk.x0_r.isna() | (chk.x0 < chk.x0_r) | (chk.y0 < chk.y0_r) | (chk.x1 > chk.x1_r + 1) | (chk.y1 > chk.y1_r + 1)
assert not bad.any(), f"{int(bad.sum())} tiles fall outside data/bmwsi/rois.csv: the notebook's tables differ from the repo's"
C = pd.read_csv(a.cells, dtype={"patient": str}); C = C[(C.source == "consensus") & (C.roi >= 0) & (~C.label.isin(EXCL))]
C = C.assign(y=C.label.isin(BLAST).astype(int))
nb = np.zeros(len(T), int); nn = np.zeros(len(T), int)
for (pat, roi), idx in T.groupby(["patient", "roi"]).groups.items():
    cc = C[(C.patient == pat) & (C.roi == roi)]
    for i in idx:
        t = T.loc[i]; inside = (cc.cx >= t.x0) & (cc.cx < t.x1) & (cc.cy >= t.y0) & (cc.cy < t.y1)
        nb[i] = int(cc.y[inside].sum()); nn[i] = int(inside.sum())
T["n_blast"] = nb; T["n_nuc"] = nn; T["crit"] = (nb > 0).astype(int)
pats = np.array(sorted(T.patient.unique())); rng = np.random.default_rng(a.seed)
fold_of = dict(zip(rng.permutation(pats), np.arange(len(pats)) % 5)); T["fold"] = T.patient.map(fold_of)
score = np.zeros(len(T))
for f in range(5):
    tr, te = (T.fold != f).values, (T.fold == f).values
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=3000, class_weight="balanced")).fit(E[tr], T.crit.values[tr])
    score[te] = clf.predict_proba(E[te])[:, 1]
T["teacher"] = score
T.to_csv(f"{a.out}/tiles_scored.csv", index=False)                  # tile labels, folds and cross-fitted teacher scores
spec = {p: g.index.values for p, g in T.groupby("patient")}


def h(N, k, n):
    return comb(N - k, n) / comb(N, n) if n <= N - k else 0.0


def raster_order(idx, r):
    g = T.loc[idx]; rois = list(g.roi.unique()); r.shuffle(rois)
    out = []
    for ro in rois:
        s = g[g.roi == ro]; key = (s.tile % 4) * 4 + s.tile // 4        # tile = i*4 + j (i column, j row) -> row-major key
        out += list(s.index.values[np.argsort(key.values)])
    return np.array([np.flatnonzero(idx == i)[0] for i in out])


def kcenter(F, r_star, r):
    D = np.sqrt(np.maximum(2 - 2 * F @ F.T, 0)); n = len(F); o = [int(r.integers(n))]; md = D[o[0]].copy(); md[o] = -1
    while len(o) < n and md.max() > r_star:
        j = int(np.argmax(md)); o.append(j); md = np.minimum(md, D[j]); md[o] = -1
    return np.array(o)


base = InfoDPPConfig(tau=None, seed=a.seed, m0_min=2); rows = []
for f in range(5):
    neg = [spec[p] for p in pats if fold_of[p] != f and T.crit.values[spec[p]].sum() == 0]
    cal = [(E[i], T.teacher.values[i]) for i in neg]
    if not cal: cal = [(E[spec[p]], T.teacher.values[spec[p]]) for p in pats if fold_of[p] != f]
    tau = calibrate_tau(cal, base, a.target_frac, rng=rng); cfg = InfoDPPConfig(**{**base.__dict__, "tau": tau})
    lo, hi = 0.0, 2.0
    for _ in range(25):
        mid = (lo + hi) / 2
        if np.median([len(kcenter(F, mid, np.random.default_rng(i))) for i, (F, _) in enumerate(cal)]) > a.target_frac * np.median([len(F) for F, _ in cal]): lo = mid
        else: hi = mid
    r_star = (lo + hi) / 2
    for p in [p for p in pats if fold_of[p] == f]:
        idx = spec[p]; F = E[idx]; crit = T.crit.values[idx]; ci = np.flatnonzero(crit); N = len(idx); k = len(ci)
        if k == 0 or k == N: continue
        for run in range(a.runs):
            o = run_infodpp(F, T.teacher.values[idx], cfg, np.random.default_rng(rng.integers(1 << 30))); n = o["k_star"]
            miss = lambda order, nn_: int(not np.isin(ci, order[:nn_]).any())
            okc = kcenter(F, r_star, np.random.default_rng(rng.integers(1 << 30)))
            rows.append(dict(patient=p, fold=f, N=N, k=k, p_s=T.n_blast.values[idx].sum() / max(T.n_nuc.values[idx].sum(), 1), run=run,
                             n=n, infodpp=miss(o["order"], n), uniform=miss(rng.permutation(N), n), raster=miss(raster_order(idx, rng), n),
                             ranked=miss(np.argsort(-T.teacher.values[idx], kind="stable"), n), kcenter=miss(okc, len(okc)),
                             kcenter_n=len(okc), h_uniform=h(N, k, n)))
D = pd.DataFrame(rows); D.to_csv(f"{a.out}/runs.csv", index=False)
rare = D[D.k / D.N <= 0.25]
summ = dict(n_tiles=int(len(T)), n_patients=int(len(pats)), teacher_auc=float(roc_auc_score(T.crit, T.teacher)),
            tile_cells_median=float(T.n_nuc.median()), crit_tile_fraction=float(T.crit.mean()),
            all_positive=D.groupby("patient")[["infodpp", "uniform", "raster", "ranked", "kcenter", "h_uniform", "n", "kcenter_n"]].mean().mean().round(4).to_dict(),
            rare_k_over_N_le_0_25=dict(n_patients=int(rare.patient.nunique()),
                                      **rare.groupby("patient")[["infodpp", "uniform", "raster", "ranked", "kcenter", "h_uniform", "n"]].mean().mean().round(4).to_dict()))
json.dump(summ, open(f"{a.out}/summary.json", "w"), indent=1); json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
