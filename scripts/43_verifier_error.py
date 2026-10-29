"""
Verifier error from real reader labels, BoneMarrowWSI (RQ4; exploratory).

The certificates assume perfect verification. Truth = consensus label; verifier = reader 0 (session 0).
(1) Reader agreement on blast status (reader 0 vs reader 1, and each vs consensus), patient bootstrap.
(2) Certification with reader 0 as verifier at theta_s = consensus p_s (every certification is false).
(3) Weakened guarantee: valid test of p_reader > theta + beta0 (1 - theta) when false blast calls are bounded
    by beta0 times the true non-blasts (and symmetrically); beta0, beta1 are plug-in estimates.
(4) Fixed margins d in {0.01, 0.02, 0.03, 0.05}.

Input: data/bmwsi_sessions/annotations_sessions.csv.gz (notebook 42), data/bmwsi/
Usage: python scripts/43_verifier_error.py --sessions data/bmwsi_sessions/annotations_sessions.csv.gz --out results/verifier_error_v1 --seed 20261034
"""
import argparse, json, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import binomtest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w

ap = argparse.ArgumentParser()
ap.add_argument("--sessions", default="data/bmwsi_sessions/annotations_sessions.csv.gz")
ap.add_argument("--out", default="results/verifier_error_v1"); ap.add_argument("--seed", type=int, default=20261034)
ap.add_argument("--n_orders", type=int, default=400); ap.add_argument("--grid", type=int, default=4)
ap.add_argument("--alpha", type=float, default=0.05); ap.add_argument("--n_boot", type=int, default=2000)
ap.add_argument("--margins", type=float, nargs="+", default=[0.01, 0.02, 0.03, 0.05])   # (4) fixed inflation of the threshold
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); t0 = time.time()
EXCL = {"artifact", "damaged_cell", "smudge_cell", "technically_unfit", "thrombocyte", "thrombocyte_aggregate", "giant_platelet",
        "erythrocyte", "spicule", "neutrophil_extracellular_trap", "no_consensus_found", "megakaryocyte", "micromegakaryocyte",
        "promegakaryocyte", "mitosis", "phagocytosis", "pseudo_gaucher_cell", "macrophage"}
BLAST = {"lymphocytic_blast", "myelocytic_blast", "monocytic_blast", "unknown_blast"}
KEY = ["patient", "x0", "y0", "x1", "y1"]

# ---- join reader labels to the consensus boxes of data/bmwsi/cells.csv (which carry the ROI index)
S = pd.read_csv(a.sessions, dtype={"patient": str})
for k in KEY[1:]: S[k] = S[k].round().astype(np.int64)
C = pd.read_csv("data/bmwsi/cells.csv", dtype={"patient": str}); C = C[C.source == "consensus"].copy()
for k in KEY[1:]: C[k] = C[k].round().astype(np.int64)
lab = S[S.source != "consensus"].pivot_table(index=KEY, columns="source", values="label", aggfunc="first")
C = C.merge(lab, left_on=KEY, right_index=True, how="left")
sess = sorted([c for c in lab.columns if c.startswith("session")])
match = {s: float(C[s].notna().mean()) for s in sess}
cons_check = S[S.source == "consensus"].merge(C[KEY + ["label"]], on=KEY, how="inner", suffixes=("", "_c"))
cons_agree = float((cons_check.label == cons_check.label_c).mean()) if len(cons_check) else float("nan")


def status(lbl):   # 1 blast, 0 counted non-blast, -1 not counted (excluded), nan unlabelled
    return np.where(pd.isna(lbl), np.nan, np.where(pd.Series(lbl).isin(BLAST), 1.0, np.where(pd.Series(lbl).isin(EXCL), -1.0, 0.0)))


C["s_cons"] = status(C.label.values)
for s in sess: C[f"s_{s}"] = status(C[s].values)

# ---- (1) agreement, ROI cells only (eligible regions as in scripts/23/30)
R = pd.read_csv("data/bmwsi/rois.csv", dtype={"patient": str})
R = R[((R.x1 - R.x0).between(1500, 3000)) & ((R.y1 - R.y0).between(1500, 3000))]
Q = C[(C.roi >= 0)].merge(R[["patient", "roi", "x0", "y0", "x1", "y1"]], on=["patient", "roi"], suffixes=("", "_r"))
rng = np.random.default_rng(a.seed)


def rates(D, x, y):
    """P(y calls blast | x calls non-blast), P(y calls non-blast | x calls blast), among boxes counted by both."""
    d = D[(D[x] >= 0) & (D[y] >= 0)]
    fp = ((d[x] == 0) & (d[y] == 1)).sum() / max((d[x] == 0).sum(), 1)
    fn = ((d[x] == 1) & (d[y] == 0)).sum() / max((d[x] == 1).sum(), 1)
    return float(fp), float(fn), int(len(d))


def boot_rates(D, x, y):
    pats = D.patient.unique(); idx = {p: g.index.values for p, g in D.groupby("patient")}; v = []
    for _ in range(a.n_boot):
        ii = np.concatenate([idx[p] for p in rng.choice(pats, len(pats))]); v.append(rates(D.loc[ii], x, y)[:2])
    v = np.array(v); return [[float(np.quantile(v[:, j], .025)), float(np.quantile(v[:, j], .975))] for j in range(2)]


agree = {}
pairs = [("s_session0", "s_session1"), ("s_cons", "s_session0"), ("s_cons", "s_session1")] if "session1" in sess else [("s_cons", "s_session0")]
for x, y in pairs:
    fp, fn, n = rates(Q, x, y); ci = boot_rates(Q, x, y)
    agree[f"{y}_given_{x}"] = dict(n_boxes=n, blast_call_given_nonblast=fp, ci95_fp=ci[0], nonblast_call_given_blast=fn, ci95_fn=ci[1])
tab = pd.crosstab(Q.s_session0, Q.s_session1) if "session1" in sess else None
# excluded-vs-counted disagreement (changes the denominator)
excl_flip = float(((Q.s_cons >= 0) & (Q.s_session0 == -1)).mean() + ((Q.s_cons == -1) & (Q.s_session0 >= 0)).mean())

# ---- (2)-(3) certification with reader 0 as verifier, real arrangements
G = a.grid
u_ = np.minimum(((Q.cx - Q.x0_r) / (Q.x1_r - Q.x0_r) * G).astype(int), G - 1)
v_ = np.minimum(((Q.cy - Q.y0_r) / (Q.y1_r - Q.y0_r) * G).astype(int), G - 1)
Q["field"] = Q.patient + "_" + Q.roi.astype(str) + "_" + (u_ * G + v_).astype(str)
b0 = agree["s_session1_given_s_session0"]["blast_call_given_nonblast"] if "s_session1_given_s_session0" in agree else 0.0
b1 = agree["s_session1_given_s_session0"]["nonblast_call_given_blast"] if "s_session1_given_s_session0" in agree else 0.0
rows = []
for pat, P in Q[Q.s_session0.notna()].groupby("patient"):
    t = P.s_cons.values; r0 = P.s_session0.values
    fields = pd.Categorical(P.field); fc = fields.codes; K = len(fields.categories)
    mT = np.bincount(fc, weights=(t >= 0), minlength=K); BT = np.bincount(fc, weights=(t == 1), minlength=K)
    mR = np.bincount(fc, weights=(r0 >= 0), minlength=K); BR = np.bincount(fc, weights=(r0 == 1), minlength=K)
    keep = mR > 0; mR, BR = mR[keep], BR[keep]
    if mT.sum() < 20 or BT.sum() == 0 or BT.sum() == mT.sum(): continue
    p_s = BT.sum() / mT.sum(); p_r = BR.sum() / mR.sum(); N = len(mR); w = mR / mR.sum(); f = BR / mR
    tests = dict(above=(f, p_s), below=(1 - f, 1 - p_s),
                 above_corr=(f, p_s + b0 * (1 - p_s)), below_corr=(1 - f, 1 - p_s * (1 - b1)),
                 **{f"above_m{d:.2f}": (f, p_s + d) for d in a.margins}, **{f"below_m{d:.2f}": (1 - f, 1 - (p_s - d)) for d in a.margins})
    for o in range(a.n_orders):
        r = dict(patient=pat, p_s=p_s, p_reader=p_r, M=int(mT.sum()), N=N, order=o)
        for nm, (ff, th) in tests.items():
            if not (0 < th < 1): r[f"{nm}_cert"] = 0; continue
            ob = betting_test_w(ff, w, np.zeros(N), sampler_uniform_w(N), th, alpha=a.alpha, rng=rng)
            r[f"{nm}_cert"] = int(ob["stop"] > 0)
        rows.append(r)
D = pd.DataFrame(rows); D.to_csv(f"{a.out}/episodes.csv", index=False)
TESTS = ["above", "below", "above_corr", "below_corr"] + [f"{d_}_m{m:.2f}" for d_ in ("above", "below") for m in a.margins]
per = D.groupby("patient").agg(p_s=("p_s", "first"), p_reader=("p_reader", "first"), M=("M", "first"),
                               **{f"{t}": (f"{t}_cert", "mean") for t in TESTS})
per.to_csv(f"{a.out}/patient_summary.csv")


def sig(v): return int(sum(binomtest(int(round(x * a.n_orders)), a.n_orders, a.alpha, alternative="greater").pvalue < 0.05 for x in v))


res = dict(sessions=sess, match_rate_to_consensus_boxes=match, consensus_label_agreement_check=cons_agree,
           agreement=agree, table_reader0_vs_reader1=None if tab is None else {f"{int(i)}|{int(j)}": int(tab.loc[i, j]) for i in tab.index for j in tab.columns},
           counted_vs_excluded_disagreement_reader0=excl_flip, beta0_plugin=b0, beta1_plugin=b1,
           n_specimens=int(len(per)), reader_minus_consensus_p=dict(mean=float((per.p_reader - per.p_s).mean()),
                                                                  max_abs=float((per.p_reader - per.p_s).abs().max())),
           false_cert={t: dict(pooled=float(D[f"{t}_cert"].mean()), max_specimen=float(per[t].max()),
                               n_specimens_gt_alpha=int((per[t] > a.alpha).sum()), n_specimens_signif=sig(per[t]))
                       for t in TESTS},
           abs_reader_bias_quantiles={q: float(np.quantile((per.p_reader - per.p_s).abs(), q)) for q in (0.5, 0.75, 0.9, 0.95, 1.0)},
           seconds=round(time.time() - t0))
json.dump(res, open(f"{a.out}/summary.json", "w"), indent=1); json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
print(json.dumps(res, indent=1))
