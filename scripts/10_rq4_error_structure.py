"""
RQ4 (exploratory): is classifier error exchangeable across cells, or does it depend on patient, morphology and
image quality? Prediction-powered and control-variate methods treat the classifier as side information whose
errors are not modelled per unit; certification remains VALID regardless (labels come from verification), but
efficiency and any method that trusts predictions for unverified units depend on this structure.

Operating point: threshold t giving 95% specificity on non-blast cells (cross-fitted scores).
  - blast recall by source patient (11 patients)      -> homogeneity test (chi-square, exact Monte-Carlo p)
  - false-positive rate by patient (blast-free)        -> dispersion vs a binomial with common rate
  - false-positive rate by true cell class (morphology)
  - within-class association of error with image-quality features (logistic regression, class fixed effects)
"""
import json, os
import numpy as np, pandas as pd
from scipy.stats import chi2_contingency
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

OUT = "results/rq4_error_v1"; os.makedirs(OUT, exist_ok=True)
df = pd.read_csv("results/embeddings_v1/cells_features.csv", dtype={"patient_id": str})
y = df.is_blast.values; p = df.probe_p_oof.values
t = float(np.quantile(p[y == 0], 0.95)); call = p >= t
res = {"threshold_95spec": t}
rng = np.random.default_rng(20260601)

# 1. blast recall heterogeneity across source patients
b = df[y == 1].assign(hit=call[y == 1])
tab = b.groupby("patient_id").hit.agg(["sum", "count"])
ct = np.c_[tab["sum"], tab["count"] - tab["sum"]]
chi2, p_asym, dof, _ = chi2_contingency(ct)
# Monte-Carlo p (small expected counts): permute hit labels across blasts
obs = chi2; sims = []
hits = b.hit.values.copy(); groups = b.patient_id.values
for _ in range(5000):
    rng.shuffle(hits)
    s = pd.DataFrame({"g": groups, "h": hits}).groupby("g").h.agg(["sum", "count"])
    sims.append(chi2_contingency(np.c_[s["sum"], s["count"] - s["sum"]])[0])
res["blast_recall_by_patient"] = {k: dict(blasts=int(v["count"]), recall=round(v["sum"] / v["count"], 3)) for k, v in tab.iterrows()}
res["blast_recall_homogeneity"] = dict(chi2=float(obs), dof=int(dof), p_asymptotic=float(p_asym),
                                       p_montecarlo=float((np.sum(np.array(sims) >= obs) + 1) / (len(sims) + 1)))

# 2. false-positive rate across blast-free patients: over-dispersion relative to a common binomial rate
blast_pat = set(df.loc[y == 1, "patient_id"])
nb = df[(y == 0) & ~df.patient_id.isin(blast_pat)].assign(fp=call[(y == 0) & ~df.patient_id.isin(blast_pat)])
fpp = nb.groupby("patient_id").fp.agg(["sum", "count"]); fpp = fpp[fpp["count"] >= 50]
rate = fpp["sum"].sum() / fpp["count"].sum()
var_obs = np.var(fpp["sum"] / fpp["count"], ddof=1); var_bin = np.mean(rate * (1 - rate) / fpp["count"])
res["fp_rate_by_patient"] = dict(patients=int(len(fpp)), pooled_rate=float(rate),
                                 min=float((fpp["sum"] / fpp["count"]).min()), max=float((fpp["sum"] / fpp["count"]).max()),
                                 q90=float((fpp["sum"] / fpp["count"]).quantile(0.9)),
                                 dispersion_ratio=float(var_obs / var_bin))

# 3. morphology: false-positive rate by true class
fpc = df[y == 0].assign(fp=call[y == 0]).groupby("cell_type").fp.agg(["mean", "count"]).sort_values("mean", ascending=False)
res["fp_rate_by_class"] = {k: dict(rate=round(v["mean"], 4), n=int(v["count"])) for k, v in fpc.iterrows()}

# 4. image quality WITHIN class: logistic regression error ~ standardized quality + class fixed effects
q = [c for c in df.columns if c.startswith("q_")]
err = (call != y.astype(bool)).astype(int)
X = pd.concat([pd.DataFrame(StandardScaler().fit_transform(df[q]), columns=q),
               pd.get_dummies(df.cell_type, prefix="cls", drop_first=True).astype(float)], axis=1)
lr = LogisticRegression(max_iter=4000, C=10.0).fit(X, err)
res["quality_within_class_logodds_per_sd"] = {c: round(float(w), 3) for c, w in zip(X.columns, lr.coef_[0]) if c.startswith("q_")}
# bootstrap by PATIENT for the quality coefficients
pats = df.patient_id.unique(); idx_by = {g: np.flatnonzero(df.patient_id.values == g) for g in pats}
boot = []
for _ in range(200):
    take = rng.choice(pats, len(pats), replace=True); ii = np.concatenate([idx_by[g] for g in take])
    try: boot.append(LogisticRegression(max_iter=4000, C=10.0).fit(X.iloc[ii], err[ii]).coef_[0][:len(q)])
    except Exception: pass
boot = np.array(boot)
res["quality_within_class_ci95"] = {c: [round(float(np.quantile(boot[:, i], .025)), 3), round(float(np.quantile(boot[:, i], .975)), 3)] for i, c in enumerate(q)}

json.dump(res, open(f"{OUT}/summary.json", "w"), indent=1)
print(json.dumps({k: res[k] for k in ["blast_recall_homogeneity", "fp_rate_by_patient", "quality_within_class_logodds_per_sd", "quality_within_class_ci95"]}, indent=1))
print("FP rate by class (top 6):", list(res["fp_rate_by_class"].items())[:6])
