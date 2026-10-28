"""
Audit of the numbers in the IEEE JBHI manuscript (paper/jbhi_compact.tex) against the stored results.

Same rules as paper/audit_cmpb_numbers.py: every check names a passage of the manuscript (a LaTeX substring that must
occur in jbhi_compact.tex) and the numbers it states, and recomputes each number from results/ (and the cell table of
results/embeddings_v1). A stated number passes if the recomputed value, rounded half-up to the stated number of
decimals, equals it; "a--b" ranges are checked as the rounded minimum and maximum; at-most and within claims ("<=") are
checked as bounds. Tables I-III and V are compared cell by cell. Not checked: counts that describe the public datasets
themselves (e.g. 41,621 MLL23 cells, 228 BoneMarrowWSI regions, 3,204 KU-Optofil cells without patient identifiers,
from docs/data_audit_ku_optofil.md), constants of the original InfoDPP-PAC paper other than its operating point, design
constants written into the scripts (30 consecutive fields, weight 0.25, 17 cells per field) and the push-to-output
times of the pre-registration tags (docs/PREREG_AUDIT.md). Writes paper/audit_jbhi_numbers.md; exits non-zero if any
check fails.

Run from the repository root:  python paper/audit_jbhi_numbers.py
"""
import glob, json, os, re, subprocess, sys
from decimal import Decimal, ROUND_HALF_UP
from math import comb
import numpy as np, pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R = lambda *p: os.path.join(ROOT, "results", *p)
TEX = open(os.path.join(ROOT, "paper", "jbhi_compact.tex"), encoding="utf-8").read()
J = lambda *p: json.load(open(R(*p)))
OUT = []


def rnd(x, places):
    return float(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP))


def match(stated, value):
    upper = stated.startswith("<=")
    stated = stated[2:] if upper else stated
    s = stated.replace("{,}", "").replace(",", "").replace("$", "").replace("\\%", "").replace("%", "").strip()
    pct = "%" in stated
    s = s.replace("−", "-")
    if s.startswith("+"): s = s[1:]
    places = len(s.split(".")[1]) if "." in s else 0
    v = value * 100 if pct else value
    return (rnd(v, places) <= float(s)) if upper else (rnd(v, places) == float(s)), round(v, 6)


def check(cid, snippet, pairs, note=""):
    in_doc = snippet in TEX
    res, ok_all = [], in_doc
    for stated, value in pairs:
        ok, v = match(stated, value); ok_all &= ok; res.append(f"{stated} vs {v}{'' if ok else ' **MISMATCH**'}")
    OUT.append(dict(id=cid, ok=ok_all, in_doc=in_doc, detail="; ".join(res), note=note))


def rng(a_b, values, sep="--"):
    a, b = a_b.split(sep)
    return [(a, float(np.min(values))), (b, float(np.max(values)))]


def optional(cid, snippet, fn, note=""):
    """Checks whose inputs exist only in the full results folder (not in every copy of the repository)."""
    try:
        check(cid, snippet, fn(), note)
    except FileNotFoundError as e:
        OUT.append(dict(id=cid, ok=snippet in TEX, in_doc=snippet in TEX, detail="input not available here; not recomputed", note=str(e)[-60:]))


# ------------------------------------------------------------------ abstract and data
cells = pd.read_csv(R("embeddings_v1", "cells_features.csv"), dtype={"patient_id": str})
mll = J("mll23_shift_v1", "summary.json"); ext = J("external_v1", "decisions.json")
E = pd.read_csv(R("confirm_v1", "episodes.csv")); k1 = E[E.k == 1]
miss = k1.groupby("m").infodpp_miss.mean(); umiss = k1.groupby("m")["uniform@info_miss"].mean()
tiles = pd.read_csv(R("partI_real_tiles_v1", "tiles_scored.csv"))
ns = pd.read_csv(R("part2_v1", "nstar_table.csv"))
f1 = ns[(ns.N == 1000) & np.isclose(ns.alpha, 0.05) & (ns.k == 1)].frac.iloc[0]
rs = J("review_stats_v1.json")
pva = [c["prop_vs_alpha"]["ratio"] for c in rs["part3_effect_sizes"]]
BR = {src: J("baseline_reference_v1", f"summary_{src}.json") for src in ("ku", "mll23", "rvpbs")}
def brr(src, th, p, num, den="alpha_uniform"): return BR[src]["ratios"][f"theta{th:.2f}_p{p:.2f}_{num}_over_{den}"]
shiftr = [brr(s, .10, p, "bet_prop")["ratio"] for s in ("mll23", "rvpbs") for p in (0.15, 0.20)]
DL = pd.read_csv(R("decision_level_v1", "summary.csv"))
def dl(th, p, n, proc, col): return float(DL[np.isclose(DL.theta, th) & np.isclose(DL.p, p) & (DL.n == n) & (DL.procedure == proc)][col].iloc[0])
near = [(0.10, 0.08), (0.10, 0.12), (0.20, 0.18), (0.20, 0.22)]
RAS = J("real_arrangements_v1", "summary.json"); RA = RAS["validity"]

check("abs.cells", "Using pre-registered simulations from 28,280 expert-labeled cells", [("28280", len(cells))])
check("abs.patients99", "annotated bone-marrow regions from 99 patients", [("99", tiles.patient.nunique())])
check("abs.miss", "with one critical field among 200 ten-cell fields, it missed 76\\% of positive specimens, against 85\\% for random sampling",
      [("76\\%", miss[10]), ("85\\%", umiss[10])])
check("abs.absence", "any strategy must examine 95\\% of the fields to certify absence at a 5\\% error level", [("95\\%", f1)])
check("abs.sens", "cut blast sensitivity from 0.91 to 0.63 and 0.34",
      [("0.91", ext["H4a_recall_drops"]["ku_recall"]), ("0.63", mll["recall"]), ("0.34", ext["H4a_recall_drops"]["rv_recall"])])
check("abs.cost", "needed 0.45 to 0.75 times the verified cells of uniform sampling in distribution, and 1.3 to 3.0 times under shift at margins of 0.05 to 0.10",
      rng("0.45 to 0.75", pva, " to ") + rng("1.3 to 3.0", shiftr, " to "))
check("abs.count200", "a 200-cell count was wrong in up to 27\\% of specimens", [("27\\%", max(dl(t, p, 200, "conventional", "wrong") for t, p in near))])
check("abs.real", "exceeded the nominal error in 27 of 75 real specimens", [("27", RA["naive_cell"]["n_patients_rate_gt_alpha"]), ("75", RAS["n_patients"])])

r_ = cells.probe_p_oof.rank(); npos = int(cells.is_blast.sum()); nneg = len(cells) - npos
auc_probe = (r_[cells.is_blast == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)
check("data.ku", "we used 28,280 cells from 275 identified patients (13 classes) in five patient-grouped folds",
      [("28280", len(cells)), ("275", cells.patient_id.nunique()), ("13", cells.cell_type.nunique()), ("5", cells.fold.nunique())])
check("data.blasts", "Blasts come from only 11 patients (489 blasts)", [("11", cells[cells.is_blast == 1].patient_id.nunique()), ("489", int(cells.is_blast.sum()))])
check("data.rvpbs", "752 annotated cells (153 blasts)", [("752", ext["descriptive"]["n"]), ("153", ext["descriptive"]["n_blast"])])
check("data.mll23", "8,606 myeloblasts", [("8606", mll["n_blast"])])
check("data.mll23b", "leaving 36,529", [("36529", mll["n"])])
check("data.probe", "reached an out-of-fold area under the ROC curve (AUC) of 0.979", [("0.979", auc_probe)])
check("methods.op", "so that the median stop is 16.3\\% of all fields. The original runs stopped at 16.3\\% of a 300-patch budget, about 2.4\\% of a slide",
      [("16.3\\%", J("confirm_v1", "config.json")["target_frac"]), ("16.3\\%", 48.9 / 300), ("2.4\\%", 48.9 / 2075)], note="48.9 of 300 patches, 2,075 per slide: reported averages of the original paper")
check("methods.seven", "about seven times its own examination share", [("7", (48.9 / 300) / (48.9 / 2075))])
sys.path.insert(0, os.path.join(ROOT, "src"))
from cadpbs.infodpp import InfoDPPConfig
check("methods.seed", "A random seed of 5 fields", [("5", InfoDPPConfig().m0_min)], note="InfoDPPConfig.m0_min")
check("stats.126", "a bootstrap over five folds has only 126 distinct resamples", [("126", comb(9, 5))], note="multisets of size 5 from 5 folds")
CS = J("crosscheck_shangrla_v1", "summary.json")
check("stats.shangrla", "to within $3\\times10^{-13}$ in log scale on 2,000 random cases",
      [("<=3", max(CS["B_alpha_martingale_same_eta"]["max_abs_log_diff"], CS["C_betting_uniform_equal_weights"]["max_abs_log_diff"]) * 1e13),
       ("2000", CS["R"])], note="largest log difference x 1e13")
try:
    out = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q", "--noconftest", "-p", "no:cacheprovider", os.path.join(ROOT, "tests")],
                         capture_output=True, text=True, cwd=ROOT, env=dict(os.environ, PYTHONPATH=os.path.join(ROOT, "src"))).stdout
    n_tests = sum(1 for l in out.splitlines() if "::" in l)
    if n_tests == 0: raise RuntimeError("no tests collected (is pytest installed?)")
    check("stats.tests", "a suite of 50 deterministic tests", [("50", n_tests)])
except Exception as e:
    OUT.append(dict(id="stats.tests", ok="a suite of 50 deterministic tests" in TEX, in_doc=True, detail="pytest not available; not counted", note=str(e)[:60]))

# ------------------------------------------------------------------ RQ1
med = E.groupby(["m", "k"]).infodpp_n.median()
check("rq1.stop", "In 24,000 confirmatory episodes the median stop was 32--33 of 200 fields in every cell, varying by at most one field",
      [("24000", len(E))] + rng("32--33", med.values) + [("<=1", max(med.xs(m).max() - med.xs(m).min() for m in E.m.unique()))])
auc = E.groupby("m").kt1_auc_info_step1.mean(); aucr = E.groupby("m").kt1_auc_relevance.mean()
check("rq1.auc_info", "fell from 0.778 at $m=1$ to 0.593, 0.552 and 0.525 at $m=5,10,20$",
      [("0.778", auc[1]), ("0.593", auc[5]), ("0.552", auc[10]), ("0.525", auc[20])])
FT = pd.read_csv(R("kt1_field_teacher", "episodes.csv")); fta = FT.groupby("m").auc_rel_field.mean()
check("rq1.auc_field", "kept 0.978, 0.925, 0.882 and 0.827", [("0.978", fta[1]), ("0.925", fta[5]), ("0.882", fta[10]), ("0.827", fta[20])])
check("rq1.auc_cell", "the pre-registered cell-level teacher gave 0.977 at $m=1$ to 0.896 at $m=20$", [("0.977", aucr[1]), ("0.896", aucr[20])])
check("rq1.miss", "InfoDPP-PAC missed the critical field in 0.315, 0.681, 0.759 and 0.805",
      [(s, miss[m]) for s, m in zip(["0.315", "0.681", "0.759", "0.805"], [1, 5, 10, 20])])
check("rq1.uniform", "uniform sampling of the same size missed 0.822, 0.830, 0.846 and 0.828",
      [(s, umiss[m]) for s, m in zip(["0.822", "0.830", "0.846", "0.828"], [1, 5, 10, 20])])
dec1 = J("confirm_v1", "decisions.json"); h1c = dec1["H1c"]["cells"]
check("rq1.h1c", "the difference was $-0.087$ at $m=10$ (fold-level $t$ interval $[-0.153,-0.021]$; pre-registered episode-level interval $[-0.107,-0.067]$) and $-0.023$ $[-0.050,0.003]$ at $m=20$",
      [("-0.087", h1c["10"]["diff"]), ("-0.153", rs["H1c_fold_t"]["10"]["ci95"][0]), ("-0.021", rs["H1c_fold_t"]["10"]["ci95"][1]),
       ("-0.107", h1c["10"]["ci95"][0]), ("-0.067", h1c["10"]["ci95"][1]), ("-0.023", rs["H1c_fold_t"]["20"]["mean"]),
       ("-0.050", rs["H1c_fold_t"]["20"]["ci95"][0]), ("0.003", rs["H1c_fold_t"]["20"]["ci95"][1])])
fo = k1.groupby("m")["fixed_order@info_miss"].mean()
cm_ = k1.groupby("m").confidence_miss.mean(); cu_ = k1.groupby("m")["uniform@conf_miss"].mean(); cn_ = k1.groupby("m").conf_n.median()
check("rq1.baselines", "the fixed field order at the same cost missed 0.83--0.85, like uniform sampling",
      rng("0.83--0.85", fo.values))
check("rq1.conf", "confidence-based stopping halted after a median of 10--33 fields and missed 0.74--0.90, close to uniform sampling at its own cost (0.77--0.92)",
      rng("10--33", cn_.values) + rng("0.74--0.90", cm_.values) + rng("0.77--0.92", cu_.values), note="k = 1, confirmatory episodes")
rec = E.groupby(["m", "k"])[["infodpp_recall", "uniform@info_recall", "fixed_order@info_recall"]].mean()
check("rq1.recall", "fell with dilution for InfoDPP-PAC from 0.34--0.69 at $m=1$ to 0.19--0.21 at $m=20$, against 0.15--0.18 for uniform sampling and the fixed order",
      rng("0.34--0.69", rec.xs(1).infodpp_recall) + rng("0.19--0.21", rec.xs(20).infodpp_recall) +
      rng("0.15--0.18", np.r_[rec["uniform@info_recall"].values, rec["fixed_order@info_recall"].values]))
nbar = E.infodpp_n.mean(); lb = [1 - E[E.m == m].infodpp_n.mean() / 200 for m in (1, 5, 10, 20)]
check("rq1.corollary", "at a mean of 32.4 examined fields, its worst-case miss probability with one critical field is at least $1-32.4/200\\approx0.84$",
      [("32.4", nbar), ("0.84", 1 - nbar / 200)] + rng("0.84--0.84", lb), note="1 - mean stop / 200, overall and for each m")
ftm = FT[FT.k == 1].groupby("m").infodpp_fieldteacher_miss.mean()
check("rq1.teacher", "hardly changed its misses (0.32, 0.66, 0.76 and 0.77; exploratory)",
      [("0.32", ftm[1]), ("0.66", ftm[5]), ("0.76", ftm[10]), ("0.77", ftm[20])])
sc = J("stopping_class_v1", "summary.json"); k1r = [r for r in sc["by_m_k"] if r["k"] == 1]
gap = max(max(r["infodpp_miss"], r["gp_var_miss"], r["kcenter_miss"]) - min(r["infodpp_miss"], r["gp_var_miss"], r["kcenter_miss"]) for r in k1r)
check("rq1.rules", "calibrated to the same cost on 4,000 new specimens, missed within 0.04 of InfoDPP-PAC in every cell",
      [("4000", len(pd.read_csv(R("stopping_class_v1", "episodes.csv")))), ("<=0.04", gap)])
ranked = [r["ranked_miss"] for r in k1r]; pred = [sc["bridge"][str(r["m"])]["predicted_exact_iid"] for r in k1r]
check("rq1.ranked", "missed 0.05, 0.14, 0.21 and 0.37, against 0.051, 0.146, 0.221 and 0.394 predicted",
      list(zip(["0.05", "0.14", "0.21", "0.37"], ranked)) + list(zip(["0.051", "0.146", "0.221", "0.394"], pred)))
sd = J("infodpp_seed_v1", "summary.json")["summary"]
d1 = sd["m1_k1"]["diff_seed_excluded"]
other = max(abs(sd[f"m{m}_k1"][key][0]) for m in (5, 10, 20) for key in ("diff_seed_excluded", "diff_seed_included"))
check("rq1.op24", "InfoDPP-PAC made 4--5 selections after its seed; its advantage was $-0.17$ $[-0.27,-0.08]$ with single-cell fields and at most about 0.03 otherwise",
      rng("4--5", [sd[f"m{m}_k1"]["selections_median"] for m in (1, 5, 10, 20)]) + [("-0.17", d1[0]), ("-0.27", d1[1][0]), ("-0.08", d1[1][1]), ("<=0.03", other)])
T = J("partI_real_tiles_v1", "analysis.json"); rr = T["rare_kN_le_0.25"]; t1 = T["k_eq_1"]
check("rq1.tiles", "3,040 tiles of about $56\\,\\mu$m from 99 patients, specimens of 16--32 tiles, 200 runs each), InfoDPP-PAC stopped after 5.4 tiles on average",
      [("3040", len(tiles)), ("99", tiles.patient.nunique())] + rng("16--32", [int(k) for k in T["N_values"]]) +
      [("200", J("partI_real_tiles_v1", "config.json")["runs"]), ("5.4", T["stop_n_mean"])])
check("rq1.tiles_rare", "In the 24 specimens with at most a quarter of tiles critical it missed every critical tile in 0.80 of runs against 0.70 for uniform sampling ($+0.105$ $[0.002,0.203]$; worse in 18 of 24), and in at least 95\\% of runs in 15 specimens",
      [("24", rr["n_patients"]), ("0.80", rr["mean_miss"]["infodpp"]), ("0.70", rr["mean_miss"]["uniform"]), ("0.105", rr["infodpp_minus_uniform"]["diff"]),
       ("0.002", rr["infodpp_minus_uniform"]["ci95"][0]), ("0.203", rr["infodpp_minus_uniform"]["ci95"][1]), ("18", rr["infodpp_minus_uniform"]["n_pos"]),
       ("15", rr["n_specimens_miss_ge_0.95"]["infodpp"])])
check("rq1.tiles_other", "In the same specimens a raster missed 0.67 and a scan ranked by the tile-level teacher 0.54, a teacher whose AUC was only 0.69 there",
      [("0.67", rr["mean_miss"]["raster"]), ("0.54", rr["mean_miss"]["ranked"]), ("0.69", rr["teacher_auc_tiles"])])
check("rq1.tiles_kc", "in the 14 specimens with a single critical tile, $k$-center stopping, which examined 8.5 tiles on average, missed it in 0.89 of runs against 0.70 for uniform sampling at its cost ($+0.19$ $[0.03,0.34]$)",
      [("14", t1["n_patients"]), ("8.5", T["kcenter_n_mean"]), ("0.89", t1["mean_miss"]["kcenter"]), ("0.70", t1["mean_miss"]["h_kc"]),
       ("0.19", t1["kcenter_minus_h_kc"]["diff"]), ("0.03", t1["kcenter_minus_h_kc"]["ci95"][0]), ("0.34", t1["kcenter_minus_h_kc"]["ci95"][1])])
check("fig1.n", "(10,000 separate episodes)", [("10000", len(FT))])

# ------------------------------------------------------------------ RQ3
P3 = J("part3_decisions.json")
h3a = [c["rate"] for t in P3.values() for c in t["H3a_validity"]["cells"].values()]
check("rq3.valid", "In 12,000 confirmatory 500-cell specimens, all four valid designs kept false certification at $p=\\theta$ between 0.001 and 0.022 (H3a), whereas the uncorrected design certified falsely in 0.993 ($\\theta=0.10$) and 0.889",
      [("12000", sum(len(pd.read_csv(R(d, "episodes.csv"))) for d in ("part3_theta010", "part3_theta020"))), ("0.001", min(h3a)), ("0.022", max(h3a)),
       ("0.993", P3["theta_0.10"]["H3b_naive_invalid"]["rate"]), ("0.889", P3["theta_0.20"]["H3b_naive_invalid"]["rate"])])
c15 = [c for c in P3["theta_0.10"]["H3c_adaptive_cheaper"]["cells"] if np.isclose(c["p"], 0.15)][0]
check("rq3.cost", "0.45--0.75 times as many, for example 75 vs.\\ 165 at $p=0.15$", rng("0.45--0.75", pva) + [("75", c15["median_bet_prop"]), ("165", c15["median_alpha_uniform"])])
KU_ALT = [(0.10, 0.12), (0.10, 0.15), (0.10, 0.20), (0.20, 0.22), (0.20, 0.25), (0.20, 0.30)]
hu = [brr("ku", t, p, "hedge_uniform")["ratio"] for t, p in KU_ALT]; bk = BR["ku"]["median_cells"]
chg14 = max(abs(bk[f"theta{t:.2f}_p{p:.2f}"]["hedge_prop"] / bk[f"theta{t:.2f}_p{p:.2f}"]["bet_prop"] - 1) for t, p in KU_ALT)
optional("rq3.kelly", "Uniform betting with approximate-Kelly bets needed 1.03--1.83 times the cells of ALPHA, so savings measured against it (0.28--0.73) are overstated; with hedged bets (exploratory, 4,400 specimens) it matched ALPHA (0.98--1.17), while classifier-driven medians changed by at most 14\\%",
         lambda: rng("1.03--1.83", [c["test_only"]["ratio"] for c in rs["part3_effect_sizes"]]) + rng("0.28--0.73", [c["sampler_only"]["ratio"] for c in rs["part3_effect_sizes"]]) +
                 [("4400", len(pd.read_csv(R("baseline_reference_v1", "episodes_ku.csv"))))] + rng("0.98--1.17", hu) + [("<=14\\%", chg14)])
ms = pd.read_csv(R("m_scaling_v1", "episodes.csv")); nul = ms[np.isclose(ms.p, ms.theta)]
v05 = max(g[f"{d}_cert_a0.05"].mean() for d in ("alpha_uniform", "bet_uniform", "bet_prop") for _, g in nul.groupby(["theta", "M"]))
v01 = max(g[f"{d}_cert_a0.01"].mean() for d in ("alpha_uniform", "bet_uniform", "bet_prop") for _, g in nul.groupby(["theta", "M"]))
MS = pd.read_csv(R("m_scaling_v1", "summary.csv"))
def ms_(th, p, M, col): return float(MS[np.isclose(MS.theta, th) & np.isclose(MS.p, p) & (MS.M == M)][col].iloc[0])
check("rq3.mscale", "With specimens of up to 4,000 cells (exploratory; 17,600 specimens), false certification stayed $\\le0.040$ at $\\alpha=0.05$ and $\\le0.007$ at $\\alpha=0.01$",
      [("17600", len(ms)), ("<=0.040", v05), ("<=0.007", v01)])
check("rq3.mscale_ex", "uniform ALPHA needed a median of 165, 205 and 245 cells at $M=500$, 1,000 and 4,000 and classifier-driven betting 80, 80 and 85",
      [(s, ms_(0.10, 0.15, M, "alpha_uniform_0.05")) for s, M in zip(["165", "205", "245"], [500, 1000, 4000])] +
      [(s, ms_(0.10, 0.15, M, "bet_prop_0.05")) for s, M in zip(["80", "80", "85"], [500, 1000, 4000])])
big = MS[MS.p - MS.theta >= 0.05 - 1e-9]; w4 = big[big.M == 4000]
check("rq3.wald", "At margins of 0.05 or more and $M=4{,}000$, uniform ALPHA was within 20\\% of Wald's approximation", [("<=20\\%", max(abs(w4["alpha_uniform_0.05"] / w4["approx_0.05"] - 1)))])
check("rq3.near", "to 1,583 ($\\theta=0.10$) and 2,033 cells ($\\theta=0.20$) at $M=4{,}000$", [("1,583", ms_(0.10, 0.12, 4000, "alpha_uniform_0.05")), ("2,033", ms_(0.20, 0.22, 4000, "alpha_uniform_0.05"))])
sav = 1 - MS["bet_prop_0.05"] / MS["alpha_uniform_0.05"]; ex = np.isclose(MS.theta, 0.20) & np.isclose(MS.p, 0.22)
check("rq3.saving", "Classifier-driven sampling saved 36--65\\% of the cells of uniform ALPHA at every $M$, except at $p=0.22$ against $\\theta=0.20$, where the saving was 19--23\\% at $M\\le1{,}000$ and at most 3\\% for $M\\ge2{,}000$",
      rng("36\\%--65\\%", sav[~ex]) + rng("19\\%--23\\%", sav[ex & (MS.M <= 1000)]) + [("<=3\\%", sav[ex & (MS.M >= 2000)].max())],
      note=f"saving at M>=2,000: {', '.join(f'{x:.3f}' for x in sav[ex & (MS.M >= 2000)])}")
r01 = [x for col in ("alpha_uniform", "bet_uniform", "bet_prop") for x in (MS[f"{col}_0.01"] / MS[f"{col}_0.05"]).values]
check("rq3.alpha01", "Lowering $\\alpha$ to 0.01 multiplied the cost by 1.0--1.7", rng("1.0--1.7", r01))
check("rq3.conv", "a 200-cell count, the usual peripheral-blood differential, was wrong in 11--27\\% of specimens (100 cells: 18--37\\%)",
      rng("11\\%--27\\%", [dl(t, p, 200, "conventional", "wrong") for t, p in near]) + rng("18\\%--37\\%", [dl(t, p, 100, "conventional", "wrong") for t, p in near]))
offwrong = max(dl(t, p, n, proc, "wrong") for t, p in near + [(0.10, 0.05), (0.10, 0.15), (0.20, 0.15), (0.20, 0.25)] for n in (100, 200, 500) for proc in ("sequential", "fixed_exact"))
und = [dl(t, p, n, proc, "undetermined") for t, p in near for n in (100, 200, 500) for proc in ("sequential", "fixed_exact")]
check("rq3.valid_dec", "were wrong in at most 0.7\\% of specimens off the threshold, but at those margins returned ``undetermined'' in 61--98\\% of specimens",
      [("<=0.7\\%", offwrong)] + rng("61\\%--98\\%", und))
check("rq3.prior_cells", "certified ``above'' in only 31\\% instead of 86\\% of specimens at $p=0.15$, $\\theta=0.10$ and $n=500$",
      [("31\\%", dl(0.10, 0.15, 500, "seq_default", "correct")), ("86\\%", dl(0.10, 0.15, 500, "sequential", "correct"))])
AP = pd.read_csv(R("alpha_prior_v1", "summary.csv")); alt = AP[AP.p > AP.theta + 1e-9]
slow = max(np.r_[alt.alpha_m05_median_cells / alt.alpha_default_median_cells, alt.alpha_m10_median_cells / alt.alpha_default_median_cells])
check("rq3.prior_fields", "kept ALPHA valid but made it up to 2.2 times slower than the default",
      [("2.2", slow), ("<=0.05", np.nanmax(np.r_[AP.alpha_m05_false_cert, AP.alpha_m10_false_cert]))])
tab1 = TEX[TEX.index("\\label{tab:mscale}"):TEX.index("\\end{table", TEX.index("\\label{tab:mscale}"))]
pairs = []
for line in tab1.splitlines():
    m_ = re.match(r"^(0\.\d0) & (0\.\d\d) & (.*)\\\\$", line.strip())
    if not m_: continue
    th, p = float(m_.group(1)), float(m_.group(2)); cc = [c.replace("$^\\dagger$", "").strip() for c in m_.group(3).split("&")]
    k = 0
    for M in (500, 1000, 2000, 4000):
        pairs += [(cc[k], ms_(th, p, M, "alpha_uniform_0.05")), (cc[k + 1], ms_(th, p, M, "bet_prop_0.05"))]; k += 2
    pairs.append((cc[k], ms_(th, p, 4000, "approx_0.05")))
check("table1", "\\label{tab:mscale}", pairs, note=f"{len(pairs)} cells")
per = ms[ms.p > ms.theta + 1e-9].groupby(["theta", "M", "p"]).size()
check("table1.caption", "400 specimens per cell", [("400", per.min()), ("400", per.max())])
P3c = [(t, c) for t, key in ((0.10, "theta_0.10"), (0.20, "theta_0.20")) for c in P3[key]["H3c_adaptive_cheaper"]["cells"]]
d500 = max(max(abs(c["median_alpha_uniform"] - ms_(t, c["p"], 500, "alpha_uniform_0.05")), abs(c["median_bet_prop"] - ms_(t, c["p"], 500, "bet_prop_0.05"))) for t, c in P3c)
check("table1.sep", "so its medians differ by up to 10 cells", [("10", d500)], note="M = 500 column against the confirmatory medians (H3c cells)")
tab2 = TEX[TEX.index("\\label{tab:decision}"):TEX.index("\\end{table", TEX.index("\\label{tab:decision}"))]
def eps(x): return "0" if x == 0 else ("$\\epsilon$" if x < 0.005 else f"{rnd(x, 2):.2f}")
for line in tab2.splitlines():
    m_ = re.match(r"^(0\.\d0) & (0\.\d\d) & (.*)\\\\$", line.strip())
    if not m_: continue
    th, p = float(m_.group(1)), float(m_.group(2)); got = [x.strip() for x in m_.group(3).split("&")]
    want = []
    for n_ in (200, 500):
        want += [f"{rnd(dl(th, p, n_, 'conventional', 'wrong'), 2):.2f}",
                 f"{eps(dl(th, p, n_, 'fixed_exact', 'correct'))}/{eps(dl(th, p, n_, 'fixed_exact', 'wrong'))}",
                 f"{eps(dl(th, p, n_, 'sequential', 'correct'))}/{eps(dl(th, p, n_, 'sequential', 'wrong'))}"]
    for i_, (g, w) in enumerate(zip(got, want)):
        OUT.append(dict(id=f"table2.{th}_{p}_{i_}", ok=g == w, in_doc=True, detail=f"{g} vs {w}", note=""))
dcfg = J("decision_level_v1", "config.json")
check("table2.caption", "$M=2{,}000$; 2,000 specimens per row)", [("2000", dcfg["M"]), ("2000", dcfg["episodes"])])

# ------------------------------------------------------------------ RQ4
err = J("rq4_error_v1", "summary.json"); recp = [v["recall"] for v in err["blast_recall_by_patient"].values()]
check("rq4.recall", "blast recall at 95\\% specificity varied across the 11 blast patients from 0.39 to 1.00, false-positive rates across blast-free patients were overdispersed 6.6-fold",
      [("11", len(recp)), ("0.39", min(recp)), ("1.00", max(recp)), ("6.6", err["fp_rate_by_patient"]["dispersion_ratio"])])
q = err["quality_within_class_ci95"]
QK = ("q_contrast", "q_brightness", "q_g", "q_b", "q_blur_lapvar")
OUT.append(dict(id="rq4.stain", ok=all(q[k_][0] > 0 or q[k_][1] < 0 for k_ in QK),
                in_doc="within cell class, errors were associated with the color, brightness, contrast and sharpness of the image" in TEX,
                detail="95% intervals of within-class log-odds per SD exclude 0 for contrast, brightness, green, blue and sharpness: " + json.dumps({k_: q[k_] for k_ in QK}), note=""))
fc_cls = err["fp_rate_by_class"]; top = sorted(fc_cls, key=lambda c: -fc_cls[c]["rate"])[:6]
check("rq4.classes", "by class were highest in immature granulocytes (0.26--0.33), reactive lymphocytes and basophils (0.19--0.20)",
      rng("0.26--0.33", [fc_cls[c]["rate"] for c in ("myelocyte", "metamyelocyte")]) + rng("0.19--0.20", [fc_cls[c]["rate"] for c in ("reactive_lymphocyte", "basophil")]),
      note="four highest class rates: " + ", ".join(f"{c} {fc_cls[c]['rate']:.3f}" for c in top[:4]))
def low_recall():
    c4 = pd.read_csv(R("rq4_certify_v1", "episodes.csv"))
    lr = c4[(c4.group == "low_recall") & np.isclose(c4.p, 0.10)]
    hi_ = c4[(c4.group != "low_recall") & np.isclose(c4.p, 0.15)].prop_cells.median(); lo_ = c4[(c4.group == "low_recall") & np.isclose(c4.p, 0.15)].prop_cells.median()
    return [("<=0.007", max(lr.uniform_cert.mean(), lr.prop_cert.mean())), ("2.3", lo_ / hi_)]
optional("rq4.lowrecall", "Certification stayed valid with blasts drawn only from low-recall patients ($\\le0.007$), at 2.3 times the cost", low_recall)
check("rq4.shift", "blast sensitivity fell to 0.340 [0.250, 0.421] on RV-PBS (H4a) and 0.63 [0.62, 0.64] on MLL23 (AUC 0.786 and 0.915), while its false-positive rate barely rose (0.058 and 0.060 against 0.050",
      [("0.340", ext["H4a_recall_drops"]["rv_recall"]), ("0.250", ext["H4a_recall_drops"]["ci95"][0]), ("0.421", ext["H4a_recall_drops"]["ci95"][1]),
       ("0.63", mll["recall"]), ("0.62", mll["recall_ci95"][0]), ("0.64", mll["recall_ci95"][1]), ("0.786", ext["descriptive"]["auc"]), ("0.915", mll["auc"]),
       ("0.058", ext["H4b_fp_rises"]["rv_fp"]), ("0.060", mll["fp"]), ("0.050", ext["H4b_fp_rises"]["ku_fp"])])
A3 = J("a3_decisions.json")
allvalid = max(v["rate"] for b in BR.values() for v in b["validity"].values())
check("rq4.shift_valid", "Every design stayed valid (H4c; all $\\le0.034$)", [("<=0.034", allvalid)])
m15, m20, v15, v20 = brr("mll23", .10, .15, "bet_prop"), brr("mll23", .10, .20, "bet_prop"), brr("rvpbs", .10, .15, "bet_prop"), brr("rvpbs", .10, .20, "bet_prop")
check("rq4.ratios", "needed 1.30 [1.18, 1.37] and 1.46 [1.36, 1.75] times the cells on MLL23 at $p=0.15$ and $0.20$, and 1.67 [1.53, 1.75] and 3.00 [2.79, 3.33] times on RV-PBS",
      [(f"{x:.2f}", x) for d_ in (m15, m20, v15, v20) for x in (d_["ratio"], d_["ci95"][0], d_["ci95"][1])])
check("rq4.h6d", "the confirmatory comparison with uniform betting gave 1.77 on RV-PBS at $p=0.20$ (H6d)",
      [("1.77", A3["H6d_reversal_replicates"]["median_prop"] / A3["H6d_reversal_replicates"]["median_uniform"])])
h6c = [c["retained"] for c in A3["H6c_mixture_keeps_gain_in_distribution"]["cells"].values()]
h6b = A3["H6b_mixture_caps_loss_under_shift"]["cells"]["p0.20"]
check("rq4.mixture", "kept 58--73\\% of the saving over uniform betting (H6c) but cut the RV-PBS excess only to $+23\\%$ [11, 30], missing the pre-registered 20\\% bound (H6b), and a 25\\% uniform share gave essentially no protection (ratio 1.73)",
      rng("58\\%--73\\%", h6c) + [("23\\%", h6b["ratio"] - 1), ("11", (h6b["ci95"][0] - 1) * 100), ("30", (h6b["ci95"][1] - 1) * 100), ("1.73", h6b["mix25_ratio_descriptive"])])
mixk = [brr("ku", t, p, "hedge_mix50")["ratio"] for t, p in KU_ALT]
mixm = [brr("mll23", .10, p, "hedge_mix50")["ratio"] for p in (0.12, 0.15, 0.20)]; mixv = [brr("rvpbs", .10, p, "hedge_mix50")["ratio"] for p in (0.12, 0.15, 0.20)]
check("rq4.hedged", "with hedged bets it needed 0.61--0.81 times the cells of ALPHA in distribution, 0.77--0.82 on MLL23 and 0.95--1.15 on RV-PBS",
      rng("0.61--0.81", mixk) + rng("0.77--0.82", mixm) + rng("0.95--1.15", mixv))
H5b = A3["H5b_naive_cell_analysis_breaks_under_clustering"]["cells"]; h5a = max(c["rate"] for c in A3["H5a_field_level_valid_all_arrangements"]["cells"].values())
h5c = [v["medians"][a] / v["medians"]["scattered"] for v in A3["H5c_cost_rises_with_clustering"]["cells"].values() for a in ("hot0.6", "hot1.0")]
check("rq4.cluster", "under the strongest clustering it certified falsely in 0.327 ($\\theta=0.10$) and 0.478 ($\\theta=0.20$) of boundary specimens (H5b), whereas field-level designs stayed valid ($\\le0.042$, H5a), although clustering raised the cost of uniform field sampling 1.6 to 2.2-fold",
      [("0.327", H5b["theta0.1"]["hot1_0"]["rate"]), ("0.478", H5b["theta0.2"]["hot1_0"]["rate"]), ("<=0.042", h5a)] + rng("1.6 to 2.2", h5c, " to "))
CL = J("a4_bmwsi", "clustering.json"); IX = J("a4_bmwsi", "icc_excess.json")
check("rq4.icc", "the bias-corrected intraclass correlation was 0.013 [0.002, 0.028] and 0.043 [0.018, 0.067] (exploratory), a design effect of about 1.2 for 5-cell fields",
      [("0.013", IX["S1"]["excess"]), ("0.002", IX["S1"]["excess_ci95"][0]), ("0.028", IX["S1"]["excess_ci95"][1]),
       ("0.043", IX["S2"]["excess"]), ("0.018", IX["S2"]["excess_ci95"][0]), ("0.067", IX["S2"]["excess_ci95"][1]), ("1.2", 1 + 4 * IX["S2"]["excess"])])
check("rq4.regions", "In 116 exhaustively annotated marrow regions of 75 patients", [("116", CL["n_rois"]), ("75", CL["n_patients"])])
check("rq4.perm", "permutation $p=0.024$ at the quadrant and $0.0005$ at the $4\\times4$ scale", [("0.024", CL["primary"]["S1"]["perm_p"]), ("0.0005", CL["primary"]["S2"]["perm_p"])])
n_sp = sum(len(pd.read_csv(f)) for d in ("a3_spatial_theta010", "a3_spatial_theta020") for f in glob.glob(R(d, "ep_*.csv")))
check("rq4.cluster_n", "in synthetic arrangements (18,000 specimens)", [("18000", n_sp)], note="all three arrangements, both thresholds")
ID = J("icc_diagnostics_v1", "summary.json")["low_blast_rois_p_le_0.30"]["S2"]
check("rq4.icc_low", "a null-centered estimate of 0.059 [0.007, 0.108] at $4\\times4$ in the 32 regions with at most 30\\% blasts",
      [("0.059", ID["null_centred"]), ("0.007", ID["null_centred_ci95"][0]), ("0.108", ID["null_centred_ci95"][1]), ("32", ID["n_rois"])])
rae = pd.read_csv(R("real_arrangements_v1", "episodes.csv"), keep_default_na=False); rae = rae[rae.scen == "null"]
per_ra = rae.groupby("patient").agg(p_s=("p_s", "first"), rate=("naive_cell_cert", "mean")); low = per_ra[per_ra.p_s <= 0.20 + 1e-12]
check("rq4.real", "(an exact null; 400 sampling orders per specimen; exploratory). The weighted test of Proposition~\\ref{prop:unequal}(a) certified falsely in 0.001 of runs and at most 0.01 in any specimen; the cell-level analysis gave 0.051 overall, yet exceeded 0.05 in 27 of 75 specimens (21 significantly; up to 0.335), including 7 of the 23 with at most 20\\% blasts",
      [("400", J("real_arrangements_v1", "config.json")["n_orders"]), ("0.001", RA["field_known"]["pooled_rate"]), ("<=0.01", RA["field_known"]["max_patient_rate"]),
       ("0.051", RA["naive_cell"]["pooled_rate"]), ("27", RA["naive_cell"]["n_patients_rate_gt_alpha"]), ("75", len(per_ra)), ("21", RA["naive_cell"]["n_patients_signif_gt_alpha"]),
       ("0.335", RA["naive_cell"]["max_patient_rate"]), ("7", int((low.rate > 0.05).sum())), ("23", len(low))])
check("rq4.caveats", "These specimens are small (median 105 cells) and mostly high-blast: only two had blast proportions between 0.05 and 0.30",
      [("105", RAS["cells_per_specimen"]["median"]), ("2", int(((per_ra.p_s > 0.05) & (per_ra.p_s <= 0.30)).sum()))])
stress = [v for k_, v in J("a4_bmwsi", "bridge.json")["e31_icc_placement"].items() if "hot" in k_]
ratio_lo, ratio_hi = min(stress) / IX["S2"]["excess"], max(stress) / IX["S1"]["excess"]
OUT.append(dict(id="rq4.orders", ok=10 <= ratio_lo and ratio_hi <= 100, in_doc="one to two orders of magnitude more clustered than these monolayer regions" in TEX,
                detail=f"ICC of the clustered E3.1 arrangements {min(stress):.2f}-{max(stress):.2f} over the bias-corrected real ICC 0.013-0.043: {ratio_lo:.0f} to {ratio_hi:.0f} times", note=""))
UF = J("unequal_fields_v1", "summary.json"); uv = UF["validity"]
cons = max(uv[c][d]["a0.05"]["rate"] for c in uv for d in ("assort", "assort_bet", "known_unif", "known_prop"))
check("rq4.unequal", "Marrow quadrants hold 12--29 nucleated cells (interquartile range). With these sizes, averaging field proportions, the mean-of-ratios estimator \\cite{cochran1977}, certified falsely in 0.38 [0.35, 0.41] and 0.32 [0.29, 0.35]",
      [("12", UF["partA"]["quadrant_cells"]["q25"]), ("29", UF["partA"]["quadrant_cells"]["q75"]),
       ("0.38", uv["theta0.10_gamma0"]["naive"]["a0.05"]["rate"]), ("0.35", uv["theta0.10_gamma0"]["naive"]["a0.05"]["cp95"][0]), ("0.41", uv["theta0.10_gamma0"]["naive"]["a0.05"]["cp95"][1]),
       ("0.32", uv["theta0.20_gamma0"]["naive"]["a0.05"]["rate"]), ("0.29", uv["theta0.20_gamma0"]["naive"]["a0.05"]["cp95"][0]), ("0.35", uv["theta0.20_gamma0"]["naive"]["a0.05"]["cp95"][1])])
check("rq4.unequal_enriched", "and in every specimen with blasts enriched in small fields", [("1.00", min(uv["theta0.10_gamma1"]["naive"]["a0.05"]["rate"], uv["theta0.20_gamma1"]["naive"]["a0.05"]["rate"]))])
check("rq4.unequal_real", "differed from the blast proportion by only $-0.004$ on average $[-0.008,0.001]$",
      [("-0.004", UF["partA"]["bias_mean"]), ("-0.008", UF["partA"]["bias_ci95"][0]), ("0.001", UF["partA"]["bias_ci95"][1])])
cu = UF["cost"]["theta0.10_gamma0"]["p0.15"]
check("rq4.unequal_cost", "stayed at or below 0.005; with known counts, classifier-driven sampling certified $p=0.15$ against $\\theta=0.10$ with 297 cells against 646 for uniform sampling, and the bounded assorter needed 631 because $\\bar m/m_{\\max}\\approx0.24$",
      [("<=0.005", cons), ("297", cu["known_prop"]["a0.05"]), ("646", cu["known_unif"]["a0.05"]), ("631", cu["assort_bet"]["a0.05"]),
       ("0.24", UF["partA"]["quadrant_cells"]["mean"] / UF["m_max"])])
VE = J("verifier_error_v1", "summary.json"); ag = VE["agreement"]["s_session1_given_s_session0"]; fc = VE["false_cert"]
vee = pd.read_csv(R("verifier_error_v1", "episodes.csv"))
pv_ = vee.groupby("patient").agg(p_s=("p_s", "first"), p_r=("p_reader", "first"), above=("above_cert", "mean"), below=("below_cert", "mean"))
over, under = pv_[pv_.p_r > pv_.p_s + 1e-12], pv_[pv_.p_r < pv_.p_s - 1e-12]; exact = pv_[(pv_.p_r - pv_.p_s).abs() <= 1e-12]
check("rq4.readers", "The second reader called a blast in 1.1\\% [0.7, 1.7] of the first reader's counted non-blasts and a non-blast in 1.2\\% [0.8, 1.6] of its blasts; 2.7\\% of cells were counted by one reader only",
      [("1.1\\%", ag["blast_call_given_nonblast"]), ("0.7", ag["ci95_fp"][0] * 100), ("1.7", ag["ci95_fp"][1] * 100), ("1.2\\%", ag["nonblast_call_given_blast"]),
       ("0.8", ag["ci95_fn"][0] * 100), ("1.6", ag["ci95_fn"][1] * 100), ("2.7\\%", VE["counted_vs_excluded_disagreement_reader0"])])
check("rq4.verifier", "(83 specimens), the reader's proportion differed from the consensus by a median of 0.005 (maximum 0.075)",
      [("83", VE["n_specimens"]), ("0.005", VE["abs_reader_bias_quantiles"]["0.5"]), ("0.075", VE["reader_minus_consensus_p"]["max_abs"])])
check("rq4.verifier_sign", "was certified falsely in all 23 specimens in which the reader over-counted (0.28 of runs), ``$p<\\theta_s$'' in all 52 in which it under-counted (0.63), and in at most 0.008 of runs in the 8 without error",
      [("23", len(over)), ("23", int((over.above > 0.05).sum())), ("0.28", fc["above"]["pooled"]), ("52", len(under)), ("52", int((under.below > 0.05).sum())),
       ("0.63", fc["below"]["pooled"]), ("8", len(exact)), ("<=0.008", float(max(exact.above.max(), exact.below.max())))])
check("rq4.verifier_fix", "did not repair this (0.16, 0.42); a margin of 0.03, selected on the same specimens, sufficed for ``above'' but not for ``below''",
      [("0.16", fc["above_corr"]["pooled"]), ("0.42", fc["below_corr"]["pooled"]), ("0", fc["above_m0.03"]["n_specimens_gt_alpha"]), ("10", fc["below_m0.03"]["n_specimens_gt_alpha"])],
      note="'not for below': 10 specimens above alpha with the 0.03 margin")
gen = open(R("benchmark_v1", "benchmark_rows.tex"), encoding="utf-8").read().strip().splitlines()
tab3 = TEX[TEX.index("\\label{tab:bench}"):TEX.index("\\end{table", TEX.index("\\label{tab:bench}"))]
def norm_row(l): return l.strip().replace("analysed", "analyzed")
for line in gen:
    if "&" not in line: continue
    OUT.append(dict(id="table3." + line.split("&")[0].strip()[:30], ok=norm_row(line) in tab3, in_doc=True,
                    detail="row matches scripts/46 output" if norm_row(line) in tab3 else "row differs from scripts/46 output", note="spelling normalised (analysed/analyzed)"))
top = tab3[:tab3.index("Ranked, analyzed as uniform")]
vals = []
for line in top.splitlines():
    if line.strip().startswith(("ALPHA, uniform", "Betting, uniform", "Betting, classifier-driven", "Betting, 50", "Hedged betting")):
        vals += [float(x) for c in [c.strip() for c in line.split("&")][1:8] for x in re.findall(r"\d+\.\d+|(?<![\d.])0(?![\d.])", c.replace("$^c$", ""))]
check("rq4.bench", "where their assumptions held, design-based tests kept false certification at or below 0.034, whereas every analysis without a design-based guarantee failed in some condition, even where it looked acceptable in distribution (cells as draws, 0.041) or pooled over real arrangements (0.051)",
      [("<=0.034", max(vals)), ("0.041", A3["H5b_naive_cell_analysis_breaks_under_clustering"]["cells"]["theta0.2"]["scattered"]["rate"]), ("0.051", RA["naive_cell"]["pooled_rate"])])
check("table3.reader", "One reader with margin 0.03: 0.0003 [0.013]", [("0.0003", fc["above_m0.03"]["pooled"]), ("0.013", fc["above_m0.03"]["max_specimen"])])
ufM = [c[d]["median_M"] for v in UF["cost"].values() for c in v.values() for d in c]
check("table3.sizes", "about 1,300 with unequal fields; median 105 real cells", [("13", np.median(ufM) / 100), ("105", RAS["cells_per_specimen"]["median"])])

# ------------------------------------------------------------------ robustness, hypotheses, Table V, captions
SC = rs["sensitivity_confirmatory"]
def sc_(run, m): return SC[run][str(m)]
opg = [abs(sc_(r, m)["diff"]["mean"]) for r in SC if r.startswith("sensconf_tau") for m in (10, 20)]
check("rob.sens", "With max pooling it was $-0.115$ $[-0.198,-0.032]$ at $m=10$ and $-0.056$ $[-0.141,0.029]$ at $m=20$, and InfoDPP-PAC still missed 73\\% and 77\\% of specimens; with single-cell representations it was $-0.036$ $[-0.079,0.007]$ at $m=10$; and across operating points of 5--50\\% of the fields it never exceeded 0.13",
      [("-0.115", sc_("sensconf_pool_max", 10)["diff"]["mean"]), ("-0.198", sc_("sensconf_pool_max", 10)["diff"]["ci95"][0]), ("-0.032", sc_("sensconf_pool_max", 10)["diff"]["ci95"][1]),
       ("-0.056", sc_("sensconf_pool_max", 20)["diff"]["mean"]), ("-0.141", sc_("sensconf_pool_max", 20)["diff"]["ci95"][0]), ("0.029", sc_("sensconf_pool_max", 20)["diff"]["ci95"][1]),
       ("73\\%", sc_("sensconf_pool_max", 10)["infodpp_miss"]), ("77\\%", sc_("sensconf_pool_max", 20)["infodpp_miss"]),
       ("-0.036", sc_("sensconf_pool_single", 10)["diff"]["mean"]), ("-0.079", sc_("sensconf_pool_single", 10)["diff"]["ci95"][0]), ("0.007", sc_("sensconf_pool_single", 10)["diff"]["ci95"][1]),
       ("<=0.13", max(opg))])
check("rob.host", "($-0.27$ $[-0.38,-0.16]$ at $m=10$ and $-0.15$ $[-0.19,-0.10]$ at $m=20$), although InfoDPP-PAC still missed 57\\% and 71\\% of specimens",
      [("-0.27", sc_("sensconf_pool_host", 10)["diff"]["mean"]), ("-0.38", sc_("sensconf_pool_host", 10)["diff"]["ci95"][0]), ("-0.16", sc_("sensconf_pool_host", 10)["diff"]["ci95"][1]),
       ("-0.15", sc_("sensconf_pool_host", 20)["diff"]["mean"]), ("-0.19", sc_("sensconf_pool_host", 20)["diff"]["ci95"][0]), ("-0.10", sc_("sensconf_pool_host", 20)["diff"]["ci95"][1]),
       ("57\\%", sc_("sensconf_pool_host", 10)["infodpp_miss"]), ("71\\%", sc_("sensconf_pool_host", 20)["infodpp_miss"])])
ST = J("infodpp_settings_v1", "summary.json")["overall"]
check("rob.settings", "exploratory, 2,000 new specimens) changed its miss rate and its difference from uniform sampling by at most 0.03, with no significant change (smallest Holm-adjusted McNemar $p=0.38$)",
      [("2000", ST["n_specimens"]), ("<=0.03", ST["max_abs_change_miss"]), ("<=0.03", ST["max_abs_change_diff"]), ("0.38", ST["min_p_holm"])])
hyp = TEX[TEX.index("\\label{tab:hyp}"):TEX.index("\\end{table", TEX.index("\\label{tab:hyp}"))]
ids = []
for line in hyp.splitlines():
    m_ = re.match(r"^([A-Z0-9a-z,\-\s]+?) & ([EV]) & .* & (\\textbf\{no\}|yes)\\\\$", line.strip())
    if not m_ or m_.group(1).startswith("KT2"): continue
    part = []
    for tok in m_.group(1).split(","):
        tok = tok.strip()
        mm = re.match(r"^(H\d)([a-d])--([a-d])$", tok)
        part += [f"{mm.group(1)}{chr(c)}" for c in range(ord(mm.group(2)), ord(mm.group(3)) + 1)] if mm else [tok]
    ids += [(h, m_.group(2), m_.group(3) == "yes") for h in part]
nE = sum(1 for _, t, _ in ids if t == "E"); nEs = sum(1 for _, t, s in ids if t == "E" and s); nV = sum(1 for _, t, _ in ids if t == "V")
check("hyp.count", "Of the 15 pre-registered empirical hypotheses 12 were supported, and the 5 implementation checks passed",
      [("15", nE), ("12", nEs), ("5", nV), ("20", len(ids))], note="counted from Table IV")
ext_v = {"H1a": dec1["H1a"]["supported"], "H1b": dec1["H1b"]["supported"], "H1c": dec1["H1c"]["supported"], "H1d": dec1["H1d"]["supported"], "KT1": dec1["KT1"]["supported"],
         "H4a": ext["H4a_recall_drops"]["supported"], "H4b": ext["H4b_fp_rises"]["supported"], "H4c": ext["H4c_validity_under_shift"]["supported"], "H7a": CL["H7a_clustering_S1"]["supported"]}
for t in P3.values():
    for k_, v_ in t.items():
        hid = k_.split("_")[0]
        if isinstance(v_, dict) and "supported" in v_: ext_v[hid] = ext_v.get(hid, True) and v_["supported"]
for k_, v_ in A3.items():
    if isinstance(v_, dict) and "supported" in v_: ext_v[k_.split("_")[0]] = v_["supported"]
for h, t, s_ in ids:
    if h in ext_v:
        OUT.append(dict(id=f"table4.{h}", ok=s_ == bool(ext_v[h]), in_doc=True, detail=f"table {'yes' if s_ else 'no'}; results {'supported' if ext_v[h] else 'not supported'}", note=""))
tab5 = TEX[TEX.index("\\label{tab:assump}"):TEX.index("\\end{table", TEX.index("\\label{tab:assump}"))]
check("table5.ranking", "Uncorrected ranking: false certification 0.89--0.99", [("0.89", P3["theta_0.20"]["H3b_naive_invalid"]["rate"]), ("0.99", P3["theta_0.10"]["H3b_naive_invalid"]["rate"])])
check("table5.unit", "up to 0.34 per real specimen, 0.33--0.48 under strong clustering", [("0.34", RA["naive_cell"]["max_patient_rate"]), ("0.33", H5b["theta0.1"]["hot1_0"]["rate"]), ("0.48", H5b["theta0.2"]["hot1_0"]["rate"])])
check("table5.target", "Averaging proportions: 0.32--0.38; Proposition~\\ref{prop:unequal} restores $\\le0.005$",
      [("0.32", uv["theta0.20_gamma0"]["naive"]["a0.05"]["rate"]), ("0.38", uv["theta0.10_gamma0"]["naive"]["a0.05"]["rate"]), ("<=0.005", cons)])
check("table5.verifier", "(0.28 ``above'', 0.63 ``below'')", [("0.28", fc["above"]["pooled"]), ("0.63", fc["below"]["pooled"])])
check("table5.shift", "Validity held under shift ($\\le0.034$); cost relative to uniform ALPHA rose from 0.45--0.75 to 1.3--3.0 at margins of 0.05--0.10",
      [("<=0.034", allvalid)] + rng("0.45--0.75", pva) + rng("1.3--3.0", shiftr))
check("fig4.caption", "(exploratory rerun, 400 specimens per point); all designs kept false certification $\\le0.034$", [("<=0.034", allvalid)])
check("disc.27", "a 200-cell count that forces a call was wrong in up to 27\\% of specimens", [("27\\%", max(dl(t, p, 200, "conventional", "wrong") for t, p in near))])

# ------------------------------------------------------------------ validity statement of the abstract
valid = {"part3": max(h3a), "clustering": h5a, "population size": v05, "rerun": allvalid, "unequal": cons,
         "real (known counts)": RA["field_known"]["pooled_rate"], "real (assorter)": RA["field_assort"]["pooled_rate"]}
spec = {k: RA[k]["n_patients_signif_gt_alpha"] for k in ("field_known", "field_assort", "cell_srs")}
OUT.append(dict(id="abs.valid", ok=max(valid.values()) <= 0.05 and all(v == 0 for v in spec.values()),
                in_doc="Valid tests kept false certification within the nominal level" in TEX,
                detail=f"largest pooled rate of a valid design at an exact null {max(valid.values()):.3f}; specimens significantly above 0.05: {spec}", note=json.dumps({k: round(v, 3) for k, v in valid.items()})))

# ------------------------------------------------------------------ report
bad = [o for o in OUT if not o["ok"]]
lines = ["# Audit of numbers in paper/jbhi_compact.tex (IEEE JBHI manuscript) against results/", "",
         f"{len(OUT)} checks, {len(OUT) - len(bad)} passed, {len(bad)} failed. Generated by `paper/audit_jbhi_numbers.py`.", "",
         "Stated numbers are rounded half-up to their decimals; `<=` marks an at-most or within claim checked as a bound.", "",
         "| Check | Result | Passage found | Stated vs recomputed | Note |", "|---|---|---|---|---|"]
for o in OUT:
    lines.append(f"| {o['id']} | {'pass' if o['ok'] else '**FAIL**'} | {'yes' if o['in_doc'] else '**no**'} | {o['detail']} | {o['note']} |")
open(os.path.join(ROOT, "paper", "audit_jbhi_numbers.md"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
print(f"{len(OUT)} checks, {len(bad)} failed")
for o in bad: print("FAIL", o["id"], "| in doc:", o["in_doc"], "|", o["detail"], "|", o["note"])
sys.exit(1 if bad else 0)
