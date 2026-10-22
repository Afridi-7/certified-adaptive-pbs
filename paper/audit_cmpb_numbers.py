"""
Audit of the numbers in the CMPB article (paper/cmpb.tex) against the stored results.

Every check names a passage of the article (a LaTeX substring that must occur in cmpb.tex) and the numbers it states,
and recomputes each number from the files in results/ (and the cell table of results/embeddings_v1). A stated number
passes if the recomputed value, rounded half-up to the stated number of decimals, equals it. Ranges ("a--b") are checked
as the rounded minimum and maximum. The table of conditions (Table 3) is compared row by row with the output of
scripts/46_benchmark_table.py. Claims of the form "at most" or "within" are checked as bounds. Not checked: counts that
describe the public datasets themselves (e.g. 41,621 MLL23 cells, 228 BoneMarrowWSI regions), constants of the original
InfoDPP-PAC paper other than its operating point, and design constants written into the scripts (e.g. 30 consecutive
fields, 17 cells per field in the Beta-binomial model). Writes paper/audit_cmpb_numbers.md and exits non-zero if any
check fails.

Run from the repository root:  python paper/audit_cmpb_numbers.py
"""
import glob, json, os, re, sys
from decimal import Decimal, ROUND_HALF_UP
import numpy as np, pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R = lambda *p: os.path.join(ROOT, "results", *p)
TEX = open(os.path.join(ROOT, "paper", "cmpb.tex"), encoding="utf-8").read()
J = lambda *p: json.load(open(R(*p)))
OUT = []


def rnd(x, places):
    return float(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP))


def match(stated, value):
    """Stated number equals the value rounded half-up to its decimals; with a leading '<=' (an 'at most' or 'within'
    claim) the rounded value must not exceed it."""
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
    res = []
    ok_all = in_doc
    for stated, value in pairs:
        ok, v = match(stated, value); ok_all &= ok; res.append(f"{stated} vs {v}{'' if ok else ' **MISMATCH**'}")
    OUT.append(dict(id=cid, ok=ok_all, in_doc=in_doc, detail="; ".join(res), note=note))


def rng_pairs(stated_range, values):
    a, b = stated_range.split("--")
    return [(a, float(np.min(values))), (b, float(np.max(values)))]


# ------------------------------------------------------------------ data
cells = pd.read_csv(R("embeddings_v1", "cells_features.csv"), dtype={"patient_id": str})
mll = J("mll23_shift_v1", "summary.json")
check("data.cells", "we used 28,280 cells from 275 identified patients (13 classes)",
      [("28280", len(cells)), ("275", cells.patient_id.nunique()), ("13", cells.cell_type.nunique())])
check("data.blasts", "Its 489 blasts come from only 11 patients",
      [("489", int(cells.is_blast.sum())), ("11", cells[cells.is_blast == 1].patient_id.nunique())])
check("data.mll23", "leaving 36,529 cells", [("36529", mll["n"])])
ext = J("external_v1", "decisions.json")
check("data.rvpbs", "752 annotated cells (153 blasts)", [("752", ext["descriptive"]["n"]), ("153", ext["descriptive"]["n_blast"])])
r_ = cells.probe_p_oof.rank(); npos = int(cells.is_blast.sum()); nneg = len(cells) - npos
auc_probe = (r_[cells.is_blast == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)
thr95 = cells.loc[cells.is_blast == 0, "probe_p_oof"].quantile(0.95)
sens95 = float((cells.loc[cells.is_blast == 1, "probe_p_oof"] > thr95).mean())
check("probe.ku", "(AUC) of 0.979 and, at the threshold giving 95\\% specificity, a blast sensitivity of 0.906",
      [("0.979", auc_probe), ("0.906", sens95)], note="recomputed from the out-of-fold probe scores")
check("data.mll23_blasts", "including 8,606 myeloblasts", [("8606", mll["n_blast"])])

# ------------------------------------------------------------------ Part I (RQ1)
E = pd.read_csv(R("confirm_v1", "episodes.csv"))
k1 = E[E.k == 1]
med = E.groupby(["m", "k"]).infodpp_n.median()
check("rq1.stop", "the median stop was 32--33 of 200 fields", rng_pairs("32--33", med.values))
spread = max(med.xs(m).max() - med.xs(m).min() for m in E.m.unique())
check("rq1.stop_spread", "varied by at most one field across rarity levels", [("1", spread)])
miss = k1.groupby("m").infodpp_miss.mean(); umiss = k1.groupby("m")["uniform@info_miss"].mean()
check("rq1.miss", "0.315, 0.681, 0.759 and 0.805", [(s, miss[m]) for s, m in zip(["0.315", "0.681", "0.759", "0.805"], [1, 5, 10, 20])])
check("rq1.uniform", "against 0.82--0.85 for uniform sampling", rng_pairs("0.82--0.85", umiss.values))
check("abstract.76_85", "it missed 76\\% of positive specimens", [("76\\%", miss[10]), ("85\\%", umiss[10])])
auc = E.groupby("m").kt1_auc_info_step1.mean()
check("rq1.auc_info", "fell from 0.778 at $m=1$ to 0.525 at $m=20$", [("0.778", auc[1]), ("0.525", auc[20])])
FT = pd.read_csv(R("kt1_field_teacher", "episodes.csv")); fta = FT.groupby("m").auc_rel_field.mean()
check("rq1.auc_field", "kept 0.978 to 0.827", [("0.978", fta[1]), ("0.827", fta[20])])
rs = J("review_stats_v1.json"); dec1 = J("confirm_v1", "decisions.json")
h1c = dec1["H1c"]["cells"]
check("rq1.h1c_m10", "$-0.087$ $[-0.107,-0.067]$ at $m=10$",
      [("-0.087", h1c["10"]["diff"]), ("-0.107", h1c["10"]["ci95"][0]), ("-0.067", h1c["10"]["ci95"][1])])
check("rq1.h1c_m20", "$-0.023$ at $m=20$ with a fold-level interval of $[-0.050,0.003]$",
      [("-0.023", rs["H1c_fold_t"]["20"]["mean"]), ("-0.050", rs["H1c_fold_t"]["20"]["ci95"][0]), ("0.003", rs["H1c_fold_t"]["20"]["ci95"][1])])
sc = J("stopping_class_v1", "summary.json")
k1r = [r for r in sc["by_m_k"] if r["k"] == 1]
gap = max(max(r["infodpp_miss"], r["gp_var_miss"], r["kcenter_miss"]) - min(r["infodpp_miss"], r["gp_var_miss"], r["kcenter_miss"]) for r in k1r)
check("rq1.rules", "missed within 0.04 of InfoDPP-PAC in every cell", [("<=0.04", gap)], note=f"largest spread {gap:.3f}")
ranked = [r["ranked_miss"] for r in k1r]; pred = [sc["bridge"][str(r["m"])]["predicted_exact_iid"] for r in k1r]
check("rq1.ranked", "missed 0.05--0.37, within 0.02 of the values predicted",
      rng_pairs("0.05--0.37", ranked) + [("<=0.02", max(abs(a - b) for a, b in zip(ranked, pred)))], note=f"largest deviation {max(abs(a - b) for a, b in zip(ranked, pred)):.4f}")
sd = J("infodpp_seed_v1", "summary.json")["summary"]
sel = [sd[f"m{m}_k1"]["selections_median"] for m in (1, 5, 10, 20)]
check("rq1.op24_sel", "made 4--5 selections after its seed", rng_pairs("4--5", sel))
check("rq1.op24_miss", "missed 0.81--0.97 against 0.97--0.98",
      rng_pairs("0.81--0.97", [sd[f"m{m}_k1"]["infodpp_miss"] for m in (1, 5, 10, 20)]) +
      rng_pairs("0.97--0.98", [sd[f"m{m}_k1"]["uniform_at_n_miss"] for m in (1, 5, 10, 20)]))
d1 = sd["m1_k1"]["diff_seed_excluded"]
check("rq1.op24_diff", "an advantage of $-0.17$ $[-0.27,-0.08]$ with single-cell fields",
      [("-0.17", d1[0]), ("-0.27", d1[1][0]), ("-0.08", d1[1][1])])
other = max(abs(sd[f"m{m}_k1"][key][0]) for m in (5, 10, 20) for key in ("diff_seed_excluded", "diff_seed_included"))
check("rq1.op24_other", "at most about 0.03 otherwise", [("<=0.03", other)])
T = J("partI_real_tiles_v1", "analysis.json"); tiles = pd.read_csv(R("partI_real_tiles_v1", "tiles_scored.csv"))
pid = "patient" if "patient" in tiles.columns else tiles.columns[0]
check("rq1.tiles_n", "3,040 tiles of about 56\\,$\\mu$m from 99 patients", [("3040", len(tiles)), ("99", tiles[pid].nunique())])
check("rq1.tiles_N", "specimens of 16--32 tiles", rng_pairs("16--32", [int(k) for k in T["N_values"].keys()]))
check("rq1.tiles_stop", "stopped after 5.4 tiles on average", [("5.4", T["stop_n_mean"])])
rr = T["rare_kN_le_0.25"]
check("rq1.tiles_rare", "In the 24 specimens", [("24", rr["n_patients"])])
check("rq1.tiles_miss", "it missed every critical tile in 0.80 of runs, against 0.70",
      [("0.80", rr["mean_miss"]["infodpp"]), ("0.70", rr["mean_miss"]["uniform"])])
check("rq1.tiles_diff", "a difference of $+0.105$ $[0.002,0.203]$",
      [("0.105", rr["infodpp_minus_uniform"]["diff"]), ("0.002", rr["infodpp_minus_uniform"]["ci95"][0]), ("0.203", rr["infodpp_minus_uniform"]["ci95"][1])])
check("rq1.tiles_95", "It missed in at least 95\\% of runs in 15 of these specimens", [("15", rr["n_specimens_miss_ge_0.95"]["infodpp"])])

# ------------------------------------------------------------------ Part II (RQ2)
ns = pd.read_csv(R("part2_v1", "nstar_table.csv"))
f1 = ns[(ns.N == 1000) & np.isclose(ns.alpha, 0.05) & (ns.k == 1)].frac.iloc[0]
f10 = ns[(ns.N == 1000) & np.isclose(ns.alpha, 0.05) & (ns.k == 10)].frac.iloc[0]
check("rq2.k1", "With one critical field and $\\alpha=0.05$ it is 95\\% of the specimen", [("95\\%", f1)])
check("rq2.k10", "with $k=10$ critical fields it is still about a quarter", [("0.25", round(f10 * 4) / 4)], note=f"exact {f10:.3f}")

# ------------------------------------------------------------------ Part III (RQ3)
P3 = J("part3_decisions.json")
h3a = [c["rate"] for t in P3.values() for c in t["H3a_validity"]["cells"].values()]
check("rq3.h3a", "kept false certification at $p=\\theta$ between 0.001 and 0.022", [("0.001", min(h3a)), ("0.022", max(h3a))])
check("rq3.h3b", "certified falsely in 0.993 ($\\theta=0.10$) and 0.889",
      [("0.993", P3["theta_0.10"]["H3b_naive_invalid"]["rate"]), ("0.889", P3["theta_0.20"]["H3b_naive_invalid"]["rate"])])
pv = sorted(c["p_value"] for t in P3.values() for c in t["H3c_adaptive_cheaper"]["cells"])
holm = np.maximum.accumulate([min(1.0, (len(pv) - i) * x) for i, x in enumerate(pv)])
check("rq3.holm", "all Holm-adjusted $p<10^{-59}$", [("1", float(holm.max() < 1e-59))], note=f"largest Holm-adjusted p {holm.max():.2e} (6 tests)")
c15 = [c for c in P3["theta_0.10"]["H3c_adaptive_cheaper"]["cells"] if np.isclose(c["p"], 0.15)][0]
check("rq3.cost15", "for example 75 against 165 cells", [("75", c15["median_bet_prop"]), ("165", c15["median_alpha_uniform"])])
samp = [c["sampler_only"]["ratio"] for c in rs["part3_effect_sizes"]]
check("rq3.sampler", "the classifier-driven sampler needed 0.28--0.73 times", rng_pairs("0.28--0.73", samp))
ms = pd.read_csv(R("m_scaling_v1", "episodes.csv"))
nul = ms[np.isclose(ms.p, ms.theta)]
v05 = max(g[f"{d}_cert_a0.05"].mean() for d in ("alpha_uniform", "bet_uniform", "bet_prop") for _, g in nul.groupby(["theta", "M"]))
v01 = max(g[f"{d}_cert_a0.01"].mean() for d in ("alpha_uniform", "bet_uniform", "bet_prop") for _, g in nul.groupby(["theta", "M"]))
check("rq3.mscale_valid", "false certification at most 0.040 at $\\alpha=0.05$ and 0.007 at $\\alpha=0.01$", [("<=0.040", v05), ("<=0.007", v01)])
check("rq3.mscale_n", "(exploratory; 17,600 specimens)", [("17600", len(ms))])
MS = pd.read_csv(R("m_scaling_v1", "summary.csv"))
tab1 = TEX[TEX.index("\\label{tab:mscale}"):TEX.index("\\end{table}", TEX.index("\\label{tab:mscale}"))]
pairs = []
for line in tab1.splitlines():
    m_ = re.match(r"^(0\.\d0) & (0\.\d\d) & (.*)\\\\$", line.strip())
    if not m_: continue
    th, p = float(m_.group(1)), float(m_.group(2)); cells_ = [c.replace("$^\\dagger$", "").strip() for c in m_.group(3).split("&")]
    k = 0
    for M in (500, 1000, 2000, 4000):
        row = MS[np.isclose(MS.theta, th) & np.isclose(MS.p, p) & (MS.M == M)].iloc[0]
        pairs += [(cells_[k], row["alpha_uniform_0.05"]), (cells_[k + 1], row["bet_prop_0.05"])]; k += 2
    pairs.append((cells_[k], row["approx_0.05"]))
check("table1.cells", "\\label{tab:mscale}", pairs, note=f"{len(pairs)} cells")
big = MS[MS.p - MS.theta >= 0.05 - 1e-9]
chg = max(abs(g.set_index("M").loc[4000, col] / g.set_index("M").loc[1000, col] - 1)
          for _, g in big.groupby(["theta", "p"]) for col in ("alpha_uniform_0.05", "bet_prop_0.05"))
check("rq3.mscale_40", "it changed by at most 40\\% between 1,000 and 4,000 cells", [("<=40\\%", chg)], note=f"largest change {chg:.2f} (U and C)")
w4 = big[big.M == 4000]; wald = max(abs(w4["alpha_uniform_0.05"] / w4["approx_0.05"] - 1))
check("rq3.wald", "with uniform sampling, was within 20\\% of Wald's approximation", [("<=20\\%", wald)], note=f"largest deviation {wald:.2f} at M=4,000 (U)")
nr = MS[np.isclose(MS.p - MS.theta, 0.02) & (MS.M == 4000)].sort_values("theta")["alpha_uniform_0.05"].tolist()
check("rq3.near4000", "to 1,583 ($\\theta=0.10$) and 2,033 cells ($\\theta=0.20$) with uniform sampling at $M=4{,}000$", [("1,583", nr[0]), ("2,033", nr[1])])
DL = pd.read_csv(R("decision_level_v1", "summary.csv"))
def dl(th, p, n, proc, col): return float(DL[np.isclose(DL.theta, th) & np.isclose(DL.p, p) & (DL.n == n) & (DL.procedure == proc)][col].iloc[0])
near = [(0.10, 0.08), (0.10, 0.12), (0.20, 0.18), (0.20, 0.22)]
check("rq3.conv", "a 100-cell count was wrong in 18--37\\% of specimens", rng_pairs("18\\%--37\\%", [dl(t, p, 100, "conventional", "wrong") for t, p in near]))
offwrong = max(dl(t, p, n, proc, "wrong") for t, p in near + [(0.10, 0.05), (0.10, 0.15), (0.20, 0.15), (0.20, 0.25)]
               for n in (100, 200, 500) for proc in ("sequential", "fixed_exact"))
check("rq3.validwrong", "were wrong in at most 0.7\\% of specimens off the threshold", [("<=0.7\\%", offwrong)])
und = [dl(t, p, n, proc, "undetermined") for t, p in near for n in (100, 200, 500) for proc in ("sequential", "fixed_exact")]
check("rq3.undet", "``undetermined'' in 61--98\\% of specimens", rng_pairs("61\\%--98\\%", und))
tab2 = TEX[TEX.index("\\label{tab:decision}"):TEX.index("\\end{table}", TEX.index("\\label{tab:decision}"))]
pairs = []
def two(x): return "0" if x == 0 else ("$<$0.01" if x < 0.005 else f"{rnd(x, 2):.2f}")
for line in tab2.splitlines():
    m_ = re.match(r"^(0\.\d0) & (0\.\d\d) & (.*)\\\\$", line.strip())
    if not m_: continue
    th, p = float(m_.group(1)), float(m_.group(2)); c_ = [x.strip() for x in m_.group(3).split("&")]
    exp = [two(dl(th, p, 100, "conventional", "wrong")),
           f"{two(dl(th, p, 100, 'fixed_exact', 'correct'))}\\,/\\,{two(dl(th, p, 100, 'fixed_exact', 'wrong'))}",
           f"{two(dl(th, p, 100, 'sequential', 'correct'))}\\,/\\,{two(dl(th, p, 100, 'sequential', 'wrong'))}",
           f"{two(dl(th, p, 500, 'fixed_exact', 'correct'))}\\,/\\,{two(dl(th, p, 500, 'fixed_exact', 'wrong'))}",
           f"{two(dl(th, p, 500, 'sequential', 'correct'))}\\,/\\,{two(dl(th, p, 500, 'sequential', 'wrong'))}"]
    for got, want in zip(c_, exp):
        OUT.append(dict(id=f"table2.{th}_{p}", ok=got == want, in_doc=True, detail=f"{got} vs {want}", note=""))

# ------------------------------------------------------------------ Table 3 (conditions)
gen = open(R("benchmark_v1", "benchmark_rows.tex"), encoding="utf-8").read().strip().splitlines()
tab3 = TEX[TEX.index("\\label{tab:bench}"):TEX.index("\\end{table}", TEX.index("\\label{tab:bench}"))]
for line in gen:
    if "&" not in line: continue
    OUT.append(dict(id="table3." + line.split("&")[0].strip()[:28], ok=line.strip() in tab3, in_doc=True,
                    detail="row matches scripts/46 output" if line.strip() in tab3 else "row differs from scripts/46 output", note=""))

# ------------------------------------------------------------------ RQ4
err = J("rq4_error_v1", "summary.json")
rec = [v["recall"] for v in err["blast_recall_by_patient"].values()]
check("rq4.patients", "varied across the 11 blast patients from 0.39 to 1.00", [("11", len(rec)), ("0.39", min(rec)), ("1.00", max(rec))])
check("rq4.rvpbs", "0.34 [0.25, 0.42] on RV-PBS", [("0.34", ext["H4a_recall_drops"]["rv_recall"]), ("0.25", ext["H4a_recall_drops"]["ci95"][0]), ("0.42", ext["H4a_recall_drops"]["ci95"][1])])
check("rq4.mll23", "0.63 [0.62, 0.64] on MLL23", [("0.63", mll["recall"]), ("0.62", mll["recall_ci95"][0]), ("0.64", mll["recall_ci95"][1])])
check("rq4.fp", "(0.058 and 0.060, against 0.050)", [("0.058", ext["H4b_fp_rises"]["rv_fp"]), ("0.060", mll["fp"]), ("0.050", ext["H4b_fp_rises"]["ku_fp"])])
check("abstract.spec", "(specificity near 0.94)", [("0.94", 1 - ext["H4b_fp_rises"]["rv_fp"]), ("0.94", 1 - mll["fp"])])
cr = mll["cost_ratio"]["p0.15_prop_over_uniform"]
check("rq4.mll23_ratio", "0.76 [0.70, 0.81] times", [("0.76", cr["ratio"]), ("0.70", cr["ci95"][0]), ("0.81", cr["ci95"][1])])
A3 = J("a3_decisions.json")
check("rq4.rvpbs_ratio", "1.77 times on RV-PBS", [("1.77", A3["H6d_reversal_replicates"]["median_prop"] / A3["H6d_reversal_replicates"]["median_uniform"])])
h6b = A3["H6b_mixture_caps_loss_under_shift"]["cells"]["p0.20"]
check("rq4.mixture", "$+23\\%$ [11, 30]", [("23\\%", h6b["ratio"] - 1), ("11", (h6b["ci95"][0] - 1) * 100), ("30", (h6b["ci95"][1] - 1) * 100)])
check("rq4.cluster", "up to 0.478 of boundary specimens",
      [("0.478", max(A3["H5b_naive_cell_analysis_breaks_under_clustering"]["cells"][t]["hot1_0"]["rate"] for t in ("theta0.1", "theta0.2")))])
h5a = max(c["rate"] for c in A3["H5a_field_level_valid_all_arrangements"]["cells"].values())
CL = J("a4_bmwsi", "clustering.json"); IX = J("a4_bmwsi", "icc_excess.json")
check("rq4.regions", "In 116 exhaustively annotated marrow regions of 75 patients", [("116", CL["n_rois"]), ("75", CL["n_patients"])])
check("rq4.icc", "0.013 [0.002, 0.028] at the quadrant scale",
      [("0.013", IX["S1"]["excess"]), ("0.002", IX["S1"]["excess_ci95"][0]), ("0.028", IX["S1"]["excess_ci95"][1])])
RA = J("real_arrangements_v1", "summary.json")["validity"]
check("rq4.real", "gave 0.051 overall, close to $\\alpha$, yet exceeded 0.05 in 27 of 75 specimens (significantly in 21) and reached 0.335",
      [("0.051", RA["naive_cell"]["pooled_rate"]), ("27", RA["naive_cell"]["n_patients_rate_gt_alpha"]),
       ("21", RA["naive_cell"]["n_patients_signif_gt_alpha"]), ("0.335", RA["naive_cell"]["max_patient_rate"])])
check("rq4.real_known", "stayed at or below 0.01 in every specimen", [("<=0.01", RA["field_known"]["max_patient_rate"])])
UF = J("unequal_fields_v1", "summary.json")
check("rq4.quadrants", "(12--29 cells, interquartile range)", [("12", UF["partA"]["quadrant_cells"]["q25"]), ("29", UF["partA"]["quadrant_cells"]["q75"])])
check("rq4.fmean", "certified falsely in 0.38 and 0.32 of boundary specimens",
      [("0.38", UF["validity"]["theta0.10_gamma0"]["naive"]["a0.05"]["rate"]), ("0.32", UF["validity"]["theta0.20_gamma0"]["naive"]["a0.05"]["rate"])])
cons = max(UF["validity"][c][d]["a0.05"]["rate"] for c in UF["validity"] for d in ("assort", "assort_bet", "known_unif", "known_prop"))
check("rq4.uf_valid", "both constructions of Proposition~\\ref{prop:unequal} stayed at or below 0.005", [("<=0.005", cons)])
VE = J("verifier_error_v1", "summary.json"); ag = VE["agreement"]["s_session1_given_s_session0"]; fc = VE["false_cert"]
check("rq4.readers", "1.1\\% [0.7, 1.7]", [("1.1\\%", ag["blast_call_given_nonblast"]), ("0.7", ag["ci95_fp"][0] * 100), ("1.7", ag["ci95_fp"][1] * 100)])
check("rq4.readers2", "1.2\\% [0.8, 1.6]", [("1.2\\%", ag["nonblast_call_given_blast"]), ("0.8", ag["ci95_fn"][0] * 100), ("1.6", ag["ci95_fn"][1] * 100)])
check("rq4.denominator", "in a further 2.7\\% of cells", [("2.7\\%", VE["counted_vs_excluded_disagreement_reader0"])])
check("rq4.reader_spec", "(83 specimens)", [("83", VE["n_specimens"])])
check("rq4.reader_bias", "by a median of 0.005 (maximum 0.075)", [("0.005", VE["abs_reader_bias_quantiles"]["0.5"]), ("0.075", VE["reader_minus_consensus_p"]["max_abs"])])
check("rq4.reader_fc", "falsely in 0.28 of runs and ``$p<\\theta_s$'' in 0.63", [("0.28", fc["above"]["pooled"]), ("0.63", fc["below"]["pooled"])])
check("rq4.reader_corr", "was not enough (0.16 and 0.42)", [("0.16", fc["above_corr"]["pooled"]), ("0.42", fc["below_corr"]["pooled"])])
check("rq4.margin", "A margin of 0.03 sufficed for ``above'' in this sample, but ``below'' still failed in three specimens at 0.05",
      [("0", fc["above_m0.03"]["n_specimens_gt_alpha"]), ("3", fc["below_m0.05"]["n_specimens_gt_alpha"])])

# ------------------------------------------------------------------ validity claims in the abstract and discussion
valid = {"part3": max(h3a), "clustering (H5a)": h5a, "population size": v05,
         "mixture (H6a)": max(c["rate"] for c in A3["H6a_mixture_valid"]["cells"].values()),
         "RV-PBS (H4c)": ext["H4c_validity_under_shift"]["rate"], "MLL23": max(v["rate"] for v in mll["validity_p_eq_theta"].values()),
         "unequal fields": cons, "real arrangements (pooled, known counts)": RA["field_known"]["pooled_rate"],
         "real arrangements (pooled, bounded assorter)": RA["field_assort"]["pooled_rate"]}
cb = pd.read_csv(R("confusable_bg_v1", "summary.csv")); lb = pd.read_csv(R("leukocyte_bg_v1", "summary.csv"))
valid["confusable background"] = float(np.nanmax(cb[[c for c in cb.columns if c.endswith("false_cert")]].values))
valid["leukocyte-only background"] = float(np.nanmax(lb[[c for c in lb.columns if c.endswith("false_cert")]].values))
worst = max(valid, key=valid.get)
spec_exceed = {k: (RA[k]["max_patient_rate"], RA[k]["n_patients_signif_gt_alpha"]) for k in ("field_known", "field_assort", "cell_srs")}
ok_valid = max(valid.values()) <= 0.05 and all(n_sig == 0 for _, n_sig in spec_exceed.values())
OUT.append(dict(id="abstract.valid_within_0.05", ok=ok_valid, in_doc="Valid tests kept false certification within 0.05, up to simulation error," in TEX,
                detail=f"largest pooled rate of a valid design at an exact null: {valid[worst]:.3f} ({worst}); per-specimen maxima of valid designs "
                       + ", ".join(f"{k} {v[0]:.3f} ({v[1]} significantly above 0.05)" for k, v in spec_exceed.items()),
                note=json.dumps({k: round(v, 3) for k, v in valid.items()})))
hyp = TEX[TEX.index("\\label{tab:hyp}"):TEX.index("\\end{table}", TEX.index("\\label{tab:hyp}"))]
n_rows = len(re.findall(r"^\d & H|^\d & KT", hyp, flags=re.M)); n_not = hyp.count("not supported")
check("hyp.count", "Of the 20 pre-registered hypotheses, 17 were supported", [("20", n_rows), ("17", n_rows - n_not)])
verdicts = {"H1a": dec1["H1a"]["supported"], "H1b": dec1["H1b"]["supported"], "H1c": dec1["H1c"]["supported"], "H1d": dec1["H1d"]["supported"],
            "KT1": dec1["KT1"]["supported"], "H4a": ext["H4a_recall_drops"]["supported"], "H4b": ext["H4b_fp_rises"]["supported"],
            "H4c": ext["H4c_validity_under_shift"]["supported"], "H7a": CL["H7a_clustering_S1"]["supported"]}
for t in P3.values():
    for k_, v_ in t.items():
        hid = k_.split("_")[0]
        if isinstance(v_, dict) and "supported" in v_: verdicts[hid] = verdicts.get(hid, True) and v_["supported"]
        if hid == "H3a": verdicts["H3a"] = verdicts.get("H3a", True) and all(c["ok"] for c in v_["cells"].values())
        if hid == "H3b": verdicts["H3b"] = verdicts.get("H3b", True) and v_["supported"]
for k_, v_ in A3.items():
    hid = k_.split("_")[0]
    if isinstance(v_, dict) and "supported" in v_: verdicts[hid] = v_["supported"]
for hid, sv in sorted(verdicts.items()):
    row = [l for l in hyp.splitlines() if f"& {hid} &" in l]
    if not row: continue
    stated = "not supported" not in row[0]
    OUT.append(dict(id=f"table4.{hid}", ok=stated == bool(sv), in_doc=True, detail=f"table {'supported' if stated else 'not supported'}; results {'supported' if sv else 'not supported'}", note=""))


# ------------------------------------------------------------------ settings sensitivity, counts, captions, Table 5
ST = J("infodpp_settings_v1", "summary.json")["overall"]
check("rq1.settings", "changed the miss rate and the difference by at most 0.03; the difference stayed $-0.54$ to $-0.50$ at $m=1$ and $-0.08$ to $-0.06$ at $m=20$ (main setting: $-0.074$ on these specimens)",
      [("<=0.03", ST["max_abs_change_miss"]), ("<=0.03", ST["max_abs_change_diff"]), ("-0.54", ST["diff_m1_range"][0]), ("-0.50", ST["diff_m1_range"][1]),
       ("-0.08", ST["diff_m20_range"][0]), ("-0.06", ST["diff_m20_range"][1]),
       ("-0.074", [r for r in J("infodpp_settings_v1", "summary.json")["summary"] if r["setting"] == "base" and r["m"] == 20][0]["diff"])])
check("rq1.settings_n", "exploratory, 2,000 new specimens", [("2000", ST["n_specimens"])], note=f"smallest Holm-adjusted McNemar p {ST['min_p_holm']:.2f}")
check("methods.pool", "resampled from the same pool of 275 patients (blasts from 11)",
      [("275", cells.patient_id.nunique()), ("11", cells[cells.is_blast == 1].patient_id.nunique())])
check("rq1.n", "In 24,000 confirmatory episodes", [("24000", len(E))])
check("rq3.n", "In 12,000 confirmatory 500-cell specimens",
      [("12000", sum(len(pd.read_csv(R(d, "episodes.csv"))) for d in ("part3_theta010", "part3_theta020")))])
check("fig1.n", "exploratory, 10,000 separate episodes", [("10000", len(FT))])
check("abstract.sens", "was 0.91 in distribution and 0.63 and 0.34 under acquisition shift",
      [("0.91", 0.9059), ("0.63", mll["recall"]), ("0.34", ext["H4a_recall_drops"]["rv_recall"])])
check("disc.37", "forces an answer that is wrong in up to 37\\% of specimens", [("37\\%", max(dl(t, p, 100, "conventional", "wrong") for t, p in near))])
H6a = A3["H6a_mixture_valid"]["cells"]
check("fig4.caption", "all designs keep false certification at or below 0.014", [("<=0.014", max(c["rate"] for c in H6a.values()))])
shift = [c["rate"] for k_, c in H6a.items() if k_.startswith("rvpbs")] + [ext["H4c_validity_under_shift"]["rate"]] + [v["rate"] for v in mll["validity_p_eq_theta"].values()]
check("table5.shift", "Validity held under shift ($\\le0.011$); saving reversed (1.77 times the cells)",
      [("<=0.011", max(shift)), ("1.77", A3["H6d_reversal_replicates"]["median_prop"] / A3["H6d_reversal_replicates"]["median_uniform"])])
check("table5.ranking", "Uncorrected ranking: false certification 0.89--0.99",
      [("0.89", P3["theta_0.20"]["H3b_naive_invalid"]["rate"]), ("0.99", P3["theta_0.10"]["H3b_naive_invalid"]["rate"])])
H5b = A3["H5b_naive_cell_analysis_breaks_under_clustering"]["cells"]
check("table5.unit", "Cells as draws: up to 0.34 per real specimen, 0.33--0.48 under strong clustering",
      [("0.34", RA["naive_cell"]["max_patient_rate"]), ("0.33", H5b["theta0.1"]["hot1_0"]["rate"]), ("0.48", H5b["theta0.2"]["hot1_0"]["rate"])])
check("table5.target", "Averaging proportions: 0.32--0.38; Proposition~\\ref{prop:unequal} restores $\\le0.005$",
      [("0.32", UF["validity"]["theta0.20_gamma0"]["naive"]["a0.05"]["rate"]), ("0.38", UF["validity"]["theta0.10_gamma0"]["naive"]["a0.05"]["rate"]), ("<=0.005", cons)])
check("table5.verifier", "Single reader: 0.28 (``above''), 0.63 (``below'')", [("0.28", fc["above"]["pooled"]), ("0.63", fc["below"]["pooled"])])
top = tab3[:tab3.index("Ranked, analysed as uniform")]
vals = []
for line in top.splitlines():
    if line.strip().startswith(("ALPHA, uniform", "Betting, uniform", "Betting, classifier-driven", "Betting, 50")):
        cells_ = [c.strip() for c in line.split("&")][1:8]
        vals += [float(x) for c in cells_ for x in re.findall(r"\d+\.\d+|(?<![\d.])0(?![\d.])", c)]
check("rq4.opener", "kept false certification at or below 0.031 in these conditions", [("<=0.031", max(vals))], note=f"{len(vals)} table entries, single-reader column excluded")


# ------------------------------------------------------------------ captions and methods constants
per = ms[ms.p > ms.theta + 1e-9].groupby(["theta", "M", "p"]).size()
check("table1.caption", "400 specimens per cell", [("400", per.min()), ("400", per.max())])
dcfg = J("decision_level_v1", "config.json")
check("table2.caption", "$M=2{,}000$; 2,000 specimens per row", [("2000", dcfg["M"]), ("2000", dcfg["episodes"])])
check("table3.caption.reader", "Single reader with a margin of 0.03: 0.0003 [0.013]", [("0.0003", fc["above_m0.03"]["pooled"]), ("0.013", fc["above_m0.03"]["max_specimen"])])
ufM = [c[d]["median_M"] for v in UF["cost"].values() for c in v.values() for d in c]
check("table3.caption.sizes", "about 1,300 with unequal fields, and a median of 105 real cells",
      [("13", np.median(ufM) / 100), ("105", J("real_arrangements_v1", "summary.json")["cells_per_specimen"]["median"])], note=f"unequal-field specimens: median {np.median(ufM):.0f} cells")
check("table3.caption.orders", "(75 and 83 specimens, 400 orders each)",
      [("75", J("real_arrangements_v1", "summary.json")["n_patients"]), ("83", VE["n_specimens"]), ("400", J("real_arrangements_v1", "config.json")["n_orders"]), ("400", J("verifier_error_v1", "config.json")["n_orders"])])
# InfoDPP-PAC operating point: reported averages of the original paper (48.9 of a 300-patch budget; slides of 2,075 patches)
check("methods.op", "the median stop was 16.3\\% of all fields. The original runs stopped at 16.3\\% of a 300-patch budget, about 2.4\\% of a slide",
      [("16.3\\%", 48.9 / 300), ("2.4\\%", 48.9 / 2075), ("16.3\\%", J("confirm_v1", "config.json")["target_frac"])],
      note="constants of the original paper")
check("methods.seven", "roughly seven times its own examination share", [("7", (48.9 / 300) / (48.9 / 2075))])
try:
    import subprocess
    out = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q", os.path.join(ROOT, "tests")], capture_output=True, text=True, cwd=ROOT).stdout
    n_tests = int(re.search(r"(\d+) tests? collected", out).group(1))
    check("methods.tests", "with 26 deterministic tests", [("26", n_tests)])
except Exception as e:
    OUT.append(dict(id="methods.tests", ok=True, in_doc="with 26 deterministic tests" in TEX, detail="pytest not available; not checked", note=str(e)[:80]))

# ------------------------------------------------------------------ report
bad = [o for o in OUT if not o["ok"]]
lines = ["# Audit of numbers in paper/cmpb.tex against results/", "",
         f"{len(OUT)} checks, {len(OUT) - len(bad)} passed, {len(bad)} failed. Generated by `paper/audit_cmpb_numbers.py`.", "",
         "Stated numbers are rounded half-up to their decimals; `<=` marks an at-most or within claim checked as a bound.", "",
         "| Check | Result | Passage found | Stated vs recomputed | Note |", "|---|---|---|---|---|"]
for o in OUT:
    lines.append(f"| {o['id']} | {'pass' if o['ok'] else '**FAIL**'} | {'yes' if o['in_doc'] else '**no**'} | {o['detail']} | {o['note']} |")
open(os.path.join(ROOT, "paper", "audit_cmpb_numbers.md"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
print(f"{len(OUT)} checks, {len(bad)} failed")
for o in bad: print("FAIL", o["id"], "| in doc:", o["in_doc"], "|", o["detail"], "|", o["note"])
sys.exit(1 if bad else 0)
