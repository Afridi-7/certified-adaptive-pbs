"""
Benchmark table across conditions, from stored results only (RQ3, RQ4).

For every design and condition: validity (certification rate of p > theta at p = theta, alpha = 0.05) and
cost (median cells to certify p = 0.15 vs theta = 0.10; uncertified specimens count all cells).
Conditions: in-distribution, clustered, confusable background, unequal fields, RV-PBS, MLL23,
real arrangements, single reader, and the paired re-run of 49_baseline_reference.py.

Usage: python scripts/46_benchmark_table.py --out results/benchmark_v1
"""
import argparse, glob, json, os
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results"); ap.add_argument("--out", default="results/benchmark_v1")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
R = a.results


def cat(pattern):
    fs = sorted(glob.glob(os.path.join(R, pattern)))
    if not fs: raise FileNotFoundError(pattern)
    return pd.concat([pd.read_csv(f, keep_default_na=False, na_values=[""]) for f in fs], ignore_index=True)


rows = []          # one row per (condition, design)


def add(cond, design, rate10=None, rate20=None, n10=None, n20=None, pooled=None, max_spec=None, n_runs=None,
        cost=None, n_cost=None, note=""):
    rates = [r for r in (rate10, rate20) if r is not None]
    rows.append(dict(condition=cond, design=design, rate_theta010=rate10, rate_theta020=rate20, n_theta010=n10,
                     n_theta020=n20, validity=(max(rates) if rates else pooled), pooled=pooled, max_specimen=max_spec,
                     n_runs=n_runs, cost_cells_p015_theta010=cost, n_cost=n_cost, note=note))


def rate(df, col): return float(df[col].astype(float).mean()), int(len(df))
def med(df, col): return float(df[col].astype(float).median()), int(len(df))


# --- in distribution (Part III confirmatory) ---------------------------------------------------------------
p10 = pd.read_csv(f"{R}/part3_theta010/episodes_p0.10.csv"); p20 = pd.read_csv(f"{R}/part3_theta020/episodes_p0.20.csv")
c15 = pd.read_csv(f"{R}/part3_theta010/episodes_p0.15.csv")
for d, col in [("alpha_u", "alpha_uniform"), ("bet_u", "bet_uniform"), ("bet_c", "bet_prop")]:
    r1, n1 = rate(p10, f"{col}_cert"); r2, n2 = rate(p20, f"{col}_cert"); m, nm = med(c15, f"{col}_cells")
    add("in_dist", d, r1, r2, n1, n2, cost=m, n_cost=nm)
r1, n1 = rate(p10, "naive_cert"); r2, n2 = rate(p20, "naive_cert"); add("in_dist", "uncorrected", r1, r2, n1, n2)
r1, n1 = rate(p10, "conv100_call_above"); r2, n2 = rate(p20, "conv100_call_above"); add("in_dist", "conv100", r1, r2, n1, n2)
s10 = cat("a3_spatial_theta010/ep_scattered_p0.10_e*.csv"); s20 = cat("a3_spatial_theta020/ep_scattered_p0.20_e*.csv")
r1, n1 = rate(s10, "naive_cell_cert"); r2, n2 = rate(s20, "naive_cell_cert")
add("in_dist", "cells_as_draws", r1, r2, n1, n2, note="scattered arrangement of E3.1")
k10 = cat("a3_mixture_ku/ep_ku_p0.10_e*.csv"); k15 = cat("a3_mixture_ku/ep_ku_p0.15_e*.csv")
r1, n1 = rate(k10, "mix50_cert"); m, nm = med(k15, "mix50_cells")
add("in_dist", "mix50", r1, None, n1, None, cost=m, n_cost=nm,
    note=f"E3.2 run; same run: betting-uniform {med(k15, 'uniform_cells')[0]:.0f}, classifier-driven {med(k15, 'prop_cells')[0]:.0f} cells")

# --- clustered (E3.1, hot fields, density 1.0) ------------------------------------------------------------------
h10 = cat("a3_spatial_theta010/ep_hot1.0_p0.10_e*.csv"); h20 = cat("a3_spatial_theta020/ep_hot1.0_p0.20_e*.csv")
h15 = cat("a3_spatial_theta010/ep_hot1.0_p0.15_e*.csv")
for d, col, has_cost in [("alpha_u", "alpha_field", 1), ("bet_u", "bet_field", 1), ("bet_c", "prop_field", 1),
                         ("cells_as_draws", "naive_cell", 0), ("conv100", None, 0)]:
    cc = f"{col}_cert" if col else "conv100_call_above"
    r1, n1 = rate(h10, cc); r2, n2 = rate(h20, cc)
    m, nm = med(h15, f"{col}_cells") if has_cost else (None, None)
    add("clustered", d, r1, r2, n1, n2, cost=m, n_cost=nm)

# --- confusable background (share 0.20) ----------------------------------------------------------------------------
cf = pd.read_csv(f"{R}/confusable_bg_v1/episodes.csv"); cf = cf[np.isclose(cf.share, 0.20)]
for d, col in [("alpha_u", "alpha_uniform"), ("bet_u", "bet_uniform"), ("bet_c", "bet_prop")]:
    r1, n1 = rate(cf[np.isclose(cf.theta, .10) & np.isclose(cf.p, .10)], f"{col}_cert")
    r2, n2 = rate(cf[np.isclose(cf.theta, .20) & np.isclose(cf.p, .20)], f"{col}_cert")
    m, nm = med(cf[np.isclose(cf.theta, .10) & np.isclose(cf.p, .15)], f"{col}_cells")
    add("confusable", d, r1, r2, n1, n2, cost=m, n_cost=nm)

# --- unequal field sizes (gamma = 0) ---------------------------------------------------------------------------------
uf = pd.read_csv(f"{R}/unequal_fields_v1/episodes.csv"); uf0 = uf[np.isclose(uf.gamma, 0)]
for d, col in [("alpha_u", "assort"), ("bet_u", "known_unif"), ("bet_c", "known_prop"), ("fmean", "naive")]:
    r1, n1 = rate(uf0[np.isclose(uf0.theta, .10) & np.isclose(uf0.p, .10)], f"{col}_cert_a0.05")
    r2, n2 = rate(uf0[np.isclose(uf0.theta, .20) & np.isclose(uf0.p, .20)], f"{col}_cert_a0.05")
    m, nm = med(uf0[np.isclose(uf0.theta, .10) & np.isclose(uf0.p, .15)], f"{col}_a0.05") if d != "fmean" else (None, None)
    g1 = uf[np.isclose(uf.gamma, 1) & np.isclose(uf.theta, .10) & np.isclose(uf.p, .10)]
    add("unequal", d, r1, r2, n1, n2, cost=m, n_cost=nm,
        note=f"gamma=1 (blasts enriched in small fields): {rate(g1, f'{col}_cert_a0.05')[0]:.3f} at theta=0.10")

# --- acquisition shift: RV-PBS (E3.2) and MLL23 ---------------------------------------------------------------------
for cond, f10, f15 in [("rvpbs", "a3_mixture_rvpbs/ep_rvpbs_p0.10_e*.csv", "a3_mixture_rvpbs/ep_rvpbs_p0.15_e*.csv"),
                       ("mll23", "mll23_shift_v1/episodes.csv", "mll23_shift_v1/episodes.csv")]:
    e10 = cat(f10); e15 = cat(f15)
    if cond == "mll23": e10 = e10[np.isclose(e10.p, .10)]; e15 = e15[np.isclose(e15.p, .15)]
    for d, col in [("bet_u", "uniform"), ("bet_c", "prop"), ("mix50", "mix50")]:
        r1, n1 = rate(e10, f"{col}_cert"); m, nm = med(e15, f"{col}_cells")
        add(cond, d, r1, None, n1, None, cost=m, n_cost=nm)

# --- real arrangements (BoneMarrowWSI, theta_s = p_s) ----------------------------------------------------------------
ra = pd.read_csv(f"{R}/real_arrangements_v1/episodes.csv", keep_default_na=False); ra = ra[ra.scen == "null"]
for d, col in [("alpha_u", "field_assort"), ("bet_u", "field_known"), ("cells_as_draws", "naive_cell"), ("fmean", "naive_fmean")]:
    per = ra.groupby("patient")[f"{col}_cert"].mean()
    add("real", d, pooled=float(ra[f"{col}_cert"].mean()), max_spec=float(per.max()), n_runs=int(len(ra)),
        note=f"{per.size} specimens, {int((per > 0.05).sum())} above 0.05")

# --- single reader as verifier -----------------------------------------------------------------------------------------
ve = pd.read_csv(f"{R}/verifier_error_v1/episodes.csv")
for col, note in [("above_cert", "reader-0 labels, threshold theta_s"), ("above_m0.03_cert", "margin 0.03")]:
    per = ve.groupby("patient")[col].mean()
    add("reader" if col == "above_cert" else "reader_margin003", "bet_u", pooled=float(ve[col].mean()),
        max_spec=float(per.max()), n_runs=int(len(ve)), note=f"{note}; {per.size} specimens, {int((per > 0.05).sum())} above 0.05")

# --- exploratory rerun on the same specimens (scripts/49) ------------------------------------------------------------
br = {src: pd.read_csv(f"{R}/baseline_reference_v1/episodes_{src}.csv") for src in ("ku", "mll23", "rvpbs")}
def bsel(src, th, p): e = br[src]; return e[np.isclose(e.theta, th) & np.isclose(e.p, p)]
for cond, src in [("rvpbs", "rvpbs"), ("mll23", "mll23")]:
    r1, n1 = rate(bsel(src, .10, .10), "alpha_uniform_cert"); m, nm = med(bsel(src, .10, .15), "alpha_uniform_cells")
    add(cond, "alpha_u", r1, None, n1, None, cost=m, n_cost=nm, note="rerun (baseline_reference_v1)")
for cond, src in [("in_dist", "ku"), ("rvpbs", "rvpbs"), ("mll23", "mll23")]:
    r1, n1 = rate(bsel(src, .10, .10), "hedge_mix50_cert")
    r2, n2 = rate(bsel(src, .20, .20), "hedge_mix50_cert") if src == "ku" else (None, None)
    m, nm = med(bsel(src, .10, .15), "hedge_mix50_cells")
    add(cond, "hedge_mix50", r1, r2, n1, n2, cost=m, n_cost=nm, note="rerun (baseline_reference_v1)")
RERUN = {("rvpbs", "alpha_u"), ("mll23", "alpha_u")}                    # cells marked $^c$ in otherwise non-rerun rows

T = pd.DataFrame(rows); T.to_csv(f"{a.out}/benchmark.csv", index=False)

# --- LaTeX body ---------------------------------------------------------------------------------------------------------
COLS = ["in_dist", "clustered", "confusable", "unequal", "rvpbs", "mll23", "real", "reader"]
VALID = [("alpha_u", "ALPHA, uniform$^{a}$"), ("bet_u", "Betting, uniform$^{a}$"), ("bet_c", "Betting, classifier-driven"),
         ("mix50", "Betting, 50\\% mixture"), ("hedge_mix50", "Hedged betting, 50\\% mixture$^c$")]
INVALID = [("uncorrected", "Ranked, analysed as uniform"), ("cells_as_draws", "Cells of fields as draws"),
           ("fmean", "Mean of field proportions$^{b}$"), ("conv100", "100-cell count, ``above''")]


def get(cond, d):
    s = T[(T.condition == cond) & (T.design == d)]
    return None if s.empty else s.iloc[0]


def fr(x):
    if x is None or (isinstance(x, float) and np.isnan(x)): return "---"
    if x == 0: return "0"
    return f"{x:.3f}" if x < 0.995 else f"{x:.2f}"


def vcell(cond, d):
    s = get(cond, d)
    if s is None: return "---"
    if cond in ("real", "reader"): return f"{fr(s.pooled)} [{fr(s.max_specimen)}]"
    return fr(s.validity) + ("$^c$" if (cond, d) in RERUN else "")


def ccell(cond, d):
    s = get(cond, d)
    if s is None or pd.isna(s.cost_cells_p015_theta010): return "---"
    v = int(np.floor(s.cost_cells_p015_theta010 + 0.5))                     # half up, as in the paper text
    return f"{v:,d}".replace(",", "{,}") + ("$^c$" if (cond, d) in RERUN else "")


lines = ["\\multicolumn{9}{@{}l}{\\emph{False certification of ``$p>\\theta$'' at $p=\\theta$}}\\\\"]
for d, lab in VALID: lines.append(lab + " & " + " & ".join(vcell(c, d) for c in COLS) + "\\\\")
lines.append("\\addlinespace")
for d, lab in INVALID: lines.append(lab + " & " + " & ".join(vcell(c, d) for c in COLS) + "\\\\")
lines.append("\\midrule")
lines.append("\\multicolumn{9}{@{}l}{\\emph{Median cells to certify $p=0.15$ against $\\theta=0.10$}}\\\\")
for d, lab in VALID: lines.append(lab + " & " + " & ".join(ccell(c, d) for c in COLS) + "\\\\")
open(f"{a.out}/benchmark_rows.tex", "w").write("\n".join(lines) + "\n")

summary = dict(rows=T.where(T.notna(), None).to_dict(orient="records"))
json.dump(summary, open(f"{a.out}/summary.json", "w"), indent=1, default=float)
json.dump(vars(a), open(f"{a.out}/config.json", "w"), indent=1)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
print(T[["condition", "design", "rate_theta010", "rate_theta020", "validity", "pooled", "max_specimen",
         "cost_cells_p015_theta010", "n_cost", "note"]].to_string(index=False))
print("\n".join(lines))
