"""
Paper figures from the FROZEN confirmatory outputs (no new analysis; every number is read from results/).
Two exploratory overlays, clearly labelled, were added on 2026-10-22 at the supervisor's request: the field-level teacher
in Fig. 1(b) (results/kt1_field_teacher) and the bias-corrected clustering estimate in Fig. 4(b)
(results/a4_bmwsi/icc_excess.json), with the raw estimate kept as a comparison.
Usage: python scripts/25_figures.py  -> figures/fig1_partI.pdf/.png ... fig4_rq4.pdf/.png + figures/captions.md
Colour roles (validated categorical slots, light mode): blue = uniform / field-level, orange = adaptive
(InfoDPP-PAC, classifier-driven) or naive, aqua = mixture; grey = reference lines. Aqua is always direct-labelled.
"""
import glob, json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import beta as Beta

BLUE, ORANGE, AQUA, GREY, INK, INK2 = "#2a78d6", "#eb6834", "#1baf7a", "#8a8985", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 7, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0",
                     "grid.linewidth": 0.6, "lines.linewidth": 2, "lines.markersize": 5, "savefig.dpi": 300,
                     "font.family": "DejaVu Sans", "pdf.fonttype": 42})
os.makedirs("figures", exist_ok=True); CAP = []


def wilson(x, n, z=1.96):
    p = x / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, c - h, c + h


def cp(x, n):
    lo = Beta.ppf(0.025, x, n - x + 1) if x > 0 else 0.0; hi = Beta.ppf(0.975, x + 1, n - x) if x < n else 1.0
    return x / n, lo, hi


def save(fig, name, caption):
    fig.savefig(f"figures/{name}.pdf", bbox_inches="tight"); fig.savefig(f"figures/{name}.png", bbox_inches="tight")
    plt.close(fig); CAP.append(f"**{name}.** {caption}\n")


# ---------------- Figure 1: Part I (RQ1) ----------------
E = pd.read_csv("results/confirm_v1/episodes.csv")
fig, ax = plt.subplots(1, 2, figsize=(7, 2.6), constrained_layout=True)
ms = [1, 5, 10, 20]; xs = np.arange(4)
for col, lab, c, mk, dx in [("infodpp_miss", "InfoDPP-PAC", ORANGE, "o", -0.08), ("uniform@info_miss", "uniform, matched cost", BLUE, "s", 0.08)]:
    pts = [wilson(E[(E.m == m) & (E.k == 1)][col].sum(), len(E[(E.m == m) & (E.k == 1)])) for m in ms]
    p = np.array([q[0] for q in pts]); lo = np.array([q[1] for q in pts]); hi = np.array([q[2] for q in pts])
    ax[0].errorbar(xs + dx, p, yerr=[p - lo, hi - p], color=c, marker=mk, capsize=0, elinewidth=1, label=lab)
nmed = E[E.k == 1].infodpp_n.median(); N = int(E.N.iloc[0])
ax[0].axhline(1 - nmed / N, color=GREY, ls="--", lw=1)
ax[0].text(-0.1, 1 - nmed / N + 0.035, f"worst case, any strategy: 1 − n/N = {1 - nmed / N:.2f}", color=INK2, ha="left", fontsize=6.5)
ax[0].set_xticks(xs, [str(m) for m in ms]); ax[0].set_xlabel("cells per field (dilution m)"); ax[0].set_ylim(0, 1)
ax[0].set_ylabel("miss rate, one critical field (k = 1)"); ax[0].set_title("(a) Premature stops at matched cost", loc="left")
ax[0].legend(frameon=False, loc="lower right")
for col, lab, c, mk in [("kt1_auc_info_step1", "information gain (stopping signal)", ORANGE, "o"),
                        ("kt1_auc_relevance", "relevance, cell-level teacher", BLUE, "s")]:
    mu, lo, hi = [], [], []
    for m in ms:
        v = E[E.m == m][col].dropna().values; rng = np.random.default_rng(0)
        bs = [rng.choice(v, len(v)).mean() for _ in range(1000)]; mu.append(v.mean()); lo.append(np.quantile(bs, .025)); hi.append(np.quantile(bs, .975))
    mu, lo, hi = map(np.array, (mu, lo, hi))
    ax[1].errorbar(xs, mu, yerr=[mu - lo, hi - mu], color=c, marker=mk, capsize=0, elinewidth=1, label=lab)
FT = pd.read_csv("results/kt1_field_teacher/episodes.csv")                 # exploratory: same pooled embeddings as the selector
mu, lo, hi = [], [], []
for m in ms:
    v = FT[FT.m == m]["auc_rel_field"].dropna().values; rng = np.random.default_rng(0)
    bs = [rng.choice(v, len(v)).mean() for _ in range(1000)]; mu.append(v.mean()); lo.append(np.quantile(bs, .025)); hi.append(np.quantile(bs, .975))
mu, lo, hi = map(np.array, (mu, lo, hi))
ax[1].errorbar(xs + 0.06, mu, yerr=[mu - lo, hi - mu], color=AQUA, marker="^", ls="--", capsize=0, elinewidth=1,
               label="relevance, field-level teacher\n(exploratory, 10,000 episodes)")
ax[1].axhline(0.5, color=GREY, ls="--", lw=1); ax[1].text(3.0, 0.508, "chance", color=INK2, fontsize=6.5, ha="right")
ax[1].set_xticks(xs, [str(m) for m in ms]); ax[1].set_xlabel("cells per field (dilution m)"); ax[1].set_ylim(0.45, 1.0)
ax[1].set_ylabel("AUC for identifying critical fields"); ax[1].set_title("(b) Kill-test 1: stopping signal vs relevance", loc="left")
ax[1].legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.33), fontsize=6.3)
save(fig, "fig1_partI",
     f"Part I (RQ1), confirmatory run (24,000 episodes; N = {N} fields; 5 patient-grouped folds). (a) Probability of stopping "
     f"without examining the single critical field, InfoDPP-PAC vs uniform sampling at the same per-episode cost "
     f"(median {int(nmed)} fields); 95% Wilson intervals. Dashed: worst-case lower bound 1 − n/N for any strategy without "
     "location information (Theorem 3). (b) AUC with which the information-gain term that drives stopping, and the "
     "cell-level teacher's relevance score, rank critical above non-critical fields; mean over episodes with 95% bootstrap "
     "intervals. Dashed (exploratory, separate 10,000 episodes): the field-level teacher, trained on the pooled embeddings the "
     "selector sees.")

# ---------------- Figure 2: Part II (RQ2) ----------------
T = pd.read_csv("results/part2_v1/nstar_table.csv")
fig, ax = plt.subplots(figsize=(3.4, 2.6), constrained_layout=True)
for al, c, mk in [(0.05, BLUE, "o"), (0.01, ORANGE, "s")]:
    g = T[(T.N == 1000) & np.isclose(T.alpha, al)].sort_values("k")
    ax.plot(g.k, g.n_star / 1000, color=c, marker=mk, label=f"exact, α = {al}")
    kk = np.logspace(0, np.log10(50), 100); ax.plot(kk, 1 - al ** (1 / kk), color=c, lw=1, ls=":")
ax.set_xscale("log"); ax.set_xticks([1, 2, 5, 10, 20, 50], ["1", "2", "5", "10", "20", "50"]); ax.set_ylim(0, 1)
ax.set_xlabel("critical fields k"); ax.set_ylabel("fraction of fields needed, n*/N")
ax.set_title("Certifying absence (N = 1000)", loc="left"); ax.legend(frameon=False, loc="upper right")
save(fig, "fig2_partII",
     "Part II (RQ2). Smallest number of fields n* that any strategy must examine so that the probability of missing all k "
     "critical fields is at most α, for N = 1000 (exact hypergeometric; solid) and the sufficient closed form "
     "1 − α^(1/k) (dotted). By Theorem 3 this is a lower bound for every strategy without location information and is "
     "attained by uniform sampling. Values are nearly independent of N (results/part2_v1/nstar_table.csv).")

# ---------------- Figure 3: Part III (RQ3) ----------------
P = pd.concat([pd.read_csv(f"results/part3_theta0{t}0/episodes.csv") for t in (1, 2)])
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7), gridspec_kw=dict(width_ratios=[1.25, 1, 1]), constrained_layout=True)
designs = [("alpha_uniform", "ALPHA\nuniform"), ("bet_uniform", "betting\nuniform"), ("bet_prop", "betting\nclassifier"),
           ("bet_prop_cv", "betting\n+ CV"), ("naive", "uncorr.\nclassifier")]
for th, c, mk, dx in [(0.10, BLUE, "o", -0.12), (0.20, ORANGE, "s", 0.12)]:
    g = P[np.isclose(P.theta, th) & np.isclose(P.p, th)]
    for i, (d, _) in enumerate(designs):
        r, lo, hi = cp(int(g[f"{d}_cert"].sum()), len(g)); r = max(r, 5e-4); lo = max(lo, 5e-4)
        ax[0].errorbar(i + dx, r, yerr=[[r - lo], [hi - r]], color=c, marker=mk, elinewidth=1, capsize=0,
                       label=(f"θ = {th:.2f} ({'ICC' if th == 0.10 else 'WHO'} 2022)" if i == 0 else None))
ax[0].axhline(0.05, color=GREY, ls="--", lw=1); ax[0].text(-0.4, 0.058, "α = 0.05", color=INK2, fontsize=6.5, ha="left")
ax[0].set_yscale("log"); ax[0].set_ylim(4e-4, 1.6); ax[0].set_xlim(-0.5, 4.5); ax[0].set_xticks(range(5), [l for _, l in designs], fontsize=6)
ax[0].set_ylabel("false certification at p = θ"); ax[0].set_title("(a) Validity at p = θ", loc="left")
ax[0].legend(frameon=False, loc="upper left", fontsize=6.5)
for j, th in enumerate((0.10, 0.20)):
    a = ax[j + 1]; g = P[np.isclose(P.theta, th) & (P.p > th + 1e-9)]
    for d, lab, c, mk in [("alpha_uniform", "ALPHA, uniform", BLUE, "o"), ("bet_prop", "betting, classifier-driven", ORANGE, "s")]:
        s = g.groupby("p")[f"{d}_cells"]; med, q1, q3 = s.median(), s.quantile(.25), s.quantile(.75)
        a.fill_between(med.index, q1, q3, color=c, alpha=0.15, lw=0); a.plot(med.index, med, color=c, marker=mk, label=lab)
    a.set_ylim(0, 520); a.axhline(500, color=GREY, lw=0.8, ls=":"); a.text(med.index.min(), 505, "census (500 cells)", color=INK2, fontsize=6)
    a.set_xlabel("true blast proportion p"); a.set_ylabel("cells examined to certify" if j == 0 else "")
    a.set_title(f"({'bc'[j]}) Cost, θ = {th:.2f} ({'ICC' if th == 0.10 else 'WHO'})", loc="left")
    if j == 0: a.legend(frameon=False, loc="upper right", bbox_to_anchor=(1, 0.92), fontsize=6.5)
save(fig, "fig3_partIII",
     "Part III (RQ3), confirmatory run (12,000 specimens of 500 KU-Optofil cells; α = 0.05; all designs paired on the same "
     "specimen). (a) Rate of certifying 'p > θ' when p = θ exactly (the hardest true null), with 95% Clopper–Pearson "
     "intervals, log scale; CV = control variate; the uncorrected design ranks fields by classifier score but analyses them as if uniformly "
     "sampled. (b, c) Median cells examined until certification (line) and interquartile range (band); specimens that "
     "never certify count as 500. Thresholds from ICC 2022 (0.10) and WHO 5th ed. 2022 (0.20) are analysed separately.")

# ---------------- Figure 4: RQ4 ----------------
X = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("results/a3_mixture_*/ep_*.csv"))])
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=True)
groups = [("ku", 0.15, "KU-Optofil\np = 0.15"), ("ku", 0.20, "KU-Optofil\np = 0.20"), ("rvpbs", 0.15, "RV-PBS (shift)\np = 0.15"), ("rvpbs", 0.20, "RV-PBS (shift)\np = 0.20")]
for d, lab, c, mk, dx in [("uniform", "uniform", BLUE, "o", -0.18), ("prop", "classifier-driven", ORANGE, "s", 0.0), ("mix50", "50/50 mixture", AQUA, "D", 0.18)]:
    med, lo, hi = [], [], []
    for src, p, _ in groups:
        v = X[(X.source == src) & np.isclose(X.p, p)][f"{d}_cells"].values; rng = np.random.default_rng(1)
        bs = [np.median(rng.choice(v, len(v))) for _ in range(1000)]; med.append(np.median(v)); lo.append(np.quantile(bs, .025)); hi.append(np.quantile(bs, .975))
    med, lo, hi = map(np.array, (med, lo, hi))
    ax[0].errorbar(np.arange(4) + dx, med, yerr=[med - lo, hi - med], color=c, marker=mk, ls="none", elinewidth=1, capsize=0, label=lab)
ax[0].set_xticks(range(4), [g[2] for g in groups], fontsize=6.3); ax[0].set_ylabel("median cells to certify (θ = 0.10)")
ax[0].set_ylim(0, 345); ax[0].axvline(1.5, color=GREY, lw=0.8); ax[0].legend(frameon=False, loc="upper center", ncol=3, fontsize=6.5, handletextpad=0.3, columnspacing=1.0)
ax[0].set_title("(a) Cost in distribution vs under shift", loc="left")
S = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("results/a3_spatial_theta0*/ep_*.csv"))])
br = json.load(open("results/a4_bmwsi/bridge.json")); cl = json.load(open("results/a4_bmwsi/clustering.json"))["primary"]["S1"]
place = br["e31_icc_placement"]
for th, c, mk in [(0.10, BLUE, "o"), (0.20, ORANGE, "s")]:
    sims = sorted([v for v in br["sims"].values() if np.isclose(v["theta"], th)], key=lambda v: v["rho"])
    ax[1].plot([max(v["rho"], 1e-3) for v in sims], [v["naive_cell_false_cert"] for v in sims], color=c, lw=1.5, marker=mk, ms=3.5)
    for arr in ["scattered", "hot0.6", "hot1.0"]:
        g = S[np.isclose(S.theta, th) & (S.arr == arr) & np.isclose(S.p, th)]
        r, lo, hi = cp(int(g.naive_cell_cert.sum()), len(g)); x = max(place[f"theta{th}_{arr}"], 1e-3)
        ax[1].errorbar(x, r, yerr=[[r - lo], [hi - r]], color=c, marker=mk, mfc="white", mew=1.5, elinewidth=1, capsize=0, ls="none")
    ax[1].plot([], [], color=c, marker=mk, label=f"θ = {th:.2f}")
ex = json.load(open("results/a4_bmwsi/icc_excess.json"))["S1"]             # bias-corrected (excess over permutation null)
ax[1].axvspan(max(cl["ci95"][0], 1e-3), cl["ci95"][1], facecolor="none", edgecolor=GREY, hatch="////", lw=0, alpha=0.5)
ax[1].axvspan(max(ex["excess_ci95"][0], 1e-3), ex["excess_ci95"][1], color=INK2, alpha=0.15, lw=0)
ax[1].axvline(ex["excess"], color=INK, lw=1.2)
ax[1].text(ex["excess"] * 0.92, 0.30, "real smears:\nbias-corrected ICC\n(95% CI)", color=INK, fontsize=6.3, ha="right")
ax[1].text(cl["ci95"][1] * 1.08, 0.13, "raw ICC\n(biased upward)", color=INK2, fontsize=6.0, ha="left")
ax[1].axhline(0.05, color=GREY, ls="--", lw=1); ax[1].text(0.9, 0.056, "α", color=INK2, fontsize=6.5)
ax[1].set_xscale("log"); ax[1].set_xlim(8e-4, 1.5); ax[1].set_ylim(0, 0.55)
ax[1].set_xlabel("within-field clustering of blasts (ICC, log scale)"); ax[1].set_ylabel("false certification at p = θ")
ax[1].set_title("(b) Cell-level analysis of field samples", loc="left"); from matplotlib.lines import Line2D
ax[1].plot([], [], color=GREY, lw=1.5, marker="o", ms=3.5, label="Beta-binomial model, 17 cells/field")
ax[1].plot([], [], color=GREY, ls="none", marker="o", mfc="white", mew=1.5, label="E3.1 composites (observed)")
ax[1].legend(frameon=False, loc="upper left", fontsize=6.3)
save(fig, "fig4_rq4",
     "RQ4. (a) Median cells to certify p > 0.10 (95% bootstrap intervals; 500 specimens per point) when the classifier is "
     "in distribution (KU-Optofil, cross-fitted) and under acquisition shift (RV-PBS, KU-trained probe). All designs keep "
     "false certification ≤ 0.014 at p = θ. (b) False certification at p = θ when cells from uniformly sampled fields are "
     "analysed as individual draws. Open markers: E3.1 composites (5 cells per field; scattered, hot0.6, hot1.0; 1,000 "
     "specimens each, 95% Clopper–Pearson). Lines: Beta-binomial model at the clustering measured in real bone-marrow "
     "monolayer ROIs (17 cells per field). Solid band and line: bias-corrected ICC of the real regions (excess over the "
     "permutation null, patient-bootstrap 95% CI; exploratory). Hatched: 95% CI of the raw estimate, shown for comparison "
     "(biased upward). Field-level analysis stayed ≤ 0.042 in every condition (not shown). ICC is plotted at 10⁻³ when 0.")

open("figures/captions.md", "w").write("# Figure captions (generated by scripts/25_figures.py)\n\n" + "\n".join(CAP))
print("written:", sorted(os.listdir("figures")))
