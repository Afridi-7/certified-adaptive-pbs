"""
Encoder robustness check with DinoBloom-S (exploratory; RQ1, RQ3, RQ4). Review of 2026-11-19.

One alternative encoder, not an encoder comparison: the question is whether the paper's conclusions about costs and
about InfoDPP-PAC's stopping survive a hematology-specific encoder. Validity of the certification tests does not
depend on the encoder (it rests on the sampling design), so only efficiency and the selector can change.

The DinoBloom-S embeddings come from notebooks/57_dinobloom_embeddings_colab.ipynb: the same cells, in the same row
order, with the same preprocessing as the DINOv2 embeddings; only the encoder differs. Every part below runs the
SAME code for both encoders, so differences are attributable to the encoder.

Parts
  classifier : patient-grouped 5-fold probe on KU-Optofil (out-of-fold AUC, sensitivity at 95% specificity); the probe
               trained on all KU-Optofil cells applied to MLL23 and RV-PBS (AUC, sensitivity and false-positive rate at
               the KU-Optofil 95%-specificity threshold).
  part3      : the paired design comparison of scripts/49_baseline_reference.py (ALPHA uniform, approximate-Kelly
               betting uniform / prop-MS, hedged betting uniform / prop-MS / 50% mixture) on N = 100 fields x 5 cells.
               Default seed family 20261101 = that of script 49, so the specimens (cells drawn) are identical to the
               stored DINOv2 run; with --encoder dinov2 this part reproduces results/baseline_reference_v1 exactly
               (checked by --check_replication).
  part1      : InfoDPP-PAC against uniform sampling at matched cost, one critical field (k = 1), N = 200 fields,
               m in {1, 5, 10, 20}, the settings of scripts/05_killtests.py (tau calibrated per fold and m on blast-free
               specimens of other folds to stop at 16.3% of fields). New seed family 20261119, run for both encoders.

Usage
  python scripts/57_dinobloom_check.py --part classifier --out results/dinobloom_check_v1
  python scripts/57_dinobloom_check.py --part part3 --encoder dinobloom --source ku --out results/dinobloom_check_v1
  python scripts/57_dinobloom_check.py --part part3 --encoder dinov2 --source ku --out /tmp/rep --check_replication
  python scripts/57_dinobloom_check.py --part part1 --encoder dinobloom --out results/dinobloom_check_v1
"""
import argparse, json, os, sys, time, platform
import numpy as np, pandas as pd
from scipy.stats import binomtest, beta, wilcoxon
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from cadpbs.alpha import alpha_martingale
from cadpbs.weighted_betting_w import betting_test_w, sampler_uniform_w, sampler_prop_w
from cadpbs.hedged_betting import hedged_betting_w

SOURCES = ["ku", "mll23", "rvpbs"]
DESIGNS = ["alpha_uniform", "bet_uniform", "bet_prop", "hedge_uniform", "hedge_prop", "hedge_mix50"]
PATHS = {
    "dinov2": dict(ku_emb="results/embeddings_v1/embeddings_dinov2_vits14.npy", ku_csv="results/embeddings_v1/cells_features.csv",
                   mll_emb="data/mll23/mll23_embeddings_dinov2_vits14.npy", mll_csv="data/mll23/mll23_cells.csv",
                   rv_emb="results/external_v1/rvpbs_embeddings_dinov2_vits14.npy", rv_csv="results/external_v1/rvpbs_features.csv"),
    "dinobloom": dict(ku_emb="results/dinobloom_v1/ku_embeddings_dinobloom_s.npy", ku_csv="results/dinobloom_v1/ku_cells_features.csv",
                      mll_emb="results/dinobloom_v1/mll23_embeddings_dinobloom_s.npy", mll_csv="results/dinobloom_v1/mll23_cells.csv",
                      rv_emb="results/dinobloom_v1/rvpbs_embeddings_dinobloom_s.npy", rv_csv="results/dinobloom_v1/rvpbs_cells.csv"),
}

ap = argparse.ArgumentParser()
ap.add_argument("--part", choices=["classifier", "part3", "part1"], required=True)
ap.add_argument("--encoder", choices=list(PATHS), default="dinobloom")
ap.add_argument("--source", choices=SOURCES, default="ku")
ap.add_argument("--out", required=True)
ap.add_argument("--seed", type=int, default=None, help="part3: 20261101 (script 49); part1: 20261119")
ap.add_argument("--n_null", type=int, default=1000); ap.add_argument("--n_alt", type=int, default=400)
ap.add_argument("--n_boot", type=int, default=4000)
ap.add_argument("--episodes", type=int, default=100, help="part1: episodes per fold and m")
ap.add_argument("--m_list", type=int, nargs="+", default=[1, 5, 10, 20])
ap.add_argument("--check_replication", action="store_true", help="part3 with --encoder dinov2: compare with results/baseline_reference_v1")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); T0 = time.time()
P = {k: os.path.join(ROOT, v) for k, v in PATHS[a.encoder].items()}


def load_ku():
    df = pd.read_csv(P["ku_csv"], dtype={"patient_id": str})
    E = np.load(P["ku_emb"]).astype(np.float32)
    assert len(df) == len(E), "embeddings and cell table differ in length"
    return E, df


def probe():
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=4000, C=0.1, class_weight="balanced"))


def threshold_95spec(df):   # as notebook 03 cell 7: 95th percentile of out-of-fold scores of non-blasts
    return float(np.quantile(df.probe_p_oof.values[df.is_blast.values == 0], 0.95))


def load_shift(source, Eku, ku):
    if source == "mll23":
        mm = pd.read_csv(P["mll_csv"]); Em = np.load(P["mll_emb"]).astype(np.float32)
        keep = ~mm.cls.str.contains("atypical|smudge|normoblast", case=False, regex=True)
        Em, mm = Em[keep.values], mm[keep].reset_index(drop=True)
        y = mm.cls.str.contains("myeloblast", case=False).values
        return y, probe().fit(Eku, ku.is_blast.values).predict_proba(Em)[:, 1]
    rv = pd.read_csv(P["rv_csv"]); Erv = np.load(P["rv_emb"]).astype(np.float32)
    return rv.is_blast.values.astype(bool), probe().fit(Eku, ku.is_blast.values).predict_proba(Erv)[:, 1]


def env():
    import numpy, scipy, sklearn
    return dict(python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__,
                pandas=pd.__version__, sklearn=sklearn.__version__)


# ----------------------------------------------------------------------------------------------------------------------
if a.part == "classifier":
    from sklearn.metrics import roc_auc_score
    out = dict(encoder={}, env=env())
    for enc in PATHS:
        P = {k: os.path.join(ROOT, v) for k, v in PATHS[enc].items()}
        if not os.path.exists(P["ku_emb"]): print("missing", P["ku_emb"], "- skipped"); continue
        Eku, ku = load_ku(); y = ku.is_blast.values.astype(bool); s = ku.probe_p_oof.values; T = threshold_95spec(ku)
        r = dict(ku_oof_auc=float(roc_auc_score(y, s)), ku_threshold_95spec=T, ku_sensitivity_95spec=float((s[y] >= T).mean()))
        for src in ["mll23", "rvpbs"]:
            if not os.path.exists(P["mll_emb" if src == "mll23" else "rv_emb"]): continue
            ys, ss = load_shift(src, Eku, ku)
            r[f"{src}_auc"] = float(roc_auc_score(ys, ss)); r[f"{src}_sensitivity_at_ku_threshold"] = float((ss[ys] >= T).mean())
            r[f"{src}_fp_at_ku_threshold"] = float((ss[~ys] >= T).mean())
        out["encoder"][enc] = r; print(enc, json.dumps(r, indent=1), flush=True)
    json.dump(out, open(f"{a.out}/classifier_summary.json", "w"), indent=1)
    sys.exit(0)

# ----------------------------------------------------------------------------------------------------------------------
if a.part == "part3":
    seed = 20261101 if a.seed is None else a.seed
    Eku, ku = load_ku(); T = threshold_95spec(ku)
    if a.source == "ku":
        y = ku.is_blast.values.astype(bool); score = ku.probe_p_oof.values
        bp = set(ku.loc[y, "patient_id"]); BL = np.flatnonzero(y); BG = np.flatnonzero(~y & ~ku.patient_id.isin(bp).values)
        grid = [(0.10, [0.10, 0.12, 0.15, 0.20]), (0.20, [0.20, 0.22, 0.25, 0.30])]
    else:
        y, score = load_shift(a.source, Eku, ku); BL = np.flatnonzero(y); BG = np.flatnonzero(~y)
        grid = [(0.10, [0.10, 0.12, 0.15, 0.20])]
    recall = float((score[y] >= T).mean()); fp = float((score[~y] >= T).mean())
    print(f"{a.encoder}/{a.source}: {len(BL)} blasts, {len(BG)} background, threshold {T:.4f}, recall {recall:.3f}, FP {fp:.3f}", flush=True)
    N, m = 100, 5; M = N * m; pi = np.full(N, 1.0 / N); si = SOURCES.index(a.source)

    def mix(S, lam):
        base = sampler_prop_w(pi, S)
        return lambda rem, t: (1 - lam) * base(rem, t) + lam / len(rem)

    rows = []
    for ti, (th, plist) in enumerate(grid):
        for p in plist:
            nb = int(round(p * M)); assert nb <= len(BL) and M - nb <= len(BG)
            n_ep = a.n_null if abs(p - th) < 1e-9 else a.n_alt
            for ep in range(n_ep):
                rng = np.random.default_rng([seed, si, ti, int(round(p * 1000)), ep])
                cells = np.r_[rng.choice(BL, nb, replace=False), rng.choice(BG, M - nb, replace=False)]
                rng.shuffle(cells); cells = cells.reshape(N, m)
                f = y[cells].mean(1).astype(float); S = score[cells].mean(1)
                r = dict(source=a.source, theta=th, p=p, ep=ep)
                _, st = alpha_martingale(f[rng.permutation(N)] / (2 * th), 1 / (2 * th), 0.5, N, alpha=0.05)
                r["alpha_uniform_stop"] = st
                r["bet_uniform_stop"] = betting_test_w(f, pi, S, sampler_uniform_w(N), th, rng=rng)["stop"]
                r["bet_prop_stop"] = betting_test_w(f, pi, S, sampler_prop_w(pi, S), th, rng=rng)["stop"]
                r["hedge_uniform_stop"] = hedged_betting_w(f, pi, sampler_uniform_w(N), th, rng=rng)["stop"]
                r["hedge_prop_stop"] = hedged_betting_w(f, pi, sampler_prop_w(pi, S), th, rng=rng)["stop"]
                r["hedge_mix50_stop"] = hedged_betting_w(f, pi, mix(S, 0.5), th, rng=rng)["stop"]
                for d in DESIGNS:
                    s = r.pop(f"{d}_stop"); r[f"{d}_cert"] = int(s > 0); r[f"{d}_cells"] = (s if s > 0 else N) * m
                rows.append(r)
            print(f"  theta={th} p={p}: {n_ep} specimens ({time.time() - T0:.0f}s)", flush=True)
    D = pd.DataFrame(rows); tag = f"{a.encoder}_{a.source}"; D.to_csv(f"{a.out}/part3_episodes_{tag}.csv", index=False)
    rb = np.random.default_rng([seed, si, 999])
    out = dict(encoder=a.encoder, source=a.source, threshold_95spec=T, recall_at_threshold=recall, fp_at_threshold=fp,
               n_blast_pool=int(len(BL)), n_background_pool=int(len(BG)), validity={}, median_cells={}, ratios={})
    for (th, p), g in D.groupby(["theta", "p"]):
        key = f"theta{th:.2f}_p{p:.2f}"
        if abs(p - th) < 1e-9:
            for d in DESIGNS:
                k = int(g[f"{d}_cert"].sum()); n = len(g)
                out["validity"][f"{key}_{d}"] = dict(rate=k / n, cp95_upper=float(beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0,
                                                     p_excess=float(binomtest(k, n, 0.05, alternative="greater").pvalue))
            continue
        out["median_cells"][key] = {d: float(g[f"{d}_cells"].median()) for d in DESIGNS}
        for num, den in [("bet_prop", "alpha_uniform"), ("bet_prop", "bet_uniform"), ("bet_uniform", "alpha_uniform"),
                         ("hedge_uniform", "alpha_uniform"), ("hedge_prop", "alpha_uniform"), ("hedge_prop", "hedge_uniform"),
                         ("hedge_mix50", "alpha_uniform")]:
            x, z = g[f"{num}_cells"].values, g[f"{den}_cells"].values; bs = []
            for _ in range(a.n_boot):
                i = rb.integers(0, len(g), len(g)); bs.append(np.median(x[i]) / np.median(z[i]))
            try: wp = float(wilcoxon(x, z).pvalue)
            except ValueError: wp = 1.0
            out["ratios"][f"{key}_{num}_over_{den}"] = dict(ratio=float(np.median(x) / np.median(z)),
                                                            ci95=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))], wilcoxon_p=wp)
    json.dump(dict(config=vars(a) | dict(seed_used=seed), env=env()), open(f"{a.out}/part3_config_{tag}.json", "w"), indent=1)
    json.dump(out, open(f"{a.out}/part3_summary_{tag}.json", "w"), indent=1)
    if a.check_replication:
        assert a.encoder == "dinov2" and seed == 20261101, "replication check needs --encoder dinov2 and the script-49 seed"
        ref = pd.read_csv(os.path.join(ROOT, f"results/baseline_reference_v1/episodes_{a.source}.csv"))
        cols = [c for c in ref.columns if c.endswith("_cells") or c.endswith("_cert")]
        same = (len(ref) == len(D)) and all((ref[c].values == D[c].values).all() for c in cols)
        print("REPLICATION of scripts/49 (episode by episode, all designs):", "IDENTICAL" if same else "DIFFERENT", flush=True)
        if not same: sys.exit(2)
    print(f"DONE part3 {tag} in {time.time() - T0:.0f}s", flush=True)
    sys.exit(0)

# ----------------------------------------------------------------------------------------------------------------------
if a.part == "part1":
    from cadpbs.composites import CellBank, make_specimen
    from cadpbs.infodpp import InfoDPPConfig, run_infodpp, calibrate_tau
    seed = 20261119 if a.seed is None else a.seed
    Eku, ku = load_ku(); bank = CellBank(Eku, ku)
    N, k, target_frac, n_calib = 200, 1, 0.163, 8          # settings of scripts/05_killtests.py
    base_cfg = InfoDPPConfig(tau=None, seed=seed)

    def calibrate(fold, m, rng):   # as scripts/05: tau from blast-free specimens of other folds
        cal = []
        for f in [x for x in np.unique(bank.fold) if x != fold][:3]:
            for _ in range(max(1, n_calib // 3)):
                s = make_specimen(bank, f, N, 0, m, rng, "pool", "mean"); cal.append((s["F"], s["teacher"]))
        return calibrate_tau(cal, base_cfg, target_frac, rng=rng)

    rows = []
    for fold in sorted(np.unique(bank.fold)):
        for m in a.m_list:
            rng = np.random.default_rng(seed + 1000 * int(fold) + m); tau = calibrate(fold, m, rng)
            cfg = InfoDPPConfig(**{**base_cfg.__dict__, "tau": tau})
            for ep in range(a.episodes):
                s = make_specimen(bank, fold, N, k, m, rng, "pool", "mean")
                res = run_infodpp(s["F"], s["teacher"], cfg, np.random.default_rng(rng.integers(1 << 30)))
                n_info = res["k_star"]; ci = np.flatnonzero(s["crit"]); uni = rng.permutation(N)
                rows.append(dict(encoder=a.encoder, fold=int(fold), m=m, ep=ep, tau=tau, infodpp_n=n_info,
                                 infodpp_miss=int(not np.isin(ci, res["order"][:n_info]).any()),
                                 uniform_miss=int(not np.isin(ci, uni[:n_info]).any())))
            print(f"  {a.encoder} fold {fold} m {m}: tau={tau:.3f} ({time.time() - T0:.0f}s)", flush=True)
    D = pd.DataFrame(rows); D.to_csv(f"{a.out}/part1_episodes_{a.encoder}.csv", index=False)
    out = dict(encoder=a.encoder, seed=seed, episodes_per_fold=a.episodes, by_m={})
    for m, g in D.groupby("m"):
        # fold-level differences with a t interval on the five folds, as in the paper
        dif = g.groupby("fold").apply(lambda h: h.infodpp_miss.mean() - h.uniform_miss.mean())
        from scipy.stats import t as tdist
        half = float(tdist.ppf(0.975, len(dif) - 1) * dif.std(ddof=1) / np.sqrt(len(dif)))
        out["by_m"][int(m)] = dict(n=int(len(g)), infodpp_miss=float(g.infodpp_miss.mean()), uniform_miss=float(g.uniform_miss.mean()),
                                   diff=float(dif.mean()), diff_ci95=[float(dif.mean() - half), float(dif.mean() + half)],
                                   median_stop_fields=float(g.infodpp_n.median()), mean_stop_fields=float(g.infodpp_n.mean()))
    json.dump(dict(config=vars(a) | dict(seed_used=seed), env=env()), open(f"{a.out}/part1_config_{a.encoder}.json", "w"), indent=1)
    json.dump(out, open(f"{a.out}/part1_summary_{a.encoder}.json", "w"), indent=1)
    print(json.dumps(out["by_m"], indent=1)); print(f"DONE part1 {a.encoder} in {time.time() - T0:.0f}s", flush=True)
