"""
Reproduce the results of the paper (paper/main.tex and paper/jbhi_compact.tex).

    python reproduce.py list [--group G] [--markdown]   what exists, which paper item it feeds, what it needs
    python reproduce.py check [--shangrla DIR]           data-free checks: validity tests, Part II exact checks,
                                                         (with --shangrla) the cross-check against the reference ALPHA
    python reproduce.py run ID [ID ...] [--dry-run]      run steps; outputs go to their results/ folders
    python reproduce.py verify ID [ID ...]               re-run a step into a temporary folder and compare its summary
                                                         files with the stored ones (numbers within 1e-9)

Every command is also listed in REPRODUCE.md, so the steps can be run by hand. Colab steps (GPU embeddings, image and
annotation downloads) cannot be run from here; they are listed with the notebook to open.
Status: P = pre-registered (confirmatory), X = exploratory, C = check, D = data preparation.
Cost: S < 5 min, M 5-60 min, L > 1 h on a laptop CPU (approximate); GPU = Colab notebook.
"""
import argparse, json, math, os, shutil, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.abspath(__file__)); PY = sys.executable


def S(*a):
    return ["{py}", *a]


STEPS = [
    # ------------------------------------------------------------------ checks without data
    dict(id="tests", group="check", status="C", cost="S", paper="Validity claims; Section 'Implementation checks'",
         needs=[], cmds=[S("-m", "pytest", "-q", "tests")], out=None, compare=[]),
    dict(id="part2", group="check", status="C", cost="S", paper="Prop. 1, Theorem 1, Corollary (exact and LP checks)",
         needs=[], cmds=[S("scripts/17_part2_verify.py", "--out", "{out}")], out="results/part2_v1",
         compare=["verify.json", "nstar_table.csv"]),
    dict(id="crosscheck", group="check", status="C", cost="S", paper="Implementation checks (reference ALPHA, SHANGRLA)",
         needs=["{shangrla}/NonnegMean.py"],
         cmds=[S("-I", "scripts/44_crosscheck_shangrla.py", "--shangrla", "{shangrla}", "--R", "2000", "--seed", "20261035", "--out", "{out}")],
         out="results/crosscheck_shangrla_v1", compare=["summary.json"]),
    # ------------------------------------------------------------------ KU-Optofil preparation
    dict(id="ku_inspect", group="data", status="D", cost="S", paper="Data (KU-Optofil counts)",
         needs=["data/ku_optofil/metadata.csv"], cmds=[S("scripts/01_inspect_kuoptofil.py", "--data", "data/ku_optofil")], out=None, compare=[]),
    dict(id="ku_folds", group="data", status="D", cost="S", paper="Patient-grouped folds",
         needs=["data/ku_optofil/metadata.csv"],
         cmds=[S("scripts/02_make_folds.py", "--data", "data/ku_optofil", "--out", "results/folds_v1", "--seed", "20260919")],
         out=None, compare=[]),
    dict(id="ku_embeddings", group="data", status="D", cost="GPU", paper="Embeddings and cross-fitted probe (AUC 0.979)",
         needs=["results/folds_v1"], colab="notebooks/03_embeddings_colab.ipynb", cmds=[], out="results/embeddings_v1", compare=[]),
    dict(id="separability", group="data", status="X", cost="S", paper="Blast separability by pooling (background for Kill-test 1)",
         needs=["results/embeddings_v1/cells_features.csv"], cmds=[S("scripts/04_separability.py", "--out", "{out}")],
         out="results/separability_v1", compare=["summary.json"]),
    # ------------------------------------------------------------------ Part I (RQ1)
    dict(id="partI_confirm", group="partI", status="P", cost="L", paper="Fig. 1; H1a-H1d; Kill-test 1",
         needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/05_killtests.py", "--out", "results/confirm_v1", "--seed", "20260401", "--episodes", "200",
                 "--cell_episodes", "10:1:600", "20:1:600"),
               S("scripts/09_confirmatory_analysis.py", "--dir", "results/confirm_v1")], out=None, compare=[]),
    dict(id="partI_sens_max", group="partI", status="P", cost="L", paper="Sensitivity: max pooling",
         needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/05_killtests.py", "--out", "results/sensconf_pool_max", "--seed", "20260402", "--episodes", "200", "--pooling", "max")],
         out=None, compare=[]),
    dict(id="partI_sens_single", group="partI", status="P", cost="L", paper="Sensitivity: single-cell representation",
         needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/05_killtests.py", "--out", "results/sensconf_pool_single", "--seed", "20260403", "--episodes", "200", "--pooling", "single")],
         out=None, compare=[]),
    dict(id="partI_sens_tau", group="partI", status="P", cost="L", paper="Sensitivity: operating point 5-50%",
         needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/05_killtests.py", "--out", f"results/sensconf_tau_{t}", "--seed", "20260404", "--episodes", "200",
                 "--k_list", "1", "--target_frac", t) for t in ("0.05", "0.10", "0.163", "0.30", "0.50")], out=None, compare=[]),
    dict(id="partI_sens_host", group="partI", status="P", cost="L", paper="Sensitivity: single-host background",
         needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/05_killtests.py", "--out", "results/sensconf_pool_host", "--seed", "20260405", "--episodes", "200", "--pool", "host")],
         out=None, compare=[]),
    dict(id="partI_op24", group="partI", status="X", cost="M", paper="InfoDPP-PAC at its own operating point (2.4%)",
         needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/05_killtests.py", "--out", "results/sens_op_paperpool", "--seed", "20261022", "--episodes", "200",
                 "--k_list", "1", "2", "5", "--target_frac", "0.024")], out=None, compare=[]),
    dict(id="partI_field_teacher", group="partI", status="X", cost="M", paper="Kill-test 1 with the field-level teacher",
         needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/26_kt1_field_teacher.py", "--out", "{out}", "--seed", "20261019")], out="results/kt1_field_teacher", compare=[]),
    dict(id="partI_rules", group="partI", status="X", cost="M", paper="Other saturation rules; ranked search (Prop. 3)",
         needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/31_stopping_class.py", "--out", "{out}", "--seed", "20261024")], out="results/stopping_class_v1",
         compare=["summary.json"]),
    dict(id="review_stats", group="partI", status="X", cost="S", paper="Fold t intervals, effect sizes, sensitivity table",
         needs=["results/confirm_v1/episodes.csv", "results/kt1_field_teacher/episodes.csv"],
         cmds=[S("scripts/29_review_stats.py")], out=None, compare=[]),
    dict(id="tiles_embeddings", group="partI", status="D", cost="GPU", paper="Real BoneMarrowWSI tiles",
         needs=["data/bmwsi/bmwsi_ann.zip (Drive)"], colab="notebooks/35_bmwsi_tiles_colab.ipynb", cmds=[], out="data/bmwsi_tiles", compare=[]),
    dict(id="partI_tiles", group="partI", status="X", cost="S", paper="Part I on real tiles",
         needs=["data/bmwsi_tiles/tiles.csv", "data/bmwsi/cells.csv"],
         cmds=[S("scripts/36_partI_real_tiles.py", "--tiles", "data/bmwsi_tiles", "--out", "results/partI_real_tiles_v1", "--seed", "20261028"),
               S("scripts/40_partI_real_tiles_analysis.py")], out=None, compare=[]),
    # ------------------------------------------------------------------ Part III (RQ3)
    dict(id="partIII_confirm", group="partIII", status="P", cost="M", paper="Fig. 3; H3a-H3d",
         needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/12_certify_confirm.py", "--out", "results/part3_theta010", "--theta", "0.10", "--p_list", "0.05", "0.08",
                 "0.10", "0.12", "0.15", "0.20", "--episodes", "1000", "--seed", "20260511"),
               S("scripts/12_certify_confirm.py", "--out", "results/part3_theta020", "--theta", "0.20", "--p_list", "0.15", "0.18",
                 "0.20", "0.22", "0.25", "0.30", "--episodes", "1000", "--seed", "20260512"),
               S("scripts/13_part3_analysis.py", "--dirs", "results/part3_theta010", "results/part3_theta020")], out=None, compare=[]),
    dict(id="m_scaling", group="partIII", status="X", cost="L", paper="Population size table",
         needs=["results/embeddings_v1/cells_features.csv"], cmds=[S("scripts/28_m_scaling.py", "--out", "{out}", "--seed", "20261021")],
         out="results/m_scaling_v1", compare=["summary.csv"]),
    dict(id="alpha_prior", group="partIII", status="X", cost="M", paper="ALPHA prior guess", needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/39_alpha_prior.py", "--out", "{out}", "--seed", "20261030")], out="results/alpha_prior_v1", compare=["summary.csv"]),
    dict(id="decisions", group="partIII", status="X", cost="M", paper="Decisions at fixed budgets (table)", needs=[],
         cmds=[S("scripts/32_decision_level.py", "--out", "{out}", "--seed", "20261025")], out="results/decision_level_v1", compare=["summary.csv"]),
    dict(id="confusable", group="partIII", status="X", cost="M", paper="Confusable backgrounds", needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/33_confusable_background.py", "--out", "{out}", "--seed", "20261026")], out="results/confusable_bg_v1", compare=["summary.csv"]),
    dict(id="leukocyte", group="partIII", status="X", cost="M", paper="Leukocyte-only background", needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/41_leukocyte_background.py", "--out", "{out}", "--seed", "20261032")], out="results/leukocyte_bg_v1", compare=["summary.csv"]),
    dict(id="detrule", group="partIII", status="X", cost="M", paper="Impact of the rule missing from the frozen ALPHA", needs=[],
         cmds=[S("scripts/45_alpha_deterministic_rule.py", "--out", "{out}", "--seed", "20261036")], out="results/alpha_detrule_v1", compare=["summary.csv"]),
    dict(id="unequal", group="partIII", status="X", cost="M", paper="Unequal field sizes; Prop. 2",
         needs=["results/embeddings_v1/cells_features.csv", "data/bmwsi/cells.csv"],
         cmds=[S("scripts/27_unequal_fields.py", "--out", "{out}", "--seed", "20261020")], out="results/unequal_fields_v1", compare=["summary.json"]),
    # ------------------------------------------------------------------ RQ4
    dict(id="ku_error", group="rq4", status="X", cost="S", paper="Classifier error is not exchangeable", needs=["results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/10_rq4_error_structure.py"), S("scripts/11_rq4_certify_by_error_group.py")], out=None, compare=[]),
    dict(id="rvpbs_crops", group="rq4", status="D", cost="S", paper="RV-PBS crops", needs=["data/rv_pbs"],
         cmds=[S("scripts/14_crop_rvpbs.py")], out=None, compare=[]),
    dict(id="rvpbs_embeddings", group="rq4", status="D", cost="GPU", paper="RV-PBS embeddings", needs=["data/external_crops.zip (Drive)"],
         colab="notebooks/15_external_embeddings_colab.ipynb", cmds=[], out="results/external_v1", compare=[]),
    dict(id="rvpbs", group="rq4", status="P", cost="M", paper="H4a-H4c (acquisition shift)", needs=["results/external_v1/rvpbs_features.csv"],
         cmds=[S("scripts/16_external_analysis.py")], out=None, compare=[]),
    dict(id="spatial_mixture", group="rq4", status="P", cost="L", paper="H5a-H5c, H6a-H6d; Fig. 4", needs=["results/external_v1/rvpbs_features.csv"],
         cmds=[S("scripts/18_spatial_confirm.py", "--out", "results/a3_spatial_theta010", "--theta", "0.10", "--p_list", "0.08", "0.10", "0.12", "0.15", "0.20", "--seed", "20260611"),
               S("scripts/18_spatial_confirm.py", "--out", "results/a3_spatial_theta020", "--theta", "0.20", "--p_list", "0.18", "0.20", "0.22", "0.25", "0.30", "--seed", "20260612"),
               S("scripts/19_mixture_confirm.py", "--out", "results/a3_mixture_rvpbs", "--source", "rvpbs", "--p_list", "0.10", "0.15", "0.20", "--seed", "20260621"),
               S("scripts/19_mixture_confirm.py", "--out", "results/a3_mixture_ku", "--source", "ku", "--p_list", "0.10", "0.12", "0.15", "0.20", "--seed", "20260622"),
               S("scripts/20_amend3_analysis.py", "--spatial", "results/a3_spatial_theta010", "results/a3_spatial_theta020",
                 "--mixture", "results/a3_mixture_rvpbs", "results/a3_mixture_ku", "--out", "results/a3_decisions.json")], out=None, compare=[]),
    dict(id="bmwsi_annotations", group="rq4", status="D", cost="GPU", paper="BoneMarrowWSI annotations (no GPU needed)", needs=[],
         colab="notebooks/21_bmwsi_annotations_colab.ipynb", cmds=[], out="data/bmwsi/bmwsi_ann.zip", compare=[]),
    dict(id="bmwsi_parse", group="rq4", status="D", cost="S", paper="cells.csv, rois.csv", needs=["data/bmwsi/ann"],
         cmds=[S("scripts/22_parse_bmwsi.py")], out=None, compare=[]),
    dict(id="bmwsi_clustering", group="rq4", status="P", cost="M", paper="H7a; Beta-binomial bridge", needs=["data/bmwsi/cells.csv"],
         cmds=[S("scripts/23_bmwsi_clustering.py", "--out", "results/a4_bmwsi/clustering.json"),
               S("scripts/24_bmwsi_bridge_sim.py", "--clustering", "results/a4_bmwsi/clustering.json", "--out", "results/a4_bmwsi/bridge.json")],
         out=None, compare=[]),
    dict(id="icc_excess", group="rq4", status="X", cost="M", paper="Excess ICC", needs=["data/bmwsi/cells.csv"],
         cmds=[S("scripts/34_icc_excess.py", "--out", "{out}/icc_excess.json", "--seed", "20261027")], out="results/a4_bmwsi", compare=["icc_excess.json"]),
    dict(id="real_arrangements", group="rq4", status="X", cost="M", paper="Observed validity on 75 real specimens", needs=["data/bmwsi/cells.csv"],
         cmds=[S("scripts/30_real_arrangement_validity.py", "--out", "{out}", "--seed", "20261023")], out="results/real_arrangements_v1", compare=["summary.json"]),
    dict(id="mll23_embeddings", group="rq4", status="D", cost="GPU", paper="MLL23 embeddings", needs=[],
         colab="notebooks/37_mll23_embeddings_colab.ipynb", cmds=[], out="data/mll23", compare=[]),
    dict(id="mll23", group="rq4", status="X", cost="S", paper="Second external source (MLL23)",
         needs=["data/mll23/mll23_embeddings_dinov2_vits14.npy", "results/embeddings_v1/cells_features.csv"],
         cmds=[S("scripts/38_mll23_shift.py", "--mll", "data/mll23", "--out", "{out}", "--seed", "20261029")], out="results/mll23_shift_v1", compare=["summary.json"]),
    dict(id="sessions_download", group="rq4", status="D", cost="GPU", paper="Reader labels (no GPU needed)", needs=[],
         colab="notebooks/42_bmwsi_sessions_colab.ipynb", cmds=[], out="data/bmwsi_sessions", compare=[]),
    dict(id="verifier", group="rq4", status="X", cost="M", paper="Imperfect verification; Remark 2",
         needs=["data/bmwsi_sessions/annotations_sessions.csv.gz", "data/bmwsi/cells.csv"],
         cmds=[S("scripts/43_verifier_error.py", "--sessions", "data/bmwsi_sessions/annotations_sessions.csv.gz", "--out", "{out}", "--seed", "20261034")],
         out="results/verifier_error_v1", compare=["summary.json"]),
    # ------------------------------------------------------------------ figures and summary tables
    dict(id="benchmark", group="paper", status="X", cost="S", paper="Benchmark table across conditions (main Table 3, compact Table IV)",
         needs=["results/part3_theta010/episodes_p0.10.csv", "results/a3_spatial_theta010", "results/a3_mixture_ku",
                "results/a3_mixture_rvpbs", "results/confusable_bg_v1/episodes.csv", "results/unequal_fields_v1/episodes.csv",
                "results/mll23_shift_v1/episodes.csv", "results/real_arrangements_v1/episodes.csv",
                "results/verifier_error_v1/episodes.csv"],
         cmds=[S("scripts/46_benchmark_table.py", "--out", "{out}")], out="results/benchmark_v1", compare=["summary.json"]),
    dict(id="figures", group="paper", status="C", cost="S", paper="Figs. 1-4 (copy figures/*.pdf to paper/figures/)",
         needs=["results/confirm_v1", "results/part3_theta010", "results/a3_decisions.json", "results/a4_bmwsi/clustering.json"],
         cmds=[S("scripts/25_figures.py")], out=None, compare=[]),
]
BY_ID = {s["id"]: s for s in STEPS}


def fmt(cmd, out, shangrla):
    return [c.replace("{py}", PY).replace("{out}", out or "").replace("{shangrla}", shangrla or "<SHANGRLA folder>") for c in cmd]


def missing_inputs(step, shangrla):
    miss = []
    for n in step["needs"]:
        if "(Drive)" in n: continue
        p = n.replace("{shangrla}", shangrla or "<none>")
        if not os.path.exists(os.path.join(ROOT, p)): miss.append(p)
    return miss


def run_step(step, out=None, shangrla=None, dry=False):
    if step.get("colab"):
        print(f"[{step['id']}] Colab step: open {step['colab']} in Google Colab and run all cells; output: {step['out']}"); return 0
    miss = missing_inputs(step, shangrla)
    if miss and not dry:
        print(f"[{step['id']}] missing inputs: {', '.join(miss)} (see REPRODUCE.md, 'Data')"); return 2
    out = out or step["out"]
    for c in step["cmds"]:
        cmd = fmt(c, out, shangrla); print(f"[{step['id']}] $ " + " ".join(cmd), flush=True)
        if dry: continue
        r = subprocess.run(cmd, cwd=ROOT)
        if r.returncode: print(f"[{step['id']}] failed with exit code {r.returncode}"); return r.returncode
    return 0


def close(a, b, tol=1e-9):
    """Numeric comparison of nested JSON-like objects or CSV tables; returns (equal, max_abs_diff)."""
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b): return False, math.inf
        res = [close(a[k], b[k], tol) for k in a]; return all(r[0] for r in res), max([r[1] for r in res] or [0.0])
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): return False, math.inf
        res = [close(x, y, tol) for x, y in zip(a, b)]; return all(r[0] for r in res), max([r[1] for r in res] or [0.0])
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if (a != a and b != b) or a == b: return True, 0.0
        d = abs(a - b); return d <= tol * max(1.0, abs(a), abs(b)), d
    return a == b, 0.0 if a == b else math.inf


def compare_file(new, old):
    if new.endswith(".json"):
        return close(json.load(open(new)), json.load(open(old)))
    import pandas as pd
    A, B = pd.read_csv(new), pd.read_csv(old)
    if list(A.columns) != list(B.columns) or len(A) != len(B): return False, math.inf
    return close(A.where(A.notna(), None).values.tolist(), B.where(B.notna(), None).values.tolist())


def verify_step(step, shangrla=None):
    if not step["compare"] or not step["out"]:
        print(f"[{step['id']}] no stored summary to compare with; use 'run' and inspect 'git diff results/'"); return 0
    tmp = tempfile.mkdtemp(prefix=f"verify_{step['id']}_")
    try:
        code = run_step(step, out=tmp, shangrla=shangrla)
        if code: return code
        ok_all = True
        for f in step["compare"]:
            new, old = os.path.join(tmp, f), os.path.join(ROOT, step["out"], f)
            if not os.path.exists(old): print(f"[{step['id']}] stored {old} not found"); ok_all = False; continue
            ok, d = compare_file(new, old); ok_all &= ok
            print(f"[{step['id']}] {f}: {'IDENTICAL' if ok and d == 0 else ('equal within 1e-9' if ok else 'DIFFERS')} (max |diff| {d:.3g})")
        return 0 if ok_all else 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("list"); p.add_argument("--group"); p.add_argument("--markdown", action="store_true")
    p = sub.add_parser("check"); p.add_argument("--shangrla")
    p = sub.add_parser("run"); p.add_argument("ids", nargs="+"); p.add_argument("--dry-run", action="store_true"); p.add_argument("--shangrla")
    p = sub.add_parser("verify"); p.add_argument("ids", nargs="+"); p.add_argument("--shangrla")
    a = ap.parse_args()
    if a.cmd == "list":
        rows = [s for s in STEPS if not a.group or s["group"] == a.group]
        if a.markdown:
            print("| ID | Status | Cost | Paper item | How |\n|---|---|---|---|---|")
            for s in rows:
                how = f"Colab: `{s['colab']}`" if s.get("colab") else "<br>".join(
                    "`" + " ".join(fmt(c, s["out"], None)).replace(PY, "python") + "`" for c in s["cmds"])
                print(f"| `{s['id']}` | {s['status']} | {s['cost']} | {s['paper']} | {how} |")
        else:
            for s in rows:
                print(f"{s['id']:<20} {s['group']:<8} {s['status']} {s['cost']:<4} {s['paper']}")
        return 0
    if a.cmd == "check":
        codes = [run_step(BY_ID["tests"]), verify_step(BY_ID["part2"])]
        if a.shangrla: codes.append(verify_step(BY_ID["crosscheck"], shangrla=a.shangrla))
        else: print("[crosscheck] skipped: pass --shangrla DIR (a folder holding SHANGRLA's shangrla/core/NonnegMean.py)")
        return max(codes)
    for i in a.ids:
        if i not in BY_ID: print(f"unknown step '{i}'; see 'python reproduce.py list'"); return 2
    codes = [(run_step(BY_ID[i], shangrla=a.shangrla, dry=a.dry_run) if a.cmd == "run" else verify_step(BY_ID[i], shangrla=a.shangrla))
             for i in a.ids]
    return max(codes)


if __name__ == "__main__":
    sys.exit(main())
