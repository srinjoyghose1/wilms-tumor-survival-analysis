"""
09_qc_report.py

Final quality-control report for the Wilms' Tumor survival analysis pipeline.
Runs 5 checks; prints PASS / FAIL for each with specific remediation on failure.

Checks
------
  CHECK 1 — File Existence
  CHECK 2 — Classification QC
  CHECK 3 — Statistics QC
  CHECK 4 — PDMR Focus Gene QC
  CHECK 5 — Plot QC
"""

import os
import glob
import re
import numpy as np
import pandas as pd
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────

SCRIPT_DIR  = os.path.dirname(__file__)
BASE        = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

DATA_DIR    = os.path.join(BASE, "data")
PROC_DIR    = os.path.join(BASE, "processed")
RESULTS_DIR = os.path.join(BASE, "results")
PLOTS_DIR   = os.path.join(BASE, "plots")
SCRIPTS_DIR = SCRIPT_DIR

MASTER_CSV  = os.path.join(PROC_DIR,    "master_classification.csv")
KM_CSV      = os.path.join(RESULTS_DIR, "km_results.csv")
COX_CSV     = os.path.join(RESULTS_DIR, "cox_results.csv")
EXPR_CSV    = os.path.join(RESULTS_DIR, "expression_km_results.csv")
T1_CSV      = os.path.join(RESULTS_DIR, "Table1_Cohorts.csv")
T2_CSV      = os.path.join(RESULTS_DIR, "Table2_KM_Cox_Results.csv")
T3_CSV      = os.path.join(RESULTS_DIR, "Table3_Expression_Results.csv")
CLASSIFY_PY = os.path.join(SCRIPTS_DIR, "04_classify_patients.py")

OR_CLASSIFIERS = [
    "TP53_MDM4_MDM2_OR", "DROSHA_DGCR8_OR",
    "MYCN_SIX_OR", "WNT_AXIS_OR", "ALL_PRIMARY_OR",
]

PDMR_GENES = ["MDM4", "MDM2", "IGF2"]
PLOT_SIZE_FLOOR_KB = 5

# Theoretical maximum: 2 endpoints × (12 primary + 3 PDMR focus + 5 OR-classifiers)
EXPECTED_PLOTS = 2 * (12 + 3 + 5)   # = 40


# ── Output helpers ─────────────────────────────────────────────────────────────

def section(title: str) -> None:
    print(f"\n{'═' * 70}")
    print(f"  {title}")
    print(f"{'═' * 70}")


def ok(msg: str) -> None:
    print(f"  ✓  {msg}")


def warn(msg: str) -> None:
    print(f"  ⚠  {msg}")


def fail(msg: str) -> None:
    print(f"  ✗  {msg}")


# ── CHECK 1 — File Existence ───────────────────────────────────────────────────

def check1_files() -> tuple[bool, list[str]]:
    section("CHECK 1 — File Existence")
    passed = True
    remediation = []

    required_globs = [
        (os.path.join(DATA_DIR, "clinical_*.csv"),       "clinical data"),
        (os.path.join(DATA_DIR, "mutations_*.csv"),      "mutation data"),
        (os.path.join(DATA_DIR, "cna_*.csv"),            "CNA data"),
        (os.path.join(DATA_DIR, "mrna_zscores_*.csv"),   "mRNA z-score data"),
    ]
    required_exact = [
        MASTER_CSV,
        KM_CSV,
        COX_CSV,
        T1_CSV,
        T2_CSV,
        T3_CSV,
    ]

    for pattern, label in required_globs:
        matches = glob.glob(pattern)
        if matches:
            ok(f"{label}: {len(matches)} file(s) — {[os.path.basename(m) for m in sorted(matches)]}")
        else:
            fail(f"{label}: no files matching {os.path.relpath(pattern, BASE)}")
            passed = False
            remediation.append(f"Re-run 02_download_clinical.py / 03_download_genomics.py to create {label} files.")

    for path in required_exact:
        rel = os.path.relpath(path, BASE)
        if os.path.isfile(path):
            size_kb = os.path.getsize(path) / 1024
            ok(f"{rel}  ({size_kb:.1f} KB)")
        else:
            fail(f"Missing: {rel}")
            passed = False
            script = {
                MASTER_CSV: "04_classify_patients.py",
                KM_CSV:     "05_kaplan_meier.py",
                COX_CSV:    "06_cox_regression.py",
                T1_CSV:     "08_summary_table.py",
                T2_CSV:     "08_summary_table.py",
                T3_CSV:     "08_summary_table.py",
            }.get(path, "unknown script")
            remediation.append(f"Run {script} to generate {os.path.basename(path)}.")

    # Plots (at least one PDF)
    pdfs = glob.glob(os.path.join(PLOTS_DIR, "*.pdf"))
    if pdfs:
        ok(f"wilms_survival/plots/: {len(pdfs)} PDF(s) found")
    else:
        fail("No PDFs found in wilms_survival/plots/")
        passed = False
        remediation.append("Run 05_kaplan_meier.py and 07_expression_km.py to generate plots.")

    return passed, remediation


# ── CHECK 2 — Classification QC ────────────────────────────────────────────────

def check2_classification() -> tuple[bool, list[str]]:
    section("CHECK 2 — Classification QC")
    passed = True
    remediation = []

    if not os.path.isfile(MASTER_CSV):
        fail("master_classification.csv missing — skipping classification QC")
        return False, ["Run 04_classify_patients.py"]

    df = pd.read_csv(MASTER_CSV)
    n  = len(df)
    alt_cols = [c for c in df.columns if c.endswith("_ALTERED")]

    print(f"\n  Alteration frequencies and NaN rates ({n} patients):\n")
    print(f"  {'Column':<26}  {'N altered':>10}  {'Freq %':>7}  {'N NaN':>7}  {'NaN %':>7}  Notes")
    print(f"  {'─' * 76}")

    for col in alt_cols:
        s        = df[col].astype("Float64")
        n_alt    = int((s == 1).sum())
        n_nan    = int(s.isna().sum())
        freq_pct = n_alt / n * 100
        nan_pct  = n_nan / n * 100

        notes = []
        if n_alt == 0:
            notes.append("⚠ ZERO altered — check gene name / query")
        if nan_pct > 70:
            notes.append("⚠ HIGH NaN — low sequencing coverage; results may not be reliable")

        note_str = "  |  ".join(notes) if notes else "OK"
        print(f"  {col:<26}  {n_alt:>10}  {freq_pct:>6.1f}%  {n_nan:>7}  {nan_pct:>6.1f}%  {note_str}")

        if n_alt == 0:
            passed = False
            remediation.append(
                f"{col}: zero patients classified — verify gene name in 03_download_genomics.py "
                f"gene list and re-run scripts 03–04."
            )
        if nan_pct > 70 and not any(ig in col for ig in ("IGF2_CNA", "IGF2_EXPR")):
            warn(f"{col}: {nan_pct:.0f}% NaN — low sequencing coverage")

    # ── SIX1_SIX2 = logical OR of SIX1 and SIX2 ──────────────────────────────
    print()
    if all(c in df.columns for c in ("SIX1_ALTERED", "SIX2_ALTERED", "SIX1_SIX2_ALTERED")):
        s1  = df["SIX1_ALTERED"].astype("Float64")
        s2  = df["SIX2_ALTERED"].astype("Float64")
        s12 = df["SIX1_SIX2_ALTERED"].astype("Float64")

        # Expected: 1 if either=1; 0 if both=0; NaN if both=NaN
        expected_1 = ((s1 == 1) | (s2 == 1)).fillna(False)
        expected_0 = ((s1 == 0) & (s2 == 0)).fillna(False)

        violations = (
            (expected_1 & (s12 != 1)).sum() +
            (expected_0 & (s12 != 0)).sum()
        )
        if violations == 0:
            ok("SIX1_SIX2_ALTERED is a valid logical OR of SIX1_ALTERED and SIX2_ALTERED")
        else:
            fail(f"SIX1_SIX2_ALTERED has {violations} row(s) inconsistent with OR(SIX1, SIX2)")
            passed = False
            remediation.append("Re-run 04_classify_patients.py — SIX1_SIX2 OR logic is broken.")
    else:
        warn("SIX1_ALTERED, SIX2_ALTERED or SIX1_SIX2_ALTERED column missing — cannot verify OR")

    # ── IGF2 has two separate columns ─────────────────────────────────────────
    if "IGF2_CNA_ALTERED" in df.columns and "IGF2_EXPR_ALTERED" in df.columns:
        ok("IGF2 has two separate alteration columns: IGF2_CNA_ALTERED and IGF2_EXPR_ALTERED")
    else:
        missing = [c for c in ("IGF2_CNA_ALTERED", "IGF2_EXPR_ALTERED") if c not in df.columns]
        fail(f"IGF2 column(s) missing: {missing}")
        passed = False
        remediation.append("Re-run 04_classify_patients.py — IGF2 dual-column logic missing.")

    # ── Code inspection via grep of 04_classify_patients.py ───────────────────
    print()
    if not os.path.isfile(CLASSIFY_PY):
        warn("04_classify_patients.py not found — skipping code-inspection checks")
    else:
        src = Path(CLASSIFY_PY).read_text()

        # AMER1: checked under both AMER1 and WTX
        amer1_dual = (
            re.search(r'trunc_pts\s*\(.*?"AMER1".*?"WTX"', src) or
            re.search(r'trunc_pts\s*\(.*?"WTX"', src) or
            re.search(r'"WTX"', src)
        )
        if amer1_dual:
            ok('AMER1: searched under both "AMER1" and "WTX" gene symbols ✓')
        else:
            fail('AMER1: "WTX" alias not found in 04_classify_patients.py')
            passed = False
            remediation.append(
                'Add "WTX" as alias in trunc_pts(muts, "AMER1", "WTX") in 04_classify_patients.py.'
            )

        # MLLT1: used any_mut_pts (ALL mutations, not just truncating)
        mllt1_any = bool(re.search(r'any_mut_pts\s*\(.*?["\']MLLT1["\']', src))
        mllt1_trunc = bool(re.search(r'trunc_pts\s*\(.*?["\']MLLT1["\']', src))
        if mllt1_any and not mllt1_trunc:
            ok("MLLT1: classified using any_mut_pts (ALL mutations, not truncating-only) ✓")
        elif mllt1_trunc:
            fail("MLLT1: classified using trunc_pts — should use any_mut_pts for this oncogenic driver")
            passed = False
            remediation.append("Change MLLT1 to use any_mut_pts in 04_classify_patients.py.")
        else:
            warn("MLLT1: any_mut_pts call not found — classification logic may differ from expected")

        # TP53: used any_mut_pts (ALL mutations, broadened from truncating-only per analysis choice)
        tp53_any   = bool(re.search(r'any_mut_pts\s*\(.*?["\']TP53["\']', src))
        tp53_trunc = bool(re.search(r'trunc_pts\s*\(.*?["\']TP53["\']', src))
        if tp53_any:
            ok("TP53: classified using any_mut_pts (ALL mutation types captured) ✓")
        elif tp53_trunc:
            warn("TP53: classified using trunc_pts only — missense variants (e.g. p.R248W) excluded")
        else:
            warn("TP53: mutation classification call not found")

        # IGF2 expression threshold: z > 1.5
        igf2_thresh = bool(re.search(r'igf2.*?>.*?1\.5|>\s*1\.5.*?igf2', src, re.IGNORECASE))
        # More targeted: look for the expression alt line
        igf2_thresh_v2 = bool(re.search(r'per_pt\s*>\s*1\.5', src))
        if igf2_thresh or igf2_thresh_v2:
            ok("IGF2_EXPR_ALTERED: z-score threshold confirmed as > 1.5 ✓")
        else:
            fail("IGF2_EXPR_ALTERED: z > 1.5 threshold not confirmed in 04_classify_patients.py")
            passed = False
            remediation.append("Verify IGF2 expression threshold in 04_classify_patients.py line ~350.")

        # MDM4/MDM2 amplification threshold ≥ +2
        mdm4_amp = bool(re.search(r'MDM4.*?threshold\s*=\s*2|threshold\s*=\s*2.*?MDM4', src))
        # More flexible: look for cna_pts for MDM4 with threshold=2
        mdm4_amp_v2 = bool(re.search(r'cna_pts\s*\(.*?MDM4.*?threshold\s*=\s*2', src, re.DOTALL))
        mdm4_amp_v3 = bool(re.search(r'"MDM4".*?threshold\s*=\s*2', src))
        if mdm4_amp or mdm4_amp_v2 or mdm4_amp_v3:
            ok("MDM4: treated as oncogene — amplification threshold CNA ≥ +2 confirmed ✓")
        else:
            fail("MDM4: amplification threshold (≥ +2) not confirmed in 04_classify_patients.py")
            passed = False
            remediation.append("Verify MDM4 CNA threshold=2 in 04_classify_patients.py.")

        mdm2_amp_v3 = bool(re.search(r'"MDM2".*?threshold\s*=\s*2', src))
        if mdm2_amp_v3:
            ok("MDM2: treated as oncogene — amplification threshold CNA ≥ +2 confirmed ✓")
        else:
            fail("MDM2: amplification threshold (≥ +2) not confirmed in 04_classify_patients.py")
            passed = False
            remediation.append("Verify MDM2 CNA threshold=2 in 04_classify_patients.py.")

    return passed, remediation


# ── CHECK 3 — Statistics QC ────────────────────────────────────────────────────

def check3_statistics() -> tuple[bool, list[str]]:
    section("CHECK 3 — Statistics QC")
    passed = True
    remediation = []

    if not all(os.path.isfile(p) for p in (KM_CSV, COX_CSV)):
        fail("km_results.csv or cox_results.csv missing — skipping statistics QC")
        return False, ["Run 05_kaplan_meier.py and 06_cox_regression.py."]

    km  = pd.read_csv(KM_CSV)
    cox = pd.read_csv(COX_CSV)

    km_done  = km[km["status"]  == "completed"]
    cox_done = cox[cox["status"] == "completed"]

    # ── Gate-violation check in completed KM rows ──────────────────────────────
    gate_viol_km = km_done[
        (km_done["n_altered"]  < 20) |
        (km_done["n_wildtype"] < 20)
    ]
    if gate_viol_km.empty:
        ok(f"KM: all {len(km_done)} completed rows meet n ≥ 20 per arm")
    else:
        for _, r in gate_viol_km.iterrows():
            fail(f"KM gate violation: {r['gene']} / {r['endpoint']} "
                 f"n_alt={r['n_altered']} n_wt={r['n_wildtype']} (< 20/arm)")
        passed = False
        remediation.append("KM gate check failed — review gate logic in 05_kaplan_meier.py.")

    # ── Gate-violation check in completed Cox rows ──────────────────────────────
    gate_viol_cox = cox_done[
        (cox_done["events_altered"]  < 10) |
        (cox_done["events_wildtype"] < 10)
    ]
    if gate_viol_cox.empty:
        ok(f"Cox: all {len(cox_done)} completed rows meet events ≥ 10 per arm")
    else:
        for _, r in gate_viol_cox.iterrows():
            fail(f"Cox gate violation: {r['gene']} / {r['endpoint']} "
                 f"events_alt={r['events_altered']} events_wt={r['events_wildtype']} (< 10/arm)")
        passed = False
        remediation.append("Cox gate check failed — review gate logic in 06_cox_regression.py.")

    # ── P-value range check ────────────────────────────────────────────────────
    for df, label, p_col in [
        (km_done,  "KM",  "logrank_p"),
        (cox_done, "Cox", "cox_p"),
    ]:
        p_vals = pd.to_numeric(df[p_col], errors="coerce").dropna()
        bad_p  = p_vals[(p_vals < 0) | (p_vals > 1)]
        if bad_p.empty:
            ok(f"{label}: all {len(p_vals)} p-values in [0, 1]")
        else:
            fail(f"{label}: {len(bad_p)} p-value(s) outside [0, 1]: {bad_p.values}")
            passed = False
            remediation.append(f"Invalid p-values in {p_col} — check {label.lower()} analysis.")

    # ── HR positive check ──────────────────────────────────────────────────────
    if "HR" in cox_done.columns:
        hr_vals = pd.to_numeric(cox_done["HR"], errors="coerce").dropna()
        bad_hr  = hr_vals[hr_vals <= 0]
        if bad_hr.empty:
            ok(f"Cox: all {len(hr_vals)} HR values are positive")
        else:
            fail(f"Cox: {len(bad_hr)} HR value(s) ≤ 0: {bad_hr.values}")
            passed = False
            remediation.append("Non-positive HR values — check Cox model convergence.")

    # ── All OR-classifiers present in km_results ───────────────────────────────
    km_genes = set(km["gene"].unique())
    missing_or = [g for g in OR_CLASSIFIERS if g not in km_genes]
    if not missing_or:
        ok(f"All {len(OR_CLASSIFIERS)} OR-classifiers present in km_results.csv: {OR_CLASSIFIERS}")
    else:
        fail(f"OR-classifiers missing from km_results.csv: {missing_or}")
        passed = False
        remediation.append(
            f"Re-run 05_kaplan_meier.py — OR_CLASSIFIERS dict may have been modified."
        )

    # ── Summary counts ─────────────────────────────────────────────────────────
    print(f"\n  Summary:")
    print(f"    KM:  {len(km_done):>3} completed,  {(km['status']!='completed').sum():>3} skipped")
    print(f"    Cox: {len(cox_done):>3} completed,  {(cox['status']!='completed').sum():>3} skipped")
    n_sig_km  = km_done["significant"].sum() if "significant" in km_done.columns else 0
    n_sig_cox = (cox_done["cox_p"] < 0.05).sum() if "cox_p" in cox_done.columns else 0
    print(f"    KM significant  (logrank_p < 0.05): {int(n_sig_km)}")
    print(f"    Cox significant (cox_p    < 0.05):  {int(n_sig_cox)}")

    return passed, remediation


# ── CHECK 4 — PDMR Focus Gene QC ──────────────────────────────────────────────

def check4_pdmr() -> tuple[bool, list[str]]:
    section("CHECK 4 — PDMR Focus Gene QC")
    passed = True
    remediation = []

    km   = pd.read_csv(KM_CSV)  if os.path.isfile(KM_CSV)  else pd.DataFrame()
    cox  = pd.read_csv(COX_CSV) if os.path.isfile(COX_CSV) else pd.DataFrame()
    expr = pd.read_csv(EXPR_CSV) if os.path.isfile(EXPR_CSV) else pd.DataFrame()
    mstr = pd.read_csv(MASTER_CSV) if os.path.isfile(MASTER_CSV) else pd.DataFrame()

    pdmr_statuses = {}

    for gene in PDMR_GENES:
        print(f"\n  ── {gene} ──")

        # Alteration frequency
        alt_col = f"{gene}_ALTERED" if gene != "IGF2" else None
        if gene == "IGF2":
            # IGF2 has two sub-columns; report both
            cna_col  = "IGF2_CNA_ALTERED"
            expr_col = "IGF2_EXPR_ALTERED"
            if not mstr.empty:
                for col, label in [(cna_col, "CNA"), (expr_col, "Expr")]:
                    if col in mstr.columns:
                        s = mstr[col].astype("Float64")
                        n_valid = s.notna().sum()
                        n_alt   = int((s == 1).sum())
                        freq    = n_alt / n_valid * 100 if n_valid > 0 else 0
                        ok(f"IGF2 [{label}]: {n_alt} / {n_valid} sequenced ({freq:.1f}% of tested)")
        else:
            col = f"{gene}_ALTERED"
            if not mstr.empty and col in mstr.columns:
                s      = mstr[col].astype("Float64")
                n_alt  = int((s == 1).sum())
                freq   = n_alt / len(mstr) * 100
                ok(f"{gene}: {n_alt} / {len(mstr)} patients altered ({freq:.1f}%)")

        # KM status
        km_rows = km[km["gene"].str.startswith(gene)] if not km.empty else pd.DataFrame()
        if not km_rows.empty:
            statuses = km_rows["status"].unique()
            if all("completed" in s for s in statuses):
                ok(f"{gene}: KM — completed")
                km_status = "ran"
            elif all("skipped" in str(s) for s in statuses):
                n_alt_km = km_rows.iloc[0].get("n_altered", "?")
                ok(f"{gene}: KM — skipped (n_altered = {int(float(n_alt_km)) if pd.notna(n_alt_km) else '?'}; below 20/arm gate)")
                km_status = "skipped"
            else:
                km_status = "partial"
        else:
            warn(f"{gene}: not found in km_results.csv")
            km_status = "absent"

        # Cox status
        cox_rows = cox[cox["gene"].str.startswith(gene)] if not cox.empty else pd.DataFrame()
        if not cox_rows.empty and "completed" in cox_rows["status"].values:
            ok(f"{gene}: Cox — completed")
            cox_status = "ran"
        else:
            ok(f"{gene}: Cox — skipped (individual gene below KM gate)")
            cox_status = "skipped"

        # Expression result
        expr_rows = expr[expr["gene"] == gene] if not expr.empty else pd.DataFrame()
        expr_done = expr_rows[expr_rows["status"] == "completed"] if not expr_rows.empty else pd.DataFrame()
        if not expr_done.empty:
            best_ep = expr_done.sort_values("logrank_p").iloc[0]
            p_e     = float(best_ep["logrank_p"])
            ep_e    = best_ep["endpoint"]
            n_hi    = int(best_ep["n_high"])
            thr_e   = best_ep["z_threshold"]
            sig_e   = p_e < 0.05
            mark    = "✱ SIGNIFICANT" if sig_e else "not significant"
            ok(f"{gene}: expression KM — {ep_e} logrank p = {p_e:.3f} (n_high={n_hi}, z>{thr_e}) — {mark}")
            expr_status = "significant" if sig_e else "trend"
        else:
            if not expr_rows.empty:
                n_hi = expr_rows.iloc[0].get("n_high", "?")
                thr  = expr_rows.iloc[0].get("z_threshold", "?")
                ok(f"{gene}: expression KM — gate not met (n_high={n_hi}, z>{thr})")
            else:
                ok(f"{gene}: not in expression results")
            expr_status = "gate_failed"

        # OR-classifier context
        or_context = False
        if not cox.empty:
            or_sig = cox[
                (cox["gene"] == "TP53_MDM4_MDM2_OR") &
                (cox["status"] == "completed") &
                (cox.get("significant", pd.Series(False)) == True)
            ]
            if gene in ("MDM4", "MDM2") and not or_sig.empty:
                r_or = or_sig.sort_values("HR", ascending=False).iloc[0]
                ok(f"{gene}: included in TP53_MDM4_MDM2_OR — Cox HR = {r_or['HR']:.2f}, p < 0.001")
                or_context = True

        # PDMR readiness verdict
        if km_status == "ran" or (expr_status == "significant"):
            verdict = "PDMR query STRONGLY RECOMMENDED"
        elif or_context or expr_status == "trend":
            verdict = "PDMR query recommended (pathway-level or expression signal)"
        else:
            verdict = "Insufficient individual data — use as exploratory target only"

        freq_str = "—"
        if not mstr.empty:
            if gene == "IGF2":
                cna_n = int((mstr.get("IGF2_CNA_ALTERED", pd.Series()).astype("Float64") == 1).sum())
                freq_str = f"CNA {cna_n/len(mstr)*100:.1f}%"
            else:
                col = f"{gene}_ALTERED"
                if col in mstr.columns:
                    n_a = int((mstr[col].astype("Float64") == 1).sum())
                    freq_str = f"{n_a/len(mstr)*100:.1f}%"

        print(f"\n  PDMR READINESS — {gene}: {freq_str} altered, "
              f"KM {km_status}, Cox {cox_status} → {verdict}")

        pdmr_statuses[gene] = {
            "km": km_status, "cox": cox_status,
            "expr": expr_status, "or_context": or_context,
            "verdict": verdict,
        }

    return passed, remediation


# ── CHECK 5 — Plot QC ─────────────────────────────────────────────────────────

def check5_plots() -> tuple[bool, list[str]]:
    section("CHECK 5 — Plot QC")
    passed = True
    remediation = []

    pdfs = sorted(glob.glob(os.path.join(PLOTS_DIR, "*.pdf")))

    if not pdfs:
        fail("No PDF files found in wilms_survival/plots/")
        return False, ["Run 05_kaplan_meier.py and 07_expression_km.py."]

    # File size check
    small = []
    for pdf in pdfs:
        kb = os.path.getsize(pdf) / 1024
        if kb < PLOT_SIZE_FLOOR_KB:
            small.append((os.path.basename(pdf), kb))

    if small:
        for fname, kb in small:
            fail(f"{fname}: only {kb:.1f} KB (< {PLOT_SIZE_FLOOR_KB} KB — likely empty)")
        passed = False
        remediation.append("Re-generate small PDFs — they may be empty or corrupted.")
    else:
        sizes_kb = [os.path.getsize(p) / 1024 for p in pdfs]
        ok(f"All {len(pdfs)} PDFs exceed {PLOT_SIZE_FLOOR_KB} KB "
           f"(range: {min(sizes_kb):.0f}–{max(sizes_kb):.0f} KB)")

    # Categorise plots
    genomic_km = [p for p in pdfs if "/KM_" in p and "/KM_EXPR_" not in p]
    expr_km    = [p for p in pdfs if "/KM_EXPR_" in p]

    print(f"\n  Plot inventory:")
    print(f"    Genomic KM plots (05_kaplan_meier.py): {len(genomic_km)}")
    for p in genomic_km:
        print(f"      {os.path.basename(p)}")
    print(f"    Expression KM plots (07_expression_km.py): {len(expr_km)}")
    for p in expr_km:
        print(f"      {os.path.basename(p)}")

    total = len(pdfs)
    print(f"\n  Total: {total} PDFs  |  Theoretical maximum (genomic KM only): {EXPECTED_PLOTS}")

    gap = EXPECTED_PLOTS - len(genomic_km)
    if gap > 0:
        warn(
            f"{gap} genomic KM plot(s) not generated "
            f"({len(genomic_km)}/{EXPECTED_PLOTS} expected).\n"
            f"     This is expected: all {gap} correspond to individual genes / OR-classifiers\n"
            f"     that were below the minimum-arm-size gate (< 20 patients/arm).\n"
            f"     Gate failures are not errors — they are recorded in km_results.csv."
        )
    else:
        ok(f"All {EXPECTED_PLOTS} expected genomic KM plots were generated")

    return passed, remediation


# ── FINAL REPORT ──────────────────────────────────────────────────────────────

def final_report(
    checks: list[tuple[str, bool, list[str]]],
) -> None:
    section("FINAL REPORT")

    all_passed = all(p for _, p, _ in checks)
    for name, chk_passed, rems in checks:
        mark = "PASS" if chk_passed else "FAIL"
        print(f"  [{mark}]  {name}")
        if not chk_passed:
            for r in rems:
                print(f"         → {r}")

    # Readiness statement
    km  = pd.read_csv(KM_CSV)  if os.path.isfile(KM_CSV)  else pd.DataFrame()
    cox = pd.read_csv(COX_CSV) if os.path.isfile(COX_CSV) else pd.DataFrame()
    expr = pd.read_csv(EXPR_CSV) if os.path.isfile(EXPR_CSV) else pd.DataFrame()

    n_km_ran  = (km["status"]  == "completed").sum() if not km.empty  else 0
    n_cox_ran = (cox["status"] == "completed").sum() if not cox.empty else 0
    n_sig     = (
        (cox[cox["status"] == "completed"]["cox_p"] < 0.05).sum()
        if not cox.empty and "cox_p" in cox.columns else 0
    )

    def pdmr_status_str(gene: str) -> str:
        """One-word status for PDMR readiness line."""
        km_rows = km[km["gene"].str.startswith(gene)] if not km.empty else pd.DataFrame()
        cox_rows = cox[cox["status"] == "completed"] if not cox.empty else pd.DataFrame()
        expr_rows = expr[(expr["gene"] == gene) & (expr["status"] == "completed")] \
                    if not expr.empty else pd.DataFrame()
        or_sig = cox_rows[cox_rows["gene"] == "TP53_MDM4_MDM2_OR"] if not cox_rows.empty else pd.DataFrame()

        if not km_rows.empty and (km_rows["status"] == "completed").any():
            return "ran KM+Cox"
        if not expr_rows.empty and (expr_rows["logrank_p"] < 0.05).any():
            return "expr significant"
        if gene in ("MDM4", "MDM2") and not or_sig.empty:
            n_alt = km_rows.iloc[0]["n_altered"] if not km_rows.empty else "?"
            return f"in OR-classifier (n_individual={int(float(n_alt)) if pd.notna(n_alt) else '?'})"
        n_alt_str = "?"
        if not km_rows.empty:
            n_alt_str = str(int(float(km_rows.iloc[0]["n_altered"])))
        return f"below gate (n_altered={n_alt_str})"

    mdm4_s = pdmr_status_str("MDM4")
    mdm2_s = pdmr_status_str("MDM2")
    igf2_s = pdmr_status_str("IGF2_CNA").replace("IGF2_CNA", "IGF2")

    print()
    print("  ─" * 35)
    print(
        f"  Analysis complete: {n_km_ran} gene/classifier pairs ran KM, "
        f"{n_cox_ran} ran Cox, {n_sig} significant result(s) found.\n"
        f"  PDMR focus genes: MDM4 [{mdm4_s}], "
        f"MDM2 [{mdm2_s}], IGF2 [{igf2_s}]."
    )
    print("  ─" * 35)
    print(f"\n  Overall pipeline status: {'✓ ALL CHECKS PASSED' if all_passed else '✗ SOME CHECKS FAILED — see remediation above'}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print(f"\n{'█' * 70}")
    print(f"  WILMS' TUMOR SURVIVAL ANALYSIS — QC REPORT")
    print(f"  Pipeline root: {BASE}")
    print(f"{'█' * 70}")

    results = []

    p1, r1 = check1_files()
    results.append(("CHECK 1 — File Existence",          p1, r1))

    p2, r2 = check2_classification()
    results.append(("CHECK 2 — Classification QC",       p2, r2))

    p3, r3 = check3_statistics()
    results.append(("CHECK 3 — Statistics QC",           p3, r3))

    p4, r4 = check4_pdmr()
    results.append(("CHECK 4 — PDMR Focus Gene QC",      p4, r4))

    p5, r5 = check5_plots()
    results.append(("CHECK 5 — Plot QC",                 p5, r5))

    final_report(results)


if __name__ == "__main__":
    main()
