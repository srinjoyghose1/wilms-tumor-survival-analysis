"""
08_summary_table.py

Compile all Wilms' Tumor survival results into publication-ready summary tables
and a plain-text results paragraph.

Outputs
-------
  wilms_survival/results/Table1_Cohorts.csv
  wilms_survival/results/Table2_KM_Cox_Results.csv
  wilms_survival/results/Table3_Expression_Results.csv
"""

import os
import numpy as np
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────

SCRIPT_DIR  = os.path.dirname(__file__)
PROC_DIR    = os.path.join(SCRIPT_DIR, "..", "processed")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "..", "results")
MASTER_CSV  = os.path.join(PROC_DIR,    "master_classification.csv")
KM_CSV      = os.path.join(RESULTS_DIR, "km_results.csv")
COX_CSV     = os.path.join(RESULTS_DIR, "cox_results.csv")
EXPR_CSV    = os.path.join(RESULTS_DIR, "expression_km_results.csv")
T1_CSV      = os.path.join(RESULTS_DIR, "Table1_Cohorts.csv")
T2_CSV      = os.path.join(RESULTS_DIR, "Table2_KM_Cox_Results.csv")
T3_CSV      = os.path.join(RESULTS_DIR, "Table3_Expression_Results.csv")

# ── Constants ─────────────────────────────────────────────────────────────────

PDMR_GENES = {"MDM4", "MDM2", "IGF2", "IGF2_CNA", "IGF2_EXPR", "TP53_MDM4_MDM2_OR"}

# Per-gene alteration logic applied in script 04
ALTERATION_TYPES = {
    "TP53":              "Truncating mutation",
    "MYCN":              "Amplification (CNA ≥ +2)",
    "SIX1":              "Hotspot mutation (Q177R)",
    "SIX2":              "Hotspot mutation (Q177R)",
    "SIX1_SIX2":         "Hotspot mutation (Q177R), SIX1 or SIX2",
    "DROSHA":            "Truncating mutation",
    "DGCR8":             "Truncating mutation",
    "CTNNB1":            "Hotspot mutation (codons 32/33/34/37/41/45)",
    "AMER1":             "Truncating mutation",
    "WT1":               "Truncating mutation / Homodeletion (CNA −2)",
    "MLLT1":             "Truncating mutation",
    "NIPBL":             "Truncating mutation",
    "MDM4":              "Amplification (CNA ≥ +2)",
    "MDM2":              "Amplification (CNA ≥ +2)",
    "IGF2_CNA":          "LOI / Gain (CNA ≥ +1)",
    "IGF2_EXPR":         "mRNA overexpression (z > 1.5)",
    "TP53_MDM4_MDM2_OR": "Composite OR — TP53 | MDM4 | MDM2",
    "DROSHA_DGCR8_OR":   "Composite OR — DROSHA | DGCR8",
    "MYCN_SIX_OR":       "Composite OR — MYCN | SIX1/SIX2",
    "WNT_AXIS_OR":       "Composite OR — CTNNB1 | AMER1 | WT1",
    "ALL_PRIMARY_OR":    "Composite OR — all primary pathway genes",
}

DISPLAY_NAMES = {
    "TP53_MDM4_MDM2_OR": "TP53/MDM4/MDM2 (OR)",
    "DROSHA_DGCR8_OR":   "DROSHA/DGCR8 (OR)",
    "MYCN_SIX_OR":       "MYCN/SIX1-SIX2 (OR)",
    "WNT_AXIS_OR":       "WNT Axis (OR)",
    "ALL_PRIMARY_OR":    "All Primary Genes (OR)",
    "IGF2_CNA":          "IGF2 (CNA)",
    "IGF2_EXPR":         "IGF2 (Expression)",
    "SIX1_SIX2":         "SIX1/SIX2",
}


# ── Helpers ────────────────────────────────────────────────────────────────────

def display(gene: str) -> str:
    return DISPLAY_NAMES.get(gene, gene)


def star(gene: str) -> str:
    return "★ " + display(gene) if gene in PDMR_GENES else display(gene)


def fmt_p(p, threshold=0.001) -> str:
    if pd.isna(p):
        return "—"
    return "< 0.001" if float(p) < threshold else f"{float(p):.3f}"


def fmt_ci(lo, hi) -> str:
    if pd.isna(lo) or pd.isna(hi):
        return "—"
    return f"[{float(lo):.2f}–{float(hi):.2f}]"


def fmt_hr(hr) -> str:
    return "—" if pd.isna(hr) else f"{float(hr):.2f}"


def fmt_med(m) -> str:
    if pd.isna(m) or str(m).strip() in ("", "nan"):
        return "NR"
    try:
        return f"{float(m):.1f}"
    except (ValueError, TypeError):
        return str(m)


def is_sig(row, col="significant") -> bool:
    v = row.get(col, False)
    return bool(v) if not pd.isna(v) else False


# ── TABLE 1: Cohort Summary ────────────────────────────────────────────────────

def build_table1(master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for study_id, grp in master.groupby("STUDY_ID"):
        n  = len(grp)
        os_valid  = grp["OS_MONTHS"].notna()  & grp["OS_STATUS"].notna()
        efs_valid = grp["EFS_MONTHS"].notna() & grp["EFS_STATUS"].notna()
        os_ev     = int(grp.loc[os_valid,  "OS_STATUS"].sum())
        efs_ev    = int(grp.loc[efs_valid, "EFS_STATUS"].sum())
        rows.append({
            "Cohort":           "Pediatric Wilms' Tumor (TARGET, 2018)",
            "Disease Type":     "Wilms' Tumor (Nephroblastoma)",
            "Study ID":         study_id,
            "N patients":       n,
            "N with EFS data":  int(efs_valid.sum()),
            "EFS events (n)":   efs_ev,
            "N with OS data":   int(os_valid.sum()),
            "OS events (n)":    os_ev,
            "Primary endpoint": "EFS",
        })
    return pd.DataFrame(rows)


# ── TABLE 2: Genomic Alteration Results ───────────────────────────────────────

def build_table2(km: pd.DataFrame, cox: pd.DataFrame) -> pd.DataFrame:
    # Bring only completed Cox rows to join
    cox_done = cox[cox["status"] == "completed"][
        ["gene", "endpoint", "HR", "CI_lower", "CI_upper", "cox_p", "significant"]
    ].rename(columns={"significant": "cox_significant"})

    merged = km.merge(cox_done, on=["gene", "endpoint"], how="left")

    rows = []
    for _, r in merged.iterrows():
        gene = r["gene"]
        n_alt = r["n_altered"]
        n_wt  = r["n_wildtype"]
        total = (n_alt + n_wt) if (pd.notna(n_alt) and pd.notna(n_wt)) else np.nan
        freq  = round(n_alt / total * 100, 1) if (pd.notna(n_alt) and total > 0) else np.nan

        km_skip   = "skipped" in str(r.get("status", ""))
        km_done   = not km_skip
        cox_sig   = bool(r["cox_significant"]) if pd.notna(r.get("cox_significant")) else False
        km_sig    = bool(r["significant"])     if pd.notna(r.get("significant"))     else False
        overall_sig = cox_sig or (km_done and km_sig)

        rows.append({
            "Gene":                  star(gene),
            "Pathway":               r.get("pathway", "—"),
            "Is PDMR Focus":         gene in PDMR_GENES,
            "Alteration Type":       ALTERATION_TYPES.get(gene, "Unknown"),
            "N altered":             int(n_alt)  if pd.notna(n_alt)  else "—",
            "N wildtype":            int(n_wt)   if pd.notna(n_wt)   else "—",
            "Freq %":                f"{freq:.1f}" if pd.notna(freq)  else "—",
            "Endpoint":              r["endpoint"],
            "Median altered (mo)":   fmt_med(r.get("median_altered_months")),
            "Median wildtype (mo)":  fmt_med(r.get("median_wildtype_months")),
            "Log-rank p":            fmt_p(r.get("logrank_p")),
            "HR (multivariate)":     fmt_hr(r.get("HR")),
            "95% CI":                fmt_ci(r.get("CI_lower"), r.get("CI_upper")),
            "Cox p":                 fmt_p(r.get("cox_p")),
            "Significant?":          "Yes" if overall_sig else ("No" if km_done else "Not tested"),
            # sort keys (dropped later)
            "_sig":    overall_sig,
            "_hr":     float(r["HR"]) if pd.notna(r.get("HR")) else 0.0,
            "_freq":   float(freq)    if pd.notna(freq)        else 0.0,
        })

    df = pd.DataFrame(rows)
    sig_df  = df[df["_sig"]].sort_values("_hr", ascending=False)
    nsig_df = df[~df["_sig"]].sort_values("_freq", ascending=False)
    result  = pd.concat([sig_df, nsig_df], ignore_index=True)
    return result.drop(columns=["_sig", "_hr", "_freq"])


# ── TABLE 3: Expression Results ────────────────────────────────────────────────

def build_table3(expr: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in expr.iterrows():
        gene      = r["gene"]
        completed = r["status"] == "completed"
        sig       = (
            bool(r["logrank_p"] < 0.05)
            if completed and pd.notna(r.get("logrank_p"))
            else False
        )
        n_hi = r.get("n_high")
        n_lo = r.get("n_low")
        total = (n_hi + n_lo) if (pd.notna(n_hi) and pd.notna(n_lo)) else np.nan
        freq_hi = round(n_hi / total * 100, 1) if pd.notna(total) and total > 0 else np.nan

        rows.append({
            "Gene":               star(gene),
            "Is PDMR Focus":      bool(r.get("is_pdmr_focus", gene in PDMR_GENES)),
            "Cohort":             r.get("cohort", "—"),
            "Endpoint":           r["endpoint"],
            "Z-score threshold":  r.get("z_threshold", "—"),
            "N high":             int(n_hi) if pd.notna(n_hi) else "—",
            "N low":              int(n_lo) if pd.notna(n_lo) else "—",
            "Freq high %":        f"{freq_hi:.1f}" if pd.notna(freq_hi) else "—",
            "Median HIGH (mo)":   fmt_med(r.get("median_high")) if completed else "—",
            "Median LOW (mo)":    fmt_med(r.get("median_low"))  if completed else "—",
            "Log-rank p":         fmt_p(r.get("logrank_p"))     if completed else "—",
            "HR (vs. LOW)":       fmt_hr(r.get("HR"))           if completed else "—",
            "95% CI":             fmt_ci(r.get("CI_lower"), r.get("CI_upper")) if completed else "—",
            "Cox p":              fmt_p(r.get("cox_p"))         if completed else "—",
            "Significant?":       "Yes" if sig else ("No" if completed else "Not tested (gate)"),
        })
    return pd.DataFrame(rows)


# ── PLAIN-TEXT PARAGRAPHS ──────────────────────────────────────────────────────

def generate_paragraphs(km: pd.DataFrame, cox: pd.DataFrame, expr: pd.DataFrame) -> None:
    # Merge to get HR alongside KM data for significant results
    cox_done = cox[cox["status"] == "completed"].copy()
    sig_rows = cox_done[cox_done["significant"] == True].sort_values("HR", ascending=False)

    # ── Paragraph 1: significant genomic alteration results ────────────────────
    print("=" * 80)
    print("RESULTS PARAGRAPH 1 — Significant Genomic Alteration Findings")
    print("=" * 80)

    if sig_rows.empty:
        print("No genomic alteration analyses reached statistical significance.\n")
    else:
        sentences = []
        for _, r in sig_rows.iterrows():
            gene     = r["gene"]
            pathway  = r["pathway"]
            ep       = r["endpoint"]
            ep_long  = "event-free survival (EFS)" if ep == "EFS" else "overall survival (OS)"
            med_alt  = fmt_med(r["median_altered_months"])
            med_wt   = fmt_med(r["median_wildtype_months"])
            hr       = float(r["HR"])
            ci_lo    = float(r["CI_lower"])
            ci_hi    = float(r["CI_upper"])
            p_cox    = float(r["cox_p"])
            p_lr     = km.loc[
                (km["gene"] == gene) & (km["endpoint"] == ep), "logrank_p"
            ].values
            p_lr_val = float(p_lr[0]) if len(p_lr) > 0 else np.nan

            direction = "worse" if hr > 1 else "better"
            dname     = display(gene)
            alt_word  = "pathway alteration" if "OR" in gene else "alteration"

            p_lr_str = "< 0.001" if (not np.isnan(p_lr_val) and p_lr_val < 0.001) else (
                f"= {p_lr_val:.3f}" if not np.isnan(p_lr_val) else "not available"
            )

            sentence = (
                f"{dname} ({pathway}) {alt_word} was associated with significantly "
                f"{direction} {ep_long} "
                f"(median {ep}: {med_alt} vs. {med_wt} months; "
                f"log-rank p {p_lr_str}; "
                f"HR = {hr:.2f} [95% CI: {ci_lo:.2f}–{ci_hi:.2f}])."
            )
            sentences.append(sentence)

        print("\n".join(sentences))

    # ── Paragraph 2: PDMR focus gene summary ──────────────────────────────────
    print()
    print("=" * 80)
    print("RESULTS PARAGRAPH 2 — PDMR Focus Gene Summary")
    print("=" * 80)

    def pdmr_gene_sentence(gene_key: str, gene_display: str) -> str:
        """Build the clause for one PDMR focus gene."""
        # Check if it has a significant Cox result as an individual gene
        indiv_sig = cox_done[
            (cox_done["gene"] == gene_key) & (cox_done["significant"] == True)
        ]
        if not indiv_sig.empty:
            r      = indiv_sig.sort_values("HR", ascending=False).iloc[0]
            ep     = r["endpoint"]
            hr     = float(r["HR"])
            ci_lo  = float(r["CI_lower"])
            ci_hi  = float(r["CI_upper"])
            p_cox  = float(r["cox_p"])
            p_str  = "< 0.001" if p_cox < 0.001 else f"= {p_cox:.3f}"
            return (
                f"{gene_display} alteration was associated with significantly worse survival "
                f"(HR = {hr:.2f} [95% CI: {ci_lo:.2f}–{ci_hi:.2f}], p {p_str})"
            )

        # Check if skipped — get n_altered from km_results
        km_rows = km[km["gene"].str.startswith(gene_key)]
        if not km_rows.empty:
            n_alt = km_rows.iloc[0]["n_altered"]
            skip_note = (
                f"did not reach statistical significance as an individual alteration "
                f"(n altered = {int(n_alt) if pd.notna(n_alt) else '?'}; "
                f"below the minimum-arm-size threshold for KM analysis)"
            )
            # Check expression result
            expr_rows = expr[
                (expr["gene"] == gene_key) & (expr["status"] == "completed")
            ]
            if not expr_rows.empty:
                best = expr_rows.sort_values("logrank_p").iloc[0]
                p_e  = float(best["logrank_p"])
                ep_e = best["endpoint"]
                n_hi = int(best["n_high"])
                hr_e = float(best["HR"]) if pd.notna(best.get("HR")) else None
                if p_e < 0.05:
                    direction_e = "worse" if (hr_e or 1) > 1 else "better"
                    skip_note += (
                        f"; however, high mRNA expression (z > {best['z_threshold']}) "
                        f"was associated with significantly {direction_e} {ep_e} "
                        f"(log-rank p = {p_e:.3f}, n high = {n_hi})"
                    )
                else:
                    skip_note += (
                        f"; mRNA expression trend (z > {best['z_threshold']}) "
                        f"was non-significant ({ep_e}, log-rank p = {p_e:.3f}, "
                        f"n high = {n_hi})"
                    )
            else:
                expr_skip = expr[expr["gene"] == gene_key]
                if not expr_skip.empty:
                    n_hi = expr_skip.iloc[0].get("n_high")
                    thr  = expr_skip.iloc[0].get("z_threshold")
                    skip_note += (
                        f"; expression analysis also below gate "
                        f"(n high = {int(n_hi) if pd.notna(n_hi) else '?'} "
                        f"with z > {thr})"
                    )
            return f"{gene_display} {skip_note}"

        return f"{gene_display}: no data available"

    mdm4_clause = pdmr_gene_sentence("MDM4", "MDM4")
    mdm2_clause = pdmr_gene_sentence("MDM2", "MDM2")
    igf2_clause = pdmr_gene_sentence("IGF2_CNA", "IGF2")  # IGF2_CNA is the primary IGF2 entry

    # Also note the OR-classifier containing MDM4/MDM2
    or_row = cox_done[
        (cox_done["gene"] == "TP53_MDM4_MDM2_OR") & (cox_done["endpoint"] == "OS")
    ]
    if not or_row.empty:
        r_or  = or_row.iloc[0]
        hr_or = float(r_or["HR"])
        ci_or = f"[{float(r_or['CI_lower']):.2f}–{float(r_or['CI_upper']):.2f}]"
        p_or  = float(r_or["cox_p"])
        p_or_str = "< 0.001" if p_or < 0.001 else f"= {p_or:.3f}"
        or_context = (
            f" Notably, the composite TP53/MDM4/MDM2 (OR) classifier "
            f"(n = {int(r_or['n_altered'])} patients) was significantly associated "
            f"with worse OS (HR = {hr_or:.2f} {ci_or}, p {p_or_str}), "
            f"suggesting pathway-level impact of p53/MDM axis alterations."
        )
    else:
        or_context = ""

    para2 = (
        f"Among the PDMR focus genes, {mdm4_clause}. {mdm2_clause}. "
        f"{igf2_clause}."
        f"{or_context}"
    )
    print(para2)


# ── Printing helpers ───────────────────────────────────────────────────────────

def print_table(df: pd.DataFrame, title: str) -> None:
    print(f"\n{'=' * 80}")
    print(f"  {title}")
    print(f"{'=' * 80}")
    with pd.option_context(
        "display.max_rows", None,
        "display.max_columns", None,
        "display.width", 200,
        "display.max_colwidth", 60,
    ):
        print(df.to_string(index=False))


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    os.makedirs(RESULTS_DIR, exist_ok=True)

    master = pd.read_csv(MASTER_CSV)
    km     = pd.read_csv(KM_CSV)
    cox    = pd.read_csv(COX_CSV)
    expr   = pd.read_csv(EXPR_CSV)

    # ── Build tables ───────────────────────────────────────────────────────────
    t1 = build_table1(master)
    t2 = build_table2(km, cox)
    t3 = build_table3(expr)

    # ── Save ───────────────────────────────────────────────────────────────────
    t1.to_csv(T1_CSV, index=False)
    t2.to_csv(T2_CSV, index=False)
    t3.to_csv(T3_CSV, index=False)
    print(f"Saved:\n  {os.path.abspath(T1_CSV)}\n  {os.path.abspath(T2_CSV)}\n  {os.path.abspath(T3_CSV)}\n")

    # ── Print results paragraphs ───────────────────────────────────────────────
    generate_paragraphs(km, cox, expr)

    # ── Print all three tables ─────────────────────────────────────────────────
    print_table(t1, "TABLE 1 — Cohort Summary")
    print_table(t2, "TABLE 2 — Genomic Alteration Results (KM + Cox)")
    print_table(t3, "TABLE 3 — Expression Analysis Results")

    # ── Quick sanity counts ────────────────────────────────────────────────────
    print(f"\n{'─' * 80}")
    print("SUMMARY")
    print(f"{'─' * 80}")
    print(f"  Table 1: {len(t1)} cohort row(s)")
    n_sig_t2 = (t2["Significant?"] == "Yes").sum()
    print(f"  Table 2: {len(t2)} rows; {n_sig_t2} significant analyses")
    n_done_t3 = (t3["Significant?"] != "Not tested (gate)").sum()
    n_sig_t3  = (t3["Significant?"] == "Yes").sum()
    print(f"  Table 3: {len(t3)} rows; {n_done_t3} completed; {n_sig_t3} significant")


if __name__ == "__main__":
    main()
