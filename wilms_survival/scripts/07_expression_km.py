"""
07_expression_km.py

Expression-based Kaplan-Meier + Cox survival analysis for Wilms' Tumor.

Stratifies patients into HIGH / LOW expression groups using mRNA z-score
thresholds, then runs log-rank test and multivariate Cox (adjusted for AGE)
for each gene × cohort × endpoint combination.

Thresholds
----------
  IGF2 : z > 1.5  (strict — LOI-driven overexpression is nearly universal)
  All others : z > 1.0

Gate rule (KM): minimum 20 patients per arm.

Cohorts processed dynamically from wilms_survival/data/:
  wt_target_2018_pub → labeled TARGET_2018  (n=125 with expression; OS + EFS)
  wt_target_gdc      → labeled TARGET_GDC   (n=125 with expression; OS only;
                        AGE omitted — all GDC clinical ages redacted to 18)

Outputs
-------
  wilms_survival/plots/KM_EXPR_{GENE}_{COHORT}_{ENDPOINT}.pdf
  wilms_survival/results/expression_km_results.csv
"""

import os
import re
import warnings
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────

SCRIPT_DIR  = os.path.dirname(__file__)
DATA_DIR    = os.path.join(SCRIPT_DIR, "..", "data")
PLOTS_DIR   = os.path.join(SCRIPT_DIR, "..", "plots")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "..", "results")
RESULTS_CSV = os.path.join(RESULTS_DIR, "expression_km_results.csv")

# ── Gene catalogue ─────────────────────────────────────────────────────────────

GENES_PRIMARY = ["TP53", "MYCN", "SIX1", "SIX2", "DROSHA", "DGCR8", "CTNNB1", "WT1", "MLLT1"]
GENES_PDMR    = ["MDM4", "MDM2", "IGF2"]
ALL_GENES     = GENES_PRIMARY + GENES_PDMR
PDMR_SET      = set(GENES_PDMR)

IGF2_THRESHOLD    = 1.5
DEFAULT_THRESHOLD = 1.0
MIN_N             = 20      # minimum patients per arm
AGE_MISS_LIMIT    = 0.30    # skip AGE if > 30% missing

# ── Plot colours ───────────────────────────────────────────────────────────────

ORANGE = "#E07B00"   # HIGH expression
TEAL   = "#1D9E75"   # LOW expression

# ── Endpoints ──────────────────────────────────────────────────────────────────

ALL_ENDPOINTS = [
    ("EFS", "EFS_MONTHS", "EFS_STATUS"),
    ("OS",  "OS_MONTHS",  "OS_STATUS"),
]

# Study ID → short display label for filenames and plot titles
COHORT_LABELS = {
    "wt_target_2018_pub": "TARGET_2018",
    "wt_target_gdc":      "TARGET_GDC",
}


# ── Helpers ────────────────────────────────────────────────────────────────────

def strip_sample_suffix(sid: str) -> str:
    """
    Convert a TARGET sample ID to patient ID by removing the trailing
    aliquot / vial suffix.
      TARGET-50-CAAAAB-01   → TARGET-50-CAAAAB
      TARGET-50-CAAAAB-01A  → TARGET-50-CAAAAB
    """
    return re.sub(r"-\d+[A-Z]*$", "", sid)


def km_median(kmf: KaplanMeierFitter) -> str:
    m = kmf.median_survival_time_
    try:
        if m is None or np.isinf(float(m)) or np.isnan(float(m)):
            return "NR"
        return f"{float(m):.1f}"
    except (TypeError, ValueError):
        return "NR"


def p_label(p: float) -> str:
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"


# ── Core analysis function ─────────────────────────────────────────────────────

def run_expr_km(
    gene:       str,
    gene_mrna:  pd.DataFrame,    # rows for this gene only; cols: patientId, value
    clin:       pd.DataFrame,    # full clinical table for this cohort
    threshold:  float,
    is_pdmr:    bool,
    ep_label:   str,
    time_col:   str,
    event_col:  str,
    cohort_id:  str,
    cohort_lbl: str,
    use_age:    bool,
) -> dict:
    """Run one gene × endpoint expression KM + Cox analysis."""

    base = dict(
        gene=gene, is_pdmr_focus=is_pdmr, cohort=cohort_lbl,
        endpoint=ep_label, z_threshold=threshold,
        n_high=None, n_low=None,
        median_high=None, median_low=None,
        events_high=None, events_low=None,
        logrank_p=None, HR=None, CI_lower=None, CI_upper=None,
        cox_p=None, status=None,
    )

    # ── Merge expression with survival ─────────────────────────────────────────
    clin_cols = ["patientId", time_col, event_col]
    if use_age and "AGE" in clin.columns:
        clin_cols.append("AGE")

    merged = (
        gene_mrna[["patientId", "value"]]
        .rename(columns={"value": "zscore"})
        .merge(clin[clin_cols], on="patientId", how="inner")
        .dropna(subset=["zscore", time_col, event_col])
    )

    if merged.empty:
        base["status"] = "skipped (no overlap after merge)"
        return base

    # ── Stratify ───────────────────────────────────────────────────────────────
    merged["group"] = np.where(merged["zscore"] > threshold, "HIGH", "LOW")
    high = merged[merged["group"] == "HIGH"]
    low  = merged[merged["group"] == "LOW"]
    n_high, n_low = len(high), len(low)
    base["n_high"] = n_high
    base["n_low"]  = n_low

    # ── Gate check ─────────────────────────────────────────────────────────────
    if n_high < MIN_N or n_low < MIN_N:
        base["status"] = (
            f"skipped (n_high={n_high}, n_low={n_low}; gate={MIN_N}/arm)"
        )
        print(f"    SKIP  {gene:<8} / {ep_label}:  "
              f"n_high={n_high:>3}  n_low={n_low:>3}")
        return base

    # ── Kaplan-Meier ───────────────────────────────────────────────────────────
    T_hi, E_hi = high[time_col].astype(float), high[event_col].astype(float)
    T_lo, E_lo = low[time_col].astype(float),  low[event_col].astype(float)

    kmf_hi = KaplanMeierFitter()
    kmf_lo = KaplanMeierFitter()
    kmf_hi.fit(T_hi, E_hi, label=f"High z-score  (n={n_high})")
    kmf_lo.fit(T_lo, E_lo, label=f"Low z-score   (n={n_low})")

    lr   = logrank_test(T_hi, T_lo, event_observed_A=E_hi, event_observed_B=E_lo)
    p_lr = float(lr.p_value)

    med_hi  = km_median(kmf_hi)
    med_lo  = km_median(kmf_lo)
    ev_hi   = int(E_hi.sum())
    ev_lo   = int(E_lo.sum())

    # ── Cox regression ─────────────────────────────────────────────────────────
    cox_df = merged[["zscore", time_col, event_col]].copy()
    cox_df["HIGH"] = (merged["zscore"] > threshold).astype(int)
    if use_age and "AGE" in merged.columns:
        cox_df["AGE"] = merged["AGE"].values
    cox_df = (
        cox_df.dropna()
        .rename(columns={time_col: "duration", event_col: "event"})
    )

    covs = ["HIGH"] + (["AGE"] if (use_age and "AGE" in cox_df.columns) else [])
    hr = ci_lo = ci_hi = p_cox = None
    try:
        cph = CoxPHFitter()
        cph.fit(
            cox_df[["duration", "event"] + covs],
            duration_col="duration", event_col="event",
            show_progress=False,
        )
        row   = cph.summary.loc["HIGH"]
        hr    = float(np.exp(row["coef"]))
        ci_lo = float(np.exp(row["coef lower 95%"]))
        ci_hi = float(np.exp(row["coef upper 95%"]))
        p_cox = float(row["p"])
    except Exception as exc:
        print(f"    COX ERROR {gene}/{ep_label}: {exc}")

    # ── KM plot ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.subplots_adjust(bottom=0.28)

    kmf_lo.plot_survival_function(
        ax=ax, color=TEAL, ci_show=True, ci_alpha=0.12, linewidth=2.2,
    )
    kmf_hi.plot_survival_function(
        ax=ax, color=ORANGE, ci_show=True, ci_alpha=0.12, linewidth=2.2,
    )

    add_at_risk_counts(kmf_lo, kmf_hi, ax=ax, fontsize=8.5,
                       rows_to_show=["At risk"])

    ep_long = "Event-Free Survival" if ep_label == "EFS" else "Overall Survival"
    ax.set_xlabel(f"{ep_long} (months)", fontsize=11)
    ax.set_ylabel("Survival probability", fontsize=11)
    ax.set_ylim(-0.02, 1.05)
    ax.spines[["top", "right"]].set_visible(False)

    # Main title
    thr_str  = f"z > {threshold}"
    ax.set_title(
        f"{gene} mRNA — {thr_str} vs. rest — {cohort_lbl} — {ep_long}",
        fontsize=12, fontweight="bold", pad=10,
    )

    # Subtitle — PDMR note or pathway label
    if is_pdmr:
        subtitle = f"PDMR Focus Gene — cross-reference with PDMR model expression"
    else:
        subtitle = f"mRNA z-score threshold: {thr_str}"

    fig.text(
        0.5, 0.96, subtitle,
        ha="center", fontsize=9.5, style="italic",
        color="#8B0000" if is_pdmr else "#555555",
        transform=fig.transFigure,
    )

    # Statistics annotation
    p_lr_str  = f"{p_lr:.3f}" if p_lr >= 0.001 else f"{p_lr:.2e}"
    sig_mk    = " ✱" if p_lr < 0.05 else ""
    hr_str    = (f"  |  Cox HR = {hr:.2f} [{ci_lo:.2f}–{ci_hi:.2f}], {p_label(p_cox)}"
                 if hr is not None else "")
    ann = (
        f"Median HIGH:  {med_hi} mo  |  Median LOW:  {med_lo} mo"
        f"\nLog-rank  p = {p_lr_str}{sig_mk}{hr_str}"
        f"\nEvents: HIGH {ev_hi}/{n_high}  •  LOW {ev_lo}/{n_low}"
    )
    ax.text(
        0.98, 0.97, ann,
        transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#FFFDE7",
                  edgecolor="#AAAAAA", alpha=0.92),
    )

    ax.legend(loc="lower left", fontsize=10, framealpha=0.85)

    os.makedirs(PLOTS_DIR, exist_ok=True)
    fname = os.path.join(PLOTS_DIR, f"KM_EXPR_{gene}_{cohort_lbl}_{ep_label}.pdf")
    fig.savefig(fname, dpi=300, bbox_inches="tight")
    plt.close(fig)

    p_lr_out = f"{p_lr:.3f}" if p_lr >= 0.001 else f"{p_lr:.2e}"
    hr_out   = f"HR={hr:.2f}" if hr is not None else "Cox-err"
    print(f"    DONE  {gene:<8} / {ep_label}:  "
          f"n_hi={n_high:>3}  n_lo={n_low:>3}  "
          f"p={p_lr_out:<10}  {hr_out}  "
          f"{'✱ SIG' if p_lr < 0.05 else ''}")

    base.update(
        n_high=n_high, n_low=n_low,
        events_high=ev_hi, events_low=ev_lo,
        median_high=med_hi, median_low=med_lo,
        logrank_p=round(p_lr, 6),
        HR=round(hr, 3) if hr is not None else None,
        CI_lower=round(ci_lo, 3) if ci_lo is not None else None,
        CI_upper=round(ci_hi, 3) if ci_hi is not None else None,
        cox_p=round(p_cox, 6) if p_cox is not None else None,
        status="completed",
    )
    return base


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    os.makedirs(PLOTS_DIR,   exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Discover cohorts from mrna_zscores_*.csv files
    mrna_files = sorted(
        f for f in os.listdir(DATA_DIR) if f.startswith("mrna_zscores_")
    )
    if not mrna_files:
        print("No mrna_zscores_*.csv files found in data/. Run 03_download_genomics.py first.")
        return

    results: list[dict] = []

    for mrna_file in mrna_files:
        study_id  = mrna_file.replace("mrna_zscores_", "").replace(".csv", "")
        clin_file = f"clinical_{study_id}.csv"
        clin_path = os.path.join(DATA_DIR, clin_file)
        mrna_path = os.path.join(DATA_DIR, mrna_file)

        print(f"\n{'='*70}")
        print(f"  Cohort: {study_id}")
        print(f"{'='*70}")

        if not os.path.exists(clin_path):
            print(f"  No matching clinical file ({clin_file}) — skipping")
            continue

        mrna = pd.read_csv(mrna_path)
        clin = pd.read_csv(clin_path)

        # Strip sample suffix → patient ID
        mrna = mrna.copy()
        mrna["patientId"] = mrna["sampleId"].apply(strip_sample_suffix)

        n_expr_pts = mrna["patientId"].nunique()
        genes_avail = sorted(mrna["hugoGeneSymbol"].unique())
        print(f"  Expression: {n_expr_pts} patients, {len(genes_avail)} genes: {genes_avail}")

        # ── Covariate decisions ────────────────────────────────────────────────
        age_col = "AGE" if "AGE" in clin.columns else None
        if age_col:
            age_std     = clin[age_col].std(skipna=True)
            age_missing = clin[age_col].isna().mean()
            use_age     = (age_std > 0.01) and (age_missing <= AGE_MISS_LIMIT)
        else:
            use_age = False

        print(f"  AGE covariate: {'INCLUDED' if use_age else 'OMITTED'}", end="")
        if age_col:
            print(f"  (std={age_std:.2f}, missing={age_missing*100:.0f}%)")
        else:
            print("  (column absent)")

        # ── Available endpoints ────────────────────────────────────────────────
        available_eps = []
        for ep_label, time_col, event_col in ALL_ENDPOINTS:
            n_valid = (
                clin[time_col].notna() & clin[event_col].notna()
            ).sum() if (time_col in clin.columns and event_col in clin.columns) else 0
            if n_valid > 0:
                available_eps.append((ep_label, time_col, event_col))
            else:
                print(f"  {ep_label}: no valid rows — endpoint skipped")

        print(f"  Endpoints: {[ep for ep,_,_ in available_eps]}")

        # ── Cohort label (short, safe for filenames) ───────────────────────────
        cohort_lbl = COHORT_LABELS.get(
            study_id,
            study_id.replace("wt_target_", "").upper().replace("_", ""),
        )

        # ── Per gene × endpoint ────────────────────────────────────────────────
        for gene in ALL_GENES:
            gene_mrna = mrna[mrna["hugoGeneSymbol"] == gene].copy()
            if gene_mrna.empty:
                print(f"  {gene}: not in expression data for this cohort — skipping")
                continue

            threshold = IGF2_THRESHOLD if gene == "IGF2" else DEFAULT_THRESHOLD
            is_pdmr   = gene in PDMR_SET

            # Quick distribution summary
            zvals  = gene_mrna["value"]
            n_high_preview = (zvals > threshold).sum()
            n_low_preview  = (zvals <= threshold).sum()
            print(f"\n  {gene:<8}  "
                  f"threshold={threshold}  "
                  f"n_high={n_high_preview}  n_low={n_low_preview}  "
                  f"median_z={zvals.median():.2f}  max_z={zvals.max():.2f}")

            for ep_label, time_col, event_col in available_eps:
                r = run_expr_km(
                    gene=gene, gene_mrna=gene_mrna, clin=clin,
                    threshold=threshold, is_pdmr=is_pdmr,
                    ep_label=ep_label, time_col=time_col, event_col=event_col,
                    cohort_id=study_id, cohort_lbl=cohort_lbl,
                    use_age=use_age,
                )
                results.append(r)

    # ── Save results ───────────────────────────────────────────────────────────
    results_df = pd.DataFrame(results)

    completed = (
        results_df[results_df["status"] == "completed"]
        .sort_values("logrank_p")
    )
    skipped = results_df[results_df["status"] != "completed"]
    out_df  = pd.concat([completed, skipped], ignore_index=True)
    out_df.to_csv(RESULTS_CSV, index=False)
    print(f"\nSaved → {os.path.abspath(RESULTS_CSV)}")

    # ── Summary table ──────────────────────────────────────────────────────────
    print("\n" + "=" * 120)
    print("EXPRESSION KM RESULTS  (sorted by log-rank p;  ★ = PDMR focus gene)")
    print("=" * 120)

    hdr = (f"{'Gene':<8} {'Coh':<12} {'EP':<4} {'Thr':>4} "
           f"{'N_hi':>5} {'N_lo':>5} {'Ev_hi':>6} {'Ev_lo':>6} "
           f"{'Med_hi':>7} {'Med_lo':>7} "
           f"{'LR_p':>10} {'HR':>6} {'CI_lower':>8} {'CI_upper':>8} {'Cox_p':>10}  Status")
    print(hdr)
    print("-" * 120)

    for _, row in out_df.iterrows():
        star = "★" if row.get("is_pdmr_focus") else " "
        if row["status"] != "completed":
            short = row["status"].split("(")[0].strip()
            print(f"{star+' '+row['gene']:<8} {str(row['cohort']):<12} "
                  f"{row['endpoint']:<4} {'—':>4} {'—':>5} {'—':>5} "
                  f"{'—':>6} {'—':>6} {'—':>7} {'—':>7} "
                  f"{'—':>10} {'—':>6} {'—':>8} {'—':>8} {'—':>10}  {short}")
        else:
            lr_str  = (f"{row['logrank_p']:.3f}" if pd.notna(row['logrank_p'])
                       and row['logrank_p'] >= 0.001 else
                       (f"{row['logrank_p']:.2e}" if pd.notna(row['logrank_p']) else "—"))
            cox_str = (f"{row['cox_p']:.3f}" if pd.notna(row['cox_p'])
                       and row['cox_p'] >= 0.001 else
                       (f"{row['cox_p']:.2e}" if pd.notna(row['cox_p']) else "—"))
            sig_mk  = " ✱" if pd.notna(row['logrank_p']) and row['logrank_p'] < 0.05 else "  "
            hr_str  = f"{row['HR']:.2f}" if pd.notna(row.get('HR')) else "—"
            cil     = f"{row['CI_lower']:.2f}" if pd.notna(row.get('CI_lower')) else "—"
            ciu     = f"{row['CI_upper']:.2f}" if pd.notna(row.get('CI_upper')) else "—"
            print(f"{star+' '+row['gene']:<8} {str(row['cohort']):<12} "
                  f"{row['endpoint']:<4} {row['z_threshold']:>4.1f} "
                  f"{int(row['n_high']):>5} {int(row['n_low']):>5} "
                  f"{int(row['events_high']):>6} {int(row['events_low']):>6} "
                  f"{str(row['median_high']):>7} {str(row['median_low']):>7} "
                  f"{lr_str:>10}{sig_mk} {hr_str:>6} {cil:>8} {ciu:>8} {cox_str:>10}  completed")

    print("=" * 120)
    n_done = (out_df["status"] == "completed").sum()
    n_skip = (out_df["status"] != "completed").sum()
    print(f"\nCompleted: {n_done}  |  Skipped (gate < {MIN_N}/arm): {n_skip}")
    if n_done > 0:
        n_sig = (
            out_df[out_df["status"] == "completed"]["logrank_p"] < 0.05
        ).sum()
        print(f"Significant (log-rank p < 0.05): {n_sig} / {n_done} completed analyses")

    # ── List generated PDFs ────────────────────────────────────────────────────
    pdfs = sorted(f for f in os.listdir(PLOTS_DIR) if f.startswith("KM_EXPR_"))
    if pdfs:
        print(f"\nGenerated plots ({len(pdfs)}):")
        for pdf in pdfs:
            print(f"  {pdf}")


if __name__ == "__main__":
    main()
