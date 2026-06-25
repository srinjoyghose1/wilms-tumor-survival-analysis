"""
05_kaplan_meier.py

Kaplan-Meier survival analysis for every pathway gene and OR-classifier in the
Wilms' Tumor (TARGET, 2018) cohort (wt_target_2018_pub).

Gate rule  : skip any gene × endpoint where EITHER arm has fewer than 20 patients.
Endpoints  : EFS (primary) and OS (secondary).
Plot colors: Altered = #C00000 (red), Wildtype = #1F4E79 (dark blue).

Expected outcome (from arm-size probe):
  All 16 individual gene analyses are SKIPPED (max n_altered = 18 for MYCN).
  OR-classifiers that PASS gate:
    TP53_MDM4_MDM2_OR   21 altered / 631 WT
    MYCN_SIX_OR         23 altered / 629 WT
    ALL_PRIMARY_OR      44 altered / 608 WT

Outputs
-------
  wilms_survival/plots/KM_{GENE}_{ENDPOINT}.pdf   (one per completed analysis)
  wilms_survival/results/km_results.csv            (all genes, status column)
"""

import os
import warnings
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")           # non-interactive backend — must precede pyplot import
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts

warnings.filterwarnings("ignore")   # suppress lifelines convergence chatter

# ── Paths ─────────────────────────────────────────────────────────────────────

SCRIPT_DIR  = os.path.dirname(__file__)
PROC_DIR    = os.path.join(SCRIPT_DIR, "..", "processed")
PLOTS_DIR   = os.path.join(SCRIPT_DIR, "..", "plots")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "..", "results")
MASTER_CSV  = os.path.join(PROC_DIR, "master_classification.csv")
RESULTS_CSV = os.path.join(RESULTS_DIR, "km_results.csv")

# ── Constants ─────────────────────────────────────────────────────────────────

RED    = "#C00000"
BLUE   = "#1F4E79"
MIN_N  = 20          # minimum patients per arm

ENDPOINTS = [
    ("EFS", "EFS_MONTHS", "EFS_STATUS"),
    ("OS",  "OS_MONTHS",  "OS_STATUS"),
]

# ── Gene catalogue ─────────────────────────────────────────────────────────────
# (display_name, column_in_master_csv, pathway_label, is_pdmr_focus_gene)

INDIVIDUAL_GENES = [
    ("TP53",      "TP53_ALTERED",       "Cell Cycle",                    False),
    ("MYCN",      "MYCN_ALTERED",       "Transcription Reg",              False),
    ("SIX1_SIX2", "SIX1_SIX2_ALTERED", "Transcription Reg (combined)",   False),
    ("SIX1",      "SIX1_ALTERED",       "Transcription Reg",              False),
    ("SIX2",      "SIX2_ALTERED",       "Transcription Reg",              False),
    ("DROSHA",    "DROSHA_ALTERED",     "RNA Processing",                 False),
    ("DGCR8",     "DGCR8_ALTERED",      "RNA Processing",                 False),
    ("CTNNB1",    "CTNNB1_ALTERED",     "WNT Signaling",                  False),
    ("AMER1",     "AMER1_ALTERED",      "WNT Signaling",                  False),
    ("WT1",       "WT1_ALTERED",        "WNT Signaling",                  False),
    ("MLLT1",     "MLLT1_ALTERED",      "Chromatin Remodeling",           False),
    ("NIPBL",     "NIPBL_ALTERED",      "PanCancer",                      False),
    ("MDM4",      "MDM4_ALTERED",       "PDMR Focus (p53 pathway)",       True),
    ("MDM2",      "MDM2_ALTERED",       "PDMR Focus (p53 pathway)",       True),
    ("IGF2_CNA",  "IGF2_CNA_ALTERED",   "PDMR Focus (imprinting/CNA)",   True),
    ("IGF2_EXPR", "IGF2_EXPR_ALTERED",  "PDMR Focus (imprinting/expr)",  True),
]

OR_CLASSIFIERS = [
    {
        "name":    "TP53_MDM4_MDM2_OR",
        "label":   "Cell Cycle + p53 Axis",
        "cols":    ["TP53_ALTERED", "MDM4_ALTERED", "MDM2_ALTERED"],
        "pdmr":    False,
    },
    {
        "name":    "DROSHA_DGCR8_OR",
        "label":   "RNA Processing Axis",
        "cols":    ["DROSHA_ALTERED", "DGCR8_ALTERED"],
        "pdmr":    False,
    },
    {
        "name":    "MYCN_SIX_OR",
        "label":   "Transcription + Progenitor Axis",
        "cols":    ["MYCN_ALTERED", "SIX1_SIX2_ALTERED"],
        "pdmr":    False,
    },
    {
        "name":    "WNT_AXIS_OR",
        "label":   "WNT Axis",
        "cols":    ["CTNNB1_ALTERED", "AMER1_ALTERED", "WT1_ALTERED"],
        "pdmr":    False,
    },
    {
        "name":    "ALL_PRIMARY_OR",
        "label":   "All Primary Pathway Genes",
        "cols":    [
            "TP53_ALTERED", "MYCN_ALTERED",
            "SIX1_ALTERED", "SIX2_ALTERED",
            "DROSHA_ALTERED", "DGCR8_ALTERED",
            "CTNNB1_ALTERED", "AMER1_ALTERED", "WT1_ALTERED",
            "MLLT1_ALTERED", "NIPBL_ALTERED",
        ],
        "pdmr":    False,
    },
]


# ── Helpers ────────────────────────────────────────────────────────────────────

def compute_or_flag(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    """
    Logical OR across multiple alteration columns, preserving NA semantics:
      1   — at least one source is 1
      0   — at least one source is 0 and none is 1
      NaN — ALL sources are NaN (patient unsequenced for every contributing gene)
    """
    result = pd.Series(pd.NA, index=df.index, dtype="Float64")
    for col in cols:
        s = df[col].astype("Float64")
        has_data = s.notna() | result.notna()
        result.loc[has_data] = result.loc[has_data].fillna(0.0)
        result.loc[(s == 1).fillna(False)] = 1.0
    return result


def median_label(kmf: KaplanMeierFitter) -> str:
    """Return median survival as 'X.X' mo string, or 'NR' if not reached."""
    m = kmf.median_survival_time_
    try:
        if m is None or np.isinf(float(m)) or np.isnan(float(m)):
            return "NR"
        return f"{float(m):.1f}"
    except (TypeError, ValueError):
        return "NR"


def endpoint_full(ep_label: str) -> str:
    return "Event-Free Survival" if ep_label == "EFS" else "Overall Survival"


# ── Core KM function ───────────────────────────────────────────────────────────

def run_km(
    gene: str,
    flag: pd.Series,
    pathway: str,
    ep_label: str,
    time_col: str,
    event_col: str,
    df: pd.DataFrame,
    is_pdmr: bool = False,
) -> dict:
    """
    Run KM analysis for one gene × endpoint.
    Always returns a result dict (status = 'completed' or 'skipped').
    Saves a PDF plot only when completed.
    """
    base = {
        "gene":                  gene,
        "pathway":               pathway,
        "pdmr_focus":            is_pdmr,
        "endpoint":              ep_label,
        "n_altered":             None,
        "n_wildtype":            None,
        "events_altered":        None,
        "events_wildtype":       None,
        "median_altered_months": None,
        "median_wildtype_months":None,
        "logrank_p":             None,
        "significant":           None,
        "status":                None,
    }

    # Merge flag with survival columns; drop rows with any missing value
    sub = df[[time_col, event_col]].copy()
    sub["flag"] = flag.values
    sub = sub.dropna(subset=[time_col, event_col, "flag"])

    alt = sub[sub["flag"] == 1]
    wt  = sub[sub["flag"] == 0]
    n_alt, n_wt = len(alt), len(wt)

    base["n_altered"]  = n_alt
    base["n_wildtype"] = n_wt

    # ── Gate check ─────────────────────────────────────────────────────────────
    if n_alt < MIN_N or n_wt < MIN_N:
        base["status"] = f"skipped (n_altered={n_alt}, n_wildtype={n_wt}; gate={MIN_N})"
        print(f"  SKIP  {gene:<22} / {ep_label}:  "
              f"n_alt={n_alt:>3}  n_wt={n_wt:>3}  — gate threshold {MIN_N}/arm")
        return base

    # ── Fit KM ─────────────────────────────────────────────────────────────────
    T_alt = alt[time_col].astype(float)
    E_alt = alt[event_col].astype(float)
    T_wt  = wt[time_col].astype(float)
    E_wt  = wt[event_col].astype(float)

    kmf_alt = KaplanMeierFitter()
    kmf_wt  = KaplanMeierFitter()
    kmf_alt.fit(T_alt, E_alt, label=f"Altered  (n={n_alt})")
    kmf_wt.fit( T_wt,  E_wt,  label=f"Wildtype (n={n_wt})")

    lr  = logrank_test(T_alt, T_wt, event_observed_A=E_alt, event_observed_B=E_wt)
    p   = float(lr.p_value)
    med_alt = median_label(kmf_alt)
    med_wt  = median_label(kmf_wt)
    ev_alt  = int(E_alt.sum())
    ev_wt   = int(E_wt.sum())

    # ── Plot ───────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.subplots_adjust(bottom=0.28)   # room for at-risk table

    kmf_wt.plot_survival_function(
        ax=ax, color=BLUE, ci_show=True, ci_alpha=0.12, linewidth=2.2,
    )
    kmf_alt.plot_survival_function(
        ax=ax, color=RED,  ci_show=True, ci_alpha=0.12, linewidth=2.2,
    )

    # At-risk counts table (WT first so colours match legend top-to-bottom)
    add_at_risk_counts(kmf_wt, kmf_alt, ax=ax, fontsize=8.5,
                       rows_to_show=["At risk"])

    # Axis labels and limits
    ep_long = endpoint_full(ep_label)
    ax.set_xlabel(f"{ep_long} (months)", fontsize=11)
    ax.set_ylabel("Survival probability", fontsize=11)
    ax.set_ylim(-0.02, 1.05)
    ax.spines[["top", "right"]].set_visible(False)

    # Title
    pdmr_tag = "  [★ PDMR Focus Gene]" if is_pdmr else ""
    ax.set_title(
        f"{gene} — {ep_long} — Wilms Tumor (TARGET 2018){pdmr_tag}",
        fontsize=13, fontweight="bold", pad=10,
    )

    # Pathway subtitle
    fig.text(
        0.5, 0.96,
        f"Pathway: {pathway}",
        ha="center", fontsize=9.5, style="italic", color="#555555",
        transform=fig.transFigure,
    )

    # Statistics annotation box
    p_str  = f"{p:.3f}" if p >= 0.001 else f"{p:.2e}"
    sig_mk = " ✱" if p < 0.05 else ""
    ann    = (
        f"Median altered:  {med_alt} mo  |  Median WT:  {med_wt} mo"
        f"\nLog-rank  p = {p_str}{sig_mk}"
        f"\nEvents: altered {ev_alt}/{n_alt}  •  WT {ev_wt}/{n_wt}"
    )
    ax.text(
        0.98, 0.97, ann,
        transform=ax.transAxes, ha="right", va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#FFFDE7",
                  edgecolor="#AAAAAA", alpha=0.92),
    )

    # Legend
    ax.legend(loc="lower left", fontsize=10, framealpha=0.85)

    # Save
    os.makedirs(PLOTS_DIR, exist_ok=True)
    fname = os.path.join(PLOTS_DIR, f"KM_{gene}_{ep_label}.pdf")
    fig.savefig(fname, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  DONE  {gene:<22} / {ep_label}:  "
          f"n_alt={n_alt:>3}  n_wt={n_wt:>3}  "
          f"p={p_str:<10}  {'★ SIGNIFICANT' if p<0.05 else ''}")

    base.update({
        "n_altered":              n_alt,
        "n_wildtype":             n_wt,
        "events_altered":         ev_alt,
        "events_wildtype":        ev_wt,
        "median_altered_months":  med_alt,
        "median_wildtype_months": med_wt,
        "logrank_p":              round(p, 6),
        "significant":            p < 0.05,
        "status":                 "completed",
    })
    return base


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    os.makedirs(PLOTS_DIR,   exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    df = pd.read_csv(MASTER_CSV)
    print(f"Loaded master_classification.csv  ({len(df)} patients)\n")

    results: list[dict] = []

    # ── Individual gene analyses ───────────────────────────────────────────────
    print("=" * 70)
    print("INDIVIDUAL GENE ANALYSES")
    print("=" * 70)

    for gene, col, pathway, is_pdmr in INDIVIDUAL_GENES:
        flag = df[col].astype("Float64")
        for ep_label, time_col, event_col in ENDPOINTS:
            r = run_km(gene, flag, pathway, ep_label, time_col, event_col, df, is_pdmr)
            results.append(r)

    # ── OR-classifier analyses ─────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("OR-CLASSIFIER ANALYSES")
    print("=" * 70)

    for cfg in OR_CLASSIFIERS:
        flag = compute_or_flag(df, cfg["cols"])
        n1 = int((flag == 1).sum())
        n0 = int((flag == 0).sum())

        # For ALL_PRIMARY_OR: report cohort capture percentage
        if cfg["name"] == "ALL_PRIMARY_OR":
            pct = round(n1 / (n1 + n0) * 100, 1)
            print(f"\n  ALL_PRIMARY_OR captures {n1}/{n1+n0} patients = {pct}% of cohort")
            print(f"  Contributing genes: {', '.join(cfg['cols'])}")

        print(f"\n  {cfg['name']}  (altered={n1}, wt={n0})")
        for ep_label, time_col, event_col in ENDPOINTS:
            r = run_km(
                cfg["name"], flag, cfg["label"],
                ep_label, time_col, event_col, df, cfg["pdmr"],
            )
            results.append(r)

    # ── Save results ───────────────────────────────────────────────────────────
    results_df = pd.DataFrame(results)

    # Sort: completed first (by p-value), skipped last
    completed = results_df[results_df["status"] == "completed"].sort_values("logrank_p")
    skipped   = results_df[results_df["status"] != "completed"]
    results_df = pd.concat([completed, skipped], ignore_index=True)
    results_df.to_csv(RESULTS_CSV, index=False)
    print(f"\nSaved → {os.path.abspath(RESULTS_CSV)}")

    # ── Print summary table ────────────────────────────────────────────────────
    print("\n" + "=" * 100)
    print("FULL KM RESULTS TABLE  (sorted by log-rank p;  ★ = PDMR focus gene)")
    print("=" * 100)

    hdr = (f"{'Gene':<24} {'EP':<4} {'N_alt':>5} {'N_wt':>5} "
           f"{'Ev_alt':>6} {'Ev_wt':>5} {'Med_alt':>8} {'Med_wt':>7} "
           f"{'p-value':>10}  {'Sig':>3}  Status")
    print(hdr)
    print("-" * 100)

    for _, row in results_df.iterrows():
        if row["status"] != "completed":
            status_short = "SKIP"
            line = (f"{'  ' + row['gene']:<24} {row['endpoint']:<4} "
                    f"{int(row['n_altered']) if pd.notna(row['n_altered']) else '?':>5} "
                    f"{int(row['n_wildtype']) if pd.notna(row['n_wildtype']) else '?':>5}  "
                    f"{'—':>6} {'—':>5} {'—':>8} {'—':>7} {'—':>10}  {'—':>3}  {status_short}")
        else:
            star    = "★" if row["pdmr_focus"] else " "
            sig_mk  = "✱" if row["significant"] else " "
            p_str   = (f"{row['logrank_p']:.3f}" if row["logrank_p"] >= 0.001
                       else f"{row['logrank_p']:.2e}")
            line = (f"{star + ' ' + row['gene']:<24} {row['endpoint']:<4} "
                    f"{int(row['n_altered']):>5} {int(row['n_wildtype']):>5} "
                    f"{int(row['events_altered']):>6} {int(row['events_wildtype']):>5} "
                    f"{row['median_altered_months']:>8} {row['median_wildtype_months']:>7} "
                    f"{p_str:>10}  {sig_mk:>3}  completed")
        print(line)

    print("=" * 100)
    print(f"\nPlots saved to: {os.path.abspath(PLOTS_DIR)}")
    n_done = (results_df["status"] == "completed").sum()
    n_skip = (results_df["status"] != "completed").sum()
    print(f"Completed: {n_done}  |  Skipped (gate < {MIN_N}/arm): {n_skip}")


if __name__ == "__main__":
    main()
