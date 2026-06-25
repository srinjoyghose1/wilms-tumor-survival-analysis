"""
06_cox_regression.py

Multivariate Cox Proportional Hazards regression for Wilms' Tumor (TARGET, 2018).

Genes analyzed
--------------
  • Every gene/OR-classifier from km_results.csv that passed the KM gate
    (status == "completed", n ≥ 20/arm).
  • All 5 OR-classifiers from script 05 (attempt regardless of KM gate).
  Union of those two sets is used.

Covariates per model
--------------------
  ALTERED      : 1 = altered, 0 = wildtype  (primary predictor)
  COHORT_CODE  : integer covariate (omitted if constant — single study)
  AGE          : patient age in years (omitted if > 30% missing)

Gate rule: minimum 10 events in BOTH altered AND wildtype arm.

Outputs
-------
  wilms_survival/results/cox_results.csv
"""

import os
import warnings
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────

SCRIPT_DIR  = os.path.dirname(__file__)
PROC_DIR    = os.path.join(SCRIPT_DIR, "..", "processed")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "..", "results")
MASTER_CSV  = os.path.join(PROC_DIR,    "master_classification.csv")
KM_CSV      = os.path.join(RESULTS_DIR, "km_results.csv")
COX_CSV     = os.path.join(RESULTS_DIR, "cox_results.csv")

# ── Constants ─────────────────────────────────────────────────────────────────

MIN_EVENTS  = 10   # minimum events per arm for Cox gate
AGE_MISSING_THRESHOLD = 0.30

ENDPOINTS = [
    ("EFS", "EFS_MONTHS", "EFS_STATUS"),
    ("OS",  "OS_MONTHS",  "OS_STATUS"),
]

PDMR_GENES = {"MDM4", "MDM2", "IGF2"}   # trigger PDMR notes

# ── OR-classifier catalogue ────────────────────────────────────────────────────

OR_CLASSIFIERS = {
    "TP53_MDM4_MDM2_OR": {
        "cols":    ["TP53_ALTERED", "MDM4_ALTERED", "MDM2_ALTERED"],
        "pathway": "Cell Cycle + p53 Axis",
    },
    "DROSHA_DGCR8_OR": {
        "cols":    ["DROSHA_ALTERED", "DGCR8_ALTERED"],
        "pathway": "RNA Processing Axis",
    },
    "MYCN_SIX_OR": {
        "cols":    ["MYCN_ALTERED", "SIX1_SIX2_ALTERED"],
        "pathway": "Transcription + Progenitor Axis",
    },
    "WNT_AXIS_OR": {
        "cols":    ["CTNNB1_ALTERED", "AMER1_ALTERED", "WT1_ALTERED"],
        "pathway": "WNT Axis",
    },
    "ALL_PRIMARY_OR": {
        "cols": [
            "TP53_ALTERED", "MYCN_ALTERED",
            "SIX1_ALTERED", "SIX2_ALTERED",
            "DROSHA_ALTERED", "DGCR8_ALTERED",
            "CTNNB1_ALTERED", "AMER1_ALTERED", "WT1_ALTERED",
            "MLLT1_ALTERED", "NIPBL_ALTERED",
        ],
        "pathway": "All Primary Pathway Genes",
    },
}


# ── Helpers ────────────────────────────────────────────────────────────────────

def compute_or_flag(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    """
    NA-aware logical OR: 1 if any source=1; 0 if any non-NaN and none=1;
    NaN if ALL sources are NaN.
    """
    result = pd.Series(pd.NA, index=df.index, dtype="Float64")
    for col in cols:
        s = df[col].astype("Float64")
        has_data = s.notna() | result.notna()
        result.loc[has_data] = result.loc[has_data].fillna(0.0)
        result.loc[(s == 1).fillna(False)] = 1.0
    return result


def km_median(T: pd.Series, E: pd.Series) -> str:
    """Return median survival as 'X.X' string, or 'NR' if not reached."""
    if len(T) == 0:
        return "NR"
    kmf = KaplanMeierFitter()
    kmf.fit(T.astype(float), E.astype(float))
    m = kmf.median_survival_time_
    try:
        if m is None or np.isinf(float(m)) or np.isnan(float(m)):
            return "NR"
        return f"{float(m):.1f}"
    except (TypeError, ValueError):
        return "NR"


def p_label(p: float) -> str:
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}"


def format_result(hr: float, ci_lo: float, ci_hi: float, p: float) -> str:
    return f"HR = {hr:.2f} [{ci_lo:.2f}–{ci_hi:.2f}], {p_label(p)}"


def is_pdmr_relevant(gene: str) -> bool:
    """True if the gene name or OR-classifier contains any PDMR focus gene."""
    return any(g in gene for g in PDMR_GENES)


# ── Core Cox function ──────────────────────────────────────────────────────────

def run_cox(
    gene: str,
    flag: pd.Series,
    pathway: str,
    ep_label: str,
    time_col: str,
    event_col: str,
    df: pd.DataFrame,
    use_age: bool,
    use_cohort: bool,
) -> dict:
    """
    Fit Cox PH for one gene × endpoint.
    Returns a result dict (status = 'completed' or 'skipped').
    """
    base = {
        "gene":                   gene,
        "pathway":                pathway,
        "pdmr_focus":             is_pdmr_relevant(gene),
        "endpoint":               ep_label,
        "n_altered":              None,
        "n_wildtype":             None,
        "events_altered":         None,
        "events_wildtype":        None,
        "median_altered_months":  None,
        "median_wildtype_months": None,
        "HR":                     None,
        "CI_lower":               None,
        "CI_upper":               None,
        "cox_p":                  None,
        "significant":            None,
        "covariates_used":        None,
        "formatted_result":       None,
        "status":                 None,
    }

    # Assemble analysis subset — drop rows missing time, event, or flag
    sub = df[[time_col, event_col]].copy()
    sub["ALTERED"] = flag.values
    if use_age:
        sub["AGE"] = df["AGE"].values
    if use_cohort:
        sub["COHORT_CODE"] = df["COHORT_CODE"].values
    sub = sub.dropna()

    alt = sub[sub["ALTERED"] == 1]
    wt  = sub[sub["ALTERED"] == 0]
    n_alt, n_wt = len(alt), len(wt)
    ev_alt = int(alt[event_col].sum())
    ev_wt  = int(wt[event_col].sum())

    base["n_altered"]  = n_alt
    base["n_wildtype"] = n_wt

    # ── Gate check ─────────────────────────────────────────────────────────────
    if ev_alt < MIN_EVENTS or ev_wt < MIN_EVENTS:
        base["status"] = (
            f"skipped (events_altered={ev_alt}, events_wildtype={ev_wt}; "
            f"gate={MIN_EVENTS} events/arm)"
        )
        print(f"  SKIP  {gene:<22} / {ep_label}:  "
              f"ev_alt={ev_alt:>3}  ev_wt={ev_wt:>3}  — gate {MIN_EVENTS} events/arm")
        return base

    # ── Median survival for result sentences ───────────────────────────────────
    med_alt = km_median(alt[time_col], alt[event_col])
    med_wt  = km_median(wt[time_col],  wt[event_col])
    base["events_altered"]         = ev_alt
    base["events_wildtype"]        = ev_wt
    base["median_altered_months"]  = med_alt
    base["median_wildtype_months"] = med_wt

    # ── Build Cox dataframe ────────────────────────────────────────────────────
    covs = ["ALTERED"]
    if use_age:
        covs.append("AGE")
    if use_cohort:
        covs.append("COHORT_CODE")
    base["covariates_used"] = "+".join(covs)

    cox_df = sub[[time_col, event_col] + covs].copy()
    cox_df = cox_df.rename(columns={time_col: "duration", event_col: "event"})

    # ── Fit ────────────────────────────────────────────────────────────────────
    cph = CoxPHFitter()
    try:
        cph.fit(cox_df, duration_col="duration", event_col="event", show_progress=False)
    except Exception as exc:
        base["status"] = f"error: {exc}"
        print(f"  ERROR {gene:<22} / {ep_label}: {exc}")
        return base

    print(f"\n{'─'*60}")
    print(f"  {gene} / {ep_label}  "
          f"(n_alt={n_alt}, n_wt={n_wt}; ev_alt={ev_alt}, ev_wt={ev_wt})")
    print(f"  Covariates: {', '.join(covs)}")
    cph.print_summary(decimals=3, style="ascii")

    # ── Extract ALTERED row ────────────────────────────────────────────────────
    summary = cph.summary
    if "ALTERED" not in summary.index:
        base["status"] = "error: ALTERED not in summary index"
        return base

    row     = summary.loc["ALTERED"]
    hr      = float(np.exp(row["coef"]))
    ci_lo   = float(np.exp(row["coef lower 95%"]))
    ci_hi   = float(np.exp(row["coef upper 95%"]))
    p       = float(row["p"])

    direction = "worse" if hr > 1 else "better"
    formatted = format_result(hr, ci_lo, ci_hi, p)

    base.update({
        "HR":               round(hr,   3),
        "CI_lower":         round(ci_lo, 3),
        "CI_upper":         round(ci_hi, 3),
        "cox_p":            round(p,    6),
        "significant":      p < 0.05,
        "formatted_result": formatted,
        "status":           "completed",
    })

    # ── Result sentence ────────────────────────────────────────────────────────
    if p < 0.05:
        ep_long   = "event-free survival" if ep_label == "EFS" else "overall survival"
        alt_type  = "alteration" if "OR" not in gene else "pathway alteration"
        med_str   = f"{med_alt} vs. {med_wt} months"
        p_str_lrk = "(see KM results)"    # logrank p not re-run here
        result_sentence = (
            f"\n  RESULT: {gene} {alt_type} was associated with significantly "
            f"{direction} {ep_long}\n"
            f"    (median {ep_label}: {med_str}; {formatted})."
        )
        print(result_sentence)

        # PDMR note
        for pdmr_gene in PDMR_GENES:
            if pdmr_gene in gene:
                print(
                    f"\n  PDMR NOTE — {pdmr_gene}: HR = {hr:.2f}. "
                    f"Consider querying PDMR for {pdmr_gene}-altered nephroblastoma\n"
                    f"  models to validate this survival signal in vitro."
                )

    return base


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    os.makedirs(RESULTS_DIR, exist_ok=True)

    df = pd.read_csv(MASTER_CSV)
    km = pd.read_csv(KM_CSV)

    print(f"Loaded master_classification.csv  ({len(df)} patients)\n")

    # ── Decide covariates ──────────────────────────────────────────────────────
    age_col     = "AGE" if "AGE" in df.columns else "AGE_YEARS"
    age_missing = df[age_col].isna().mean() if age_col in df.columns else 1.0
    cohort_var  = df["COHORT_CODE"].nunique() > 1

    use_age    = age_col in df.columns and age_missing <= AGE_MISSING_THRESHOLD
    use_cohort = cohort_var

    print("=== COVARIATE DECISIONS ===")
    if use_age:
        print(f"  AGE          : INCLUDED  "
              f"({age_missing*100:.1f}% missing, threshold ≤{AGE_MISSING_THRESHOLD*100:.0f}%)")
    else:
        reason = (f"{age_missing*100:.1f}% missing > {AGE_MISSING_THRESHOLD*100:.0f}% threshold"
                  if age_col in df.columns else "column not found")
        print(f"  AGE          : OMITTED   ({reason})")

    if use_cohort:
        print(f"  COHORT_CODE  : INCLUDED  ({df['COHORT_CODE'].nunique()} unique values)")
    else:
        print(f"  COHORT_CODE  : OMITTED   (constant = {df['COHORT_CODE'].unique()[0]}; "
              f"single-study cohort — zero variance)")

    # ── Build gene catalogue ───────────────────────────────────────────────────

    # KM-passed individual genes (columns already in master_classification.csv)
    km_genes_passed = (
        km[km["status"] == "completed"]["gene"]
        .unique()
        .tolist()
    )
    # Remove OR-classifiers from the individual-gene set (handled separately)
    indiv_genes_passed = [g for g in km_genes_passed if g not in OR_CLASSIFIERS]

    # Map gene display name → master_classification column name
    # Derived by appending _ALTERED; if column doesn't exist, skip
    indiv_catalogue: list[tuple[str, str, str]] = []   # (name, col, pathway)
    km_pathway = km.drop_duplicates("gene").set_index("gene")["pathway"].to_dict()
    for g in indiv_genes_passed:
        col = f"{g}_ALTERED"
        if col not in df.columns:
            print(f"  WARNING: column '{col}' not found in master_classification.csv — skipping {g}")
            continue
        indiv_catalogue.append((g, col, km_pathway.get(g, "Unknown")))

    print(f"\n=== GENE CATALOGUE ===")
    print(f"  Individual genes from KM: {len(indiv_catalogue)}  "
          f"{[g for g,_,_ in indiv_catalogue]}")
    print(f"  OR-classifiers          : {len(OR_CLASSIFIERS)}  {list(OR_CLASSIFIERS.keys())}")

    # ── Analysis loop ──────────────────────────────────────────────────────────
    results: list[dict] = []

    # --- Individual genes ---
    if indiv_catalogue:
        print("\n" + "=" * 70)
        print("INDIVIDUAL GENE ANALYSES")
        print("=" * 70)
        for gene, col, pathway in indiv_catalogue:
            flag = df[col].astype("Float64")
            for ep_label, time_col, event_col in ENDPOINTS:
                r = run_cox(gene, flag, pathway, ep_label,
                            time_col, event_col, df, use_age, use_cohort)
                results.append(r)
    else:
        print("\n  (No individual genes passed KM gate — all skipped at KM stage.)\n")

    # --- OR-classifiers ---
    print("\n" + "=" * 70)
    print("OR-CLASSIFIER ANALYSES")
    print("=" * 70)
    for or_name, cfg in OR_CLASSIFIERS.items():
        flag = compute_or_flag(df, cfg["cols"])
        n1   = int((flag == 1).sum())
        n0   = int((flag == 0).sum())
        print(f"\n  {or_name}  (altered={n1}, wt={n0})")
        for ep_label, time_col, event_col in ENDPOINTS:
            r = run_cox(or_name, flag, cfg["pathway"], ep_label,
                        time_col, event_col, df, use_age, use_cohort)
            results.append(r)

    # ── Save results ───────────────────────────────────────────────────────────
    results_df = pd.DataFrame(results)

    # Sort: completed by HR desc, skipped last
    completed = results_df[results_df["status"] == "completed"].sort_values("HR", ascending=False)
    skipped   = results_df[results_df["status"] != "completed"]
    results_df = pd.concat([completed, skipped], ignore_index=True)
    results_df.to_csv(COX_CSV, index=False)
    print(f"\nSaved → {os.path.abspath(COX_CSV)}")

    # ── Print summary table ────────────────────────────────────────────────────
    print("\n" + "=" * 110)
    print("COX RESULTS TABLE  (sorted by HR descending;  ★ = PDMR focus gene)")
    print("=" * 110)

    hdr = (f"{'Gene':<24} {'EP':<4} {'Ev_alt':>6} {'Ev_wt':>6} "
           f"{'Med_alt':>8} {'Med_wt':>7} "
           f"{'HR':>6} {'CI_lower':>8} {'CI_upper':>8} {'p-value':>10}  Sig  Status")
    print(hdr)
    print("-" * 110)

    for _, row in results_df.iterrows():
        star = "★" if row["pdmr_focus"] else " "
        if row["status"] != "completed":
            status_short = row["status"].split("(")[0].strip()
            ev_alt = f"{int(row['n_altered']):>3}" if pd.notna(row["n_altered"]) else "?"
            ev_wt  = f"{int(row['n_wildtype']):>3}" if pd.notna(row["n_wildtype"]) else "?"
            line = (
                f"{star+' '+row['gene']:<24} {row['endpoint']:<4} "
                f"{'—':>6} {'—':>6} {'—':>8} {'—':>7} "
                f"{'—':>6} {'—':>8} {'—':>8} {'—':>10}  {'—':>3}  {status_short}"
            )
        else:
            sig_mk = "✱" if row["significant"] else " "
            p_str  = (f"{row['cox_p']:.3f}" if row["cox_p"] >= 0.001
                      else f"{row['cox_p']:.2e}")
            line = (
                f"{star+' '+row['gene']:<24} {row['endpoint']:<4} "
                f"{int(row['events_altered']):>6} {int(row['events_wildtype']):>6} "
                f"{row['median_altered_months']:>8} {row['median_wildtype_months']:>7} "
                f"{row['HR']:>6.2f} {row['CI_lower']:>8.2f} {row['CI_upper']:>8.2f} "
                f"{p_str:>10}  {sig_mk:>3}  {row['formatted_result']}"
            )
        print(line)

    print("=" * 110)
    n_done = (results_df["status"] == "completed").sum()
    n_skip = (results_df["status"] != "completed").sum()
    print(f"\nCompleted: {n_done}  |  Skipped (gate < {MIN_EVENTS} events/arm): {n_skip}")
    print(f"Significant (cox_p < 0.05): "
          f"{(results_df['significant'] == True).sum()} analyses")


if __name__ == "__main__":
    main()
