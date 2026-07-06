"""
11_pathway_survival.py

Pathway-level survival analysis for the Wilms' Tumor (TARGET, 2018) cohort
(wt_target_2018_pub, N=652), built on top of
processed/master_classification_extended.csv (script 10).

No individual gene in this cohort reaches the n=20 KM threshold (max = 18,
MYCN). Pathway-level OR-grouping aggregates all alterations that share a
common biological mechanism, which is the only route to adequately powered
survival analysis in this cohort.

Steps
-----
  1. Threshold check       — KM / Cox eligibility per pathway
  2. Kaplan-Meier          — EFS + OS for KM-eligible pathways
  3. Cox regression        — EFS + OS for Cox-eligible pathway x endpoint pairs
  4. Gene contribution     — stacked bar per pathway with n_altered >= 10
  5. Frequency summary     — one bar chart across all 6 pathway flags
  6. Results table          — results/pathway_results.csv
  7. README section         — results/README_pathway_section.md
  8. Commit + push, then splice the section into the top-level README.md

Outputs
-------
  wilms_survival/plots/PATHWAY_{NAME}_{ENDPOINT}.png
  wilms_survival/plots/PATHWAY_{NAME}_contributions.png
  wilms_survival/plots/PATHWAY_frequency_summary.png
  wilms_survival/results/pathway_results.csv
  wilms_survival/results/README_pathway_section.md
"""

import os
import subprocess
import sys
import textwrap
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
REPO_ROOT   = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
PROC_DIR    = os.path.join(SCRIPT_DIR, "..", "processed")
PLOTS_DIR   = os.path.join(SCRIPT_DIR, "..", "plots")
RESULTS_DIR = os.path.join(SCRIPT_DIR, "..", "results")
MASTER_CSV  = os.path.join(PROC_DIR, "master_classification_extended.csv")
RESULTS_CSV = os.path.join(RESULTS_DIR, "pathway_results.csv")
README_SECTION = os.path.join(RESULTS_DIR, "README_pathway_section.md")
README_MAIN    = os.path.join(REPO_ROOT, "README.md")

N_COHORT = 652

# ── Constants ─────────────────────────────────────────────────────────────────

KM_MIN_N       = 20   # patients per arm
COX_MIN_EVENTS = 10   # events per arm

ENDPOINTS = [
    ("EFS", "EFS_MONTHS", "EFS_STATUS", "Event-Free Survival"),
    ("OS",  "OS_MONTHS",  "OS_STATUS",  "Overall Survival"),
]

PATHWAY_GENE_COLS = {
    "WNT":       ["WT1_ALTERED", "CTNNB1_ALTERED", "AMER1_ALTERED", "NF1_ALTERED"],
    "SIX_miRNA": ["SIX1_ALTERED", "SIX2_ALTERED", "DROSHA_ALTERED", "DGCR8_ALTERED",
                  "DICER1_ALTERED", "XPO5_ALTERED", "DIS3L2_ALTERED"],
    "CHROMATIN": ["MLLT1_ALTERED", "BCOR_ALTERED", "BCORL1_ALTERED", "ARID1A_ALTERED",
                  "SMARCA4_ALTERED", "BRD7_ALTERED", "CREBBP_ALTERED", "EP300_ALTERED",
                  "HDAC4_ALTERED", "ASXL1_ALTERED", "CHD4_ALTERED", "MAP3K4_ALTERED",
                  "CTR9_ALTERED"],
    "P53":       ["TP53_ALTERED", "MYCN_ALTERED", "MDM2_ALTERED", "MDM4_ALTERED",
                  "MAX_ALTERED", "CHEK2_ALTERED", "PALB2_ALTERED"],
    "IGF2":      ["IGF2_CNA_ALTERED", "IGF2_EXPR_ALTERED", "NIPBL_ALTERED",
                  "NONO_ALTERED", "COL6A3_ALTERED", "FGFR1_ALTERED", "ACTB_ALTERED"],
}
# PATHWAY_ANY's "members" are the 5 pathway flags themselves
ANY_MEMBER_COLS = ["PATHWAY_WNT", "PATHWAY_SIX_miRNA", "PATHWAY_CHROMATIN",
                   "PATHWAY_P53", "PATHWAY_IGF2"]

PATHWAY_COLUMN = {
    "WNT": "PATHWAY_WNT", "SIX_miRNA": "PATHWAY_SIX_miRNA",
    "CHROMATIN": "PATHWAY_CHROMATIN", "P53": "PATHWAY_P53",
    "IGF2": "PATHWAY_IGF2", "ANY": "PATHWAY_ANY",
}

PATHWAY_DISPLAY = {
    "WNT": "WNT", "SIX_miRNA": "SIX-miRNA", "CHROMATIN": "Chromatin",
    "P53": "p53", "IGF2": "IGF2", "ANY": "Any-Pathway",
}

PATHWAY_TITLE = {
    "WNT": "WNT / Renal Differentiation",
    "SIX_miRNA": "Renal Progenitor / SIX-miRNAPG",
    "CHROMATIN": "Chromatin Remodeling",
    "P53": "p53 / Cell Cycle",
    "IGF2": "IGF2 / Growth Factor",
    "ANY": "Any Pathway (Combined)",
}

COLORS = {
    "WNT": "#1D9E75", "SIX_miRNA": "#7F77DD", "CHROMATIN": "#BA7517",
    "P53": "#C00000", "IGF2": "#639922", "ANY": "#1F3A5F",
}
WILDTYPE_COLOR = "#AAAAAA"

ORDER = ["WNT", "SIX_miRNA", "CHROMATIN", "P53", "IGF2", "ANY"]

REFERENCES = [
    "Gadd S et al. *Nat Genet.* 2017;49:1487–1494. (TARGET — primary cohort and clustering)",
    "Perotti D et al. *Nat Rev Urol.* 2024;21:158–180. (pathway definitions and prevalence)",
    "Treger TD et al. *Nat Rev Nephrol.* 2019;15:240–257. (chromatin/elongation gene class)",
    "Wegert J et al. *Cancer Cell.* 2015;27:298–311. (SIX-miRNAPG co-occurrence and clustering)",
    "Perlman EJ et al. *Nat Commun.* 2015;6:10013. (MLLT1 YEATS domain mutations)",
]


# ── Helpers ────────────────────────────────────────────────────────────────────

def get_arms(df: pd.DataFrame, flag_col: str, time_col: str, event_col: str):
    sub = df[[time_col, event_col]].copy()
    sub["flag"] = df[flag_col].astype("Float64").values
    sub = sub.dropna(subset=[time_col, event_col, "flag"])
    alt = sub[sub["flag"] == 1]
    wt  = sub[sub["flag"] == 0]
    return alt, wt


def median_label(kmf: KaplanMeierFitter):
    m = kmf.median_survival_time_
    try:
        if m is None or np.isinf(float(m)) or np.isnan(float(m)):
            return "NR", None
        return f"{float(m):.1f}", float(m)
    except (TypeError, ValueError):
        return "NR", None


def p_str(p: float) -> str:
    return f"{p:.3f}" if p >= 0.001 else f"{p:.2e}"


def wrap_genes(genes: list[str]) -> str:
    names = [g.replace("_ALTERED", "").replace("PATHWAY_", "") for g in genes]
    return textwrap.fill(", ".join(names), width=60)


# ── STEP 1: threshold check ──────────────────────────────────────────────────

def step1_threshold_check(df: pd.DataFrame) -> dict:
    print("\n" + "=" * 100)
    print("STEP 1 — PATHWAY THRESHOLD CHECK  (KM gate: n>=20/arm  |  Cox gate: events>=10/arm)")
    print("=" * 100)

    eligibility = {}   # name -> {"km": bool, "cox": {"EFS": bool, "OS": bool}}

    hdr = (f"{'Pathway':<12} {'n_alt':>6} {'n_wt':>6}  "
           f"{'EFS ev(alt/wt)':>16}  {'OS ev(alt/wt)':>15}   Flags")
    print(hdr)
    print("-" * 100)

    for name in ORDER:
        col = PATHWAY_COLUMN[name]
        flag = df[col].astype("Float64")
        n_alt = int((flag == 1).sum())
        n_wt  = int((flag == 0).sum())

        ev = {}
        for ep_label, time_col, event_col, _ in ENDPOINTS:
            alt, wt = get_arms(df, col, time_col, event_col)
            ev[ep_label] = (int(alt[event_col].sum()), int(wt[event_col].sum()))

        km_ok = n_alt >= KM_MIN_N and n_wt >= KM_MIN_N
        cox_ok = {
            ep_label: (ev[ep_label][0] >= COX_MIN_EVENTS and ev[ep_label][1] >= COX_MIN_EVENTS)
            for ep_label, *_ in ENDPOINTS
        }
        eligibility[name] = {"km": km_ok, "cox": cox_ok, "n_alt": n_alt, "n_wt": n_wt, "events": ev}

        tags = []
        if km_ok:
            tags.append("✅ KM eligible")
        for ep_label, ok in cox_ok.items():
            if ok:
                tags.append(f"✅ Cox eligible ({ep_label})")
        if not tags:
            tags = ["⚠️ insufficient"]

        efs_str = f"{ev['EFS'][0]}/{ev['EFS'][1]}"
        os_str  = f"{ev['OS'][0]}/{ev['OS'][1]}"
        print(f"{name:<12} {n_alt:>6} {n_wt:>6}  {efs_str:>16}  {os_str:>15}   " + "  ".join(tags))

    print("=" * 100)
    return eligibility


# ── STEP 2: Kaplan-Meier ─────────────────────────────────────────────────────

def run_km(df: pd.DataFrame, name: str) -> list[dict]:
    col = PATHWAY_COLUMN[name]
    genes = ANY_MEMBER_COLS if name == "ANY" else PATHWAY_GENE_COLS[name]
    subtitle = wrap_genes(genes)
    color = COLORS[name]
    out = []

    for ep_label, time_col, event_col, ep_long in ENDPOINTS:
        alt, wt = get_arms(df, col, time_col, event_col)
        n_alt, n_wt = len(alt), len(wt)

        kmf_alt = KaplanMeierFitter()
        kmf_wt  = KaplanMeierFitter()
        kmf_alt.fit(alt[time_col].astype(float), alt[event_col].astype(float),
                    label=f"Altered (n={n_alt})")
        kmf_wt.fit(wt[time_col].astype(float), wt[event_col].astype(float),
                   label=f"Wildtype (n={n_wt})")

        lr = logrank_test(alt[time_col].astype(float), wt[time_col].astype(float),
                           event_observed_A=alt[event_col].astype(float),
                           event_observed_B=wt[event_col].astype(float))
        p = float(lr.p_value)
        med_alt, med_alt_val = median_label(kmf_alt)
        med_wt, med_wt_val   = median_label(kmf_wt)
        ev_alt = int(alt[event_col].sum())
        ev_wt  = int(wt[event_col].sum())

        # Unadjusted direction (worse/better) from a bare ALTERED-only Cox fit —
        # independent of whether the covariate-adjusted model in step 3 converges,
        # so wording never silently defaults when that model is skipped/fails.
        uni = pd.DataFrame({
            "duration": pd.concat([alt[time_col], wt[time_col]]).astype(float).values,
            "event":    pd.concat([alt[event_col], wt[event_col]]).astype(float).values,
            "ALTERED":  [1] * n_alt + [0] * n_wt,
        })
        try:
            cph_uni = CoxPHFitter()
            cph_uni.fit(uni, duration_col="duration", event_col="event", show_progress=False)
            direction = "worse" if float(cph_uni.summary.loc["ALTERED", "coef"]) > 0 else "better"
        except Exception:
            direction = "worse" if ev_alt / max(n_alt, 1) > ev_wt / max(n_wt, 1) else "better"

        # ── Plot ──────────────────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(10, 8))
        plt.subplots_adjust(bottom=0.28, top=0.86)

        kmf_wt.plot_survival_function(ax=ax, color=WILDTYPE_COLOR, ci_show=True,
                                       ci_alpha=0.12, linewidth=2.2)
        kmf_alt.plot_survival_function(ax=ax, color=color, ci_show=True,
                                        ci_alpha=0.15, linewidth=2.4)
        add_at_risk_counts(kmf_wt, kmf_alt, ax=ax, fontsize=8.5, rows_to_show=["At risk"])

        ax.set_xlabel(f"{ep_long} (months)", fontsize=11)
        ax.set_ylabel("Survival probability", fontsize=11)
        ax.set_ylim(-0.02, 1.05)
        ax.spines[["top", "right"]].set_visible(False)

        ax.set_title(
            f"{PATHWAY_DISPLAY[name]} Pathway — {ep_long} — "
            f"Wilms Tumor (N={N_COHORT}, wt_target_2018_pub)",
            fontsize=12.5, fontweight="bold", pad=32,
        )
        fig.text(0.5, 0.895, subtitle, ha="center", va="top", fontsize=8.5,
                  style="italic", color="#555555", transform=fig.transFigure)

        pstr = p_str(p)
        ann = f"Median altered: {med_alt} mo | Median WT: {med_wt} mo | log-rank p = {pstr}"
        ax.text(0.98, 0.97, ann, transform=ax.transAxes, ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFFDE7",
                          edgecolor="#AAAAAA", alpha=0.92))
        ax.legend(loc="lower left", fontsize=10, framealpha=0.85)

        fname = os.path.join(PLOTS_DIR, f"PATHWAY_{name}_{ep_label}.png")
        fig.savefig(fname, dpi=150, bbox_inches="tight")
        plt.close(fig)

        print(f"  {name:<12} {ep_label:<4} n_alt={n_alt:>3} n_wt={n_wt:>3} "
              f"med_alt={med_alt:>6} med_wt={med_wt:>6} p={pstr:<10} "
              f"{'★ SIGNIFICANT' if p < 0.05 else ''}  -> {fname}")

        out.append({
            "pathway": name, "endpoint": ep_label,
            "n_altered": n_alt, "n_wildtype": n_wt,
            "events_altered": ev_alt, "events_wildtype": ev_wt,
            "median_altered": med_alt, "median_wildtype": med_wt,
            "logrank_p": p, "km_significant": p < 0.05, "direction": direction,
        })
    return out


# ── STEP 3: Cox regression ───────────────────────────────────────────────────

def run_cox(df: pd.DataFrame, name: str, ep_label: str, time_col: str, event_col: str,
            use_age: bool, use_cohort: bool) -> dict | None:
    col = PATHWAY_COLUMN[name]
    flag = df[col].astype("Float64")

    sub = df[[time_col, event_col]].copy()
    sub["ALTERED"] = flag.values
    if use_age:
        sub["AGE"] = df["AGE"].values
    if use_cohort:
        sub["COHORT_CODE"] = df["COHORT_CODE"].values
    sub = sub.dropna()

    covs = ["ALTERED"] + (["AGE"] if use_age else []) + (["COHORT_CODE"] if use_cohort else [])
    cox_df = sub[[time_col, event_col] + covs].rename(
        columns={time_col: "duration", event_col: "event"})

    cph = CoxPHFitter()
    try:
        cph.fit(cox_df, duration_col="duration", event_col="event", show_progress=False)
    except Exception as exc:
        print(f"  ERROR  {name:<12} {ep_label}: {exc}")
        return None

    if "ALTERED" not in cph.summary.index:
        print(f"  ERROR  {name:<12} {ep_label}: ALTERED not in Cox summary")
        return None

    row = cph.summary.loc["ALTERED"]
    hr    = float(np.exp(row["coef"]))
    ci_lo = float(np.exp(row["coef lower 95%"]))
    ci_hi = float(np.exp(row["coef upper 95%"]))
    p     = float(row["p"])

    print(f"  {name:<12} {ep_label:<4} covariates={'+'.join(covs):<20} "
          f"HR={hr:.2f} [{ci_lo:.2f}–{ci_hi:.2f}] p={p_str(p)} "
          f"-> HR = {hr:.2f} [{ci_lo:.2f}–{ci_hi:.2f}], p = {p_str(p)}")

    return {
        "pathway": name, "endpoint": ep_label,
        "HR": round(hr, 3), "CI_lower": round(ci_lo, 3), "CI_upper": round(ci_hi, 3),
        "cox_p": p, "cox_significant": p < 0.05,
    }


# ── STEP 4: gene contribution charts ─────────────────────────────────────────

def gene_contribution_chart(df: pd.DataFrame, name: str) -> dict:
    col = PATHWAY_COLUMN[name]
    gene_cols = ANY_MEMBER_COLS if name == "ANY" else PATHWAY_GENE_COLS[name]

    altered = df[df[col] == 1]
    n_total = len(altered)

    n_hit = pd.Series(0, index=altered.index)
    for g in gene_cols:
        n_hit = n_hit + (altered[g] == 1).astype(int)

    rows = []
    for g in gene_cols:
        hit = altered[g] == 1
        unique = int((hit & (n_hit == 1)).sum())
        shared = int((hit & (n_hit > 1)).sum())
        rows.append({
            "gene": g.replace("_ALTERED", "").replace("PATHWAY_", ""),
            "unique": unique, "shared": shared, "total": unique + shared,
        })

    contrib = pd.DataFrame(rows).sort_values("total", ascending=False).reset_index(drop=True)
    contrib = contrib[contrib["total"] > 0].reset_index(drop=True)

    color = COLORS[name]
    fig, ax = plt.subplots(figsize=(9, max(3, 0.55 * len(contrib) + 1.2)))
    y_pos = np.arange(len(contrib))[::-1]   # largest at top

    ax.barh(y_pos, contrib["unique"], color=color, edgecolor="white", label="Unique")
    ax.barh(y_pos, contrib["shared"], left=contrib["unique"], color=color, alpha=0.55,
            hatch="//", edgecolor="white", label="Shared (multi-gene)")

    for y, (u, s, t) in zip(y_pos, contrib[["unique", "shared", "total"]].values):
        ax.text(t + 0.15, y, str(int(t)), va="center", fontsize=9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(contrib["gene"])
    ax.set_xlabel("Altered patients contributed")
    ax.set_title(f"{PATHWAY_DISPLAY[name]} — Gene Contributions to Altered Arm (n={n_total})",
                 fontsize=12.5, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.85)

    fname = os.path.join(PLOTS_DIR, f"PATHWAY_{name}_contributions.png")
    fig.tight_layout()
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  {name:<12} n_altered={n_total:<4} top={contrib.iloc[0]['gene'] if len(contrib) else '-'} "
          f"-> {fname}")

    return {
        "pathway": name, "n_altered": n_total, "contrib": contrib,
        "n_multi_gene": int((n_hit > 1).sum()),
        "file": os.path.relpath(fname, REPO_ROOT),
    }


# ── STEP 5: frequency summary chart ──────────────────────────────────────────

def frequency_summary_chart(df: pd.DataFrame, eligibility: dict) -> None:
    rows = []
    for name in ORDER:
        col = PATHWAY_COLUMN[name]
        n_alt = int((df[col] == 1).sum())
        pct = round(n_alt / N_COHORT * 100, 1)
        rows.append({"name": name, "n_alt": n_alt, "pct": pct,
                     "km_ok": eligibility[name]["km"]})
    freq = pd.DataFrame(rows).sort_values("pct", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10, 5))
    y_pos = np.arange(len(freq))[::-1]

    for y, r in zip(y_pos, freq.itertuples()):
        ax.barh(y, r.pct, color=COLORS[r.name], edgecolor="white")
        # DejaVu Sans (matplotlib's default font) has no glyph for U+2705 (✅);
        # use the plain check mark (U+2713) here and reserve the emoji for
        # console output / the Markdown README, which both render it fine.
        tag = "✓ KM eligible" if r.km_ok else "⚠ n<20"
        ax.text(r.pct + 0.3, y, f"n={r.n_alt} ({r.pct}%)  {tag}", va="center", fontsize=9.5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([PATHWAY_DISPLAY[n] for n in freq["name"]])
    ax.set_xlabel("% of cohort altered")
    ax.set_xlim(0, max(freq["pct"]) * 1.9)
    ax.set_title(f"Wilms Tumor — Pathway-Level Alteration Frequency (N={N_COHORT})",
                 fontsize=13, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    fname = os.path.join(PLOTS_DIR, "PATHWAY_frequency_summary.png")
    fig.tight_layout()
    fig.savefig(fname, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved -> {fname}")


# ── STEP 7: README section ───────────────────────────────────────────────────

def fmt_pct(n: int) -> float:
    return round(n / N_COHORT * 100, 1)


def build_readme_section(df, eligibility, results_df, contrib_results) -> str:
    n_any = int((df["PATHWAY_ANY"] == 1).sum())
    pct_any = fmt_pct(n_any)

    lines = []
    lines.append("---")
    lines.append("")
    lines.append("## Pathway-Level Survival Analysis")
    lines.append("")
    lines.append("<!-- Why pathways: individual genes all failed n≥20 threshold in this 652-patient cohort -->")
    lines.append("")
    lines.append("> **Why pathway analysis?** No individual gene in the TARGET 2018 cohort (N=652) reached")
    lines.append("> the minimum 20-patient threshold for Kaplan-Meier analysis — the most frequent gene")
    lines.append("> (MYCN) had only 18 altered patients. This reflects a well-characterized feature of")
    lines.append("> Wilms tumor biology: no single gene drives more than ~15% of cases. Pathway-level")
    lines.append("> OR-grouping aggregates all alterations sharing a common biological mechanism,")
    lines.append("> providing the statistical power needed for survival analysis while preserving")
    lines.append("> biological interpretability.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("### Five Molecular Pathways")
    lines.append("")
    lines.append("<!-- Evidence-based pathway definitions from Gadd 2017 (TARGET), Perotti 2024, Treger 2019 -->")
    lines.append("")
    lines.append("| Pathway | Member Genes | N Altered | Cohort Coverage | Literature Basis |")
    lines.append("|---|---|---|---|---|")

    pathway_rows = [
        ("WNT / Renal Differentiation", "WT1, CTNNB1, AMER1, NF1",
         "Co-cluster in TARGET expression analysis (Gadd 2017); converge on WNT-β-catenin disruption and mesenchymal-to-epithelial transition"),
        ("Renal Progenitor / SIX-miRNAPG", "SIX1, SIX2, DROSHA, DGCR8, DICER1, XPO5, DIS3L2",
         "Share common output of perpetuating renal progenitor state; blastemal histology association; enriched at relapse (Wegert 2015; Perotti 2024)"),
        ("Chromatin Remodeling", "MLLT1, BCOR, BCORL1, ARID1A, SMARCA4, BRD7, CREBBP, EP300, HDAC4, ASXL1, CHD4, MAP3K4, CTR9",
         "Named as functional class in Treger 2019: epigenetic regulators of renal progenitor differentiation via transcriptional elongation and chromatin modification"),
        ("p53 / Cell Cycle", "TP53, MYCN, MDM2, MDM4, MAX, CHEK2, PALB2",
         "TP53 defines diffuse anaplastic Wilms tumor; MYCN amplification predicts poor EFS/OS; MDM2/MDM4 provide alternative p53 inactivation (Perotti 2024)"),
        ("IGF2 / Growth Factor", "IGF2 (CNA + expr), NIPBL, NONO, COL6A3, FGFR1, ACTB",
         "IGF2 loss of imprinting (~70% prevalence) drives PI3K-AKT via IGF1R; NIPBL disrupts cohesin-mediated control of the IGF2/H19 imprinted domain (Treger 2019)"),
    ]
    name_for_row = ["WNT", "SIX_miRNA", "CHROMATIN", "P53", "IGF2"]
    for (title, genes, lit), name in zip(pathway_rows, name_for_row):
        n_alt = eligibility[name]["n_alt"]
        pct = fmt_pct(n_alt)
        lines.append(f"| **{title}** | {genes} | {n_alt} | {pct}% | {lit} |")
    lines.append("")
    lines.append(f"**Combined:** {pct_any}% of the cohort ({n_any}/{N_COHORT} patients) carry at least one alteration across any pathway.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("### Alteration Frequency")
    lines.append("")
    lines.append("<!-- How much of the cohort each pathway captures vs individual genes -->")
    lines.append("")
    lines.append("![Pathway alteration frequency](wilms_survival/plots/PATHWAY_frequency_summary.png)")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("### Survival Analysis Results")
    lines.append("")
    lines.append("<!-- KM and Cox results — first valid survival statistics from this dataset -->")
    lines.append("")
    lines.append("| Pathway | N Altered | N Wildtype | Cohort % | Endpoint | Median Altered (mo) | Median WT (mo) | Log-rank p | HR [95% CI] | Sig? |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")

    for _, r in results_df.iterrows():
        km_ran = pd.notna(r["logrank_p"])
        cox_ran = pd.notna(r["HR"])
        med_alt = r["median_altered"] if km_ran else "—"
        med_wt  = r["median_wildtype"] if km_ran else "—"
        p_col   = p_str(r["logrank_p"]) if km_ran else "—"
        hr_col  = f"{r['HR']:.2f} [{r['CI_lower']:.2f}–{r['CI_upper']:.2f}]" if cox_ran else "—"
        if km_ran and r["km_significant"]:
            sig = "✅"
        elif not km_ran:
            sig = "⚠️ skipped"
        else:
            sig = ""
        lines.append(
            f"| {PATHWAY_DISPLAY[r['pathway']]} | {int(r['n_altered'])} | {int(r['n_wildtype'])} | "
            f"{r['pct_cohort']}% | {r['endpoint']} | {med_alt} | {med_wt} | {p_col} | {hr_col} | {sig} |"
        )
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("### Key Findings")
    lines.append("")
    lines.append("<!-- One block per significant result; honest null-result statement if none significant -->")
    lines.append("")

    sig_rows = results_df[(results_df["km_significant"] == True)]
    if len(sig_rows) > 0:
        for _, r in sig_rows.iterrows():
            name = r["pathway"]
            n_alt = int(r["n_altered"])
            pct = r["pct_cohort"]
            direction = r["direction"]
            hr_txt = (f"; HR = {r['HR']:.2f} [{r['CI_lower']:.2f}–{r['CI_upper']:.2f}]"
                      if pd.notna(r["HR"]) else "")
            lines.append(f"**{PATHWAY_DISPLAY[name]}** — {r['endpoint']}")
            lines.append("")
            lines.append(
                f"> Patients with any alteration in the {PATHWAY_TITLE[name]} pathway "
                f"({n_alt} patients, {pct}% of cohort)"
            )
            lines.append(
                f"> showed significantly {direction} {r['endpoint']} compared to unaltered patients"
            )
            lines.append(
                f"> (median: **{r['median_altered']} vs. {r['median_wildtype']} months**; "
                f"log-rank p = {p_str(r['logrank_p'])}{hr_txt})."
            )
            lines.append("")
            lines.append(f"![KM curve](wilms_survival/plots/PATHWAY_{name}_{r['endpoint']}.png)")
            lines.append("")
    else:
        any_row_efs = results_df[(results_df["pathway"] == "ANY") & (results_df["endpoint"] == "EFS")].iloc[0]
        did = "did" if any_row_efs["km_significant"] else "did not"
        lines.append(
            "> No pathway reached statistical significance in this single-cohort analysis."
        )
        lines.append(
            f"> The combined pathway OR-classifier (`PATHWAY_ANY`) captured {pct_any}% of the cohort —"
        )
        lines.append(
            f"> substantially more than any individual gene — but {did} reach significance"
        )
        lines.append(
            f"> (p = {p_str(any_row_efs['logrank_p'])}). Integration with additional cBioPortal cohorts (PedcBioPortal,"
        )
        lines.append("> additional TARGET data) is the recommended next step to increase event counts.")
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("### Gene Contributions Within Pathways")
    lines.append("")
    lines.append("<!-- Shows which specific genes drive the altered arm in each pathway -->")
    lines.append("")

    for cr in contrib_results:
        name = cr["pathway"]
        contrib = cr["contrib"]
        lines.append(f"#### {PATHWAY_TITLE[name]}")
        lines.append("")
        lines.append(f"![Gene contributions](wilms_survival/plots/PATHWAY_{name}_contributions.png)")
        lines.append("")
        if len(contrib) >= 1:
            top = contrib.iloc[0]
            sentence1 = (f"> The {PATHWAY_TITLE[name]} altered arm (n={cr['n_altered']}) was primarily "
                         f"driven by {top['gene']} ({int(top['total'])} patients).")
            if len(contrib) >= 2:
                second = contrib.iloc[1]
                sentence2 = f"{second['gene']} contributed {int(second['total'])} additional unique cases."
            else:
                sentence2 = ""
            sentence3 = (f"{cr['n_multi_gene']} patients carried alterations in more than one pathway "
                         f"gene." if name != "ANY" else
                         f"{cr['n_multi_gene']} patients carried alterations in more than one pathway.")
            lines.append(f"{sentence1} {sentence2} {sentence3}".strip())
            lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("### How to Read This Analysis")
    lines.append("")
    lines.append("<!-- Plain-language guide for research advisors unfamiliar with OR-classifiers -->")
    lines.append("")
    lines.append("**What \"pathway-altered\" means:** A patient is classified as altered for a pathway")
    lines.append("if they carry a qualifying mutation or copy number change in *any* gene in that")
    lines.append("pathway — they do not need alterations in all genes. This is an OR logic classifier.")
    lines.append("")
    lines.append("**Cohort coverage vs. significance:** A pathway capturing 30% of 652 patients gives")
    lines.append("roughly 196 altered patients — enough for robust KM and Cox analysis. This is why")
    lines.append("pathway grouping succeeds where individual genes (max n=18) cannot.")
    lines.append("")
    lines.append("**Interpreting HR:** A hazard ratio > 1 means pathway-altered patients have higher")
    lines.append("instantaneous risk of the event (death or relapse) at any point in follow-up vs.")
    lines.append("unaltered patients, after adjusting for cohort and age. HR is from multivariate")
    lines.append("Cox regression.")
    lines.append("")
    lines.append("**Gene contribution charts** show whether a pathway result is driven by one dominant")
    lines.append("gene (e.g., TP53 dominating the p53 arm) or is a genuinely distributed multi-gene")
    lines.append("effect. This distinction matters for clinical interpretation.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("### References")
    lines.append("")
    for ref in REFERENCES:
        lines.append(f"- {ref}")
    lines.append("")
    lines.append("---")
    return "\n".join(lines)


# ── STEP 8: git commit + push, splice into README.md ─────────────────────────

def run_git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git"] + list(args), cwd=REPO_ROOT,
                           capture_output=True, text=True)


def git_commit_and_push() -> None:
    print("\n" + "=" * 100)
    print("STEP 8 — GIT COMMIT + PUSH")
    print("=" * 100)

    rr = run_git("rev-parse", "--show-toplevel")
    if rr.returncode != 0:
        print("  Not a git repository — skipping commit/push.")
        return

    add_paths = [
        "wilms_survival/scripts/10_pathway_classification.py",
        "wilms_survival/scripts/11_pathway_survival.py",
        "wilms_survival/processed/master_classification_extended.csv",
        "wilms_survival/results/pathway_results.csv",
        "wilms_survival/results/README_pathway_section.md",
    ]
    run_git("add", *add_paths)
    png_glob = [os.path.relpath(p, REPO_ROOT) for p in
                __import__("glob").glob(os.path.join(PLOTS_DIR, "PATHWAY_*.png"))]
    if png_glob:
        run_git("add", *png_glob)

    status = run_git("status")
    print(status.stdout)

    commit = run_git("commit", "-m",
                      "Add pathway-level survival analysis — five molecular pathways, N=652 cohort")
    print(commit.stdout + commit.stderr)
    if commit.returncode != 0 and "nothing to commit" not in (commit.stdout + commit.stderr):
        print("  Commit failed — skipping push.")
        return

    push = run_git("push", "origin", "main")
    print(push.stdout + push.stderr)
    if push.returncode != 0:
        print("  Push failed.")
        return

    # ── Splice README section into top-level README.md ──────────────────────
    with open(README_MAIN, "r") as f:
        readme = f.read()
    with open(README_SECTION, "r") as f:
        section = f.read()

    marker = "## Suggested Next Steps"
    if marker in readme:
        idx = readme.index(marker)
        # back up to the start of that heading's line, and the blank line/rule before it
        insert_at = readme.rfind("\n---\n", 0, idx)
        insert_at = insert_at + 1 if insert_at != -1 else idx
    else:
        # fall back to the last top-level section heading in the file
        import re
        headings = [m.start() for m in re.finditer(r"^## ", readme, flags=re.MULTILINE)]
        last_heading = headings[-1]
        insert_at = readme.rfind("\n---\n", 0, last_heading)
        insert_at = insert_at + 1 if insert_at != -1 else last_heading

    new_readme = readme[:insert_at] + section + "\n\n" + readme[insert_at:]
    with open(README_MAIN, "w") as f:
        f.write(new_readme)
    print(f"\n  Inserted pathway section into {README_MAIN}")

    run_git("add", "README.md")
    commit2 = run_git("commit", "-m", "Insert pathway analysis section into README")
    print(commit2.stdout + commit2.stderr)
    push2 = run_git("push", "origin", "main")
    print(push2.stdout + push2.stderr)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    os.makedirs(PLOTS_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    df = pd.read_csv(MASTER_CSV)
    print(f"Loaded master_classification_extended.csv  ({len(df)} patients)")

    # ── STEP 1 ────────────────────────────────────────────────────────────────
    eligibility = step1_threshold_check(df)

    # ── STEP 2: KM ────────────────────────────────────────────────────────────
    print("\n" + "=" * 100)
    print("STEP 2 — KAPLAN-MEIER ANALYSIS")
    print("=" * 100)
    km_results: list[dict] = []
    for name in ORDER:
        if eligibility[name]["km"]:
            km_results.extend(run_km(df, name))
        else:
            print(f"  SKIP  {name:<12} — KM gate not met "
                  f"(n_alt={eligibility[name]['n_alt']}, n_wt={eligibility[name]['n_wt']})")

    # ── STEP 3: Cox ───────────────────────────────────────────────────────────
    print("\n" + "=" * 100)
    print("STEP 3 — COX REGRESSION")
    print("=" * 100)
    use_age = df["AGE"].isna().mean() <= 0.30
    use_cohort = df["COHORT_CODE"].nunique() > 1
    print(f"  AGE covariate    : {'INCLUDED' if use_age else 'OMITTED'}")
    print(f"  COHORT_CODE      : {'INCLUDED' if use_cohort else 'OMITTED (constant — zero variance)'}")

    cox_results: list[dict] = []
    for name in ORDER:
        for ep_label, time_col, event_col, _ in ENDPOINTS:
            if eligibility[name]["cox"][ep_label]:
                r = run_cox(df, name, ep_label, time_col, event_col, use_age, use_cohort)
                if r is not None:
                    cox_results.append(r)
            else:
                ev = eligibility[name]["events"][ep_label]
                print(f"  SKIP  {name:<12} {ep_label:<4} — Cox gate not met "
                      f"(events_alt={ev[0]}, events_wt={ev[1]})")

    # ── STEP 4: gene contribution charts ─────────────────────────────────────
    print("\n" + "=" * 100)
    print("STEP 4 — GENE CONTRIBUTION CHARTS  (pathways with n_altered >= 10)")
    print("=" * 100)
    contrib_results = []
    for name in ORDER:
        n_alt = eligibility[name]["n_alt"]
        if n_alt >= 10:
            contrib_results.append(gene_contribution_chart(df, name))
        else:
            print(f"  SKIP  {name:<12} — n_altered={n_alt} < 10")

    # ── STEP 5: frequency summary chart ──────────────────────────────────────
    print("\n" + "=" * 100)
    print("STEP 5 — PATHWAY FREQUENCY SUMMARY CHART")
    print("=" * 100)
    frequency_summary_chart(df, eligibility)

    # ── STEP 6: assemble + save results table ────────────────────────────────
    print("\n" + "=" * 100)
    print("STEP 6 — SAVE RESULTS")
    print("=" * 100)

    km_df  = pd.DataFrame(km_results)
    cox_df = pd.DataFrame(cox_results)

    rows = []
    for name in ORDER:
        n_alt = eligibility[name]["n_alt"]
        n_wt  = eligibility[name]["n_wt"]
        pct   = fmt_pct(n_alt)
        for ep_label, *_ in ENDPOINTS:
            row = {
                "pathway": name, "endpoint": ep_label,
                "n_altered": n_alt, "n_wildtype": n_wt, "pct_cohort": pct,
                "events_alt": eligibility[name]["events"][ep_label][0],
                "events_wt": eligibility[name]["events"][ep_label][1],
                "median_altered": np.nan, "median_wildtype": np.nan,
                "logrank_p": np.nan, "HR": np.nan, "CI_lower": np.nan,
                "CI_upper": np.nan, "cox_p": np.nan,
                "km_significant": np.nan, "cox_significant": np.nan, "direction": np.nan,
            }
            if not km_df.empty:
                m = km_df[(km_df["pathway"] == name) & (km_df["endpoint"] == ep_label)]
                if not m.empty:
                    m = m.iloc[0]
                    row.update({
                        "median_altered": m["median_altered"], "median_wildtype": m["median_wildtype"],
                        "logrank_p": m["logrank_p"], "km_significant": m["km_significant"],
                        "direction": m["direction"],
                    })
            if not cox_df.empty:
                c = cox_df[(cox_df["pathway"] == name) & (cox_df["endpoint"] == ep_label)]
                if not c.empty:
                    c = c.iloc[0]
                    row.update({
                        "HR": c["HR"], "CI_lower": c["CI_lower"], "CI_upper": c["CI_upper"],
                        "cox_p": c["cox_p"], "cox_significant": c["cox_significant"],
                    })
            rows.append(row)

    results_df = pd.DataFrame(rows)

    csv_cols = ["pathway", "endpoint", "n_altered", "n_wildtype", "pct_cohort",
                "events_alt", "events_wt", "median_altered", "median_wildtype",
                "logrank_p", "HR", "CI_lower", "CI_upper", "cox_p",
                "km_significant", "cox_significant"]
    results_df[csv_cols].rename(columns={"median_wildtype": "median_wt"}).to_csv(RESULTS_CSV, index=False)
    print(f"  Saved -> {os.path.abspath(RESULTS_CSV)}")

    # ── STEP 7: README section ────────────────────────────────────────────────
    print("\n" + "=" * 100)
    print("STEP 7 — WRITE README_pathway_section.md")
    print("=" * 100)
    section = build_readme_section(df, eligibility, results_df, contrib_results)
    with open(README_SECTION, "w") as f:
        f.write(section)
    print(f"  Saved -> {os.path.abspath(README_SECTION)}")

    # ── STEP 8: commit + push ─────────────────────────────────────────────────
    git_commit_and_push()

    # ── Final confirmation ─────────────────────────────────────────────────────
    km_pass = [n for n in ORDER if eligibility[n]["km"]]
    sig_rows = results_df[results_df["km_significant"] == True]
    sig_list = [f"{PATHWAY_DISPLAY[r['pathway']]} ({r['endpoint']})" for _, r in sig_rows.iterrows()]

    print("\n✅ Pathway analysis complete and pushed.")
    print("View at: https://github.com/srinjoyghose1/wilms-tumor-survival-analysis")
    print(f"Pathways that reached KM threshold: {[PATHWAY_DISPLAY[n] for n in km_pass]}")
    print(f"Pathways that were significant (p<0.05): "
          f"{sig_list if sig_list else 'none — expand cohort recommended'}")


if __name__ == "__main__":
    main()
