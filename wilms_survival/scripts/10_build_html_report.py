"""
10_build_html_report.py
=======================
Build a self-contained mentor-ready HTML report of Wilms tumor survival findings.
Re-renders KM plots as PNG (via matplotlib + lifelines, no pdf2image dependency).
"""

import matplotlib
matplotlib.use("Agg")  # Must be before any pyplot import

import io
import re
import html
import base64
import traceback
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# ── Paths ────────────────────────────────────────────────────────────────────
BASE        = Path("/Users/srinjoy/Documents/Ravi Research/WIlms Tumor/wilms_survival")
DATA_DIR    = BASE / "data"
PROC_DIR    = BASE / "processed"
RESULTS_DIR = BASE / "results"
PLOTS_DIR   = BASE / "plots"
OUT_HTML    = RESULTS_DIR / "WilmsTumor_cBioPortal_Findings.html"
MASTER_CSV  = PROC_DIR / "master_classification.csv"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ── OR-classifier definitions ─────────────────────────────────────────────────
OR_DEFS = {
    "TP53_MDM4_MDM2_OR": ["TP53_ALTERED", "MDM4_ALTERED", "MDM2_ALTERED"],
    "MYCN_SIX_OR":       ["MYCN_ALTERED", "SIX1_SIX2_ALTERED"],
    "ALL_PRIMARY_OR":    ["TP53_ALTERED", "MYCN_ALTERED", "SIX1_ALTERED", "SIX2_ALTERED",
                          "DROSHA_ALTERED", "DGCR8_ALTERED", "CTNNB1_ALTERED",
                          "AMER1_ALTERED", "WT1_ALTERED", "MLLT1_ALTERED", "NIPBL_ALTERED"],
}
OR_GENE_LABELS = {
    "TP53_MDM4_MDM2_OR": "TP53, MDM4, MDM2",
    "DROSHA_DGCR8_OR":   "DROSHA, DGCR8",
    "MYCN_SIX_OR":       "MYCN, SIX1/SIX2",
    "WNT_AXIS_OR":       "CTNNB1, AMER1, WT1",
    "ALL_PRIMARY_OR":    "TP53, MYCN, SIX1, SIX2, DROSHA, DGCR8, CTNNB1, AMER1, WT1, MLLT1, NIPBL",
}

# Pathway → color mapping
PATHWAY_COLORS = {
    "p53 / MDM Axis":                "#C00000",
    "Transcription + Progenitor Axis":"#1F4E79",
    "WNT Signaling Axis":             "#375623",
    "miRNA Processing":               "#7030A0",
    "Epigenetic / Cohesin":           "#833C00",
    "IGF2 / Growth":                  "#375623",
    "All Primary Pathway Genes":      "#1F4E79",
}

GENE_PATHWAYS = {
    "TP53_ALTERED":      "p53 / MDM Axis",
    "MDM4_ALTERED":      "p53 / MDM Axis",
    "MDM2_ALTERED":      "p53 / MDM Axis",
    "MYCN_ALTERED":      "Transcription + Progenitor Axis",
    "SIX1_ALTERED":      "Transcription + Progenitor Axis",
    "SIX2_ALTERED":      "Transcription + Progenitor Axis",
    "SIX1_SIX2_ALTERED": "Transcription + Progenitor Axis",
    "CTNNB1_ALTERED":    "WNT Signaling Axis",
    "AMER1_ALTERED":     "WNT Signaling Axis",
    "WT1_ALTERED":       "WNT Signaling Axis",
    "DROSHA_ALTERED":    "miRNA Processing",
    "DGCR8_ALTERED":     "miRNA Processing",
    "MLLT1_ALTERED":     "Epigenetic / Cohesin",
    "NIPBL_ALTERED":     "Epigenetic / Cohesin",
    "IGF2_CNA_ALTERED":  "IGF2 / Growth",
    "IGF2_EXPR_ALTERED": "IGF2 / Growth",
}

GENE_DISPLAY = {
    "TP53_ALTERED":      "TP53",
    "MDM4_ALTERED":      "MDM4",
    "MDM2_ALTERED":      "MDM2",
    "MYCN_ALTERED":      "MYCN",
    "SIX1_ALTERED":      "SIX1",
    "SIX2_ALTERED":      "SIX2",
    "SIX1_SIX2_ALTERED": "SIX1/SIX2",
    "CTNNB1_ALTERED":    "CTNNB1",
    "AMER1_ALTERED":     "AMER1",
    "WT1_ALTERED":       "WT1",
    "DROSHA_ALTERED":    "DROSHA",
    "DGCR8_ALTERED":     "DGCR8",
    "MLLT1_ALTERED":     "MLLT1",
    "NIPBL_ALTERED":     "NIPBL",
    "IGF2_CNA_ALTERED":  "IGF2 (CNA)",
    "IGF2_EXPR_ALTERED": "IGF2 (Expr)",
}

COHORT_MAP = {
    "TARGET_2018": "wt_target_2018_pub",
    "TARGET_GDC":  "wt_target_gdc",
}

# ── Build log ─────────────────────────────────────────────────────────────────
build_log = []

def safe_section(name, fn, *args, **kwargs):
    try:
        result = fn(*args, **kwargs)
        build_log.append(("✓", name))
        return result
    except Exception as e:
        tb = traceback.format_exc()
        build_log.append(("✗", name, str(e), tb))
        return (
            f'<div class="container"><div class="error-notice">'
            f'⚠ Section "{html.escape(name)}" failed to build: '
            f'{html.escape(str(e))}</div></div>'
        )

# ── Helper formatters ─────────────────────────────────────────────────────────
def fmt_p(p):
    if pd.isna(p): return "—"
    p = float(p)
    if p < 0.001: return "&lt;0.001"
    return f"{p:.3f}"

def fmt_med(m):
    if pd.isna(m) or str(m).strip() in ("", "nan", "NaN", "NR"):
        return "NR"
    try:
        return f"{float(m):.1f}"
    except Exception:
        return str(m)

def fmt_ci(lo, hi):
    if pd.isna(lo) or pd.isna(hi): return "—"
    return f"[{float(lo):.2f}–{float(hi):.2f}]"

def fmt_hr_ci(hr, lo, hi):
    if pd.isna(hr): return "—"
    ci = fmt_ci(lo, hi)
    return f"{float(hr):.2f} {ci}"

def esc(s):
    if pd.isna(s):
        return ""
    return html.escape(str(s))

def _p_is_sig(p_str):
    """Return True if the (possibly fmt_p-formatted) p-value string is significant."""
    clean = str(p_str).replace("&lt;", "<").replace("&gt;", ">").strip()
    if clean.startswith("<"):
        try:
            return float(clean[1:]) <= 0.05
        except Exception:
            return False
    try:
        return float(clean) < 0.05
    except Exception:
        return False

def img_to_b64(fig, dpi=150):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return b64

def pathway_badge(pathway):
    color = PATHWAY_COLORS.get(pathway, "#666")
    return f'<span class="pathway-badge" style="background:{color}">{esc(pathway)}</span>'

# ── OR flag computation ───────────────────────────────────────────────────────
def compute_or_flag(df, cols):
    result = pd.Series(pd.NA, index=df.index, dtype="Float64")
    for col in cols:
        if col not in df.columns:
            continue
        s = df[col].astype("Float64")
        has_data = s.notna() | result.notna()
        result.loc[has_data] = result.loc[has_data].fillna(0.0)
        result.loc[(s == 1).fillna(False)] = 1.0
    return result

# ── Parse stem ───────────────────────────────────────────────────────────────
def parse_stem(stem):
    if stem.startswith("KM_EXPR_"):
        rest = stem[8:]
        ep = rest.rsplit("_", 1)[1]
        rest2 = rest.rsplit("_", 1)[0]
        for cohort in ["TARGET_2018", "TARGET_GDC"]:
            if rest2.endswith("_" + cohort):
                gene = rest2[:-(len(cohort) + 1)]
                return {"type": "expr", "gene": gene, "cohort": cohort, "ep": ep}
        return None
    elif stem.startswith("KM_"):
        rest = stem[3:]
        ep = rest.rsplit("_", 1)[1]
        gene = rest.rsplit("_", 1)[0]
        return {"type": "genomic", "gene": gene, "ep": ep}
    return None

# ── Plot renderers ────────────────────────────────────────────────────────────
def render_genomic_km(gene, ep, master_df):
    """Render a genomic KM plot and return base64 PNG."""
    df = master_df.copy()

    # Compute flag
    if gene in OR_DEFS:
        df["_FLAG"] = compute_or_flag(df, OR_DEFS[gene])
        gene_label = OR_GENE_LABELS.get(gene, gene)
    else:
        col = f"{gene}_ALTERED"
        if col not in df.columns:
            raise ValueError(f"Column {col} not found in master_classification.csv")
        df["_FLAG"] = df[col].astype("Float64")
        gene_label = gene.replace("_ALTERED", "").replace("_", "/")

    time_col = f"{ep}_MONTHS"
    event_col = f"{ep}_STATUS"
    if time_col not in df.columns or event_col not in df.columns:
        raise ValueError(f"Missing columns {time_col} or {event_col}")

    sub = df[[time_col, event_col, "_FLAG"]].dropna()
    sub = sub[sub["_FLAG"].isin([0.0, 1.0])]

    alt = sub[sub["_FLAG"] == 1.0]
    wt  = sub[sub["_FLAG"] == 0.0]

    if len(alt) < 2 or len(wt) < 2:
        raise ValueError(f"Too few samples: altered={len(alt)}, wildtype={len(wt)}")

    # Log-rank
    lr = logrank_test(
        alt[time_col], wt[time_col],
        event_observed_A=alt[event_col],
        event_observed_B=wt[event_col]
    )
    p_val = lr.p_value

    fig, ax = plt.subplots(figsize=(8, 5.5))
    fig.subplots_adjust(bottom=0.22)

    kmf_alt = KaplanMeierFitter()
    kmf_wt  = KaplanMeierFitter()

    kmf_alt.fit(alt[time_col], alt[event_col], label=f"Altered (n={len(alt)})")
    kmf_wt.fit(wt[time_col],  wt[event_col],  label=f"Wildtype (n={len(wt)})")

    kmf_alt.plot_survival_function(ax=ax, color="#C00000", ci_show=True, ci_alpha=0.12, linewidth=2)
    kmf_wt.plot_survival_function(ax=ax,  color="#1F4E79", ci_show=True, ci_alpha=0.12, linewidth=2)

    from lifelines.plotting import add_at_risk_counts
    add_at_risk_counts(kmf_wt, kmf_alt, ax=ax, fontsize=7, rows_to_show=["At risk"])

    ep_label = "Event-Free Survival" if ep == "EFS" else "Overall Survival"
    ax.set_title(f"{gene_label} — {ep_label}", fontsize=13, fontweight="bold", pad=10)
    ax.set_xlabel("Time (months)", fontsize=11)
    ax.set_ylabel("Survival Probability", fontsize=11)
    ax.set_ylim(-0.02, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    p_str = f"p < 0.001" if p_val < 0.001 else f"p = {p_val:.3f}"
    ax.text(0.98, 0.97, p_str, transform=ax.transAxes,
            ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#CCC", alpha=0.9))

    ax.legend(loc="lower left", fontsize=9, frameon=True)
    fig.tight_layout()

    return img_to_b64(fig)


def render_expr_km(gene, cohort, ep):
    """Render an expression KM plot and return base64 PNG."""
    study_id = COHORT_MAP.get(cohort)
    if study_id is None:
        raise ValueError(f"Unknown cohort: {cohort}")

    mrna_path = DATA_DIR / f"mrna_zscores_{study_id}.csv"
    clin_path  = DATA_DIR / f"clinical_{study_id}.csv"

    mrna = pd.read_csv(mrna_path)
    clin = pd.read_csv(clin_path)

    # Extract patient ID from sample ID
    mrna = mrna[mrna["hugoGeneSymbol"] == gene].copy()
    if mrna.empty:
        raise ValueError(f"No expression data for gene {gene} in {mrna_path.name}")

    mrna["patientId"] = mrna["sampleId"].apply(
        lambda s: re.sub(r"-\d+[A-Z]*$", "", str(s))
    )

    merged = mrna.merge(clin, on="patientId", how="inner")

    time_col  = f"{ep}_MONTHS"
    event_col = f"{ep}_STATUS"

    if time_col not in merged.columns or event_col not in merged.columns:
        raise ValueError(f"Missing {time_col}/{event_col} in clinical data")

    sub = merged[["value", time_col, event_col]].dropna()

    threshold = 1.5 if gene == "IGF2" else 1.0

    high = sub[sub["value"] >  threshold]
    low  = sub[sub["value"] <= threshold]

    if len(high) < 2 or len(low) < 2:
        raise ValueError(f"Too few samples: high={len(high)}, low={len(low)}")

    lr = logrank_test(
        high[time_col], low[time_col],
        event_observed_A=high[event_col],
        event_observed_B=low[event_col]
    )
    p_val = lr.p_value

    fig, ax = plt.subplots(figsize=(8, 5.5))
    fig.subplots_adjust(bottom=0.22)

    kmf_high = KaplanMeierFitter()
    kmf_low  = KaplanMeierFitter()

    kmf_high.fit(high[time_col], high[event_col], label=f"High expr (n={len(high)})")
    kmf_low.fit(low[time_col],   low[event_col],  label=f"Low expr (n={len(low)})")

    kmf_high.plot_survival_function(ax=ax, color="#E07B00", ci_show=True, ci_alpha=0.12, linewidth=2)
    kmf_low.plot_survival_function(ax=ax,  color="#1D9E75", ci_show=True, ci_alpha=0.12, linewidth=2)

    from lifelines.plotting import add_at_risk_counts
    add_at_risk_counts(kmf_low, kmf_high, ax=ax, fontsize=7, rows_to_show=["At risk"])

    ep_label = "Event-Free Survival" if ep == "EFS" else "Overall Survival"
    cohort_disp = cohort.replace("_", " ")
    ax.set_title(f"{gene} Expression — {cohort_disp} {ep_label}", fontsize=13, fontweight="bold", pad=10)
    ax.set_xlabel("Time (months)", fontsize=11)
    ax.set_ylabel("Survival Probability", fontsize=11)
    ax.set_ylim(-0.02, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    p_str = f"p < 0.001" if p_val < 0.001 else f"p = {p_val:.3f}"
    ax.text(0.98, 0.97, p_str, transform=ax.transAxes,
            ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#CCC", alpha=0.9))

    ax.legend(loc="lower left", fontsize=9, frameon=True)
    fig.tight_layout()

    return img_to_b64(fig)


# ── CSS ───────────────────────────────────────────────────────────────────────
CSS = """
:root {
  --navy: #1F3A5F; --blue: #2E75B6; --red: #C00000;
  --teal: #1D9E75; --lgray: #F5F7FA; --text: #2C2C2C;
}
body { font-family: 'Inter', -apple-system, sans-serif; color: var(--text); background: #fff; margin: 0; }
.container { max-width: 1200px; margin: 0 auto; padding: 0 24px; }
nav { position: sticky; top: 0; background: var(--navy); padding: 8px 24px; z-index: 100; display: flex; gap: 4px; flex-wrap: wrap; }
nav a { color: #C8D8F0; text-decoration: none; padding: 6px 12px; border-radius: 4px; font-size: 0.82rem; font-weight: 500; }
nav a:hover { background: rgba(255,255,255,0.15); color: #fff; }
.download-btn { margin-left: auto; background: var(--teal); color: #fff !important; padding: 6px 16px !important; border-radius: 4px; font-weight: 600 !important; }
.section-heading { color: var(--navy); border-left: 4px solid var(--blue); padding-left: 14px; margin: 2.5rem 0 1.2rem; font-size: 1.5rem; }
.data-table { width: 100%; border-collapse: collapse; font-size: 0.88rem; margin: 1rem 0; }
.data-table th { background: var(--navy); color: #fff; padding: 10px 12px; text-align: left; font-weight: 600; font-size: 0.82rem; }
.data-table td { padding: 8px 12px; border-bottom: 1px solid #E5E8ED; vertical-align: top; }
.data-table tr:nth-child(even) td { background: var(--lgray); }
.sig-row td { background: #FFF0F0 !important; }
.p-sig { color: var(--red); font-weight: 700; }
.p-ns { color: #888; }
.hr-bad { color: var(--red); font-weight: 600; }
.hr-good { color: var(--teal); font-weight: 600; }
.badge { display: inline-block; padding: 2px 7px; border-radius: 10px; font-size: 0.72rem; font-weight: 600; }
.badge-sig { background: #FFEBEB; color: var(--red); }
.badge-ns { background: #F0F0F0; color: #666; }
.badge-efs { background: #E3EDF9; color: var(--blue); }
.badge-os { background: #E5F5EE; color: var(--teal); }
.badge-skip { background: #F0F0F0; color: #888; font-size: 0.7rem; }
.card { background: #fff; border: 1px solid #D8DFE8; border-radius: 8px; padding: 20px 24px; margin: 16px 0; box-shadow: 0 2px 6px rgba(0,0,0,0.06); }
.finding-card { border-left: 5px solid var(--red); }
.summary-card { background: var(--lgray); border-radius: 8px; padding: 20px 24px; margin: 16px 0; }
.next-steps { background: #EAF4FB; border: 1px solid #BDD9EC; border-radius: 8px; padding: 20px 24px; margin: 16px 0; }
.method-note { background: #FFFDE7; border: 1px solid #F0D060; border-radius: 6px; padding: 14px 18px; margin: 12px 0; font-size: 0.87rem; }
.plot-grid { display: grid; grid-template-columns: repeat(2, 1fr); gap: 20px; margin: 1rem 0; }
.plot-card { border: 1px solid #D8DFE8; border-radius: 8px; overflow: hidden; }
.plot-card img { width: 100%; display: block; }
.plot-card-info { padding: 10px 14px; background: var(--lgray); font-size: 0.82rem; }
.placeholder-plot { background: #F5F5F5; border: 2px dashed #CCC; border-radius: 8px; padding: 32px 20px; text-align: center; color: #999; }
.bar-chart { margin: 1rem 0; }
.bar-row { display: flex; align-items: center; margin: 5px 0; }
.bar-label { width: 140px; text-align: right; padding-right: 10px; font-size: 0.82rem; color: #444; flex-shrink: 0; }
.bar-track { flex: 1; height: 22px; background: #F0F0F0; border-radius: 3px; position: relative; overflow: hidden; }
.bar-fill { height: 100%; border-radius: 3px; display: flex; align-items: center; padding-left: 6px; }
.bar-val { font-size: 0.74rem; color: #333; margin-left: 6px; white-space: nowrap; }
.legend { display: flex; flex-wrap: wrap; gap: 12px; margin: 12px 0; }
.legend-item { display: flex; align-items: center; gap: 5px; font-size: 0.8rem; }
.legend-dot { width: 12px; height: 12px; border-radius: 50%; flex-shrink: 0; }
.pathway-badge { display: inline-block; padding: 2px 8px; border-radius: 10px; font-size: 0.72rem; color: #fff; font-weight: 600; margin-right: 6px; }
.error-notice { background: #FFF0F0; border: 1px solid #F8BBBB; border-radius: 6px; padding: 12px 16px; color: var(--red); margin: 12px 0; }
.header-banner { background: var(--navy); color: #fff; padding: 36px 24px 30px; margin-bottom: 0; }
.header-banner h1 { font-size: 2.2rem; font-weight: 700; margin-bottom: 6px; color: #fff; }
.header-banner .subtitle { font-size: 1.2rem; color: #9DC3E6; margin-bottom: 12px; }
.header-banner .meta { font-size: 0.85rem; color: #8EB4D4; }
.header-banner .genes { font-size: 0.85rem; color: #A8C8E0; margin-top: 8px; }
.star-note { font-size: 0.75rem; color: #7FB4D0; }
.cohort-stat { text-align: center; padding: 12px; }
.cohort-stat .stat-num { font-size: 2rem; font-weight: 700; color: var(--navy); }
.cohort-stat .stat-label { font-size: 0.8rem; color: #666; margin-top: 2px; }
.stats-grid { display: grid; grid-template-columns: repeat(4, 1fr); gap: 0; border: 1px solid #D8DFE8; border-radius: 8px; overflow: hidden; margin: 16px 0; }
.stats-grid .cohort-stat { border-right: 1px solid #D8DFE8; }
.stats-grid .cohort-stat:last-child { border-right: none; }
.info-box { background: #EAF4FB; border: 1px solid #BDD9EC; border-radius: 8px; padding: 16px 20px; margin: 12px 0; font-size: 0.9rem; color: #1F3A5F; }
@media (max-width: 800px) {
  .plot-grid { grid-template-columns: 1fr; }
  .stats-grid { grid-template-columns: repeat(2,1fr); }
}
@media print {
  nav, .download-btn, .no-print { display: none !important; }
  .section-heading { page-break-before: always; }
  * { print-color-adjust: exact; -webkit-print-color-adjust: exact; }
  .container { max-width: 100%; padding: 0 12px; }
}
"""

# ── Section builders ──────────────────────────────────────────────────────────

def build_header():
    today = date.today().strftime("%B %d, %Y")
    return f"""
<div class="header-banner">
  <div class="container">
    <h1>Wilms Tumor Survival Analysis</h1>
    <div class="subtitle">Kaplan-Meier &amp; Cox Proportional Hazards — cBioPortal Cohort</div>
    <div class="meta">Generated: {today} &nbsp;|&nbsp; Data source: cBioPortal (Pediatric Wilms' Tumor, TARGET 2018)</div>
    <div class="genes">Genes analyzed: TP53, MDM4, MDM2, MYCN, SIX1, SIX2, DROSHA, DGCR8, CTNNB1, AMER1, WT1, MLLT1, NIPBL, IGF2 + OR classifiers</div>
    <div class="star-note" style="margin-top:10px">★ = PDMR focus gene</div>
  </div>
</div>
"""


def build_nav():
    return """
<nav>
  <a href="#cohorts">Cohorts</a>
  <a href="#landscape">Alteration Landscape</a>
  <a href="#km">KM Results</a>
  <a href="#or-classifiers">OR Classifiers</a>
  <a href="#expression">Expression</a>
  <a href="#plots">KM Plots</a>
  <a href="#summary">Summary</a>
  <a href="#methods">Methods</a>
  <a href="#" class="download-btn" onclick="window.print();return false;">⬇ Download PDF</a>
</nav>
"""


def build_cohorts_section():
    t1 = pd.read_csv(RESULTS_DIR / "Table1_Cohorts.csv")
    rows = ""
    for _, r in t1.iterrows():
        cells = "".join(f"<td>{esc(v)}</td>" for v in r)
        rows += f"<tr>{cells}</tr>\n"
    headers = "".join(f"<th>{esc(c)}</th>" for c in t1.columns)

    return f"""
<div class="container" id="cohorts">
  <h2 class="section-heading">1. Study Cohorts</h2>
  <table class="data-table">
    <thead><tr>{headers}</tr></thead>
    <tbody>{rows}</tbody>
  </table>

  <div class="stats-grid">
    <div class="cohort-stat">
      <div class="stat-num">652</div>
      <div class="stat-label">Total Patients</div>
    </div>
    <div class="cohort-stat">
      <div class="stat-num">652</div>
      <div class="stat-label">Patients with EFS</div>
    </div>
    <div class="cohort-stat">
      <div class="stat-num">211</div>
      <div class="stat-label">EFS Events</div>
    </div>
    <div class="cohort-stat">
      <div class="stat-num">114</div>
      <div class="stat-label">OS Events</div>
    </div>
  </div>
</div>
"""


def build_landscape_section():
    master = pd.read_csv(MASTER_CSV)
    total = len(master)

    alt_cols = [c for c in master.columns if c.endswith("_ALTERED")]
    rows_data = []
    for col in alt_cols:
        s = master[col].astype("Float64")
        n_alt = int((s == 1.0).sum())
        n_wt  = int((s == 0.0).sum())
        n_uns = total - n_alt - n_wt
        freq  = n_alt / total * 100
        display = GENE_DISPLAY.get(col, col.replace("_ALTERED", ""))
        pathway = GENE_PATHWAYS.get(col, "—")
        rows_data.append({
            "col": col, "display": display, "pathway": pathway,
            "n_alt": n_alt, "n_wt": n_wt, "n_uns": n_uns, "freq": freq
        })

    rows_data.sort(key=lambda x: x["n_alt"], reverse=True)
    max_n = max(r["n_alt"] for r in rows_data) if rows_data else 1

    # Bar chart
    bar_rows = ""
    for r in rows_data:
        color = PATHWAY_COLORS.get(r["pathway"], "#999")
        pct   = r["n_alt"] / max_n * 100 if max_n > 0 else 0
        bar_rows += f"""
    <div class="bar-row">
      <div class="bar-label">{esc(r['display'])}</div>
      <div class="bar-track">
        <div class="bar-fill" style="width:{pct:.1f}%;background:{color}">
          <span class="bar-val">{r['n_alt']} ({r['freq']:.1f}%)</span>
        </div>
      </div>
    </div>"""

    # Legend
    seen = {}
    for r in rows_data:
        p = r["pathway"]
        if p not in seen:
            seen[p] = PATHWAY_COLORS.get(p, "#999")
    legend_html = '<div class="legend">'
    for pathway, color in seen.items():
        legend_html += f'<div class="legend-item"><div class="legend-dot" style="background:{color}"></div>{esc(pathway)}</div>'
    legend_html += "</div>"

    # Table
    table_rows = ""
    for i, r in enumerate(rows_data):
        note = ""
        if r["n_alt"] == 0:
            note = '<span style="color:#888"> (0 — not detected in cohort)</span>'
        elif r["n_alt"] < 20:
            note = '<span style="color:#888"> (below KM gate)</span>'
        tr_class = "sig-row" if i % 2 == 0 else ""
        table_rows += f"""
    <tr class="{tr_class}">
      <td>{esc(r['display'])}{note}</td>
      <td>{pathway_badge(r['pathway'])}</td>
      <td>{r['n_alt']}</td>
      <td>{r['n_wt']}</td>
      <td>{r['n_uns']}</td>
      <td>{r['freq']:.1f}%</td>
    </tr>"""

    return f"""
<div class="container" id="landscape">
  <h2 class="section-heading">2. Alteration Landscape</h2>
  {legend_html}
  <div class="bar-chart">{bar_rows}</div>

  <table class="data-table" style="margin-top:2rem">
    <thead>
      <tr>
        <th>Gene</th><th>Pathway</th><th>N Altered</th>
        <th>N Wildtype</th><th>N Unsequenced</th><th>Freq %</th>
      </tr>
    </thead>
    <tbody>{table_rows}</tbody>
  </table>
</div>
"""


def build_km_section():
    t2 = pd.read_csv(RESULTS_DIR / "Table2_KM_Cox_Results.csv")
    cox = pd.read_csv(RESULTS_DIR / "cox_results.csv")

    # Main results table
    table_rows = ""
    for _, r in t2.iterrows():
        sig_val = str(r.get("Significant?", "")).strip()
        is_sig = sig_val.lower() in ("yes", "true", "1")
        row_cls = "sig-row" if is_sig else ""

        gene_val  = esc(r.get("Gene", ""))
        path_val  = esc(r.get("Pathway", ""))
        ep_val    = str(r.get("Endpoint", ""))
        ep_badge  = f'<span class="badge badge-efs">EFS</span>' if "EFS" in ep_val else f'<span class="badge badge-os">OS</span>'
        n_alt     = esc(r.get("N altered", ""))
        n_wt      = esc(r.get("N wildtype", ""))
        med_alt   = fmt_med(r.get("Median altered (mo)", ""))
        med_wt    = fmt_med(r.get("Median wildtype (mo)", ""))
        lr_p_raw  = r.get("Log-rank p", "")
        lr_p      = esc(lr_p_raw) if not pd.isna(lr_p_raw) else "—"
        if str(lr_p_raw).strip() in ("<0.001", "< 0.001"):
            lr_p_disp = f'<span class="p-sig">&lt;0.001</span>'
        elif "0.001" in str(lr_p_raw) and float(str(lr_p_raw).replace("<","").replace(">","").strip() or "1") < 0.05:
            lr_p_disp = f'<span class="p-sig">{lr_p}</span>'
        else:
            lr_p_disp = f'<span class="p-ns">{lr_p}</span>'

        hr_raw  = r.get("HR (multivariate)", "")
        ci_raw  = r.get("95% CI", "")
        cox_p_raw = r.get("Cox p", "")
        hr_str = "—"
        if not pd.isna(hr_raw) and str(hr_raw).strip() not in ("", "—"):
            ci_str = esc(ci_raw) if not pd.isna(ci_raw) else "—"
            try:
                hr_f = float(str(hr_raw).strip())
                hr_color = "hr-bad" if hr_f > 1 else "hr-good"
                hr_str = f'<span class="{hr_color}">{hr_f:.2f}</span> {ci_str}'
            except Exception:
                hr_str = f"{esc(hr_raw)} {esc(ci_raw)}"

        sig_badge = '<span class="badge badge-sig">★ Sig</span>' if is_sig else '<span class="badge badge-ns">ns</span>'

        stat_val = str(r.get("Alteration Type", "")).lower()
        if "skip" in stat_val or "below" in stat_val:
            sig_badge = '<span class="badge badge-skip">skipped</span>'

        table_rows += f"""
    <tr class="{row_cls}">
      <td><strong>{gene_val}</strong></td>
      <td>{path_val}</td>
      <td>{ep_badge}</td>
      <td>{n_alt}</td><td>{n_wt}</td>
      <td>{med_alt}</td><td>{med_wt}</td>
      <td>{lr_p_disp}</td>
      <td>{hr_str}</td>
      <td>{sig_badge}</td>
    </tr>"""

    # Callout cards for significant findings
    cards = ""
    sig_cox = cox[(cox.get("status", cox.get("status", "")) == "completed") &
                  (cox["significant"] == True)].copy()
    if "HR" in sig_cox.columns:
        sig_cox = sig_cox.sort_values("HR", ascending=False)

    for _, r in sig_cox.iterrows():
        gene    = str(r["gene"])
        pathway = str(r.get("pathway", ""))
        ep      = str(r["endpoint"])
        ep_badge_cls = "badge-efs" if ep == "EFS" else "badge-os"
        med_alt = fmt_med(r.get("median_altered_months", ""))
        med_wt  = fmt_med(r.get("median_wildtype_months", ""))
        p_str   = fmt_p(r.get("cox_p", r.get("logrank_p", "")))
        hr_val  = r.get("HR", "")
        ci_lo   = r.get("CI_lower", "")
        ci_hi   = r.get("CI_upper", "")
        hr_str2 = fmt_hr_ci(hr_val, ci_lo, ci_hi)

        ep_full = "event-free survival" if ep == "EFS" else "overall survival"
        sentence = (
            f"<strong>{esc(gene)}</strong> pathway alteration was associated with significantly worse "
            f"{ep_full} (median {ep}: {esc(med_alt)} vs. {esc(med_wt)} months; "
            f"log-rank p {esc(p_str)}; HR = {esc(hr_str2)})."
        )

        # Median bar viz
        try:
            ma = float(med_alt) if med_alt != "NR" else None
            mw = float(med_wt) if med_wt != "NR" else None
        except Exception:
            ma, mw = None, None

        bar_viz = ""
        if ma is not None or mw is not None:
            ref = max(v for v in [ma, mw] if v is not None)
            if ref > 0:
                pct_a = (ma / ref * 100) if ma is not None else 0
                pct_w = (mw / ref * 100) if mw is not None else 0
                bar_viz = f"""
        <div style="margin-top:12px">
          <div class="bar-row">
            <div class="bar-label">Altered</div>
            <div class="bar-track">
              <div class="bar-fill" style="width:{pct_a:.1f}%;background:#C00000">
                <span class="bar-val">{esc(med_alt)} mo</span>
              </div>
            </div>
          </div>
          <div class="bar-row">
            <div class="bar-label">Wildtype</div>
            <div class="bar-track">
              <div class="bar-fill" style="width:{pct_w:.1f}%;background:#1F4E79">
                <span class="bar-val">{esc(med_wt)} mo</span>
              </div>
            </div>
          </div>
        </div>"""

        cards += f"""
    <div class="card finding-card">
      <div style="font-size:1.4rem;font-weight:700;color:var(--navy);margin-bottom:6px">{esc(gene)}</div>
      {pathway_badge(pathway)}
      <span class="badge {ep_badge_cls}">{esc(ep)}</span>
      <p style="margin:12px 0 6px">{sentence}</p>
      {bar_viz}
    </div>"""

    return f"""
<div class="container" id="km">
  <h2 class="section-heading">3. KM &amp; Cox Results</h2>
  <div style="overflow-x:auto">
  <table class="data-table">
    <thead>
      <tr>
        <th>Gene</th><th>Pathway</th><th>Endpoint</th>
        <th>N Alt</th><th>N WT</th>
        <th>Median Alt (mo)</th><th>Median WT (mo)</th>
        <th>Log-rank p</th><th>HR [95%CI]</th><th>Sig?</th>
      </tr>
    </thead>
    <tbody>{table_rows}</tbody>
  </table>
  </div>

  <h3 style="color:var(--navy);margin-top:2rem">Significant Findings</h3>
  {cards if cards else '<div class="info-box">No significant individual-gene findings.</div>'}
</div>
"""


def build_or_section():
    km = pd.read_csv(RESULTS_DIR / "km_results.csv")
    cox = pd.read_csv(RESULTS_DIR / "cox_results.csv")

    or_km  = km[km["gene"].str.endswith("_OR", na=False)].copy()
    or_cox = cox[cox["gene"].str.endswith("_OR", na=False)].copy()

    total = 652

    table_rows = ""
    for _, r in or_km.iterrows():
        gene  = str(r["gene"])
        genes_combined = OR_GENE_LABELS.get(gene, gene)
        n_alt = r.get("n_altered", "")
        n_wt  = r.get("n_wildtype", "")
        try:
            pct_cap = f"{float(n_alt) / total * 100:.1f}%"
        except Exception:
            pct_cap = "—"
        med_alt = fmt_med(r.get("median_altered_months", ""))
        med_wt  = fmt_med(r.get("median_wildtype_months", ""))
        lr_p    = fmt_p(r.get("logrank_p", ""))
        ep      = str(r.get("endpoint", ""))

        # get HR from cox
        cox_row = or_cox[(or_cox["gene"] == gene) & (or_cox["endpoint"] == ep)]
        if len(cox_row) > 0:
            crow = cox_row.iloc[0]
            hr_str = fmt_hr_ci(crow.get("HR", ""), crow.get("CI_lower", ""), crow.get("CI_upper", ""))
        else:
            hr_str = "—"

        sig = str(r.get("significant", "")).lower() in ("true", "yes", "1")
        sig_badge = '<span class="badge badge-sig">★ Sig</span>' if sig else '<span class="badge badge-ns">ns</span>'
        row_cls = "sig-row" if sig else ""
        ep_badge_cls = "badge-efs" if ep == "EFS" else "badge-os"

        table_rows += f"""
    <tr class="{row_cls}">
      <td><strong>{esc(gene)}</strong></td>
      <td><em>{esc(genes_combined)}</em></td>
      <td><span class="badge {ep_badge_cls}">{esc(ep)}</span></td>
      <td>{esc(str(int(float(n_alt)) if not pd.isna(n_alt) else '—'))}</td>
      <td>{esc(str(int(float(n_wt)) if not pd.isna(n_wt) else '—'))}</td>
      <td>{esc(pct_cap)}</td>
      <td>{esc(med_alt)}</td><td>{esc(med_wt)}</td>
      <td><span class="{'p-sig' if lr_p not in ('—','') and _p_is_sig(lr_p) else 'p-ns'}">{esc(lr_p)}</span></td>
      <td>{esc(hr_str)}</td>
      <td>{sig_badge}</td>
    </tr>"""

    return f"""
<div class="container" id="or-classifiers">
  <h2 class="section-heading">4. OR-Classifier Analysis</h2>
  <div class="method-note">
    OR-classifiers aggregate all patients with <em>any</em> alteration within the pathway.
    A significant OR-classifier result indicates pathway-level prognostic impact even when
    individual genes are too rare to analyze independently.
  </div>
  <div style="overflow-x:auto">
  <table class="data-table">
    <thead>
      <tr>
        <th>Classifier</th><th>Genes Combined</th><th>Endpoint</th>
        <th>N Altered</th><th>N Wildtype</th><th>% Cohort Captured</th>
        <th>Median Alt (mo)</th><th>Median WT (mo)</th>
        <th>LR p</th><th>HR [95%CI]</th><th>Sig?</th>
      </tr>
    </thead>
    <tbody>{table_rows}</tbody>
  </table>
  </div>
</div>
"""


def build_expression_section():
    t3 = pd.read_csv(RESULTS_DIR / "Table3_Expression_Results.csv")

    # Sort: completed first by p-value, then skipped
    def sort_key(r):
        sig_raw = str(r.get("Significant?", "")).lower()
        p_raw   = r.get("Log-rank p", 1.0)
        try:
            p_f = float(str(p_raw).strip())
        except Exception:
            p_f = 1.0
        skipped = "skip" in str(r.get("N high", "")).lower() or p_f == 1.0
        return (1 if skipped else 0, p_f)

    t3_sorted = t3.copy()
    try:
        t3_sorted["_sort"] = t3_sorted.apply(sort_key, axis=1)
        t3_sorted = t3_sorted.sort_values("_sort").drop(columns=["_sort"])
    except Exception:
        pass

    table_rows = ""
    for _, r in t3_sorted.iterrows():
        sig_val = str(r.get("Significant?", "")).lower()
        is_sig  = sig_val in ("yes", "true", "1")
        row_cls = "sig-row" if is_sig else ""
        gene    = esc(r.get("Gene", ""))
        cohort  = esc(r.get("Cohort", ""))
        ep      = str(r.get("Endpoint", ""))
        ep_badge_cls = "badge-efs" if ep == "EFS" else "badge-os"
        n_high  = esc(r.get("N high", ""))
        n_low   = esc(r.get("N low", ""))
        med_h   = fmt_med(r.get("Median HIGH (mo)", ""))
        med_l   = fmt_med(r.get("Median LOW (mo)", ""))
        lr_p    = esc(r.get("Log-rank p", "—"))
        hr_val  = r.get("HR (vs. LOW)", "")
        ci_val  = esc(r.get("95% CI", ""))
        hr_str  = f"{esc(str(hr_val))} {ci_val}" if not pd.isna(hr_val) else "—"
        sig_badge = '<span class="badge badge-sig">★ Sig</span>' if is_sig else '<span class="badge badge-ns">ns</span>'

        table_rows += f"""
    <tr class="{row_cls}">
      <td><strong>{gene}</strong></td>
      <td>{cohort}</td>
      <td><span class="badge {ep_badge_cls}">{esc(ep)}</span></td>
      <td>{n_high}</td><td>{n_low}</td>
      <td>{med_h}</td><td>{med_l}</td>
      <td>{lr_p}</td><td>{hr_str}</td>
      <td>{sig_badge}</td>
    </tr>"""

    no_sig_note = ""
    if t3_sorted["Significant?"].astype(str).str.lower().isin(["yes", "true", "1"]).sum() == 0:
        no_sig_note = """
    <div class="info-box" style="margin-top:12px">
      No expression analyses reached statistical significance (p&lt;0.05). Notable trend:
      MDM4 high expression showed a borderline EFS association (p=0.059, HR=1.40 [0.94–2.07]).
    </div>"""

    return f"""
<div class="container" id="expression">
  <h2 class="section-heading">5. Expression Analysis</h2>
  <div class="method-note">
    mRNA z-scores from cBioPortal were used to stratify patients into HIGH (z&gt;threshold)
    and LOW (z≤threshold) expression groups. Default threshold z=1.0; IGF2 uses z=1.5.
    Analyses performed only for PDMR-focus genes with sufficient sample sizes.
  </div>
  {no_sig_note}
  <div style="overflow-x:auto">
  <table class="data-table">
    <thead>
      <tr>
        <th>Gene</th><th>Cohort</th><th>Endpoint</th>
        <th>N High</th><th>N Low</th>
        <th>Median High (mo)</th><th>Median Low (mo)</th>
        <th>Log-rank p</th><th>HR [95%CI]</th><th>Sig?</th>
      </tr>
    </thead>
    <tbody>{table_rows}</tbody>
  </table>
  </div>
</div>
"""


def build_plots_section(master_df):
    # Collect all stems
    stems = []
    for pdf in sorted(PLOTS_DIR.glob("KM_*.pdf")):
        stem = pdf.stem
        parsed = parse_stem(stem)
        if parsed:
            stems.append((stem, parsed))

    # Determine significance order
    km  = pd.read_csv(RESULTS_DIR / "km_results.csv")
    expr = pd.read_csv(RESULTS_DIR / "expression_km_results.csv")

    sig_genomic = set(
        km[(km["significant"] == True) & (~km["gene"].str.endswith("_OR", na=False))]["gene"].tolist()
    )
    sig_or = set(
        km[(km["significant"] == True) & (km["gene"].str.endswith("_OR", na=False))]["gene"].tolist()
    )

    def sort_key_plot(item):
        stem, p = item
        if p["type"] == "genomic":
            gene = p["gene"]
            if gene in sig_or or gene in sig_genomic:
                return (0, stem)
            return (1, stem)
        else:
            return (2, stem)

    stems.sort(key=sort_key_plot)

    grid_items = ""
    skipped_rows = ""
    n_rendered = 0
    n_failed = 0

    for stem, parsed in stems:
        ptype = parsed["type"]
        gene  = parsed["gene"]
        ep    = parsed.get("ep", "")
        cohort = parsed.get("cohort", "")

        ep_badge_cls = "badge-efs" if ep == "EFS" else "badge-os"
        ep_badge_html = f'<span class="badge {ep_badge_cls}">{ep}</span>'

        try:
            if ptype == "genomic":
                b64 = render_genomic_km(gene, ep, master_df)
                title = f"{gene.replace('_OR','').replace('_',' ')} — {ep}"
            else:
                b64 = render_expr_km(gene, cohort, ep)
                title = f"{gene} Expression ({cohort.replace('_',' ')}) — {ep}"

            grid_items += f"""
      <div class="plot-card">
        <img src="data:image/png;base64,{b64}" alt="{esc(stem)}">
        <div class="plot-card-info">
          <strong>{esc(title)}</strong> {ep_badge_html}
          <span style="float:right;color:#888;font-size:0.75rem">{esc(stem)}</span>
        </div>
      </div>"""
            n_rendered += 1

        except Exception as e:
            n_failed += 1
            skipped_rows += f"<tr><td>{esc(stem)}</td><td>{esc(str(e))}</td></tr>"
            grid_items += f"""
      <div class="placeholder-plot">
        <div><strong>{esc(stem)}</strong></div>
        <div style="margin-top:8px;font-size:0.8rem;color:#B00">{esc(str(e))}</div>
      </div>"""

    skipped_table = ""
    if skipped_rows:
        skipped_table = f"""
    <h3 style="color:var(--navy);margin-top:2rem">Skipped / Failed Plots</h3>
    <table class="data-table">
      <thead><tr><th>Stem</th><th>Reason</th></tr></thead>
      <tbody>{skipped_rows}</tbody>
    </table>"""

    return f"""
<div class="container" id="plots">
  <h2 class="section-heading">6. Kaplan-Meier Plots</h2>
  <p style="color:#666;font-size:0.88rem">
    {n_rendered} plots rendered successfully; {n_failed} skipped.
    Significant OR-classifier plots shown first.
  </p>
  <div class="plot-grid">{grid_items}</div>
  {skipped_table}
</div>
"""


def build_summary_section():
    cox  = pd.read_csv(RESULTS_DIR / "cox_results.csv")
    expr = pd.read_csv(RESULTS_DIR / "expression_km_results.csv")

    sig_cox = cox[cox["significant"] == True].copy()
    or_sig  = sig_cox[sig_cox["gene"].str.endswith("_OR", na=False)]
    ind_sig = sig_cox[~sig_cox["gene"].str.endswith("_OR", na=False)]

    # OR classifier summary
    or_bullets = ""
    for _, r in or_sig.iterrows():
        gene    = r["gene"]
        ep      = r["endpoint"]
        hr      = fmt_hr_ci(r.get("HR", ""), r.get("CI_lower", ""), r.get("CI_upper", ""))
        p       = fmt_p(r.get("cox_p", ""))
        med_alt = fmt_med(r.get("median_altered_months", ""))
        med_wt  = fmt_med(r.get("median_wildtype_months", ""))
        genes_combined = OR_GENE_LABELS.get(gene, gene)
        or_bullets += f"""
      <li><strong>{esc(gene)}</strong> ({esc(genes_combined)}): {esc(ep)} median {esc(med_alt)} vs. {esc(med_wt)} months; HR = {esc(hr)}, p {esc(p)}</li>"""

    # Individual gene summary
    ind_bullets = ""
    for _, r in ind_sig.iterrows():
        gene = r["gene"]
        ep   = r["endpoint"]
        hr   = fmt_hr_ci(r.get("HR", ""), r.get("CI_lower", ""), r.get("CI_upper", ""))
        p    = fmt_p(r.get("cox_p", ""))
        ind_bullets += f"<li><strong>{esc(gene)}</strong> — {esc(ep)}, HR = {esc(hr)}, p {esc(p)}</li>"

    # Non-significant
    ns_cox = cox[(cox.get("status", cox.get("status","")) == "completed") &
                 (cox["significant"] != True)]
    ns_genes = ", ".join(esc(str(g)) for g in ns_cox["gene"].unique() if not str(g).endswith("_OR"))

    # Expression summary
    expr_trends = expr[expr["status"] == "completed"].copy()
    if "logrank_p" in expr_trends.columns:
        expr_trends = expr_trends.sort_values("logrank_p")

    expr_bullets = ""
    for _, r in expr_trends.iterrows():
        g  = r["gene"]
        ep = r["endpoint"]
        p  = fmt_p(r.get("logrank_p", ""))
        hr = fmt_hr_ci(r.get("HR", ""), r.get("CI_lower", ""), r.get("CI_upper", ""))
        expr_bullets += f"<li><strong>{esc(g)}</strong> ({esc(ep)}): p = {esc(p)}, HR = {esc(hr)}</li>"

    n_sig = len(sig_cox)
    today = date.today().strftime("%B %d, %Y")

    # Next steps
    next_steps = """
    <ul>
      <li>Validate TP53/MDM4/MDM2 pathway findings in an independent cohort (e.g., TARGET GDC or SIOP datasets)</li>
      <li>Investigate MDM4 expression trend (p=0.059) with larger sample sizes or alternative z-score thresholds</li>
      <li>Perform multivariable analysis controlling for stage and histology alongside molecular markers</li>
      <li>Examine interaction between OR-classifier status and treatment protocol (NWTS-4 vs. NTWS-4)</li>
      <li>Consider functional validation of MYCN/SIX1-SIX2 axis findings in Wilms tumor preclinical models</li>
      <li>Explore IGF2 expression stratification with a CNA-aware approach</li>
    </ul>"""

    return f"""
<div class="container" id="summary">
  <h2 class="section-heading">7. Summary of Findings</h2>

  <div class="summary-card">
    <h3 style="margin-top:0;color:var(--navy)">Significant Associations ({n_sig} results)</h3>
    <p>Across 1 Wilms tumor cohort from cBioPortal (N = 652, Pediatric Wilms' Tumor TARGET 2018),
    Kaplan-Meier and Cox proportional hazards analysis identified the following statistically
    significant survival associations:</p>

    <h4>OR-Classifier Findings</h4>
    <ul>{or_bullets if or_bullets else '<li>None identified.</li>'}</ul>

    <h4>Individual Gene Findings</h4>
    <ul>{ind_bullets if ind_bullets else '<li>No individual genes reached significance at the n≥20 threshold.</li>'}</ul>
  </div>

  <div class="summary-card">
    <h3 style="margin-top:0;color:var(--navy)">Non-Significant Genes</h3>
    <p>The following genes were analyzed but did not reach statistical significance (p≥0.05)
    or were below the minimum-events gate: {ns_genes if ns_genes else 'None below threshold.'}.</p>
    <p>Note: Some genes were excluded from KM analysis due to insufficient events (n_altered &lt; 20
    required by the pre-specified gate).</p>
  </div>

  <div class="summary-card">
    <h3 style="margin-top:0;color:var(--navy)">Expression Analysis Summary</h3>
    <p>mRNA z-score–based survival analysis did not identify statistically significant associations
    after multiple testing consideration. Noteworthy trends:</p>
    <ul>{expr_bullets if expr_bullets else '<li>No expression data available.</li>'}</ul>
    <p><em>MDM4 high expression showed a borderline EFS association (p≈0.059) that warrants
    follow-up in a larger expression cohort.</em></p>
  </div>

  <div class="next-steps">
    <h3 style="margin-top:0;color:#1F3A5F">Recommended Next Steps</h3>
    {next_steps}
  </div>
</div>
"""


def build_methods_section():
    today = date.today().strftime("%B %d, %Y")
    return f"""
<div class="container" id="methods">
  <h2 class="section-heading">8. Methods</h2>
  <div class="card" style="font-size:0.85rem;line-height:1.7">
    <p><strong>Data source:</strong> Clinical and genomic data were retrieved from cBioPortal
    (cbioportal.org) for the Pediatric Wilms' Tumor (TARGET, 2018) study
    (study ID: wt_target_2018_pub; N=652 patients).</p>

    <p><strong>Genomic alteration classification:</strong> Mutations and copy-number alterations
    were combined into a binary "altered/wildtype" flag per gene. OR-classifiers were constructed
    by flagging any patient with at least one alteration in the specified gene set.</p>

    <p><strong>Minimum events gate:</strong> Individual gene analyses required ≥20 patients in the
    altered group to ensure statistical power. Genes below this threshold were excluded from KM
    and Cox analyses.</p>

    <p><strong>Survival analysis:</strong> Kaplan-Meier curves were estimated using the lifelines
    library (v0.27+). Log-rank tests were used for univariate comparison. Cox proportional hazards
    models were fit using age as a covariate (ALTERED + AGE). Primary endpoint was event-free
    survival (EFS); overall survival (OS) was a secondary endpoint.</p>

    <p><strong>Expression analysis:</strong> mRNA z-scores (relative to diploid samples) were
    obtained from cBioPortal. Patients were stratified into HIGH (z&gt;1.0; z&gt;1.5 for IGF2)
    and LOW (z≤threshold) groups. Log-rank and Cox analyses were performed as above.</p>

    <p><strong>Report generated:</strong> {today} using Python (pandas, lifelines, matplotlib).
    All analyses were conducted using the cBioPortal REST API and publicly available data.</p>
  </div>
</div>
"""


def build_build_log_section():
    rows = ""
    for entry in build_log:
        icon  = entry[0]
        name  = entry[1]
        error = entry[2] if len(entry) > 2 else ""
        color = "var(--teal)" if icon == "✓" else "var(--red)"
        rows += f"""
    <tr>
      <td style="color:{color};font-weight:700">{esc(icon)}</td>
      <td>{esc(name)}</td>
      <td style="color:#888;font-size:0.8rem">{esc(error)}</td>
    </tr>"""

    return f"""
<div class="container" style="margin-top:2rem;margin-bottom:3rem">
  <h2 class="section-heading">Build Log</h2>
  <table class="data-table">
    <thead><tr><th>Status</th><th>Section</th><th>Error</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
</div>
"""


# ── Main assembler ────────────────────────────────────────────────────────────
def main():
    print("Loading master_classification.csv ...")
    master_df = pd.read_csv(MASTER_CSV)
    print(f"  → {len(master_df)} rows loaded")

    print("Building HTML sections ...")
    header_html  = safe_section("Header Banner",       build_header)
    cohorts_html = safe_section("Cohorts",             build_cohorts_section)
    land_html    = safe_section("Alteration Landscape",build_landscape_section)
    km_html      = safe_section("KM & Cox Results",    build_km_section)
    or_html      = safe_section("OR Classifiers",      build_or_section)
    expr_html    = safe_section("Expression Analysis", build_expression_section)

    print("Rendering KM plots (this may take a minute) ...")
    plots_html   = safe_section("KM Plots",            build_plots_section, master_df)

    summary_html = safe_section("Summary",             build_summary_section)
    methods_html = safe_section("Methods",             build_methods_section)
    log_html     = build_build_log_section()

    print("Assembling final HTML ...")
    today = date.today().strftime("%B %d, %Y")
    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Wilms Tumor cBioPortal Analysis — {today}</title>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap" rel="stylesheet">
<style>
{CSS}
</style>
</head>
<body>
{build_nav()}
{header_html}
{cohorts_html}
{land_html}
{km_html}
{or_html}
{expr_html}
{plots_html}
{summary_html}
{methods_html}
{log_html}
</body>
</html>
"""

    OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    OUT_HTML.write_text(html_doc, encoding="utf-8")

    size_kb = OUT_HTML.stat().st_size / 1024
    print(f"\n✓ HTML report written to:")
    print(f"  {OUT_HTML}")
    print(f"  Size: {size_kb:.1f} KB")
    print("\nBuild log:")
    for entry in build_log:
        icon = entry[0]
        name = entry[1]
        err  = entry[2] if len(entry) > 2 else ""
        print(f"  {icon} {name}" + (f" — {err}" if err else ""))


if __name__ == "__main__":
    main()
