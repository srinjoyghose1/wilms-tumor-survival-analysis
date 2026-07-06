"""
10_pathway_classification.py

Extends the master gene-level classification (04_classify_patients.py) with 23
additional genes and rolls everything up into five pathway-level OR-flags.

STEP 1 — download mutations + GISTIC CNA for the extended gene list from
         wt_target_2018_pub via cBioPortal's GET endpoints.
STEP 2 — classify each extended gene ALTERED / WILDTYPE / NaN using the same
         wildtype-if-sequenced / NaN-if-absent logic as 04_classify_patients.py.
STEP 3 — merge with processed/master_classification.csv and build five
         pathway OR-flags (WNT, SIX/miRNA, CHROMATIN, P53, IGF2) plus
         PATHWAY_ANY.
STEP 4 — print a per-pathway and per-gene alteration summary.

API quirks discovered in probe
-------------------------------
  GET /molecular-profiles/{id}/mutations           entrezGeneId is REQUIRED
      -> must be called once per gene (23 calls).
  GET /molecular-profiles/{id}/discrete-copy-number has NO entrezGeneId
      parameter at all (confirmed against the cBioPortal OpenAPI spec) -> it
      always returns every gene for the sample list. Since data/cna_{STUDY_ID}.csv
      already holds that exact full-genome pull (same profile, same sample
      list, discreteCopyNumberEventType=ALL — see 03_download_genomics.py),
      we reuse it instead of re-downloading several hundred MB of JSON; the
      live GET call is used only as a fallback if that file is missing.

Outputs
-------
  wilms_survival/data/mutations_extended.csv
  wilms_survival/data/cna_extended.csv
  wilms_survival/processed/master_classification_extended.csv
"""

import os
import re
import sys
import requests
import pandas as pd

BASE_URL   = "https://www.cbioportal.org/api"
SCRIPT_DIR = os.path.dirname(__file__)
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "data")
PROC_DIR   = os.path.join(SCRIPT_DIR, "..", "processed")

STUDY_ID = "wt_target_2018_pub"

EXTENDED_GENES = [
    "BCOR", "BCORL1", "NF1", "NONO", "COL6A3", "ARID1A", "MAX", "MAP3K4",
    "ASXL1", "BRD7", "FGFR1", "CHD4", "HDAC4", "ACTB", "CREBBP", "EP300",
    "DICER1", "XPO5", "DIS3L2", "SMARCA4", "CHEK2", "PALB2", "CTR9",
]

TRUNCATING = frozenset([
    "Nonsense_Mutation",
    "Frame_Shift_Ins",
    "Frame_Shift_Del",
    "Splice_Site",
    "Splice_Region",
    "Translation_Start_Site",
    "Nonstop_Mutation",
])

# DICER1 RNase IIIb metal-binding hotspot codons (D1709, E1705, D1810, E1813)
DICER1_HOTSPOT_CODONS = frozenset([1705, 1709, 1810, 1813])

# Gene groupings per the classification rules
TUMOR_SUPPRESSORS_TRUNC = [
    "BCOR", "BCORL1", "NF1", "ARID1A", "SMARCA4", "BRD7", "CHD4", "CREBBP",
    "EP300", "ASXL1", "HDAC4", "CTR9", "CHEK2", "PALB2", "XPO5", "DIS3L2",
]
ALL_MUT_TYPES_WITH_DEL = ["MAP3K4", "NONO"]
ONCOGENES_ANY_MUT_AMP  = ["MAX", "FGFR1"]
# DICER1 handled separately (domain-restricted truncating, no CNA)
# COL6A3 handled separately (CNA only, both directions)
# ACTB handled separately (mutation only, no CNA)

PATHWAYS = {
    "PATHWAY_WNT": [
        "WT1_ALTERED", "CTNNB1_ALTERED", "AMER1_ALTERED", "NF1_ALTERED",
    ],
    "PATHWAY_SIX_miRNA": [
        "SIX1_ALTERED", "SIX2_ALTERED", "DROSHA_ALTERED", "DGCR8_ALTERED",
        "DICER1_ALTERED", "XPO5_ALTERED", "DIS3L2_ALTERED",
    ],
    "PATHWAY_CHROMATIN": [
        "MLLT1_ALTERED", "BCOR_ALTERED", "BCORL1_ALTERED", "ARID1A_ALTERED",
        "SMARCA4_ALTERED", "BRD7_ALTERED", "CREBBP_ALTERED", "EP300_ALTERED",
        "HDAC4_ALTERED", "ASXL1_ALTERED", "CHD4_ALTERED", "MAP3K4_ALTERED",
        "CTR9_ALTERED",
    ],
    "PATHWAY_P53": [
        "TP53_ALTERED", "MYCN_ALTERED", "MDM2_ALTERED", "MDM4_ALTERED",
        "MAX_ALTERED", "CHEK2_ALTERED", "PALB2_ALTERED",
    ],
    "PATHWAY_IGF2": [
        "IGF2_CNA_ALTERED", "IGF2_EXPR_ALTERED", "NIPBL_ALTERED",
        "NONO_ALTERED", "COL6A3_ALTERED", "FGFR1_ALTERED", "ACTB_ALTERED",
    ],
}


# ── HTTP helpers ────────────────────────────────────────────────────────────

def _get(path: str, params: dict | None = None) -> list | dict | None:
    url = f"{BASE_URL}/{path.lstrip('/')}"
    try:
        r = requests.get(url, params=params, timeout=120)
        r.raise_for_status()
        return r.json()
    except requests.exceptions.HTTPError as e:
        if e.response.status_code != 404:
            print(f"  HTTP {e.response.status_code} — {url}", file=sys.stderr)
        return None
    except requests.exceptions.RequestException as e:
        print(f"  Request error — {e}", file=sys.stderr)
        return None


def get_sample_lists(study_id: str) -> list[dict]:
    data = _get(f"studies/{study_id}/sample-lists")
    return data if isinstance(data, list) else []


def pick_sample_list(sample_lists: list[dict], priority_categories: list[str]) -> str | None:
    by_cat = {sl["category"]: sl["sampleListId"] for sl in sample_lists}
    for cat in priority_categories:
        if cat in by_cat:
            return by_cat[cat]
    return by_cat.get("all_cases_in_study") or (sample_lists[0]["sampleListId"] if sample_lists else None)


def resolve_entrez_ids(gene_symbols: list[str]) -> dict[str, int]:
    sym_to_entrez: dict[str, int] = {}
    for sym in gene_symbols:
        data = _get(f"genes/{sym}")
        if data and "entrezGeneId" in data:
            sym_to_entrez[sym] = data["entrezGeneId"]
        else:
            print(f"  WARNING: could not resolve gene '{sym}' — skipping")
    return sym_to_entrez


def extract_hugo(df: pd.DataFrame) -> pd.DataFrame:
    if "gene" not in df.columns:
        return df
    sample = df["gene"].dropna()
    if sample.empty or not isinstance(sample.iloc[0], dict):
        return df
    df = df.copy()
    df["hugoGeneSymbol"] = df["gene"].apply(
        lambda g: g.get("hugoGeneSymbol", "") if isinstance(g, dict) else ""
    )
    return df


# ── STEP 1: download ────────────────────────────────────────────────────────

def download_mutations_extended(sym_to_entrez: dict[str, int], sample_list_id: str) -> pd.DataFrame:
    """GET /molecular-profiles/{id}/mutations — entrezGeneId required, so one call per gene."""
    profile_id = f"{STUDY_ID}_mutations"
    frames = []
    for sym, eid in sym_to_entrez.items():
        data = _get(
            f"molecular-profiles/{profile_id}/mutations",
            params={"sampleListId": sample_list_id, "entrezGeneId": eid, "projection": "DETAILED"},
        )
        n = len(data) if data else 0
        print(f"    {sym:<8} entrezGeneId={eid:<8} → {n} mutation records")
        if data:
            frames.append(pd.DataFrame(data))

    if not frames:
        return pd.DataFrame()

    df = pd.concat(frames, ignore_index=True)
    df = extract_hugo(df)
    if "chr" in df.columns and "chromosome" not in df.columns:
        df = df.rename(columns={"chr": "chromosome"})

    keep = ["patientId", "sampleId", "hugoGeneSymbol", "mutationType",
            "proteinChange", "variantType", "chromosome", "startPosition"]
    df = df[[c for c in keep if c in df.columns]].copy()
    df.insert(0, "studyId", STUDY_ID)
    df = df[df["hugoGeneSymbol"].isin(EXTENDED_GENES)].reset_index(drop=True)
    return df


def download_cna_extended(sample_list_id: str) -> pd.DataFrame:
    """
    GET /molecular-profiles/{id}/discrete-copy-number has no gene filter — it
    always returns the whole profile. data/cna_{STUDY_ID}.csv already holds
    that exact pull (same profile/sample-list/eventType=ALL), so reuse it
    instead of re-downloading. Falls back to a live call if absent.
    """
    cached = os.path.join(DATA_DIR, f"cna_{STUDY_ID}.csv")
    if os.path.isfile(cached):
        print(f"    Reusing already-downloaded full CNA pull → {cached}")
        df = pd.read_csv(cached)
        return df[df["hugoGeneSymbol"].isin(EXTENDED_GENES)].reset_index(drop=True)

    print("    No cached full CNA file found — calling GET discrete-copy-number "
          "(whole profile, will filter client-side)")
    profile_id = f"{STUDY_ID}_gistic"
    data = _get(
        f"molecular-profiles/{profile_id}/discrete-copy-number",
        params={"sampleListId": sample_list_id, "discreteCopyNumberEventType": "ALL",
                "projection": "DETAILED"},
    )
    if not data:
        return pd.DataFrame()

    df = pd.DataFrame(data)
    df = extract_hugo(df)
    if "alteration" in df.columns:
        df = df.rename(columns={"alteration": "value"})
    keep = ["patientId", "sampleId", "hugoGeneSymbol", "value"]
    df = df[[c for c in keep if c in df.columns]].copy()
    df.insert(0, "studyId", STUDY_ID)
    return df[df["hugoGeneSymbol"].isin(EXTENDED_GENES)].reset_index(drop=True)


# ── STEP 2: classification helpers (same semantics as 04_classify_patients.py) ─

def empty_flag(idx: pd.Index) -> pd.Series:
    return pd.Series(pd.NA, index=idx, dtype="Float64")


def set_wt(flags: pd.Series, universe: set) -> pd.Series:
    flags = flags.copy()
    in_u = flags.index.isin(universe)
    flags.loc[in_u] = flags.loc[in_u].fillna(0.0)
    return flags


def set_alt(flags: pd.Series, altered: set) -> pd.Series:
    flags = flags.copy()
    flags.loc[flags.index.isin(altered)] = 1.0
    return flags


def or_flags(*series_list: pd.Series) -> pd.Series:
    """
    Logical OR across any number of flag series, NA semantics:
      1   — at least one source is 1
      0   — at least one source is 0 and none is 1
      NaN — every source is NaN for that patient
    """
    idx = series_list[0].index
    result = empty_flag(idx)
    has_data = pd.Series(False, index=idx)
    is_1 = pd.Series(False, index=idx)
    for s in series_list:
        has_data = has_data | s.notna()
        is_1 = is_1 | (s == 1).fillna(False)
    result.loc[has_data] = 0.0
    result.loc[is_1] = 1.0
    return result


def muts_for(muts: pd.DataFrame, *symbols: str) -> pd.DataFrame:
    return muts[muts["hugoGeneSymbol"].isin(symbols)]


def any_mut_pts(muts: pd.DataFrame, *symbols: str) -> set:
    return set(muts_for(muts, *symbols)["patientId"])


def trunc_pts(muts: pd.DataFrame, *symbols: str) -> set:
    df = muts_for(muts, *symbols)
    return set(df.loc[df["mutationType"].isin(TRUNCATING), "patientId"])


def cna_pts(cna: pd.DataFrame, *symbols: str, op: str, threshold: int | float) -> set:
    df = cna[cna["hugoGeneSymbol"].isin(symbols)]
    if op == "le":
        return set(df.loc[df["value"] <= threshold, "patientId"])
    if op == "ge":
        return set(df.loc[df["value"] >= threshold, "patientId"])
    raise ValueError(f"op must be 'le' or 'ge', got {op!r}")


def build_flag(idx, universe, altered_pts) -> pd.Series:
    return set_alt(set_wt(empty_flag(idx), universe), altered_pts)


def _codon_nums(pc) -> list:
    if pd.isna(pc):
        return []
    return [int(n) for n in re.findall(r"\d+", str(pc))]


def is_dicer1_hotspot(pc) -> bool:
    return any(c in DICER1_HOTSPOT_CODONS for c in _codon_nums(pc))


def classify_extended(muts: pd.DataFrame, cna_df: pd.DataFrame, idx: pd.Index,
                       mut_u: set, cna_u: set) -> dict[str, pd.Series]:
    flags: dict[str, pd.Series] = {}

    # ── Tumor suppressors: truncating OR GISTIC ≤ -2 ────────────────────────
    for gene in TUMOR_SUPPRESSORS_TRUNC:
        m = build_flag(idx, mut_u, trunc_pts(muts, gene))
        c = build_flag(idx, cna_u, cna_pts(cna_df, gene, op="le", threshold=-2))
        flags[f"{gene}_ALTERED"] = or_flags(m, c)
        print(f"  {gene:<8} trunc={int((m==1).sum()):>3}  cna_del={int((c==1).sum()):>3}"
              f"  → combined {int((flags[f'{gene}_ALTERED']==1).sum()):>3}")

    # ── All mutation types OR GISTIC ≤ -2 ───────────────────────────────────
    for gene in ALL_MUT_TYPES_WITH_DEL:
        m = build_flag(idx, mut_u, any_mut_pts(muts, gene))
        c = build_flag(idx, cna_u, cna_pts(cna_df, gene, op="le", threshold=-2))
        flags[f"{gene}_ALTERED"] = or_flags(m, c)
        print(f"  {gene:<8} mut={int((m==1).sum()):>3}  cna_del={int((c==1).sum()):>3}"
              f"  → combined {int((flags[f'{gene}_ALTERED']==1).sum()):>3}")

    # ── Oncogenes: any mutation OR GISTIC ≥ +2 ──────────────────────────────
    for gene in ONCOGENES_ANY_MUT_AMP:
        m = build_flag(idx, mut_u, any_mut_pts(muts, gene))
        c = build_flag(idx, cna_u, cna_pts(cna_df, gene, op="ge", threshold=2))
        flags[f"{gene}_ALTERED"] = or_flags(m, c)
        print(f"  {gene:<8} mut={int((m==1).sum()):>3}  cna_amp={int((c==1).sum()):>3}"
              f"  → combined {int((flags[f'{gene}_ALTERED']==1).sum()):>3}")

    # ── DICER1: RNase IIIb hotspot mutations only; fall back to all truncating ─
    dicer1_all = muts_for(muts, "DICER1")
    dicer1_hot = set(dicer1_all.loc[dicer1_all["proteinChange"].apply(is_dicer1_hotspot), "patientId"])
    if not dicer1_hot:
        print("  DICER1: no RNase IIIb hotspot found; falling back to all truncating mutations")
        dicer1_hot = trunc_pts(muts, "DICER1")
    flags["DICER1_ALTERED"] = build_flag(idx, mut_u, dicer1_hot)
    print(f"  DICER1   hotspot/trunc altered={int((flags['DICER1_ALTERED']==1).sum()):>3}")

    # ── COL6A3: CNA only, both directions ───────────────────────────────────
    col6a3_alt = cna_pts(cna_df, "COL6A3", op="le", threshold=-2) | cna_pts(cna_df, "COL6A3", op="ge", threshold=2)
    flags["COL6A3_ALTERED"] = build_flag(idx, cna_u, col6a3_alt)
    print(f"  COL6A3   cna_del_or_amp altered={int((flags['COL6A3_ALTERED']==1).sum()):>3}")

    # ── ACTB: any mutation only, no CNA ──────────────────────────────────────
    flags["ACTB_ALTERED"] = build_flag(idx, mut_u, any_mut_pts(muts, "ACTB"))
    print(f"  ACTB     mut altered={int((flags['ACTB_ALTERED']==1).sum()):>3}")

    return flags


# ── STEP 3: pathway OR-flags ─────────────────────────────────────────────────

def build_pathways(master: pd.DataFrame) -> pd.DataFrame:
    idx = master.index
    pathway_cols: dict[str, pd.Series] = {}
    for pathway, cols in PATHWAYS.items():
        series_list = [master[c].astype("Float64") for c in cols if c in master.columns]
        pathway_cols[pathway] = or_flags(*series_list)
    pathway_df = pd.DataFrame(pathway_cols, index=idx)
    pathway_df["PATHWAY_ANY"] = or_flags(*[pathway_df[p] for p in PATHWAYS])
    return pathway_df


# ── STEP 4: summary ──────────────────────────────────────────────────────────

def print_pathway_summary(master: pd.DataFrame, n_cohort: int) -> None:
    width = 100
    print("\n" + "=" * width)
    print("PATHWAY ALTERATION SUMMARY")
    print("=" * width)
    for pathway, cols in list(PATHWAYS.items()) + [("PATHWAY_ANY", list(PATHWAYS.keys()))]:
        s = master[pathway].astype("Float64")
        n_alt = int((s == 1).sum())
        n_wt  = int((s == 0).sum())
        n_nan = int(s.isna().sum())
        pct   = round(n_alt / n_cohort * 100, 1) if n_cohort else float("nan")
        n_genes = len(cols)
        print(f"\n{pathway} | n_genes={n_genes} | n_altered={n_alt} | "
              f"n_wildtype={n_wt} | n_NaN={n_nan} | pct_cohort_altered={pct}%")
        for c in cols:
            if c not in master.columns:
                print(f"    {c:<22} (not available)")
                continue
            gs = master[c].astype("Float64")
            gene_name = c.replace("_ALTERED", "")
            print(f"    {gene_name:<22} n_altered={int((gs==1).sum())}")
    print("\n" + "=" * width)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(PROC_DIR, exist_ok=True)

    # ── Load existing master classification ─────────────────────────────────
    master_path = os.path.join(PROC_DIR, "master_classification.csv")
    print(f"Loading existing classification → {master_path}")
    master = pd.read_csv(master_path)
    idx = pd.Index(master["PATIENT_ID"].unique(), name="PATIENT_ID")
    print(f"  patients already classified : {len(idx)}  (study={master['STUDY_ID'].iloc[0]!r})")

    # ── STEP 1: download extended gene set ───────────────────────────────────
    print(f"\n{'='*70}\nSTEP 1 — download extended gene set ({len(EXTENDED_GENES)} genes)\n{'='*70}")
    sample_lists = get_sample_lists(STUDY_ID)
    mut_sl = pick_sample_list(sample_lists, ["all_cases_with_mutation_data", "all_cases_in_study"])
    cna_sl = pick_sample_list(sample_lists, [
        "all_cases_with_cna_data", "all_cases_with_mutation_and_cna_data", "all_cases_in_study",
    ])
    print(f"  mutation sample list : {mut_sl}")
    print(f"  CNA sample list      : {cna_sl}")

    print("\nResolving Entrez IDs...")
    sym_to_entrez = resolve_entrez_ids(EXTENDED_GENES)

    print("\n[MUTATIONS] — one GET call per gene (entrezGeneId is a required param)")
    muts = download_mutations_extended(sym_to_entrez, mut_sl)
    mut_out = os.path.join(DATA_DIR, "mutations_extended.csv")
    muts.to_csv(mut_out, index=False)
    print(f"  Saved {len(muts):,} rows → {mut_out}")

    print("\n[CNA]")
    cna_df = download_cna_extended(cna_sl)
    cna_out = os.path.join(DATA_DIR, "cna_extended.csv")
    cna_df.to_csv(cna_out, index=False)
    print(f"  Saved {len(cna_df):,} rows → {cna_out}")

    # ── STEP 2: classify extended genes ──────────────────────────────────────
    print(f"\n{'='*70}\nSTEP 2 — classify extended genes\n{'='*70}")
    mut_u = set(idx)                              # WES study — all patients sequenced
    cna_u = set(cna_df["patientId"].unique()) if not cna_df.empty else set()
    print(f"  mutation universe : {len(mut_u)}")
    print(f"  CNA universe      : {len(cna_u)}\n")

    ext_flags = classify_extended(muts, cna_df, idx, mut_u, cna_u)

    ext_df = pd.DataFrame(ext_flags, index=idx)
    for col in ext_df.columns:
        ext_df[col] = ext_df[col].astype("Float64")

    # ── STEP 3: merge + build pathways ───────────────────────────────────────
    print(f"\n{'='*70}\nSTEP 3 — merge and build pathway OR-flags\n{'='*70}")
    merged = master.set_index("PATIENT_ID").join(ext_df, how="left")

    pathway_df = build_pathways(merged)
    merged = merged.join(pathway_df)

    # Convert all *_ALTERED and PATHWAY_* columns to nullable Int8 for clean CSV output
    flag_cols = [c for c in merged.columns if c.endswith("_ALTERED") or c.startswith("PATHWAY_")]
    for col in flag_cols:
        merged[col] = merged[col].astype("Float64").astype("Int8")

    merged = merged.reset_index()

    # ── STEP 4: summary ───────────────────────────────────────────────────────
    print_pathway_summary(merged, n_cohort=len(idx))

    # ── Save ───────────────────────────────────────────────────────────────────
    out = os.path.join(PROC_DIR, "master_classification_extended.csv")
    merged.to_csv(out, index=False)
    n_pathway_cols = len(PATHWAYS) + 1  # + PATHWAY_ANY
    print(f"\nSaved master_classification_extended.csv — {len(merged)} patients, "
          f"{n_pathway_cols} pathway columns")


if __name__ == "__main__":
    main()
