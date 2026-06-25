"""
04_classify_patients.py

Classifies each Wilms' Tumor patient as ALTERED (1) or WILDTYPE (0) for every
pathway gene, following the rules from the prostate cancer manuscript methodology.

WILDTYPE RULE: a patient is WILDTYPE only if they are in the sequencing universe
for that gene AND no qualifying alteration is detected. Un-sequenced → NaN.

Sequencing universes (wt_target_2018_pub)
------------------------------------------
  Mutation   : all 652 patients (WES, allSampleCount == sequencedSampleCount)
  CNA (GISTIC): 124 patients with array CNA data
  mRNA z-score: 130 samples (RNA-Seq)

Key data quirks discovered in probe
-------------------------------------
  AMER1/WTX : completely absent from GISTIC CNA (X-chromosome excluded by GISTIC).
              CNA component will be all NaN → only mutation flag is used.
  DROSHA    : dominant-negative hotspot E1147K/D1151G are Missense → excluded by
              truncating-only rule. Only Q46* and R414* (Nonsense) qualify.
  DGCR8     : all three mutations are E518K (Missense) → 0 qualify under truncating-only.
              DGCR8_ALTERED will show 0 from mutations; CNA -2 also absent.

Outputs
-------
  wilms_survival/processed/master_classification.csv
"""

import os
import re
import sys
import pandas as pd

SCRIPT_DIR = os.path.dirname(__file__)
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "data")
PROC_DIR   = os.path.join(SCRIPT_DIR, "..", "processed")

STUDY_ID = "wt_target_2018_pub"

TRUNCATING = frozenset([
    "Nonsense_Mutation",
    "Frame_Shift_Ins",
    "Frame_Shift_Del",
    "Splice_Site",
    "Splice_Region",
    "Translation_Start_Site",
    "Nonstop_Mutation",
])

# CTNNB1 exon 3 phospho-site codons (GSK3β / CK1 sites)
CTNNB1_HOTSPOT = frozenset([32, 33, 34, 37, 41, 45])

PATHWAY = {
    "TP53_ALTERED":         "Cell Cycle",
    "MYCN_ALTERED":         "Transcription Reg",
    "SIX1_ALTERED":         "Transcription Reg",
    "SIX2_ALTERED":         "Transcription Reg",
    "SIX1_SIX2_ALTERED":    "Transcription Reg (combined)",
    "DROSHA_ALTERED":       "RNA Processing",
    "DGCR8_ALTERED":        "RNA Processing",
    "CTNNB1_ALTERED":       "WNT Signaling",
    "AMER1_ALTERED":        "WNT Signaling",
    "WT1_ALTERED":          "WNT Signaling",
    "MLLT1_ALTERED":        "Chromatin Remodeling",
    "NIPBL_ALTERED":        "PanCancer",
    "MDM4_ALTERED":         "PDMR Focus (p53 pathway)",
    "MDM2_ALTERED":         "PDMR Focus (p53 pathway)",
    "IGF2_CNA_ALTERED":     "PDMR Focus (imprinting/CNA)",
    "IGF2_EXPR_ALTERED":    "PDMR Focus (imprinting/expr)",
}


# ── Core flag helpers ──────────────────────────────────────────────────────────

def empty_flag(idx: pd.Index) -> pd.Series:
    """All-NaN series indexed by patientId (unclassified = unsequenced)."""
    return pd.Series(pd.NA, index=idx, dtype="Float64")


def set_wt(flags: pd.Series, universe: set) -> pd.Series:
    """Set WILDTYPE (0) for all patients in the sequencing universe who are still NaN."""
    flags = flags.copy()
    in_u = flags.index.isin(universe)
    flags.loc[in_u] = flags.loc[in_u].fillna(0.0)
    return flags


def set_alt(flags: pd.Series, altered: set) -> pd.Series:
    """Override to ALTERED (1) for specified patients."""
    flags = flags.copy()
    flags.loc[flags.index.isin(altered)] = 1.0
    return flags


def or_flags(s1: pd.Series, s2: pd.Series) -> pd.Series:
    """
    Logical OR with NA semantics:
      1   — at least one source is 1
      0   — at least one source is 0 and none is 1
      NaN — all sources are NaN (patient unsequenced for this gene entirely)
    """
    result = empty_flag(s1.index)
    has_data = s1.notna() | s2.notna()
    result.loc[has_data] = 0.0
    is_1 = (s1 == 1).fillna(False) | (s2 == 1).fillna(False)
    result.loc[is_1] = 1.0
    return result


# ── Query helpers ──────────────────────────────────────────────────────────────

def muts_for(muts: pd.DataFrame, *symbols: str) -> pd.DataFrame:
    """Filter mutation DataFrame to one or more gene symbols."""
    return muts[muts["hugoGeneSymbol"].isin(symbols)]


def any_mut_pts(muts: pd.DataFrame, *symbols: str) -> set:
    """patientIds with ANY mutation in the given genes."""
    return set(muts_for(muts, *symbols)["patientId"])


def trunc_pts(muts: pd.DataFrame, *symbols: str) -> set:
    """patientIds with a truncating mutation (TRUNCATING set) in the given genes."""
    df = muts_for(muts, *symbols)
    return set(df.loc[df["mutationType"].isin(TRUNCATING), "patientId"])


def cna_pts(cna: pd.DataFrame, *symbols: str, op: str, threshold: int | float) -> set:
    """
    patientIds where GISTIC value meets the threshold.
    op='le' for deletion (≤ threshold), op='ge' for amplification (≥ threshold).
    """
    df = cna[cna["hugoGeneSymbol"].isin(symbols)]
    if op == "le":
        return set(df.loc[df["value"] <= threshold, "patientId"])
    if op == "ge":
        return set(df.loc[df["value"] >= threshold, "patientId"])
    raise ValueError(f"op must be 'le' or 'ge', got {op!r}")


def build_flag(idx, universe, altered_pts) -> pd.Series:
    """Convenience: empty → set_wt → set_alt."""
    return set_alt(set_wt(empty_flag(idx), universe), altered_pts)


# ── CTNNB1 helpers ─────────────────────────────────────────────────────────────

def _codon_nums(pc) -> list:
    if pd.isna(pc):
        return []
    return [int(n) for n in re.findall(r'\d+', str(pc))]


def is_ctnnb1_hotspot(pc) -> bool:
    """True if any codon number in proteinChange falls in CTNNB1_HOTSPOT codons."""
    return any(c in CTNNB1_HOTSPOT for c in _codon_nums(pc))


# ── Summary ────────────────────────────────────────────────────────────────────

def print_summary(master: pd.DataFrame) -> None:
    rows = []
    for col, pathway in PATHWAY.items():
        if col not in master.columns:
            continue
        s = master[col].astype("Float64")
        n_alt = int((s == 1).sum())
        n_wt  = int((s == 0).sum())
        n_nan = int(s.isna().sum())
        n_seq = n_alt + n_wt
        freq  = round(n_alt / n_seq * 100, 1) if n_seq > 0 else float("nan")
        rows.append(dict(
            gene=col.replace("_ALTERED", ""),
            pathway=pathway,
            n_altered=n_alt,
            n_wildtype=n_wt,
            n_unsequenced=n_nan,
            alteration_pct=freq,
        ))
    df = (
        pd.DataFrame(rows)
        .sort_values("alteration_pct", ascending=False)
        .reset_index(drop=True)
    )
    width = 90
    print("\n" + "=" * width)
    print("ALTERATION FREQUENCY SUMMARY  (wt_target_2018_pub)")
    print("=" * width)
    print(df.to_string(index=False))
    print("=" * width)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    os.makedirs(PROC_DIR, exist_ok=True)

    # ── Load data ──────────────────────────────────────────────────────────────
    print(f"Loading {STUDY_ID} ...")
    clinical = pd.read_csv(os.path.join(DATA_DIR, f"clinical_{STUDY_ID}.csv"))
    muts     = pd.read_csv(os.path.join(DATA_DIR, f"mutations_{STUDY_ID}.csv"))
    cna_df   = pd.read_csv(os.path.join(DATA_DIR, f"cna_{STUDY_ID}.csv"))
    mrna_df  = pd.read_csv(os.path.join(DATA_DIR, f"mrna_zscores_{STUDY_ID}.csv"))

    # mRNA: add patientId (TARGET-50-CXXXXX-01 → TARGET-50-CXXXXX)
    mrna_df = mrna_df.copy()
    mrna_df["patientId"] = mrna_df["sampleId"].str.rsplit("-", n=1).str[0]

    # ── Sequencing universes ───────────────────────────────────────────────────
    idx   = pd.Index(clinical["patientId"].unique(), name="patientId")
    mut_u = set(idx)                              # all patients — WES study
    cna_u = set(cna_df["patientId"].unique())     # subset with GISTIC array CNA
    rna_u = set(mrna_df["patientId"].unique())    # subset with RNA-Seq

    print(f"  patients in study : {len(idx)}")
    print(f"  mutation universe : {len(mut_u)}")
    print(f"  CNA universe      : {len(cna_u)}")
    print(f"  mRNA universe     : {len(rna_u)}")

    # ── Classify genes ─────────────────────────────────────────────────────────
    print("\nClassifying ...")
    flags: dict[str, pd.Series] = {}

    # ── TP53 ── all mutations OR CNA ≤ -2 ──────────────────────────────────────
    tp53_m = build_flag(idx, mut_u, any_mut_pts(muts, "TP53"))
    tp53_c = build_flag(idx, cna_u, cna_pts(cna_df, "TP53", op="le", threshold=-2))
    flags["TP53_ALTERED"] = or_flags(tp53_m, tp53_c)
    print(f"  TP53        mut={int((tp53_m==1).sum()):>3}  cna_del={int((tp53_c==1).sum()):>3}"
          f"  → combined {int((flags['TP53_ALTERED']==1).sum()):>3}")

    # ── MYCN ── any mutation OR CNA ≥ +2 ───────────────────────────────────────
    mycn_m = build_flag(idx, mut_u, any_mut_pts(muts, "MYCN"))
    mycn_c = build_flag(idx, cna_u, cna_pts(cna_df, "MYCN", op="ge", threshold=2))
    flags["MYCN_ALTERED"] = or_flags(mycn_m, mycn_c)
    print(f"  MYCN        mut={int((mycn_m==1).sum()):>3}  cna_amp={int((mycn_c==1).sum()):>3}"
          f"  → combined {int((flags['MYCN_ALTERED']==1).sum()):>3}")

    # ── SIX1 ── Q177R hotspot only; no CNA ─────────────────────────────────────
    six1_hot = set(
        muts_for(muts, "SIX1")
        .loc[muts_for(muts, "SIX1")["proteinChange"].str.contains("Q177R", na=False),
             "patientId"]
    )
    flags["SIX1_ALTERED"] = build_flag(idx, mut_u, six1_hot)
    print(f"  SIX1 Q177R  altered={int((flags['SIX1_ALTERED']==1).sum()):>3}")

    # ── SIX2 ── Q177R hotspot only; no CNA ─────────────────────────────────────
    six2_hot = set(
        muts_for(muts, "SIX2")
        .loc[muts_for(muts, "SIX2")["proteinChange"].str.contains("Q177R", na=False),
             "patientId"]
    )
    flags["SIX2_ALTERED"] = build_flag(idx, mut_u, six2_hot)
    print(f"  SIX2 Q177R  altered={int((flags['SIX2_ALTERED']==1).sum()):>3}")

    # ── SIX1_SIX2 combined ──────────────────────────────────────────────────────
    flags["SIX1_SIX2_ALTERED"] = or_flags(flags["SIX1_ALTERED"], flags["SIX2_ALTERED"])
    print(f"  SIX1|SIX2   altered={int((flags['SIX1_SIX2_ALTERED']==1).sum()):>3}"
          f"  (SIX1 Q177R OR SIX2 Q177R)")

    # ── DROSHA ── truncating only OR CNA ≤ -2 ──────────────────────────────────
    # NOTE: The dominant-negative hotspot missense mutations E1147K (×3) and
    # D1151G/D1151A are excluded by the truncating-only rule.
    # Only Q46* and R414* (Nonsense_Mutation) qualify.
    drosha_m = build_flag(idx, mut_u, trunc_pts(muts, "DROSHA"))
    drosha_c = build_flag(idx, cna_u, cna_pts(cna_df, "DROSHA", op="le", threshold=-2))
    flags["DROSHA_ALTERED"] = or_flags(drosha_m, drosha_c)
    print(f"  DROSHA      trunc={int((drosha_m==1).sum()):>3}  cna_del={int((drosha_c==1).sum()):>3}"
          f"  → combined {int((flags['DROSHA_ALTERED']==1).sum()):>3}"
          f"  (NOTE: E1147K missense excluded by truncating-only rule)")

    # ── DGCR8 ── truncating only OR CNA ≤ -2 ───────────────────────────────────
    # NOTE: All three DGCR8 mutations are E518K (Missense_Mutation) — none qualify
    # under the truncating-only rule. GISTIC -2 also absent for DGCR8.
    # DGCR8_ALTERED will be 0 for all sequenced patients.
    dgcr8_m = build_flag(idx, mut_u, trunc_pts(muts, "DGCR8"))
    dgcr8_c = build_flag(idx, cna_u, cna_pts(cna_df, "DGCR8", op="le", threshold=-2))
    flags["DGCR8_ALTERED"] = or_flags(dgcr8_m, dgcr8_c)
    print(f"  DGCR8       trunc={int((dgcr8_m==1).sum()):>3}  cna_del={int((dgcr8_c==1).sum()):>3}"
          f"  → combined {int((flags['DGCR8_ALTERED']==1).sum()):>3}"
          f"  (NOTE: E518K missense hotspot excluded by truncating-only rule)")

    # ── CTNNB1 ── exon 3 phospho-site hotspot codons 32/33/34/37/41/45 ─────────
    ctnnb1_all = muts_for(muts, "CTNNB1")
    ctnnb1_hot = ctnnb1_all[ctnnb1_all["proteinChange"].apply(is_ctnnb1_hotspot)]
    if ctnnb1_hot.empty and not ctnnb1_all.empty:
        print("  CTNNB1: no exon-3 hotspot found; falling back to all CTNNB1 mutations")
        ctnnb1_hot = ctnnb1_all
    hotspot_detail = ctnnb1_hot["proteinChange"].value_counts().to_dict()
    flags["CTNNB1_ALTERED"] = build_flag(idx, mut_u, set(ctnnb1_hot["patientId"]))
    print(f"  CTNNB1      hotspot_altered={int((flags['CTNNB1_ALTERED']==1).sum()):>3}"
          f"  detail={hotspot_detail}")

    # ── AMER1 / WTX ── truncating OR CNA ≤ -2 ──────────────────────────────────
    # NOTE: AMER1 (X-linked) is absent from the GISTIC CNA array output entirely.
    # The CNA component will be all-NaN for CNA-universe patients but since all
    # patients are also in the mutation universe, the combined flag collapses to
    # the mutation flag only (no patient will be NaN).
    amer1_m = build_flag(idx, mut_u, trunc_pts(muts, "AMER1", "WTX"))
    amer1_c = build_flag(idx, cna_u, cna_pts(cna_df, "AMER1", "WTX", op="le", threshold=-2))
    flags["AMER1_ALTERED"] = or_flags(amer1_m, amer1_c)
    print(f"  AMER1/WTX   trunc={int((amer1_m==1).sum()):>3}  cna_del={int((amer1_c==1).sum()):>3}"
          f"  → combined {int((flags['AMER1_ALTERED']==1).sum()):>3}"
          f"  (CNA: AMER1 absent from GISTIC — X-chromosome excluded)")

    # ── WT1 ── truncating OR CNA ≤ -2 ──────────────────────────────────────────
    wt1_m = build_flag(idx, mut_u, trunc_pts(muts, "WT1"))
    wt1_c = build_flag(idx, cna_u, cna_pts(cna_df, "WT1", op="le", threshold=-2))
    flags["WT1_ALTERED"] = or_flags(wt1_m, wt1_c)
    print(f"  WT1         trunc={int((wt1_m==1).sum()):>3}  cna_del={int((wt1_c==1).sum()):>3}"
          f"  → combined {int((flags['WT1_ALTERED']==1).sum()):>3}")

    # ── MLLT1 ── ALL mutations (incl. in-frame) OR CNA ≤ -2 ────────────────────
    mllt1_m = build_flag(idx, mut_u, any_mut_pts(muts, "MLLT1"))
    mllt1_c = build_flag(idx, cna_u, cna_pts(cna_df, "MLLT1", op="le", threshold=-2))
    flags["MLLT1_ALTERED"] = or_flags(mllt1_m, mllt1_c)
    print(f"  MLLT1       mut={int((mllt1_m==1).sum()):>3}  cna_del={int((mllt1_c==1).sum()):>3}"
          f"  → combined {int((flags['MLLT1_ALTERED']==1).sum()):>3}")

    # ── NIPBL ── truncating OR CNA ≤ -2 ────────────────────────────────────────
    nipbl_m = build_flag(idx, mut_u, trunc_pts(muts, "NIPBL"))
    nipbl_c = build_flag(idx, cna_u, cna_pts(cna_df, "NIPBL", op="le", threshold=-2))
    flags["NIPBL_ALTERED"] = or_flags(nipbl_m, nipbl_c)
    print(f"  NIPBL       trunc={int((nipbl_m==1).sum()):>3}  cna_del={int((nipbl_c==1).sum()):>3}"
          f"  → combined {int((flags['NIPBL_ALTERED']==1).sum()):>3}")

    # ── MDM4 ── any mutation OR CNA ≥ +2 ───────────────────────────────────────
    mdm4_m = build_flag(idx, mut_u, any_mut_pts(muts, "MDM4"))
    mdm4_c = build_flag(idx, cna_u, cna_pts(cna_df, "MDM4", op="ge", threshold=2))
    flags["MDM4_ALTERED"] = or_flags(mdm4_m, mdm4_c)
    print(f"  MDM4        mut={int((mdm4_m==1).sum()):>3}  cna_amp={int((mdm4_c==1).sum()):>3}"
          f"  → combined {int((flags['MDM4_ALTERED']==1).sum()):>3}")

    # ── MDM2 ── any mutation OR CNA ≥ +2 ───────────────────────────────────────
    mdm2_m = build_flag(idx, mut_u, any_mut_pts(muts, "MDM2"))
    mdm2_c = build_flag(idx, cna_u, cna_pts(cna_df, "MDM2", op="ge", threshold=2))
    flags["MDM2_ALTERED"] = or_flags(mdm2_m, mdm2_c)
    print(f"  MDM2        mut={int((mdm2_m==1).sum()):>3}  cna_amp={int((mdm2_c==1).sum()):>3}"
          f"  → combined {int((flags['MDM2_ALTERED']==1).sum()):>3}")

    # ── IGF2 ── CNA ≥ +1 (LOI gain) AND expression z > 1.5 — separate columns ─
    # CNA flag: threshold ≥ +1 because LOI produces moderate gain, not focal amplification
    igf2_cna_f = build_flag(idx, cna_u, cna_pts(cna_df, "IGF2", op="ge", threshold=1))
    flags["IGF2_CNA_ALTERED"] = igf2_cna_f

    # Expression flag: aggregate samples per patient by max z-score, threshold > 1.5
    igf2_expr = mrna_df[mrna_df["hugoGeneSymbol"] == "IGF2"].copy()
    igf2_per_pt = igf2_expr.groupby("patientId")["value"].max()
    igf2_expr_alt = set(igf2_per_pt[igf2_per_pt > 1.5].index)
    igf2_expr_f = build_flag(idx, rna_u, igf2_expr_alt)
    flags["IGF2_EXPR_ALTERED"] = igf2_expr_f

    print(f"  IGF2_CNA    cna_gain(≥+1)={int((igf2_cna_f==1).sum()):>3}"
          f"  NaN={int(igf2_cna_f.isna().sum()):>3}")
    print(f"  IGF2_EXPR   z>1.5={int((igf2_expr_f==1).sum()):>3}"
          f"  NaN={int(igf2_expr_f.isna().sum()):>3}")

    # ── Build master DataFrame ─────────────────────────────────────────────────
    clin = clinical.set_index("patientId")

    master = pd.DataFrame(index=idx)
    master.index.name = "PATIENT_ID"

    # Clinical columns
    for src, dst in [
        ("STUDY_ID",    "STUDY_ID"),
        ("COHORT_CODE", "COHORT_CODE"),
        ("OS_MONTHS",   "OS_MONTHS"),
        ("OS_STATUS",   "OS_STATUS"),
        ("EFS_MONTHS",  "EFS_MONTHS"),
        ("EFS_STATUS",  "EFS_STATUS"),
        ("AGE_YEARS",   "AGE"),
    ]:
        master[dst] = clin[src] if src in clin.columns else pd.NA

    # Gene flags — in the column order specified
    ORDERED_FLAGS = [
        "TP53_ALTERED",
        "MYCN_ALTERED",
        "SIX1_ALTERED", "SIX2_ALTERED", "SIX1_SIX2_ALTERED",
        "DROSHA_ALTERED", "DGCR8_ALTERED",
        "CTNNB1_ALTERED", "AMER1_ALTERED", "WT1_ALTERED",
        "MLLT1_ALTERED",
        "NIPBL_ALTERED",
        "MDM4_ALTERED", "MDM2_ALTERED",
        "IGF2_CNA_ALTERED", "IGF2_EXPR_ALTERED",
    ]
    for col in ORDERED_FLAGS:
        master[col] = flags[col]

    master = master.reset_index()

    # Convert flag columns to nullable Int8 (clean 0 / 1 / <NA> in CSV)
    for col in ORDERED_FLAGS:
        master[col] = master[col].astype("Int8")

    # ── Save ───────────────────────────────────────────────────────────────────
    out = os.path.join(PROC_DIR, "master_classification.csv")
    master.to_csv(out, index=False)
    print(f"\nSaved → {os.path.abspath(out)}")
    print(f"Shape  : {master.shape[0]} rows × {master.shape[1]} columns")
    print(f"Preview:")
    print(master[["PATIENT_ID", "OS_MONTHS", "OS_STATUS", "TP53_ALTERED",
                  "MYCN_ALTERED", "SIX1_SIX2_ALTERED", "IGF2_CNA_ALTERED",
                  "IGF2_EXPR_ALTERED"]].head(8).to_string(index=False))

    # ── Summary table ──────────────────────────────────────────────────────────
    print_summary(master)


if __name__ == "__main__":
    main()
