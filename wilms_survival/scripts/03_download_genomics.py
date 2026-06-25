"""
03_download_genomics.py

Download mutation, copy-number, and mRNA z-score data from cBioPortal for all
Wilms tumor studies that have survival data (listed in available_studies.csv).

All three data types use POST endpoints — the GET /mutations endpoint requires
a mandatory entrezGeneId parameter and cannot return all-gene data.

Outputs (wilms_survival/data/)
------------------------------
  mutations_{study_id}.csv   — one row per somatic mutation
  cna_{study_id}.csv         — one row per gene × sample GISTIC call
  mrna_zscores_{study_id}.csv — one row per gene × sample z-score (target genes only)
"""

import os
import sys
import requests
import pandas as pd

BASE_URL   = "https://www.cbioportal.org/api"
SCRIPT_DIR = os.path.dirname(__file__)
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "data")
CSV_INDEX  = os.path.join(DATA_DIR, "available_studies.csv")

GENES_OF_INTEREST = [
    "TP53", "MYCN", "SIX1", "SIX2", "DROSHA", "DGCR8",
    "CTNNB1", "AMER1", "WT1", "MLLT1", "NIPBL",
    "MDM4", "MDM2", "IGF2",
]

# ── HTTP helpers ───────────────────────────────────────────────────────────────

def _get(path: str, params: dict | None = None) -> list | dict | None:
    url = f"{BASE_URL}/{path.lstrip('/')}"
    try:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        return r.json()
    except requests.exceptions.HTTPError as e:
        if e.response.status_code != 404:
            print(f"  HTTP {e.response.status_code} — {url}", file=sys.stderr)
        return None
    except requests.exceptions.RequestException as e:
        print(f"  Request error — {e}", file=sys.stderr)
        return None


def _post(path: str, payload: dict, params: dict | None = None) -> list | dict | None:
    url = f"{BASE_URL}/{path.lstrip('/')}"
    try:
        r = requests.post(url, json=payload, params=params, timeout=300)
        r.raise_for_status()
        return r.json()
    except requests.exceptions.HTTPError as e:
        print(f"  HTTP {e.response.status_code} — {url}", file=sys.stderr)
        return None
    except requests.exceptions.RequestException as e:
        print(f"  Request error — {e}", file=sys.stderr)
        return None


# ── Study metadata helpers ─────────────────────────────────────────────────────

def get_profiles(study_id: str) -> list[dict]:
    data = _get(f"studies/{study_id}/molecular-profiles")
    return data if isinstance(data, list) else []


def get_sample_lists(study_id: str) -> list[dict]:
    data = _get(f"studies/{study_id}/sample-lists")
    return data if isinstance(data, list) else []


def pick_sample_list(sample_lists: list[dict], priority_categories: list[str]) -> str | None:
    """Return the sampleListId of the first list matching a category in priority order."""
    by_cat = {sl["category"]: sl["sampleListId"] for sl in sample_lists}
    for cat in priority_categories:
        if cat in by_cat:
            return by_cat[cat]
    # fallback to the all-samples list
    return by_cat.get("all_cases_in_study") or (sample_lists[0]["sampleListId"] if sample_lists else None)


# ── Gene resolution ────────────────────────────────────────────────────────────

def resolve_entrez_ids(gene_symbols: list[str]) -> tuple[dict[str, int], dict[int, str]]:
    """
    Return (symbol→entrezId, entrezId→symbol) for each symbol.
    Prints a warning for any symbol that cannot be resolved.
    """
    sym_to_entrez: dict[str, int] = {}
    for sym in gene_symbols:
        data = _get(f"genes/{sym}")
        if data and "entrezGeneId" in data:
            sym_to_entrez[sym] = data["entrezGeneId"]
        else:
            print(f"  WARNING: could not resolve gene '{sym}' — skipping")
    entrez_to_sym = {v: k for k, v in sym_to_entrez.items()}
    return sym_to_entrez, entrez_to_sym


# ── Column helpers ─────────────────────────────────────────────────────────────

def extract_hugo(df: pd.DataFrame) -> pd.DataFrame:
    """
    If 'gene' column contains nested dicts (DETAILED projection),
    extract gene.hugoGeneSymbol into a flat column.
    """
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


def log_summary(study_id: str, dtype: str, df: pd.DataFrame) -> None:
    n_rec = len(df)
    if "patientId" in df.columns:
        n_pts = df["patientId"].nunique()
        label = "patients"
    elif "sampleId" in df.columns:
        n_pts = df["sampleId"].nunique()
        label = "samples"
    else:
        n_pts, label = "?", ""
    print(f"    → {n_rec:>9,} records   {n_pts:>4} unique {label}")


# ── Profile selectors ──────────────────────────────────────────────────────────

def find_mutation_profile(profiles: list[dict]) -> dict | None:
    return next(
        (p for p in profiles if p.get("molecularAlterationType") == "MUTATION_EXTENDED"),
        None,
    )


def find_cna_profile(profiles: list[dict]) -> dict | None:
    return next(
        (p for p in profiles
         if p.get("molecularAlterationType") == "COPY_NUMBER_ALTERATION"
         and p.get("datatype") == "DISCRETE"),
        None,
    )


def find_zscore_profile(profiles: list[dict]) -> dict | None:
    """
    Pick the best mRNA z-score profile.
    Priority: RNA-Seq TPM > RNA-Seq FPKM > any other RNA-Seq > microarray.
    """
    zscores = [
        p for p in profiles
        if p.get("molecularAlterationType") == "MRNA_EXPRESSION"
        and p.get("datatype") == "Z-SCORE"
    ]
    if not zscores:
        return None
    for keyword in ("tpm", "rpkm", "rna_seq", "rnaseq", "rna"):
        for p in zscores:
            if keyword in p.get("molecularProfileId", "").lower():
                return p
    return zscores[0]   # microarray fallback


# ── Download functions ─────────────────────────────────────────────────────────

def download_mutations(
    study_id: str,
    profiles: list[dict],
    sample_lists: list[dict],
) -> pd.DataFrame:
    """
    Endpoint: POST /molecular-profiles/{id}/mutations/fetch
    Projection DETAILED includes the nested 'gene' object with hugoGeneSymbol.
    """
    profile = find_mutation_profile(profiles)
    if profile is None:
        print("    No MUTATION_EXTENDED profile — skipping")
        return pd.DataFrame()

    pid    = profile["molecularProfileId"]
    sl_id  = pick_sample_list(sample_lists, [
        "all_cases_with_mutation_data", "all_cases_in_study",
    ])
    print(f"    profile    : {pid}")
    print(f"    sample list: {sl_id}")

    data = _post(
        f"molecular-profiles/{pid}/mutations/fetch",
        payload={"sampleListId": sl_id},
        params={"projection": "DETAILED", "pageSize": 100_000},
    )
    if not data:
        print("    Fetch returned empty result")
        return pd.DataFrame()

    df = pd.DataFrame(data)
    df = extract_hugo(df)

    keep = ["patientId", "sampleId", "hugoGeneSymbol", "mutationType",
            "proteinChange", "variantType", "chromosome", "startPosition"]
    df = df[[c for c in keep if c in df.columns]].copy()
    df.insert(0, "studyId", study_id)
    return df


def download_cna(
    study_id: str,
    profiles: list[dict],
    sample_lists: list[dict],
) -> pd.DataFrame:
    """
    Endpoint: POST /molecular-profiles/{id}/discrete-copy-number/fetch
    DETAILED projection includes patientId and nested 'gene' object.
    Field 'alteration' holds the GISTIC score (-2 to 2); renamed to 'value'.

    Note: discreteCopyNumberEventType=ALL downloads all genes including
    neutral (0). This produces ~2.6M rows for 129 samples — the save step
    may take a moment.
    """
    profile = find_cna_profile(profiles)
    if profile is None:
        print("    No DISCRETE CNA (GISTIC) profile — skipping")
        return pd.DataFrame()

    pid   = profile["molecularProfileId"]
    sl_id = pick_sample_list(sample_lists, [
        "all_cases_with_cna_data",
        "all_cases_with_mutation_and_cna_data",
        "all_cases_in_study",
    ])
    print(f"    profile    : {pid}")
    print(f"    sample list: {sl_id}")

    data = _post(
        f"molecular-profiles/{pid}/discrete-copy-number/fetch",
        payload={"sampleListId": sl_id},
        params={"projection": "DETAILED", "discreteCopyNumberEventType": "ALL"},
    )
    if not data:
        print("    Fetch returned empty result")
        return pd.DataFrame()

    df = pd.DataFrame(data)
    df = extract_hugo(df)

    # Rename GISTIC field to 'value' per spec
    if "alteration" in df.columns:
        df = df.rename(columns={"alteration": "value"})

    keep = ["patientId", "sampleId", "hugoGeneSymbol", "value"]
    df = df[[c for c in keep if c in df.columns]].copy()
    df.insert(0, "studyId", study_id)
    return df


def download_expression(
    study_id: str,
    profiles: list[dict],
    sample_lists: list[dict],
    entrez_to_sym: dict[int, str],
) -> pd.DataFrame:
    """
    Endpoint: POST /molecular-profiles/{id}/molecular-data/fetch
    Downloads z-scores for GENES_OF_INTEREST only.
    Response does not include patientId — only sampleId is kept.
    """
    profile = find_zscore_profile(profiles)
    if profile is None:
        print("    No mRNA z-score profile — skipping")
        return pd.DataFrame()

    pid   = profile["molecularProfileId"]
    sl_id = pick_sample_list(sample_lists, [
        "all_cases_with_mrna_rnaseq_data",
        "all_cases_with_mrna_array_data",
        "all_cases_in_study",
    ])
    print(f"    profile    : {pid}")
    print(f"    sample list: {sl_id}")
    print(f"    genes      : {list(entrez_to_sym.values())}")

    data = _post(
        f"molecular-profiles/{pid}/molecular-data/fetch",
        payload={
            "sampleListId": sl_id,
            "entrezGeneIds": list(entrez_to_sym.keys()),
        },
        params={"projection": "SUMMARY"},
    )
    if not data:
        print("    Fetch returned empty result")
        return pd.DataFrame()

    df = pd.DataFrame(data)
    df["hugoGeneSymbol"] = df["entrezGeneId"].map(entrez_to_sym)

    keep = ["sampleId", "hugoGeneSymbol", "value"]
    df = df[[c for c in keep if c in df.columns]].copy()
    df.insert(0, "studyId", study_id)
    return df


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)

    # Load study list
    df_idx = pd.read_csv(CSV_INDEX)
    study_ids = df_idx.loc[df_idx["has_survival_data"] == True, "study_id"].tolist()
    print(f"Studies to process : {study_ids}")

    # Resolve gene Entrez IDs once (shared across studies)
    print(f"\nResolving Entrez IDs for {len(GENES_OF_INTEREST)} genes...")
    sym_to_entrez, entrez_to_sym = resolve_entrez_ids(GENES_OF_INTEREST)
    print("  " + "  ".join(f"{s}={e}" for s, e in sym_to_entrez.items()))

    # ── Per-study loop ─────────────────────────────────────────────────────────
    for study_id in study_ids:
        print(f"\n{'=' * 70}")
        print(f"  {study_id}")
        print(f"{'=' * 70}")

        profiles     = get_profiles(study_id)
        sample_lists = get_sample_lists(study_id)

        if not profiles:
            print("  Could not retrieve molecular profiles — skipping study")
            continue

        print(f"  Profiles available ({len(profiles)}):")
        for p in profiles:
            print(f"    {p['molecularAlterationType']:<28} {p['datatype']:<12} {p['molecularProfileId']}")

        # ── Mutations ──────────────────────────────────────────────────────────
        print(f"\n  [MUTATIONS]")
        mut_df = download_mutations(study_id, profiles, sample_lists)
        if not mut_df.empty:
            out = os.path.join(DATA_DIR, f"mutations_{study_id}.csv")
            mut_df.to_csv(out, index=False)
            log_summary(study_id, "mutations", mut_df)
            print(f"    Saved → {os.path.abspath(out)}")
            # Gene breakdown
            top = (
                mut_df.groupby("hugoGeneSymbol")["sampleId"]
                .nunique()
                .sort_values(ascending=False)
                .head(10)
            )
            print(f"    Top mutated genes: " +
                  "  ".join(f"{g}({n})" for g, n in top.items()))

        # ── CNA ────────────────────────────────────────────────────────────────
        print(f"\n  [COPY NUMBER ALTERATIONS]")
        cna_df = download_cna(study_id, profiles, sample_lists)
        if not cna_df.empty:
            out = os.path.join(DATA_DIR, f"cna_{study_id}.csv")
            print(f"    Saving {len(cna_df):,} rows (may take a moment)...")
            cna_df.to_csv(out, index=False)
            log_summary(study_id, "CNA", cna_df)
            val_dist = cna_df["value"].value_counts().sort_index().to_dict()
            print(f"    GISTIC distribution: {val_dist}")
            print(f"    Saved → {os.path.abspath(out)}")

        # ── mRNA z-scores ──────────────────────────────────────────────────────
        print(f"\n  [mRNA Z-SCORES]")
        expr_df = download_expression(study_id, profiles, sample_lists, entrez_to_sym)
        if not expr_df.empty:
            out = os.path.join(DATA_DIR, f"mrna_zscores_{study_id}.csv")
            expr_df.to_csv(out, index=False)
            log_summary(study_id, "mRNA z-scores", expr_df)
            # Pivot preview: median z-score per gene
            pivot = (
                expr_df.groupby("hugoGeneSymbol")["value"]
                .agg(["median", "count"])
                .rename(columns={"median": "median_z", "count": "n_samples"})
                .sort_values("median_z", ascending=False)
            )
            print(f"    Median z-score per gene:\n{pivot.to_string()}")
            print(f"    Saved → {os.path.abspath(out)}")

    # ── Final recap ────────────────────────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print("DOWNLOAD SUMMARY")
    print(f"{'=' * 70}")
    print(f"  {'Study':<30}  {'File':<25}  {'Rows':>10}")
    print(f"  {'-' * 68}")
    for study_id in study_ids:
        for dtype in ("mutations", "cna", "mrna_zscores"):
            path = os.path.join(DATA_DIR, f"{dtype}_{study_id}.csv")
            if os.path.isfile(path):
                n = sum(1 for _ in open(path)) - 1   # row count excl. header
                print(f"  {study_id:<30}  {dtype:<25}  {n:>10,}")
            else:
                print(f"  {study_id:<30}  {dtype:<25}  {'(not downloaded)':>10}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
