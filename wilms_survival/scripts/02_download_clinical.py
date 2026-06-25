"""
02_download_clinical.py

Download and standardize patient-level clinical / survival data from cBioPortal
for all Wilms tumor studies that have survival data.

Study IDs attempted
-------------------
  target_wt_2013     (user-specified primary; will warn and skip if not found)
  + all studies in wilms_survival/data/available_studies.csv with has_survival_data=True

Outputs
-------
  wilms_survival/data/clinical_{study_id}.csv   one per study
"""

import os
import sys
import re
import requests
import pandas as pd
import numpy as np

BASE_URL   = "https://www.cbioportal.org/api"
SCRIPT_DIR = os.path.dirname(__file__)
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "data")
CSV_INDEX  = os.path.join(DATA_DIR, "available_studies.csv")

# User-requested primary study (may not exist on cBioPortal)
PRIMARY_STUDY_ID = "target_wt_2013"

# ── API helpers ───────────────────────────────────────────────────────────────

def get(path: str, params: dict | None = None) -> list | dict | None:
    url = f"{BASE_URL}/{path.lstrip('/')}"
    try:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        return r.json()
    except requests.exceptions.HTTPError as e:
        code = e.response.status_code
        if code == 404:
            print(f"  WARNING: study not found (404) — {url}")
        else:
            print(f"  HTTP error {code} — {url}", file=sys.stderr)
        return None
    except requests.exceptions.RequestException as e:
        print(f"  Request failed — {e}", file=sys.stderr)
        return None


def fetch_patient_clinical(study_id: str) -> pd.DataFrame:
    """
    Download all PATIENT-level clinical attributes and pivot to wide format.
    Returns an empty DataFrame on failure.
    """
    raw = get(
        f"studies/{study_id}/clinical-data",
        params={"clinicalDataType": "PATIENT", "pageSize": 100_000},
    )
    if not raw:
        return pd.DataFrame()
    long = pd.DataFrame(raw)
    wide = (
        long.pivot_table(
            index="patientId",
            columns="clinicalAttributeId",
            values="value",
            aggfunc="first",
        )
        .reset_index()
    )
    wide.columns.name = None
    return wide


# ── column detection ──────────────────────────────────────────────────────────

# Priority-ordered candidates for each concept.
# First hit wins; conversion function optionally applied.
OS_TIME_CANDIDATES  = ["OS_MONTHS", "OVERALL_SURVIVAL", "OS_DAYS"]
OS_STATUS_CANDIDATES = ["OS_STATUS", "VITAL_STATUS", "PATIENT_STATUS"]

EFS_TIME_CANDIDATES   = ["EFS_MONTHS", "EFS_EVENT_FREE_SURVIVAL_MONTHS",
                          "EVENT_FREE_SURVIVAL", "DFS_MONTHS", "DAYS_TO_EVENT"]
EFS_STATUS_CANDIDATES = ["EFS_STATUS", "DFS_STATUS", "RFS_STATUS", "EVENT_TYPE"]

AGE_CANDIDATES    = ["AGE", "AGE_AT_DIAGNOSIS", "AGE_IN_DAYS"]
STAGE_CANDIDATES  = ["CLINICAL_STAGE", "STAGE", "AJCC_PATHOLOGIC_TUMOR_STAGE",
                     "TUMOR_STAGE", "PATH_T_STAGE"]
HIST_CANDIDATES   = ["HISTOLOGY", "CANCER_TYPE_DETAILED", "PRIMARY_DIAGNOSIS",
                     "HISTOLOGY_CLASSIFICATION_IN_PRIMARY_TUMOR"]


def _first_present(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for col in candidates:
        if col in df.columns:
            return col
    return None


# ── status recoding ───────────────────────────────────────────────────────────

def _recode_os_status(series: pd.Series) -> pd.Series:
    """
    Map any cBioPortal OS_STATUS / VITAL_STATUS variant to 1 (event) / 0 (censored).

    Observed values across TARGET Wilms studies:
      '1:DECEASED', '0:LIVING'   (standard cBioPortal format)
      'DECEASED', 'LIVING'
      'Dead', 'Alive'
      '1', '0'
    """
    s = series.astype(str).str.strip().str.upper()
    event_patterns = re.compile(r"^(1|DECEASED|DEAD|1:DECEASED)$")
    censor_patterns = re.compile(r"^(0|LIVING|ALIVE|0:LIVING)$")
    out = pd.Series(np.nan, index=series.index, dtype="Float64")
    out[s.str.match(event_patterns)]  = 1.0
    out[s.str.match(censor_patterns)] = 0.0
    return out


def _recode_efs_status(series: pd.Series) -> pd.Series:
    """
    Map EFS_STATUS / EVENT_TYPE variants to 1 (event) / 0 (censored).

    Observed in wt_target_2018_pub EVENT_TYPE:
      'Relapse', 'Progression' → event = 1
      'None'                   → censored = 0

    Standard cBioPortal EFS_STATUS format (if present):
      '1:Recurred/Progressed', '1:RELAPSED', '1:EVENT'  → 1
      '0:DiseaseFree', '0:CENSORED'                      → 0
    """
    s = series.astype(str).str.strip().str.upper()
    event_pats  = re.compile(r"^(1|RELAPSE|RELAPSED|PROGRESSION|PROGRESSED|"
                              r"1:RECURRED.*|1:RELAPSED|1:EVENT|1:DECEASED)$")
    censor_pats = re.compile(r"^(0|NONE|DISEASEFREE|CENSORED|0:DISEASEFREE|0:CENSORED)$")
    out = pd.Series(np.nan, index=series.index, dtype="Float64")
    out[s.str.match(event_pats)]  = 1.0
    out[s.str.match(censor_pats)] = 0.0
    return out


def _to_months(series: pd.Series, source_col: str) -> pd.Series:
    """Convert to numeric months. Days → months if source column name suggests days."""
    numeric = pd.to_numeric(series, errors="coerce")
    if "DAY" in source_col.upper():
        numeric = numeric / 30.4375
    return numeric


def _to_years(series: pd.Series, source_col: str) -> pd.Series:
    """Convert age to years. Days → years if source column name suggests days."""
    numeric = pd.to_numeric(series, errors="coerce")
    if "DAY" in source_col.upper():
        numeric = numeric / 365.25
    return numeric


# ── per-study standardization ─────────────────────────────────────────────────

def standardize(df: pd.DataFrame, study_id: str, cohort_code: int) -> pd.DataFrame:
    """
    Add canonical survival columns, recode status variables, keep all original
    columns, and prepend STUDY_ID / COHORT_CODE.
    """
    out = df.copy()

    # OS time
    os_src = _first_present(out, OS_TIME_CANDIDATES)
    if os_src:
        out["OS_MONTHS"] = _to_months(out[os_src], os_src)
        if os_src != "OS_MONTHS":
            print(f"    OS time sourced from '{os_src}', converted to months")
    else:
        out["OS_MONTHS"] = np.nan
        print("    WARNING: no OS time column found")

    # OS status
    os_st_src = _first_present(out, OS_STATUS_CANDIDATES)
    if os_st_src:
        out["OS_STATUS"] = _recode_os_status(out[os_st_src])
        if os_st_src != "OS_STATUS":
            print(f"    OS status sourced from '{os_st_src}'")
    else:
        out["OS_STATUS"] = np.nan
        print("    WARNING: no OS status column found")

    # EFS time
    efs_src = _first_present(out, EFS_TIME_CANDIDATES)
    if efs_src:
        out["EFS_MONTHS"] = _to_months(out[efs_src], efs_src)
        if efs_src != "EFS_MONTHS":
            print(f"    EFS time sourced from '{efs_src}', converted to months")
    else:
        out["EFS_MONTHS"] = np.nan
        print("    INFO: no EFS time column found")

    # EFS status
    efs_st_src = _first_present(out, EFS_STATUS_CANDIDATES)
    if efs_st_src:
        out["EFS_STATUS"] = _recode_efs_status(out[efs_st_src])
        if efs_st_src != "EFS_STATUS":
            print(f"    EFS status sourced from '{efs_st_src}'")
    else:
        out["EFS_STATUS"] = np.nan
        print("    INFO: no EFS status column found")

    # Age
    age_src = _first_present(out, AGE_CANDIDATES)
    if age_src:
        out["AGE_YEARS"] = _to_years(out[age_src], age_src)
    else:
        out["AGE_YEARS"] = np.nan

    # Stage (keep raw string; one canonical column)
    stage_src = _first_present(out, STAGE_CANDIDATES)
    out["STAGE"] = out[stage_src].astype(str) if stage_src else np.nan

    # Histology
    hist_src = _first_present(out, HIST_CANDIDATES)
    out["HISTOLOGY"] = out[hist_src].astype(str) if hist_src else np.nan

    # Housekeeping
    out.insert(0, "COHORT_CODE", cohort_code)
    out.insert(0, "STUDY_ID",    study_id)

    return out


# ── summary printer ───────────────────────────────────────────────────────────

def print_summary(df: pd.DataFrame, study_id: str) -> None:
    total = len(df)
    os_ok   = df["OS_MONTHS"].notna() & df["OS_STATUS"].notna()
    efs_ok  = df["EFS_MONTHS"].notna() & df["EFS_STATUS"].notna()

    os_events  = int(df.loc[os_ok,  "OS_STATUS"].sum())
    efs_events = int(df.loc[efs_ok, "EFS_STATUS"].sum())

    print(f"\n  {'Patients total':<35}: {total}")
    print(f"  {'With OS data (time + status)':<35}: {os_ok.sum()}  "
          f"  [{os_events} events / {os_ok.sum() - os_events} censored]")
    print(f"  {'With EFS data (time + status)':<35}: {efs_ok.sum()}"
          f"  [{efs_events} events / {efs_ok.sum() - efs_events} censored]")

    age_col = "AGE_YEARS"
    if age_col in df.columns and df[age_col].notna().any():
        a = df[age_col].dropna()
        print(f"  {'Age (years) median [range]':<35}: "
              f"{a.median():.1f}  [{a.min():.1f} – {a.max():.1f}]")

    if "SEX" in df.columns:
        sex_counts = df["SEX"].value_counts().to_dict()
        print(f"  {'Sex':<35}: {sex_counts}")

    if "STAGE" in df.columns:
        stage_counts = (
            df["STAGE"]
            .replace("nan", np.nan)
            .dropna()
            .value_counts()
            .to_dict()
        )
        if stage_counts:
            print(f"  {'Stage distribution':<35}: {stage_counts}")

    if "HISTOLOGY" in df.columns:
        hist_counts = (
            df["HISTOLOGY"]
            .replace("nan", np.nan)
            .dropna()
            .value_counts()
            .head(5)
            .to_dict()
        )
        if hist_counts:
            print(f"  {'Histology (top 5)':<35}: {hist_counts}")


# ── main ──────────────────────────────────────────────────────────────────────

def build_study_list() -> list[str]:
    """
    Combine the user-specified primary study with those from the discovery CSV.
    Preserves order; primary study first.
    """
    ids: list[str] = [PRIMARY_STUDY_ID]

    if os.path.isfile(CSV_INDEX):
        df_idx = pd.read_csv(CSV_INDEX)
        survivors = df_idx.loc[df_idx["has_survival_data"] == True, "study_id"].tolist()
        for sid in survivors:
            if sid not in ids:
                ids.append(sid)
    else:
        print(f"WARNING: {CSV_INDEX} not found — only attempting {PRIMARY_STUDY_ID}")

    return ids


def main() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    study_ids = build_study_list()
    print(f"Studies to process: {study_ids}\n")
    print("=" * 70)

    cohort_code = 0
    processed: list[str] = []

    for study_id in study_ids:
        print(f"\n[{cohort_code}] {study_id}")
        print("-" * 50)

        # Validate study exists
        meta = get(f"studies/{study_id}")
        if meta is None:
            print(f"  Skipping — study '{study_id}' could not be retrieved.")
            print(f"  (It may not exist on cBioPortal. Studies found in discovery: "
                  f"wt_target_2018_pub, wt_target_gdc)")
            continue

        print(f"  Name   : {meta.get('name')}")
        print(f"  Samples: {meta.get('allSampleCount')}")

        # Download and pivot clinical data
        print("  Downloading clinical data...")
        raw = fetch_patient_clinical(study_id)
        if raw.empty:
            print("  Skipping — no clinical data returned.")
            continue

        print(f"  Raw columns ({len(raw.columns)}): {sorted(raw.columns.tolist())}")

        # Standardize
        std = standardize(raw, study_id, cohort_code)
        print_summary(std, study_id)

        # Save
        out_path = os.path.join(DATA_DIR, f"clinical_{study_id}.csv")
        std.to_csv(out_path, index=False)
        print(f"\n  Saved → {os.path.abspath(out_path)}  ({len(std)} rows × {len(std.columns)} cols)")

        processed.append(study_id)
        cohort_code += 1

    # Final recap
    print("\n" + "=" * 70)
    print(f"DONE. Successfully processed {len(processed)} / {len(study_ids)} studies:")
    for i, sid in enumerate(processed):
        path = os.path.join(DATA_DIR, f"clinical_{sid}.csv")
        df = pd.read_csv(path)
        os_n = (df["OS_MONTHS"].notna() & df["OS_STATUS"].notna()).sum()
        print(f"  [{i}] {sid:<30}  {len(df):>4} patients  {os_n:>4} with OS data")
    print("=" * 70)


if __name__ == "__main__":
    main()
