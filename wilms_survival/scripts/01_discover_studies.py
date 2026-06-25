"""
01_discover_studies.py

Search cBioPortal for all Wilms tumor / nephroblastoma studies, inventory
their molecular profiles and survival data, and write a summary CSV.

Output: wilms_survival/data/available_studies.csv
"""

import sys
import os
import requests
import pandas as pd

BASE_URL = "https://www.cbioportal.org/api"
OUT_CSV  = os.path.join(os.path.dirname(__file__), "..", "data", "available_studies.csv")
KEYWORDS = ["wilms", "nephroblastoma"]


# ── helpers ───────────────────────────────────────────────────────────────────

def get(path: str, params: dict | None = None) -> list | dict | None:
    url = f"{BASE_URL}/{path.lstrip('/')}"
    try:
        r = requests.get(url, params=params, timeout=30)
        r.raise_for_status()
        return r.json()
    except requests.exceptions.HTTPError as e:
        print(f"  HTTP error {e.response.status_code} — {url}", file=sys.stderr)
    except requests.exceptions.RequestException as e:
        print(f"  Request failed — {e}", file=sys.stderr)
    return None


def search_studies(keyword: str) -> list[dict]:
    data = get("studies", params={"keyword": keyword, "pageSize": 500})
    return data if isinstance(data, list) else []


def get_molecular_profiles(study_id: str) -> list[dict]:
    data = get(f"studies/{study_id}/molecular-profiles")
    return data if isinstance(data, list) else []


def get_clinical_attributes(study_id: str) -> list[dict]:
    data = get(f"studies/{study_id}/clinical-attributes")
    return data if isinstance(data, list) else []


# ── profile classification ─────────────────────────────────────────────────

def classify_profiles(profiles: list[dict]) -> dict[str, bool]:
    """
    Scan a list of molecular profile objects and return boolean flags
    for the data types we care about.
    """
    flags = {"has_mutations": False, "has_cna": False, "has_mrna_zscores": False}
    for p in profiles:
        alt = p.get("molecularAlterationType", "")
        dtype = p.get("datatype", "")
        if alt == "MUTATION_EXTENDED":
            flags["has_mutations"] = True
        if alt == "COPY_NUMBER_ALTERATION" and dtype == "DISCRETE":
            flags["has_cna"] = True
        if alt == "MRNA_EXPRESSION" and dtype == "Z-SCORE":
            flags["has_mrna_zscores"] = True
    return flags


def has_survival_data(clinical_attrs: list[dict]) -> bool:
    """Return True if OS_MONTHS (or OS_STATUS) is present as a clinical attribute."""
    ids = {a.get("clinicalAttributeId", "") for a in clinical_attrs}
    return bool(ids & {"OS_MONTHS", "OS_STATUS", "OS_DAYS", "OVERALL_SURVIVAL"})


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    # 1. Collect studies from all keywords, de-duplicate by studyId
    print("Searching cBioPortal for Wilms tumor / nephroblastoma studies...")
    seen: set[str] = set()
    studies: list[dict] = []
    for kw in KEYWORDS:
        found = search_studies(kw)
        print(f"  keyword '{kw}' → {len(found)} result(s)")
        for s in found:
            sid = s.get("studyId", "")
            if sid and sid not in seen:
                seen.add(sid)
                studies.append(s)

    if not studies:
        print("No studies found. Check network / API availability.")
        sys.exit(1)

    print(f"\nTotal unique studies found: {len(studies)}\n")
    print("=" * 70)

    # 2. For each study, retrieve profiles and clinical attributes
    rows: list[dict] = []
    for s in studies:
        sid   = s.get("studyId", "N/A")
        name  = s.get("name", "N/A")
        desc  = s.get("description", "")
        n_pts = s.get("allSampleCount", 0)          # closest proxy to patient count
        cancer = (s.get("cancerType") or {}).get("name", s.get("cancerTypeId", "N/A"))

        print(f"Study ID   : {sid}")
        print(f"Name       : {name}")
        print(f"Cancer type: {cancer}")
        print(f"Samples    : {n_pts}")
        # strip HTML tags from description for readability
        import re
        clean_desc = re.sub(r"<[^>]+>", "", desc).strip()
        if clean_desc:
            # wrap at ~80 chars
            words, line, lines = clean_desc.split(), "", []
            for w in words:
                if len(line) + len(w) + 1 > 80:
                    lines.append(line)
                    line = w
                else:
                    line = (line + " " + w).strip()
            if line:
                lines.append(line)
            print(f"Description: {lines[0]}")
            for l in lines[1:]:
                print(f"             {l}")

        # molecular profiles
        profiles = get_molecular_profiles(sid)
        pf_flags = classify_profiles(profiles)
        profile_names = [p.get("name", "") for p in profiles]
        print(f"Profiles   : {len(profiles)} found")
        for pname in profile_names:
            print(f"             • {pname}")
        print(f"  has_mutations   : {pf_flags['has_mutations']}")
        print(f"  has_cna         : {pf_flags['has_cna']}")
        print(f"  has_mrna_zscores: {pf_flags['has_mrna_zscores']}")

        # survival data
        clinical_attrs = get_clinical_attributes(sid)
        surv = has_survival_data(clinical_attrs)
        print(f"  has_survival    : {surv}")
        print("-" * 70)

        rows.append({
            "study_id":          sid,
            "name":              name,
            "cancer_type":       cancer,
            "n_patients":        n_pts,
            "has_mutations":     pf_flags["has_mutations"],
            "has_cna":           pf_flags["has_cna"],
            "has_mrna_zscores":  pf_flags["has_mrna_zscores"],
            "has_survival_data": surv,
        })

    # 3. Build summary DataFrame and save
    df = pd.DataFrame(rows).sort_values("n_patients", ascending=False).reset_index(drop=True)

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved → {os.path.abspath(OUT_CSV)}\n")

    # 4. Print summary table
    print("SUMMARY TABLE")
    print("=" * 70)
    col_widths = {"study_id": 30, "name": 36, "n_patients": 9,
                  "has_mutations": 4, "has_cna": 4, "has_mrna_zscores": 4, "has_survival_data": 4}

    header = (
        f"{'Study ID':<30}  {'Name':<36}  {'N':>6}  "
        f"{'Mut':>4}  {'CNA':>4}  {'Expr':>4}  {'Surv':>4}"
    )
    print(header)
    print("-" * len(header))
    for _, row in df.iterrows():
        mut  = "yes" if row["has_mutations"]     else "no"
        cna  = "yes" if row["has_cna"]           else "no"
        expr = "yes" if row["has_mrna_zscores"]  else "no"
        surv = "yes" if row["has_survival_data"] else "no"
        print(
            f"{row['study_id']:<30}  {row['name']:<36.36}  "
            f"{int(row['n_patients']):>6}  {mut:>4}  {cna:>4}  {expr:>4}  {surv:>4}"
        )
    print("=" * 70)


if __name__ == "__main__":
    main()
