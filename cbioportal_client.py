"""
cBioPortal REST API client.

Base URL: https://www.cbioportal.org/api
Docs:     https://docs.cbioportal.org/web-api-and-clients/
Swagger:  https://www.cbioportal.org/api/swagger-ui/index.html

No authentication is required for the public instance.
"""

import requests
import pandas as pd
from typing import Optional


BASE_URL = "https://www.cbioportal.org/api"


class CBioPortalClient:
    """
    Thin wrapper around the cBioPortal REST API.

    All methods return pandas DataFrames for easy downstream analysis.
    Raw JSON is always available via the corresponding _raw() pattern
    (just call the underlying requests session directly if needed).
    """

    def __init__(self, base_url: str = BASE_URL):
        self.base_url = base_url.rstrip("/")
        self.session = requests.Session()
        self.session.headers.update({"Accept": "application/json"})

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _get(self, path: str, params: Optional[dict] = None) -> list | dict:
        url = f"{self.base_url}/{path.lstrip('/')}"
        response = self.session.get(url, params=params)
        response.raise_for_status()
        return response.json()

    def _post(self, path: str, payload: dict, params: Optional[dict] = None) -> list | dict:
        url = f"{self.base_url}/{path.lstrip('/')}"
        response = self.session.post(url, json=payload, params=params)
        response.raise_for_status()
        return response.json()

    # ------------------------------------------------------------------
    # Studies
    # ------------------------------------------------------------------

    def get_study(self, study_id: str) -> dict:
        """Return metadata for a single study."""
        return self._get(f"studies/{study_id}")

    def get_all_studies(self) -> pd.DataFrame:
        """Return all public studies as a DataFrame."""
        return pd.DataFrame(self._get("studies", params={"pageSize": 10_000}))

    # ------------------------------------------------------------------
    # Clinical data
    # ------------------------------------------------------------------

    def get_clinical_attributes(self, study_id: str) -> pd.DataFrame:
        data = self._get(f"studies/{study_id}/clinical-attributes")
        return pd.DataFrame(data)

    def get_patient_clinical_data(self, study_id: str) -> pd.DataFrame:
        """
        Return all patient-level clinical data for a study as a wide DataFrame
        (one row per patient, one column per attribute).
        """
        raw = self._get(
            f"studies/{study_id}/clinical-data",
            params={"clinicalDataType": "PATIENT", "pageSize": 100_000},
        )
        if not raw:
            return pd.DataFrame()
        long = pd.DataFrame(raw)
        return long.pivot_table(
            index="patientId", columns="clinicalAttributeId", values="value", aggfunc="first"
        ).reset_index()

    def get_sample_clinical_data(self, study_id: str) -> pd.DataFrame:
        """
        Return all sample-level clinical data for a study as a wide DataFrame
        (one row per sample, one column per attribute).
        """
        raw = self._get(
            f"studies/{study_id}/clinical-data",
            params={"clinicalDataType": "SAMPLE", "pageSize": 100_000},
        )
        if not raw:
            return pd.DataFrame()
        long = pd.DataFrame(raw)
        return long.pivot_table(
            index="sampleId", columns="clinicalAttributeId", values="value", aggfunc="first"
        ).reset_index()

    # ------------------------------------------------------------------
    # Samples & patients
    # ------------------------------------------------------------------

    def get_samples(self, study_id: str) -> pd.DataFrame:
        data = self._get(f"studies/{study_id}/samples", params={"pageSize": 100_000})
        return pd.DataFrame(data)

    def get_patients(self, study_id: str) -> pd.DataFrame:
        data = self._get(f"studies/{study_id}/patients", params={"pageSize": 100_000})
        return pd.DataFrame(data)

    def get_sample_lists(self, study_id: str) -> pd.DataFrame:
        data = self._get(f"studies/{study_id}/sample-lists")
        return pd.DataFrame(data)

    def get_sample_list_ids(self, study_id: str, sample_list_id: str) -> list[str]:
        """Return the raw list of sample IDs in a sample list."""
        return self._get(f"sample-lists/{sample_list_id}/sample-ids")

    # ------------------------------------------------------------------
    # Molecular profiles
    # ------------------------------------------------------------------

    def get_molecular_profiles(self, study_id: str) -> pd.DataFrame:
        data = self._get(f"studies/{study_id}/molecular-profiles")
        return pd.DataFrame(data)

    # ------------------------------------------------------------------
    # Mutations
    # ------------------------------------------------------------------

    def get_mutations(
        self,
        molecular_profile_id: str,
        sample_list_id: str,
        projection: str = "DETAILED",
    ) -> pd.DataFrame:
        """
        Fetch somatic mutations for all samples in a sample list.

        projection options: ID | SUMMARY | DETAILED | META
        Endpoint: POST /molecular-profiles/{molecularProfileId}/mutations/fetch
        Total mutations for this study: ~750 (from WES of 657 samples).
        """
        payload = {"sampleListId": sample_list_id}
        data = self._post(
            f"molecular-profiles/{molecular_profile_id}/mutations/fetch",
            payload=payload,
            params={"projection": projection, "pageSize": 10_000},
        )
        return pd.DataFrame(data)

    def get_mutations_for_gene(
        self,
        molecular_profile_id: str,
        sample_list_id: str,
        entrez_gene_id: int,
        projection: str = "DETAILED",
    ) -> pd.DataFrame:
        """
        Fetch mutations for a specific gene (by Entrez ID).
        Endpoint: POST /molecular-profiles/{molecularProfileId}/mutations/fetch
        """
        payload = {
            "sampleListId": sample_list_id,
            "entrezGeneIds": [entrez_gene_id],
        }
        data = self._post(
            f"molecular-profiles/{molecular_profile_id}/mutations/fetch",
            payload=payload,
            params={"projection": projection, "pageSize": 10_000},
        )
        return pd.DataFrame(data)

    # ------------------------------------------------------------------
    # Copy-number alterations
    # ------------------------------------------------------------------

    def get_discrete_cna(
        self,
        molecular_profile_id: str,
        sample_list_id: str,
        projection: str = "DETAILED",
    ) -> pd.DataFrame:
        """
        Fetch discrete CNA calls (GISTIC values: -2, -1, 0, 1, 2).
        Endpoint: POST /molecular-profiles/{molecularProfileId}/discrete-copy-number/fetch
        """
        payload = {"sampleListId": sample_list_id}
        data = self._post(
            f"molecular-profiles/{molecular_profile_id}/discrete-copy-number/fetch",
            payload=payload,
            params={"projection": projection, "discreteCopyNumberEventType": "ALL"},
        )
        return pd.DataFrame(data)

    # ------------------------------------------------------------------
    # Gene expression / molecular data
    # ------------------------------------------------------------------

    def get_molecular_data(
        self,
        molecular_profile_id: str,
        sample_list_id: str,
        entrez_gene_ids: list[int],
    ) -> pd.DataFrame:
        """
        Fetch molecular data (expression, methylation, etc.) for a list of genes.
        Endpoint: POST /molecular-profiles/{molecularProfileId}/molecular-data/fetch
        """
        payload = {
            "sampleListId": sample_list_id,
            "entrezGeneIds": entrez_gene_ids,
        }
        data = self._post(
            f"molecular-profiles/{molecular_profile_id}/molecular-data/fetch",
            payload=payload,
            params={"projection": "SUMMARY"},
        )
        return pd.DataFrame(data)

    # ------------------------------------------------------------------
    # Gene lookups
    # ------------------------------------------------------------------

    def get_gene(self, gene_symbol: str) -> dict:
        """Resolve a HGNC gene symbol to Entrez ID and other metadata."""
        data = self._get(f"genes/{gene_symbol}")
        return data

    def search_genes(self, query: str) -> pd.DataFrame:
        data = self._get("genes", params={"keyword": query, "pageSize": 50})
        return pd.DataFrame(data)
