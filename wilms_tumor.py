"""
Pediatric Wilms' Tumor (TARGET, 2018) — cBioPortal dataset accessor.

Study ID  : wt_target_2018_pub
Reference : https://www.cbioportal.org/study?id=wt_target_2018_pub
Cohort    : 657 WES samples (NCI-TARGET consortium, hg19)

Data layers available
---------------------
Mutations       wt_target_2018_pub_mutations          657 samples  (MAF / WES)
CNA (GISTIC)    wt_target_2018_pub_gistic              129 samples
RNA-Seq RPKM    wt_target_2018_pub_rna_seq_mrna        130 samples
mRNA microarray wt_target_2018_pub_mrna                127 samples
Methylation     wt_target_2018_pub_methylation_hm450   126 samples
miRNA           wt_target_2018_pub_mirna               (continuous)

Sample lists
------------
wt_target_2018_pub_all            657  all samples
wt_target_2018_pub_sequenced      657  samples with mutations
wt_target_2018_pub_cna            129  samples with CNA
wt_target_2018_pub_rna_seq_mrna   130  samples with RNA-Seq
wt_target_2018_pub_mrna           127  samples with microarray
wt_target_2018_pub_3way_complete  101  samples with mutation + CNA + expression
wt_target_2018_pub_discovery      133  DISCOVERY cohort
wt_target_2018_pub_validation     524  VALIDATION cohort
"""

from cbioportal_client import CBioPortalClient
import pandas as pd

STUDY_ID = "wt_target_2018_pub"

# Molecular profile IDs
PROFILES = {
    "mutations":          "wt_target_2018_pub_mutations",
    "cna_gistic":         "wt_target_2018_pub_gistic",
    "rna_seq_rpkm":       "wt_target_2018_pub_rna_seq_mrna",
    "rna_seq_zscore":     "wt_target_2018_pub_rna_seq_mrna_median_Zscores",
    "rna_seq_zscore_all": "wt_target_2018_pub_rna_seq_mrna_median_all_sample_Zscores",
    "mrna_array":         "wt_target_2018_pub_mrna",
    "mrna_zscore":        "wt_target_2018_pub_mrna_median_Zscores",
    "mrna_zscore_all":    "wt_target_2018_pub_mrna_median_all_sample_Zscores",
    "methylation_hm450":  "wt_target_2018_pub_methylation_hm450",
    "mirna":              "wt_target_2018_pub_mirna",
}

# Sample list IDs
SAMPLE_LISTS = {
    "all":             "wt_target_2018_pub_all",
    "sequenced":       "wt_target_2018_pub_sequenced",
    "cna":             "wt_target_2018_pub_cna",
    "rna_seq":         "wt_target_2018_pub_rna_seq_mrna",
    "mrna":            "wt_target_2018_pub_mrna",
    "complete":        "wt_target_2018_pub_3way_complete",
    "discovery":       "wt_target_2018_pub_discovery",
    "validation":      "wt_target_2018_pub_validation",
    "methylation":     "wt_target_2018_pub_methylation_hm450",
}


class WilmsTumorDataset:
    """
    High-level accessor for the Pediatric Wilms' Tumor (TARGET, 2018) dataset.

    Example usage
    -------------
    >>> from wilms_tumor import WilmsTumorDataset
    >>> wt = WilmsTumorDataset()
    >>> print(wt.study_info())
    >>> clinical = wt.clinical_data()
    >>> mutations = wt.mutations()
    >>> tp53_expr = wt.expression_for_genes(["TP53", "WT1", "CTNNB1"])
    """

    def __init__(self):
        self.client = CBioPortalClient()
        self.study_id = STUDY_ID
        self.profiles = PROFILES
        self.sample_lists = SAMPLE_LISTS

    # ------------------------------------------------------------------
    # Study overview
    # ------------------------------------------------------------------

    def study_info(self) -> dict:
        """Return study-level metadata."""
        return self.client.get_study(self.study_id)

    def molecular_profiles(self) -> pd.DataFrame:
        return self.client.get_molecular_profiles(self.study_id)

    def available_sample_lists(self) -> pd.DataFrame:
        return self.client.get_sample_lists(self.study_id)

    def clinical_attributes(self) -> pd.DataFrame:
        return self.client.get_clinical_attributes(self.study_id)

    # ------------------------------------------------------------------
    # Clinical data
    # ------------------------------------------------------------------

    def clinical_data(self) -> pd.DataFrame:
        """
        Patient-level clinical data in wide format.

        Key columns include: SEX, AGE, CANCER_TYPE_DETAILED,
        HISTOLOGY_CLASSIFICATION_IN_PRIMARY_TUMOR, OS_MONTHS, OS_STATUS,
        ANALYSIS_COHORT, SOMATIC_STATUS, MUTATION_COUNT,
        FRACTION_GENOME_ALTERED, TMB_NONSYNONYMOUS.
        """
        return self.client.get_patient_clinical_data(self.study_id)

    def sample_clinical_data(self) -> pd.DataFrame:
        """Sample-level clinical annotations in wide format."""
        return self.client.get_sample_clinical_data(self.study_id)

    def samples(self) -> pd.DataFrame:
        """All 657 samples with patientId and sampleId."""
        return self.client.get_samples(self.study_id)

    def patients(self) -> pd.DataFrame:
        """All patients in the study."""
        return self.client.get_patients(self.study_id)

    # ------------------------------------------------------------------
    # Mutations
    # ------------------------------------------------------------------

    def mutations(self, cohort: str = "all") -> pd.DataFrame:
        """
        Somatic mutations from whole-exome sequencing (657 samples).

        Parameters
        ----------
        cohort : one of the keys in SAMPLE_LISTS (default 'all')

        Returns a DataFrame with columns including:
        sampleId, patientId, gene.hugoGeneSymbol, mutationType,
        proteinChange, chr, startPosition, endPosition,
        referenceAllele, variantAllele, variantType,
        functionalImpactScore, fisValue, linkXvar, etc.
        """
        sample_list_id = self.sample_lists.get(cohort, self.sample_lists["all"])
        return self.client.get_mutations(
            molecular_profile_id=self.profiles["mutations"],
            sample_list_id=sample_list_id,
        )

    def mutations_for_gene(self, gene_symbol: str, cohort: str = "all") -> pd.DataFrame:
        """
        Fetch mutations for a single gene (e.g. 'WT1', 'TP53', 'CTNNB1').
        """
        gene = self.client.get_gene(gene_symbol)
        entrez_id = gene["entrezGeneId"]
        sample_list_id = self.sample_lists.get(cohort, self.sample_lists["all"])
        return self.client.get_mutations_for_gene(
            molecular_profile_id=self.profiles["mutations"],
            sample_list_id=sample_list_id,
            entrez_gene_id=entrez_id,
        )

    def mutation_frequency(self, cohort: str = "all", top_n: int = 20) -> pd.Series:
        """
        Return top_n most frequently mutated genes as a Series (gene → sample count).
        """
        df = self.mutations(cohort=cohort)
        if df.empty:
            return pd.Series(dtype=int)
        # Each row is one mutation; count unique samples per gene
        gene_col = "gene" if "gene" in df.columns else None
        if gene_col and isinstance(df[gene_col].iloc[0], dict):
            df["hugoSymbol"] = df["gene"].apply(lambda g: g.get("hugoGeneSymbol", ""))
        elif "hugoGeneSymbol" in df.columns:
            df["hugoSymbol"] = df["hugoGeneSymbol"]
        else:
            df["hugoSymbol"] = df.get("hugoGeneSymbol", "Unknown")
        freq = (
            df.groupby("hugoSymbol")["sampleId"]
            .nunique()
            .sort_values(ascending=False)
            .head(top_n)
        )
        return freq

    # ------------------------------------------------------------------
    # Copy-number alterations
    # ------------------------------------------------------------------

    def cna(self, cohort: str = "cna") -> pd.DataFrame:
        """
        Discrete CNA calls from GISTIC 2.0 (129 samples).

        Values: -2 homozygous deletion | -1 hemizygous deletion |
                 0 neutral | 1 gain | 2 high-level amplification
        """
        sample_list_id = self.sample_lists.get(cohort, self.sample_lists["cna"])
        return self.client.get_discrete_cna(
            molecular_profile_id=self.profiles["cna_gistic"],
            sample_list_id=sample_list_id,
        )

    # ------------------------------------------------------------------
    # Gene expression
    # ------------------------------------------------------------------

    def expression_for_genes(
        self,
        gene_symbols: list[str],
        data_type: str = "rna_seq_rpkm",
        cohort: str = "rna_seq",
    ) -> pd.DataFrame:
        """
        Return expression values for a list of gene symbols.

        Parameters
        ----------
        gene_symbols : list of HGNC symbols, e.g. ['WT1', 'TP53', 'CTNNB1']
        data_type    : key in PROFILES — 'rna_seq_rpkm' (default),
                       'rna_seq_zscore', 'rna_seq_zscore_all',
                       'mrna_array', 'mrna_zscore', 'mrna_zscore_all'
        cohort       : key in SAMPLE_LISTS — 'rna_seq' (default) or 'mrna'

        Returns a wide DataFrame: rows = samples, columns = gene symbols.
        """
        profile_id = self.profiles.get(data_type, self.profiles["rna_seq_rpkm"])
        sample_list_id = self.sample_lists.get(cohort, self.sample_lists["rna_seq"])

        entrez_ids = []
        symbol_map = {}
        for sym in gene_symbols:
            try:
                gene = self.client.get_gene(sym)
                entrez_ids.append(gene["entrezGeneId"])
                symbol_map[gene["entrezGeneId"]] = gene["hugoGeneSymbol"]
            except Exception as e:
                print(f"Warning: could not resolve gene '{sym}': {e}")

        if not entrez_ids:
            return pd.DataFrame()

        df = self.client.get_molecular_data(
            molecular_profile_id=profile_id,
            sample_list_id=sample_list_id,
            entrez_gene_ids=entrez_ids,
        )
        if df.empty:
            return df

        df["hugoSymbol"] = df["entrezGeneId"].map(symbol_map)
        wide = df.pivot_table(
            index="sampleId", columns="hugoSymbol", values="value", aggfunc="first"
        ).reset_index()
        return wide

    # ------------------------------------------------------------------
    # Convenience summary
    # ------------------------------------------------------------------

    def summary(self) -> None:
        """Print a human-readable summary of the dataset."""
        info = self.study_info()
        print("=" * 60)
        print(f"Study : {info['name']}")
        print(f"ID    : {info['studyId']}")
        print(f"Genome: {info['referenceGenome']}")
        print(f"Desc  : Whole-exome sequencing of {info['allSampleCount']} samples")
        print("-" * 60)
        print(f"Total samples          : {info['allSampleCount']}")
        print(f"Sequenced (WES)        : {info['sequencedSampleCount']}")
        print(f"CNA (GISTIC)           : {info['cnaSampleCount']}")
        print(f"mRNA (RNA-Seq)         : {info['mrnaRnaSeqSampleCount']}")
        print(f"mRNA (microarray)      : {info['mrnaMicroarraySampleCount']}")
        print(f"Complete (mut+CNA+expr): {info['completeSampleCount']}")
        print("=" * 60)
        print("\nMolecular profile IDs:")
        for key, pid in PROFILES.items():
            print(f"  {key:<22} {pid}")
        print("\nSample list IDs:")
        for key, slid in SAMPLE_LISTS.items():
            print(f"  {key:<14} {slid}")
