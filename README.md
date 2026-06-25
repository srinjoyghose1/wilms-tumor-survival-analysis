# Wilms Tumor Survival Analysis — cBioPortal (TARGET 2018)

Kaplan-Meier and Cox proportional hazards survival analysis of genomic alterations and mRNA expression data for pediatric Wilms tumor, using publicly available data from the [cBioPortal for Cancer Genomics](https://www.cbioportal.org/).

---

## Dataset

| Field | Value |
|---|---|
| Study name | Pediatric Wilms' Tumor (TARGET, 2018) |
| cBioPortal study ID | `wt_target_2018_pub` |
| Patients | 652 |
| Disease | Wilms' Tumor (Nephroblastoma) |
| Primary endpoint | Event-Free Survival (EFS) |
| Secondary endpoint | Overall Survival (OS) |
| EFS events | 211 / 652 |
| OS events | 114 / 652 |
| mRNA expression subset | 125 / 652 (RNA-seq z-scores) |

Data were retrieved programmatically via the cBioPortal REST API v3. No local files beyond the API outputs are required to reproduce the analysis.

---

## Genes Analyzed

| Gene | Alteration Rule | Pathway |
|---|---|---|
| TP53 | Any mutation | Cell Cycle / p53 Axis |
| MYCN | Amplification (CNA ≥ +2) | Transcription Regulation |
| SIX1 | Hotspot Q177R | Transcription / Progenitor |
| SIX2 | Hotspot Q177R | Transcription / Progenitor |
| DROSHA | Truncating mutation | RNA Processing |
| DGCR8 | Truncating mutation | RNA Processing |
| CTNNB1 | Hotspot codons 32–45 | WNT Signaling |
| AMER1 (WTX) | Truncating mutation | WNT Signaling |
| WT1 | Truncating mutation or homodeletion (CNA ≤ −2) | WNT Signaling |
| MLLT1 | Truncating mutation | Chromatin Remodeling |
| NIPBL | Truncating mutation | Pan-Cancer |
| MDM4 ★ | Amplification (CNA ≥ +2) | p53 Axis (PDMR focus) |
| MDM2 ★ | Amplification (CNA ≥ +2) | p53 Axis (PDMR focus) |
| IGF2 ★ | CNA ≥ +1 or mRNA z-score > 1.5 | Imprinting (PDMR focus) |

★ = Pediatric Drug-Resistant / Patient-Derived Models Repository (PDMR) focus genes

Because individual gene mutation rates are low in this cohort (a known biological feature of pediatric Wilms tumor), pathway-level OR-classifiers were used to aggregate patients with any alteration within a pathway.

---

## Results Summary

### Pathway-level (OR-classifier) analyses — all significant

| Classifier | Genes Combined | N Altered | N Wildtype | Endpoint | Median Altered (mo) | Median Wildtype (mo) | Log-rank p | HR (multivariate) | 95% CI |
|---|---|---|---|---|---|---|---|---|---|
| All Primary Genes (OR) | TP53, MYCN, SIX1, SIX2, DROSHA, DGCR8, CTNNB1, AMER1, WT1, MLLT1, NIPBL | 44 | 608 | EFS | 9.9 | NR | < 0.001 | 4.21 | [2.90–6.12] |
| MYCN/SIX1-SIX2 (OR) | MYCN, SIX1/SIX2 | 23 | 629 | EFS | 10.5 | NR | < 0.001 | 3.78 | [2.31–6.19] |
| TP53/MDM4/MDM2 (OR) ★ | TP53, MDM4, MDM2 | 21 | 631 | OS | NR | NR | < 0.001 | 3.55 | [1.85–6.83] |
| All Primary Genes (OR) | TP53, MYCN, SIX1, SIX2, DROSHA, DGCR8, CTNNB1, AMER1, WT1, MLLT1, NIPBL | 44 | 608 | OS | NR | NR | < 0.001 | 3.38 | [2.06–5.54] |
| MYCN/SIX1-SIX2 (OR) | MYCN, SIX1/SIX2 | 23 | 629 | OS | NR | NR | < 0.001 | 3.27 | [1.70–6.31] |
| TP53/MDM4/MDM2 (OR) ★ | TP53, MDM4, MDM2 | 21 | 631 | EFS | 14.2 | NR | 0.001 | 2.34 | [1.33–4.11] |

HR = hazard ratio from multivariate Cox PH regression, adjusted for patient age. NR = Not Reached.

### Individual gene analyses
All 16 individual genes fell below the minimum sample size gate (n ≥ 20 altered patients per arm) required for Kaplan-Meier analysis. This is consistent with the known low somatic mutation rate in Wilms tumor. No individual gene reached statistical significance in this cohort.

### mRNA expression analyses (n = 125 patients with RNA-seq data)
None of the 12 expression-based analyses (6 genes × 2 endpoints) reached p < 0.05. MDM4 high expression showed a trend toward shorter EFS (p = 0.059, uncorrected). All analyses were underpowered due to the small expression subset.

---

## Pipeline Scripts

Scripts are run sequentially. All data are fetched live from the cBioPortal API — no manual downloads required.

| Script | Description | Output |
|---|---|---|
| `01_discover_studies.py` | Discover Wilms tumor studies and molecular profiles via cBioPortal API | `results/study_summary.csv` |
| `02_download_clinical.py` | Download and standardize clinical data (EFS, OS, age) | `data/clinical_*.csv` |
| `03_download_genomics.py` | Download mutation and copy-number data for target genes | `data/mutations_*.csv`, `data/cna_*.csv` |
| `04_classify_patients.py` | Classify each patient as ALTERED / WILDTYPE per gene-specific rules | `processed/master_classification.csv` |
| `05_kaplan_meier.py` | Kaplan-Meier survival analysis for all genes and OR-classifiers | `results/km_results.csv`, `plots/*.pdf` |
| `06_cox_regression.py` | Multivariate Cox PH regression (gene flag + age covariates) | `results/cox_results.csv` |
| `07_expression_km.py` | mRNA expression-based KM and Cox analysis (z-score stratification) | `results/expression_km_results.csv`, `plots/KM_EXPR_*.pdf` |
| `08_summary_table.py` | Compile publication-ready summary tables and narrative paragraphs | `results/Table*.csv`, `results/findings_paragraphs.txt` |
| `09_qc_report.py` | Quality control checks across all pipeline outputs | Printed QC report (5 checks) |
| `10_build_html_report.py` | Build self-contained HTML presentation report with embedded plots | `results/WilmsTumor_cBioPortal_Findings.html` |

---

## Reproduction

### Requirements

```
python >= 3.9
requests
pandas
numpy
matplotlib
lifelines==0.30.3
scipy
seaborn
```

Install:
```bash
pip install requests pandas numpy matplotlib lifelines==0.30.3 scipy seaborn
```

### Run

```bash
cd wilms_survival/scripts
python 01_discover_studies.py
python 02_download_clinical.py
python 03_download_genomics.py
python 04_classify_patients.py
python 05_kaplan_meier.py
python 06_cox_regression.py
python 07_expression_km.py
python 08_summary_table.py
python 09_qc_report.py
python 10_build_html_report.py
```

The final HTML report is written to `wilms_survival/results/WilmsTumor_cBioPortal_Findings.html`.

---

## Project Structure

```
wilms_survival/
├── scripts/           # Analysis scripts 01–10
├── data/              # Raw API downloads (clinical, mutation, CNA, mRNA)
├── processed/         # master_classification.csv (652 × 24 columns)
├── results/           # Summary tables, KM/Cox CSVs, HTML report
└── plots/             # KM PDF plots (18 files)
```

---

## Data Source and Limitations

- All data are sourced exclusively from cBioPortal (`wt_target_2018_pub`). No external databases, manual curation, or private data are used.
- Alteration frequencies are low (< 7% for any single classifier), consistent with published Wilms tumor genomics literature.
- This is a retrospective, observational analysis of publicly available data. Results should be interpreted in that context and are not intended to establish clinical recommendations.
- mRNA expression analysis is limited by the small sequenced subset (125/652 patients) and should be considered exploratory.
- All AGE values in the cohort are integers (pediatric patients); the AGE covariate in Cox models serves as a continuous adjustment variable.

---

## API Reference

[cBioPortal Web API Documentation](https://docs.cbioportal.org/web-api-and-clients/)

Base URL used: `https://www.cbioportal.org/api/v3`
