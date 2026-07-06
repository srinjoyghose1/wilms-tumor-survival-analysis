---

## Pathway-Level Survival Analysis

<!-- Why pathways: individual genes all failed n≥20 threshold in this 652-patient cohort -->

> **Why pathway analysis?** No individual gene in the TARGET 2018 cohort (N=652) reached
> the minimum 20-patient threshold for Kaplan-Meier analysis — the most frequent gene
> (MYCN) had only 18 altered patients. This reflects a well-characterized feature of
> Wilms tumor biology: no single gene drives more than ~15% of cases. Pathway-level
> OR-grouping aggregates all alterations sharing a common biological mechanism,
> providing the statistical power needed for survival analysis while preserving
> biological interpretability.

---

### Five Molecular Pathways

<!-- Evidence-based pathway definitions from Gadd 2017 (TARGET), Perotti 2024, Treger 2019 -->

| Pathway | Member Genes | N Altered | Cohort Coverage | Literature Basis |
|---|---|---|---|---|
| **WNT / Renal Differentiation** | WT1, CTNNB1, AMER1, NF1 | 10 | 1.5% | Co-cluster in TARGET expression analysis (Gadd 2017); converge on WNT-β-catenin disruption and mesenchymal-to-epithelial transition |
| **Renal Progenitor / SIX-miRNAPG** | SIX1, SIX2, DROSHA, DGCR8, DICER1, XPO5, DIS3L2 | 10 | 1.5% | Share common output of perpetuating renal progenitor state; blastemal histology association; enriched at relapse (Wegert 2015; Perotti 2024) |
| **Chromatin Remodeling** | MLLT1, BCOR, BCORL1, ARID1A, SMARCA4, BRD7, CREBBP, EP300, HDAC4, ASXL1, CHD4, MAP3K4, CTR9 | 9 | 1.4% | Named as functional class in Treger 2019: epigenetic regulators of renal progenitor differentiation via transcriptional elongation and chromatin modification |
| **p53 / Cell Cycle** | TP53, MYCN, MDM2, MDM4, MAX, CHEK2, PALB2 | 35 | 5.4% | TP53 defines diffuse anaplastic Wilms tumor; MYCN amplification predicts poor EFS/OS; MDM2/MDM4 provide alternative p53 inactivation (Perotti 2024) |
| **IGF2 / Growth Factor** | IGF2 (CNA + expr), NIPBL, NONO, COL6A3, FGFR1, ACTB | 21 | 3.2% | IGF2 loss of imprinting (~70% prevalence) drives PI3K-AKT via IGF1R; NIPBL disrupts cohesin-mediated control of the IGF2/H19 imprinted domain (Treger 2019) |

**Combined:** 10.6% of the cohort (69/652 patients) carry at least one alteration across any pathway.

---

### Alteration Frequency

<!-- How much of the cohort each pathway captures vs individual genes -->

![Pathway alteration frequency](wilms_survival/plots/PATHWAY_frequency_summary.png)

---

### Survival Analysis Results

<!-- KM and Cox results — first valid survival statistics from this dataset -->

| Pathway | N Altered | N Wildtype | Cohort % | Endpoint | Median Altered (mo) | Median WT (mo) | Log-rank p | HR [95% CI] | Sig? |
|---|---|---|---|---|---|---|---|---|---|
| WNT | 10 | 642 | 1.5% | EFS | — | — | — | — | ⚠️ skipped |
| WNT | 10 | 642 | 1.5% | OS | — | — | — | — | ⚠️ skipped |
| SIX-miRNA | 10 | 642 | 1.5% | EFS | — | — | — | — | ⚠️ skipped |
| SIX-miRNA | 10 | 642 | 1.5% | OS | — | — | — | — | ⚠️ skipped |
| Chromatin | 9 | 643 | 1.4% | EFS | — | — | — | — | ⚠️ skipped |
| Chromatin | 9 | 643 | 1.4% | OS | — | — | — | — | ⚠️ skipped |
| p53 | 35 | 617 | 5.4% | EFS | 12.2 | NR | 6.67e-10 | 3.26 [2.15–4.95] | ✅ |
| p53 | 35 | 617 | 5.4% | OS | 25.0 | NR | 1.16e-11 | 4.53 [2.72–7.55] | ✅ |
| IGF2 | 21 | 631 | 3.2% | EFS | 10.8 | NR | 6.10e-08 | 3.48 [2.06–5.88] | ✅ |
| IGF2 | 21 | 631 | 3.2% | OS | NR | NR | 6.71e-04 | — | ✅ |
| Any-Pathway | 69 | 583 | 10.6% | EFS | 10.8 | NR | 4.23e-26 | 4.56 [3.32–6.26] | ✅ |
| Any-Pathway | 69 | 583 | 10.6% | OS | NR | NR | 4.64e-12 | 3.71 [2.44–5.66] | ✅ |

---

### Key Findings

<!-- One block per significant result; honest null-result statement if none significant -->

**p53** — EFS

> Patients with any alteration in the p53 / Cell Cycle pathway (35 patients, 5.4% of cohort)
> showed significantly worse EFS compared to unaltered patients
> (median: **12.2 vs. NR months**; log-rank p = 6.67e-10; HR = 3.26 [2.15–4.95]).

![KM curve](wilms_survival/plots/PATHWAY_P53_EFS.png)

**p53** — OS

> Patients with any alteration in the p53 / Cell Cycle pathway (35 patients, 5.4% of cohort)
> showed significantly worse OS compared to unaltered patients
> (median: **25.0 vs. NR months**; log-rank p = 1.16e-11; HR = 4.53 [2.72–7.55]).

![KM curve](wilms_survival/plots/PATHWAY_P53_OS.png)

**IGF2** — EFS

> Patients with any alteration in the IGF2 / Growth Factor pathway (21 patients, 3.2% of cohort)
> showed significantly worse EFS compared to unaltered patients
> (median: **10.8 vs. NR months**; log-rank p = 6.10e-08; HR = 3.48 [2.06–5.88]).

![KM curve](wilms_survival/plots/PATHWAY_IGF2_EFS.png)

**IGF2** — OS

> Patients with any alteration in the IGF2 / Growth Factor pathway (21 patients, 3.2% of cohort)
> showed significantly worse OS compared to unaltered patients
> (median: **NR vs. NR months**; log-rank p = 6.71e-04).

![KM curve](wilms_survival/plots/PATHWAY_IGF2_OS.png)

**Any-Pathway** — EFS

> Patients with any alteration in the Any Pathway (Combined) pathway (69 patients, 10.6% of cohort)
> showed significantly worse EFS compared to unaltered patients
> (median: **10.8 vs. NR months**; log-rank p = 4.23e-26; HR = 4.56 [3.32–6.26]).

![KM curve](wilms_survival/plots/PATHWAY_ANY_EFS.png)

**Any-Pathway** — OS

> Patients with any alteration in the Any Pathway (Combined) pathway (69 patients, 10.6% of cohort)
> showed significantly worse OS compared to unaltered patients
> (median: **NR vs. NR months**; log-rank p = 4.64e-12; HR = 3.71 [2.44–5.66]).

![KM curve](wilms_survival/plots/PATHWAY_ANY_OS.png)

---

### Gene Contributions Within Pathways

<!-- Shows which specific genes drive the altered arm in each pathway -->

#### WNT / Renal Differentiation

![Gene contributions](wilms_survival/plots/PATHWAY_WNT_contributions.png)

> The WNT / Renal Differentiation altered arm (n=10) was primarily driven by WT1 (5 patients). CTNNB1 contributed 3 additional unique cases. 0 patients carried alterations in more than one pathway gene.

#### Renal Progenitor / SIX-miRNAPG

![Gene contributions](wilms_survival/plots/PATHWAY_SIX_miRNA_contributions.png)

> The Renal Progenitor / SIX-miRNAPG altered arm (n=10) was primarily driven by SIX1 (4 patients). XPO5 contributed 2 additional unique cases. 0 patients carried alterations in more than one pathway gene.

#### p53 / Cell Cycle

![Gene contributions](wilms_survival/plots/PATHWAY_P53_contributions.png)

> The p53 / Cell Cycle altered arm (n=35) was primarily driven by MYCN (18 patients). TP53 contributed 15 additional unique cases. 6 patients carried alterations in more than one pathway gene.

#### IGF2 / Growth Factor

![Gene contributions](wilms_survival/plots/PATHWAY_IGF2_contributions.png)

> The IGF2 / Growth Factor altered arm (n=21) was primarily driven by IGF2_EXPR (10 patients). IGF2_CNA contributed 5 additional unique cases. 1 patients carried alterations in more than one pathway gene.

#### Any Pathway (Combined)

![Gene contributions](wilms_survival/plots/PATHWAY_ANY_contributions.png)

> The Any Pathway (Combined) altered arm (n=69) was primarily driven by P53 (35 patients). IGF2 contributed 21 additional unique cases. 14 patients carried alterations in more than one pathway.

---

### How to Read This Analysis

<!-- Plain-language guide for research advisors unfamiliar with OR-classifiers -->

**What "pathway-altered" means:** A patient is classified as altered for a pathway
if they carry a qualifying mutation or copy number change in *any* gene in that
pathway — they do not need alterations in all genes. This is an OR logic classifier.

**Cohort coverage vs. significance:** A pathway capturing 30% of 652 patients gives
roughly 196 altered patients — enough for robust KM and Cox analysis. This is why
pathway grouping succeeds where individual genes (max n=18) cannot.

**Interpreting HR:** A hazard ratio > 1 means pathway-altered patients have higher
instantaneous risk of the event (death or relapse) at any point in follow-up vs.
unaltered patients, after adjusting for cohort and age. HR is from multivariate
Cox regression.

**Gene contribution charts** show whether a pathway result is driven by one dominant
gene (e.g., TP53 dominating the p53 arm) or is a genuinely distributed multi-gene
effect. This distinction matters for clinical interpretation.

---

### References

- Gadd S et al. *Nat Genet.* 2017;49:1487–1494. (TARGET — primary cohort and clustering)
- Perotti D et al. *Nat Rev Urol.* 2024;21:158–180. (pathway definitions and prevalence)
- Treger TD et al. *Nat Rev Nephrol.* 2019;15:240–257. (chromatin/elongation gene class)
- Wegert J et al. *Cancer Cell.* 2015;27:298–311. (SIX-miRNAPG co-occurrence and clustering)
- Perlman EJ et al. *Nat Commun.* 2015;6:10013. (MLLT1 YEATS domain mutations)

---