# TCGA-LUAD-DGE-Analysis
This project looks at what goes different in Lung Adenocarcinoma at the genomic level (LUAD) using the RNA- sequencing data from The Cancer Genome Atlas (TCGA). Several approaches were integrated to understand what is occuring at the molecular level: differential gene expression analysis to identify genes demonstrating different behavior in tumours, Gene ontology (GO) enrichment to understand the gene function,Kaplan-Meier Analysis to observe how they relate to patient outcomes and protein -protein interaction (PPI) networks to observe the associations between the proteins.
The analysis was performed in R using Bioconductor packages , with STRING and Cytoscape used for protein interaction analysis and network visualization.

## RESEARCH QUESTION
- 1. Which genes are differentially expressed between LUAD tumour and normal tissue?
- 2. Which biological processes are found to be enriched among the differentially expressed genes?
- 3. Does the analyzed patient cohort's overall survival correlate with gene expression?
- 4. Which patterns of interaction can be observed between the indicated upregulated genes?
     
## DATA

RNA-sequencing gene expression data from the TCGA-LUAD project, downloaded with TCGAbiolinks from the Genomic Data Commons: 539 primary tumor and 59 solid tissue normal samples, analyzed as raw counts. The survival analysis used clinical data and included 504 patients with matched tumor expression and valid survival information.


## TOOLS USED:

- R and RStudio - statistical computing and analysis

- TCGAbiolinks - TCGA data acquisition 

- DESeq2 - differential gene expression analysis  

- clusterProfiler and org.Hs.eg.db - functional enrichment analysis 

- EnhancedVolcano -volcano plot visualization  

- pheatmap - data visualization and heatmap generation

- survival - Kaplan–Meier and log-rank survival analysis 

- STRING - protein–protein interaction analysis 

- Cytoscape - network visualization and layout editing

## ANALYSIS WORKFLOW:

### 1. Data acquisition and preparation

TCGA-LUAD RNA-sequencing data were retrieved with TCGAbiolinks and organized into a count matrix with tumor and normal samples labeled for comparison. Genes with a total read count below 10 across all samples were removed.

### 2.Differential Gene Expression Analysis

Differential Gene Expression Analysis was asssessed using DESeq2 for identification of genes whose expression differed significantly between tumour and normal samples. Genes were classified as differentially expressed when they met the specfic criteria:-

- Adjusted p-value < 0.05 (Benjamini–Hochberg correction)

- Absolute log2 fold change > 1

Results were visualized using a volcano plot , which shows the difference between statistical significance and effect size and a heatmap that shows how the expression of the top genes varies across samples.

### 3. Gene Ontology Enrichment Analysis

With the use of clusterProfiler tool, the GO biological process enrichment was conducted in order to determine biological processes that are significantly represented among the differentially expressed genes in LUAD.Ensembl IDs were converted to Entrez IDs, and processes with an adjusted p-value < 0.05 (Benjamini–Hochberg) were considered enriched.


### 4. Survival Analysis


The Kaplan-Meier analysis was implemented to evaluate the correlation between SACK1A (FAM83A) expression and overall survival in TCGA-LUAD patients.

- Gene- SACK1A(FAM83A), EnsembleID-ENSG00000147689.16

- Expression source- DESeq2 normalized counts,primary tumor samples, one per patient

- Groups-  Patients were split into low expression and high expression (median split)

- Statistical Test- Log rank test used to compare survival curves


### 5. Protein -Protein Interaction Network Analysis

The selected upregulated genes were examined in STRING  to explore the known and anticipated protein connections. The obtained netwok was further visualized and manually arranged in Cytoscape.

- Input- Selected upregulated genes (100 upregulated genes with lowest adjusted p-values)

- Database- STRING ( known and predicted genes)

- Confidence Threshold- Minimum combined interaction score of 0.700

- Visualization- Cytoscape with manual layout



## KEY FINDINGS:

### 1. Differential Gene Expression

The analysis results in identification of 15, 156 differentially expressed genes that followed the specific criterion (Adjusted p-value < 0.05
and log2 fold change > 1)

Upregulated genes- 11, 868 (78.3%)

Downregulated genes- 3,288 (21.7%)

The top 5 most significant genes are -

| Gene | Ensembl ID | Direction | log2 fold change | Adjusted p-value |
|---|---|---|---|---|
| FAM83A (SACK1A) | ENSG00000147689.16 | Up | 6.79 | 4.3 × 10⁻¹⁷⁶ |
| PYCR1 | ENSG00000183010.17 | Up | 3.69 | 3.0 × 10⁻¹⁶⁵ |
| OTUD1 | ENSG00000165312.6 | Down | -2.12 | 9.4 × 10⁻¹⁴⁶ |
| AFAP1-AS1 | ENSG00000272620.2 | Up | 6.31 | 6.3 × 10⁻¹²⁸ |
| EPAS1 | ENSG00000116016.14 | Down | -2.72 | 4.2 × 10⁻¹²⁸ |

### 2. Functional Enrichment Analysis

GO Biological Process enrichment of the upregulated genes identified 237 significant terms (processes). The five terms with the highest gene ratios in the dot plot were:

- Organelle fission

- Nuclear divison

- Skeletal system development

- Meiotic cell cycle

- Nuclear Chromosome segregation

These findings demonstrate the enrichment of genes associated with cell division and chromosomal organization. However, enrichment does not imply that these processes are active,rather, it shows an overrepresentation of genes in these processes.


### 3. SURVIVAL ANALYSIS

 The given table shows the Kaplan Meier analysis of overall survival in TCGA-LUAD patients, divided into high expression and low expression groups.

  | Measure | Result |
|---|---|
| Total patients | 504 |
| Low-expression group | 252 |
| High-expression group | 252 |
| Observed deaths, low-expression group | 76 (30.2%) |
| Observed deaths, high-expression group | 106 (42.1%) |
| Log-rank chi-square | 11.65 |
| Log-rank p-value | 0.000642 |

The high expression group showed lower survival probability for most of the follow up period. As it is an association only, clinical variable such as tumor stage, age and treatment were not adjusted for, so independent prognostic value is not established.


NOTE-The survival script designates Ensembl ID ENSG00000147689.16 as SACK1A, whereas STRING displays FAM83A. These seem to be the same gene (SACK1A is another name), however the identification should be checked to the dataset's reference annotation before making genespecific biological claims.

### 4. Protein-Protein Interaction Network Analysis
Tools such as STRING and Cytoscape were used to assess protein associations among the selected upregulated genes
(minimum confidence score 0.700)

 Main cluster- DLGAP5, BUB1B, TOP2A, CCNB1, KIF2C, CDC20

 Small separate groups- COL10A1, COL11A1, and MMP11; and FUT2 and B3GNT3

 Biological process in association- cell cycle regulation, mitotic progression , chromosome segregation

 These genes have potential for further investigation, however network connectivity does not determine biological importance,functional activity or causality.
 Note- FAM83A (SACK1A)  was recognized by STRING results and it formed no interactions at minimum combined score of 0.700 in the network visualized. However connections appeared at medium confidence threshold (>0.400), although they rest on weaker evidence and are treated here as exploratory.This interaction was visible in the Cytoscape network.


 ## LIMITATIONS
 
  
- Results depend on the TCGA samples, the tumor/normal comparison (539 vs. 59 samples, an unbalanced design), preprocessing choices, and statistical thresholds.
 
- Batch effects and RNA-seq technical variation may influence expression measurements, and bulk RNA-seq mixes signals from tumor, immune, and stromal cells.
 
- GO enrichment shows overrepresentation of genes in a process, not that the process is active, and it was run on upregulated genes only.
  
- The survival analysis used a median split and did not adjust for clinical variables; independent prognostic value has not been tested with Cox regression.
  
- The PPI network used only the top 100 upregulated genes, so downregulated genes such as OTUD1 and EPAS1 were not examined, and STRING associations may reflect experimental data, curated knowledge, co-expression, or predictions.
  
- Independent datasets and experimental studies are needed to validate the candidate genes and interpretations.






## Repository Contents
- TCGA-LUAD Gene Expression Analysis Script - R scripts for each analysis step
- TCGA -LUAD Analysis Figures– volcano plot, heatmap, GO dot plot,  SACK1A survival curve, and  STRING / Cytoscape network images
- TCGA-LUAD Results – differentially expressed gene tables ,GO BP enrichment table, survival results (SACK1A survival log rank and survival summary), Top upregulated gene table and STRING( top 100 gene id, upregulated gene id, string edge list)
- PDF Report
 
