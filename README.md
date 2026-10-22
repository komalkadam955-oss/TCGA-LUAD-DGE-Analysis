# TCGA-LUAD-DGE-Analysis
This project looks at what goes different in Lung Adenocarcinoma at the genomic level (LUAD) using the RNA- sequencing data from The Cancer Genome Atlas (TCGA). Several approaches were integrated to understand what is occuring at the molecular level: differential gene expression analysis to identify genes demonstrating different behavior in tumours, Gene ontology (GO) enrichment to understand the gene function,Kaplan-Meier Analysis to observe how they relate to patient outcomes and protein -protein interaction (PPI) networks to observe the associations between the proteins.
The analysis was performed in R using Bioconductor pages , with STRING and Cytoscape used for protein interaction analysis and network visualization.

OBJECTIVES :

-Identify the genes which are differentially expressed in TCGA-LUAD

-Examine biological processes involved in different gene expression

-Investigate the relationship between gene expression and overall survival

-Analyze the protein interaction patterns of chosen elevated genes

-Identify potential genes and pathways for further inquiry.

TOOLS USED:

R and RStudio - statistical computing and analysis

TCGAbiolinks - TCGA data acquisition 

DESeq2 - differential gene expression analysis  

clusterProfiler and org.Hs.eg.db - functional enrichment analysis 

EnhancedVolcano -volcano plot visualization  

pheatmap - data visualization and heatmap generation

survival - Kaplan–Meier and log-rank survival analysis 

STRING - protein–protein interaction analysis 

Cytoscape - network visualization and layout editing

ANALYSIS WORKFLOW:

1. Data acquisition and preparation

DESeq2 was used for retrieving RNA sequencing data for the TCGA-LUAD cohort and prepared it for analysis, which includes organizing the count data and aligning the tumor and normal sample for comparison.

2.Differential Gene Expression Analysis

Differential Gene Expression Analysis was asssessed using DESeq2 for identification of genes whose expression differed significantly between tumour and normal samples. Genes were classified as differentially expressed when they met the specfic criteria:-

-Adjusted p-value < 0.05

Absolute log2 fold change > 1

Results were visualized using a volcano plot , which shows the difference between statistical significance and effect size and a heatmap that shows how the expression of the top genes varies across samples.

3. Gene Ontology Enrichment Analysis

With the use of clusterProfiler tool, the GO biological process enrichment was conducted in order to determine biological processes that are significantly represented among the differentially expressed genes in LUAD.


4. Survival Analysis


The Kaplan-Meier analysis was implemented to evaluate the correlation between SACK1A (FAM83A) expression and overall survival in TCGA-LUAD patients.

Gene- SACK1A(FAM83A), EnsembleID-ENSG00000147689.16

Expression source- DESeq2 normalized counts 

Groups-  Patients were split into low expression and high expression

Statistical Test- Log rank test used to compare survival curves


5. Protein -Protein Interaction Network Analysis

The selected upregulated genes were examined in STRING  to explore the known and anticipated protein connections. The obtained netwok was further visualized and manually arranged in Cytoscape.

Input- Selected upregulated genes

Database- STRING ( known and predicted genes)

Confidence Threshold- Minimum combined interaction score of 0.700

Visualization- Cytoscape with manual layout

Network genes-DLGAP5, BUB1B, TOP2A, CCNB1, KIF2C, CDC20

Note- FAM83A (SACK1A)  was present in the STRING results but appears to be connected to SPP1  at a minimum combined score of 0.700 in the network visualized. This interaction was visible in the Cytoscape network, however its corresponding score should be confirmed in the exported interaction table.


KEY FINDINGS:

1. Differential Gene Expression

The analysis results in identification of 15, 156 differentially expressed genes that followed the specific criterion (Adjusted p-value < 0.05
and log2 fold change > 1)

Upregulated genes- 11, 868

Downregulated genes- 3,288

2. Functional Enrichment Analysis

The GO biological Enrichment analysis of differentially expressed genes accentuated the top five terms appearing on the dot plot.

-Organelle fission

-Nuclear divison

-Skeletal system development

-Meiotic cell cycle

-Chromosome segregation

These findings demonstrate the enrichment of genes associated with cell division, chromosomal organization, and skeletal system development. However, enrichment does not imply that these processes are active,rather, it shows an overrepresentation of genes in these processes.


3. SURVIVAL ANALYSIS

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

A statistically significant difference between the survival curves was revealed by the log-rank test, indicating a correlation between expression group and overall survival. Neither causality nor independent prognostic value are established by this.

NOTE-The survival script designates Ensembl ID ENSG00000147689.16 as SACK1A, whereas STRING displays FAM83A. These seem to be the same gene (SACK1A is another name), however the identification should be checked to the dataset's reference annotation before making genespecific biological claims.

4. Protein-Protein Interaction Network Analysis
Tools such as STRING and Cytoscape were used to assess protein associations among the selected upregulated genes
(minimum confidence score 0.700)

 Genes in the network- DLGAP5, BUB1B, TOP2A, CCNB1, KIF2C, CDC20

 Biological process inn association- cell cycle regulation, mitotic progression , chromosome segregation

 These genes have potential for further investigation, however network connectivity does not determine biological importance,functional activity or causality.
