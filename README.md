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

3.
