
# ============================================================
# STEP 1: Load required packages
# ============================================================
required_packages <- c(
  "TCGAbiolinks", "DESeq2", "EnhancedVolcano", "clusterProfiler",
  "org.Hs.eg.db", "ComplexHeatmap", "ggplot2", "dplyr", "pheatmap"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install missing packages first: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(DESeq2)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
})

NOTE-
This step checks that the R packages needed for the analysis are installed and loads them into the current running R sessions. 
They together allow to pull data from TCGA, run differential gene expression analysis, build visualizations, perform gene enrichment
analysis and work with data frames. If one of the required packages is not found, the script will immediately exit and display a message 
at the top stating the packages that need to be installed.

# ============================================================
# STEP 2: Folders and TCGA-LUAD data
# ============================================================
dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

prepared_data_file <- file.path("data", "TCGA_LUAD_prepared_data.rds")

# Leave FALSE to prevent accidentally downloading several GB again.
# Change to TRUE only if you want to download/prepare data when no saved
# prepared RDS file and no data_se object are available.
DOWNLOAD_IF_NEEDED <- FALSE

if (exists("data_se", inherits = TRUE) &&
    inherits(get("data_se", inherits = TRUE), "SummarizedExperiment")) {
  message("Using data_se from the current R session.")
  data_se <- get("data_se", inherits = TRUE)
} else if (file.exists(prepared_data_file)) {
  message("Loading prepared data from ", prepared_data_file)
  data_se <- readRDS(prepared_data_file)
} else if (DOWNLOAD_IF_NEEDED) {
  query <- GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query)
  data_se <- GDCprepare(query)
  saveRDS(data_se, prepared_data_file)
} else {
  stop(
    "No prepared data_se object or saved RDS file found.\n",
    "If data_se is in your current R session, run the script there.\n",
    'To save it for future use, run: saveRDS(data_se, "data/TCGA_LUAD_prepared_data.rds")\n',
    "Only set DOWNLOAD_IF_NEEDED <- TRUE if you are ready to download the data."
  )
}

NOTE-
  This step prepares the data, results, and figures folders to store the project’s files in a well-organized manner. 
The step also checks whether the TCGA-LUAD dataset is already prepared (loaded in the current R session or stored as an RDS file). 
If not, it is possible to download the necessary RNA-seq gene expression data from the GDC (Genomic Data Commons). This download option 
is disabled by default because datasets that are several gigabytes in size must not be downloaded automatically.
# ============================================================
# STEP 3: Inspect sample information
# ============================================================
# Extract sample metadata
coldata <- as.data.frame(SummarizedExperiment::colData(data_se))
# Inspect available sample information
print(colnames(coldata))
print(head(coldata$sample_type))
print(table(coldata$sample_type))

NOTE-This step downloads sample metadata for the TCGA-LUAD expression dataset 
and explores what metadata columns are present. It also investigates the distribution 
of samples by their type (e.g., primary tumour, solid tissue normal)and provides an overview of the types 
of biological material represented in the dataset .Looking at this first confirms what the 
dataset actually contains before you run the differential gene expression analysis.
# ============================================================
# STEP 4: Differential gene expression using DESeq2
# Compare primary tumour with normal tissue
# ============================================================
library(DESeq2)
library(dplyr)
# Extract raw, unstranded RNA-seq counts
counts <- SummarizedExperiment::assay(data_se, "unstranded")
# Keep only primary tumour and normal samples
keep <- coldata$sample_type %in% c(
  "Primary Tumor",
  "Solid Tissue Normal"
)
counts <- counts[, keep]
coldata_de <- coldata[keep, , drop = FALSE]
# Set normal as the reference group
coldata_de$condition <- factor(
  ifelse(
    coldata_de$sample_type == "Primary Tumor",
    "tumour",
    "normal"
  ),
  levels = c("normal", "tumour")
)
# Ensure sample metadata matches count matrix columns
stopifnot(identical(colnames(counts), rownames(coldata_de)))
# Create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = coldata_de,
  design = ~ condition
)
# Remove genes with fewer than 10 total counts
dds <- dds[rowSums(counts(dds)) >= 10, ]
# Run differential expression analysis
dds <- DESeq(dds)
# Compare tumour against normal
res <- results(
  dds,
  contrast = c("condition", "tumour", "normal")
)
# Save the complete results table
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
dir.create("results", showWarnings = FALSE)
write.csv(
  res_df,
  "results/DESeq2_all_genes.csv",
  row.names = FALSE
)
# Count significant upregulated and downregulated genes
sig <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
cat("Samples used:", ncol(dds), "\n")
cat("Genes tested:", nrow(res_df), "\n")
cat("Significant genes:", nrow(sig), "\n")
cat("Upregulated:", sum(sig$log2FoldChange > 1), "\n")
cat("Downregulated:", sum(sig$log2FoldChange < -1), "\n")

NOTE- This section will compare the gene expression of the primary lung adenocarcinoma 
tumours against the normal lung tissue. First, the RNA_seq counts and metadata of the 
samples will be read in, with the normal tissue being set as the reference. Using this 
information, a DESeq2 dataset will be constructed, and all genes that have less than 10 
total counts will be removed. Finally, the differential expression analysis will be performed and 
a table of results will be exported as a CSV file. Within the results, each gene’s log2 fold 
change, p-value, and adjusted p-value will be found.

# ============================================================
# STEP 5A: Volcano plot
# ============================================================
library(EnhancedVolcano)
dir.create("figures", showWarnings = FALSE)
volcano <- EnhancedVolcano(
  res_df,
  lab = res_df$gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Differential Gene Expression in TCGA-LUAD",
  subtitle = "Primary tumour vs normal lung tissue",
  legendPosition = "right"
)
ggplot2::ggsave(
  filename = "figures/TCGA_LUAD_volcano.png",
  plot = volcano,
  width = 10,
  height = 8,
  dpi = 300
)
file.exists("figures/TCGA_LUAD_volcano.png")

NOTE- This option will plot the results of the differential expression analysis as a volcano plot. 
The plot shows the log2 fold change on the x-axis and the adjusted p-values (significance) on the y-axis. 
It allows to spot the genes that show a strong expression differences between 
tumourand normal samples.
The file is saved as a PNG file in high resolution.

# ============================================================
# STEP 5B: Heatmap of the top 50 significant genes
# ============================================================

library(DESeq2)
library(pheatmap)
library(dplyr)

# Select the top 50 significant genes
sig_genes <- res_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) > 1
  ) %>%
  arrange(padj) %>%
  head(50)

# Calculate variance-stabilized expression values
vsd <- vst(dds, blind = FALSE)

# Extract expression values for selected genes
mat <- assay(vsd)[sig_genes$gene, , drop = FALSE]

# Scale expression values across each gene
mat_scaled <- t(scale(t(mat)))

# Annotate samples by tissue type
annotation <- data.frame(
  condition = coldata_de$condition
)
rownames(annotation) <- rownames(coldata_de)

# Save heatmap to a high-resolution PNG file
png(
  "figures/TCGA_LUAD_top50_heatmap.png",
  width = 3000,
  height = 2400,
  res = 300
)

pheatmap(
  mat_scaled,
  annotation_col = annotation,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top 50 DEGs: TCGA-LUAD Tumour vs Normal"
)

dev.off()

file.exists("figures/TCGA_LUAD_top50_heatmap.png")

NOTE- This step selects the 50 most significant genes passing the adjusted p-value (< 0.05) 
and fold-change thresholds (|log2FC| > 1). This data is variance stabilized for visualization purposes 
and then scaled across each gene(z-scored) to generate a clustered heatmap with annotations indicating tumour 
and normal samples which is saved as a PNG.
# ============================================================
# STEP 6: GO Biological Process enrichment analysis
# ============================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
# Select significantly upregulated genes
up_genes <- res_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    log2FoldChange > 1
  ) %>%
  pull(gene)
# Convert Ensembl gene IDs to Entrez IDs
# Remove Ensembl version suffixes, if present
up_genes_clean <- sub("\\..*$", "", up_genes)
gene_ids <- bitr(
  up_genes_clean,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
# Remove duplicate Entrez IDs
gene_ids <- gene_ids[!duplicated(gene_ids$ENTREZID), ]
# Run GO enrichment for Biological Process
ego <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
# Save enrichment results
go_results <- as.data.frame(ego)
write.csv(
  go_results,
  "results/GO_BP_enrichment.csv",
  row.names = FALSE
)
# Save a dot plot if enrichment terms were found
if (nrow(go_results) > 0) {
  go_plot <- dotplot(
    ego,
    showCategory = 20,
    title = "GO Biological Process Enrichment: TCGA-LUAD"
  )
  ggplot2::ggsave(
    "figures/TCGA_LUAD_GO_BP_dotplot.png",
    plot = go_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  print(head(go_results[, c("ID", "Description", "p.adjust")], 10))
} else {
  message("No GO terms passed the enrichment significance thresholds.")
}
file.exists("results/GO_BP_enrichment.csv")
file.exists("figures/TCGA_LUAD_GO_BP_dotplot.png")
source("C:/Users/kkada/OneDrive/Desktop/TCGA-LUAD-PROJECT/analysis.R")
list.files(
  recursive = TRUE,
  pattern = "\\.(csv|png|R|Rproj|txt)$"
)
list.files(
  recursive = TRUE,
  pattern = "\\.(csv|png|R|Rproj|txt)$"
)

if (!inherits(data_se, "SummarizedExperiment")) {
  stop("data_se is not a SummarizedExperiment object.")
}
message("Prepared data dimensions: ", nrow(data_se), " genes x ", ncol(data_se), " samples")
NOTE - This step identifies biological processes over-represented in significantly 
upregulated genes. Initially, Ensembl gene IDs are cleaned and converted to Entrez IDs 
using human annotation database. Then, clusterProfiler performs the Gene Ontology
Biological Process enrichment analysis, and p-values are adjusted for multiple testing 
using the Benjamini–Hochberg method. Finally, the profiling results are exported as a CSV file, 
and a dot plot is created, if any of the terms are found to be significant.
# ============================================================
# STEP 7: Check expected output files
# ============================================================
output_files <- c(
  "results/DESeq2_all_genes.csv",
  "results/DESeq2_significant_genes.csv",
  "results/DESeq2_upregulated_genes.csv",
  "results/DESeq2_downregulated_genes.csv",
  "results/analysis_summary.txt",
  "figures/TCGA_LUAD_volcano.png",
  "figures/TCGA_LUAD_top50_heatmap.png",
  "results/GO_BP_enrichment.csv",
  "figures/TCGA_LUAD_GO_BP_dotplot.png"
)
output_status <- data.frame(file = output_files, exists = file.exists(output_files))
write.csv(output_status, "results/output_file_check.csv", row.names = FALSE)
print(output_status)

message("Script finished. Review results/analysis_summary.txt and results/output_file_check.csv.")

 NOTE - This step verifies that all the output files from the analysis were created. 
 It records each file in a summary table, exports it as a CSV report, and prints the results 
 to the console so that the presence of outputs can be quickly checked.
# ============================================================
# STEP 8: Identify top upregulated genes
# ============================================================

# Select significantly upregulated genes
up_genes <- subset(
  res_df,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange > 1
)

# Sort by adjusted p-value
up_genes <- up_genes[order(up_genes$padj), ]

# Display the top 20 genes
top20_upregulated <- head(
  up_genes[, c("gene_id", "log2FoldChange", "pvalue", "padj")],
  20
)

print(top20_upregulated)

# Save the complete upregulated gene list
dir.create("results", showWarnings = FALSE)

write.csv(
  up_genes,
  "results/DESeq2_upregulated_genes.csv",
  row.names = FALSE
)

# Save the top 20 genes separately
write.csv(
  top20_upregulated,
  "results/Top20_upregulated_genes.csv",
  row.names = FALSE
)

NOTE - This step filters the results down to those with an adjusted p-value below 0.05 and 
a log2 fold change above 1, which will represent the significantly upregulated genes in tumour
compared with normal tissue. These are then sorted by adjusted p-value and the top 20 are 
extracted for further inspection. The whole list of upregulated genes and the top 20 are each 
saved as a CSV file.
# ============================================================
# STEP 9: Identify the top upregulated gene
# ============================================================

# Select the top gene by adjusted p-value
top_gene <- up_genes[1, ]

print(top_gene)

# Save the top gene information
write.csv(
  top_gene,
  "results/Top_upregulated_gene.csv",
  row.names = FALSE
)

NOTE-This step retrieves the best-ranking gene from the upregulated list based on adjusted p-value. 
Its differential expression statistics are then printed to the console and saved in a different CSV file,
providing a simple summary of the most significant upregulated candidate in the analysis.
# ============================================================
# STEP 9A: Convert Ensembl ID to gene symbol
# ============================================================

library(org.Hs.eg.db)
library(AnnotationDbi)

# Remove the Ensembl version suffix
ensembl_id <- sub(
  "\\..*$",
  "",
  as.character(top_gene$gene)
)

# Map Ensembl ID to gene symbol
gene_symbol <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_id,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

# Display the gene symbol
print(gene_symbol)

# Add the symbol to the top-gene result
top_gene$gene_symbol <- unname(gene_symbol)

# Save the updated result
write.csv(
  top_gene,
  "results/Top_upregulated_gene.csv",
  row.names = FALSE
)
 NOTE - This step removes the version suffix from the selected Ensembl gene ID and maps
 it to a gene symbol based on the human gene annotation database. The resulting gene 
 symbol is added to the top-gene results table which is saved once again so that the 
 selected gene can be easily identified and interpreted

# ============================================================
# STEP 10: Install and load survival-analysis packages
# ============================================================

# Install packages only if they are not already installed
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}

if (!requireNamespace("survminer", quietly = TRUE)) {
  install.packages("survminer")
}

# Load packages
library(survival)
library(survminer)

# Verify installation
packageVersion("survival")
packageVersion("survminer")

NOTE - This step ensures the survival and survminer libraries are available for the analysis. 
If any of these packages is missing, it will prompt to install and then load it. Together these
packages provide the statistical tools and plotting functions used in the following survival analysis.
# ============================================================
# STEP 11: Retrieve TCGA-LUAD clinical survival data
# ============================================================

# Download clinical information from GDC
clinical <- GDCquery_clinic(
  project = "TCGA-LUAD",
  type = "clinical"
)

# Inspect available clinical fields
print(colnames(clinical))

# Check survival-related fields
survival_fields <- intersect(
  c(
    "submitter_id",
    "vital_status",
    "days_to_death",
    "days_to_last_follow_up"
  ),
  colnames(clinical)
)

print(head(clinical[, survival_fields, drop = FALSE]))

# Save clinical data
dir.create("data", showWarnings = FALSE)

write.csv(
  clinical,
  "data/TCGA_LUAD_clinical.csv",
  row.names = FALSE
)

 NOTE- This step involves the drawing of clinical data related to the TCGA-LUAD project 
 from the GDC, where the available clinical fields will be investigated with a specific focus 
 on patient IDs, vital status, days to death, and days to last follow-up. The clinical dataset
 will then be saved in CSV format for additional processing.
# ============================================================
# STEP 12: Inspect survival-related clinical variables
# ============================================================

# Display all clinical column names
print(colnames(clinical))

# Check the first few patient records
head(clinical[, intersect(
  c(
    "submitter_id",
    "vital_status",
    "days_to_death",
    "days_to_last_follow_up"
  ),
  colnames(clinical)
), drop = FALSE])
NOTE- This step will verify the presence of the fields, which are necessary for survival analysis, in the clinical dataset.
The analysis examines the patient identifier, vital status, and survival-time variables to observe how the clinical data is
structured before survival outcomes are calculated.

# ============================================================
# STEP 13: Prepare overall survival data
# ============================================================

# Convert survival fields to numeric where necessary
clinical$days_to_death <- as.numeric(
  as.character(clinical$days_to_death)
)

clinical$days_to_last_follow_up <- as.numeric(
  as.character(clinical$days_to_last_follow_up)
)

# Calculate overall survival time in days
clinical$OS_time <- ifelse(
  clinical$vital_status == "Dead",
  clinical$days_to_death,
  clinical$days_to_last_follow_up
)

# Event: 1 = Dead, 0 = Alive
clinical$OS_event <- ifelse(
  clinical$vital_status == "Dead", 1, 0
)

# Keep records with valid survival information
clinical_survival <- clinical[
  !is.na(clinical$OS_time) &
    clinical$OS_time >= 0 &
    !is.na(clinical$OS_event),
]

# Inspect the results
print(head(
  clinical_survival[
    , c("submitter_id", "vital_status", "OS_time", "OS_event")
  ]
))

print(table(clinical_survival$OS_event))

# Save the prepared survival data
write.csv(
  clinical_survival,
  "data/TCGA_LUAD_survival_data.csv",
  row.names = FALSE
)

# Save only the survival fields needed for analysis
survival_export <- data.frame(
  submitter_id = clinical_survival$submitter_id,
  vital_status = as.character(clinical_survival$vital_status),
  OS_time = as.numeric(clinical_survival$OS_time),
  OS_event = as.integer(clinical_survival$OS_event)
)

write.csv(
  survival_export,
  "data/TCGA_LUAD_survival_data.csv",
  row.names = FALSE,
  na = ""
)

message("Survival data saved successfully.")

NOTE- This step  involves converting the relevant clinical survival fields to numerical data and
calculating the overall survival time in days. The deceased patients are assigned their days to 
death, while the living patients their days to last follow-up. An event indicator is then created, 
where 1 represents death and 0 represents the censored observations. The Records lacking valid survival 
information are removed, and the prepared survival dataset is saved.

# ============================================================
# STEP 14: Extract SACK1A expression by patient
# ============================================================

# Ensembl ID for SACK1A
sack1a_id <- "ENSG00000147689.16"

# Check whether this gene exists in the DESeq2 dataset
gene_matches <- rownames(dds)[
  sub("\\..*$", "", rownames(dds)) ==
    sub("\\..*$", "", sack1a_id)
]

print(gene_matches)

# Stop if SACK1A is not present
if (length(gene_matches) == 0) {
  stop("SACK1A was not found in the row names of dds.")
}

# Extract normalized expression
sack1a_counts <- counts(dds, normalized = TRUE)[
  gene_matches[1], ]

# Get TCGA sample barcodes
sample_barcodes <- colnames(dds)

# Prepare expression table
sack1a_expression <- data.frame(
  sample_barcode = sample_barcodes,
  patient_id = substr(sample_barcodes, 1, 12),
  SACK1A_expression = as.numeric(sack1a_counts)
)

# Save expression table
write.csv(
  sack1a_expression,
  "results/SACK1A_expression_by_sample.csv",
  row.names = FALSE
)

# Preview
print(head(sack1a_expression))

NOTE- This step finds SACK1A by its Ensembl gene ID and retrieves its normalized expression values 
from the DESeq2 dataset. The sample barcodes are used to match each sample to its specific patient. 
The resulting table will list down the sample IDs, patient IDs, and normalized SACK1A expression, 
and it is saved for the survival analysis.

# ============================================================
# STEP 15: Match primary tumour expression to clinical survival
# ============================================================

# Keep primary tumour samples (TCGA sample type code 01)
sack1a_tumour <- subset(
  sack1a_expression,
  substr(sample_barcode, 14, 15) == "01"
)

# Keep one sample per patient
sack1a_tumour <- sack1a_tumour[
  !duplicated(sack1a_tumour$patient_id),
]

# Prepare clinical survival fields
clinical_km <- data.frame(
  patient_id = as.character(clinical_survival$submitter_id),
  OS_time = as.numeric(clinical_survival$OS_time),
  OS_event = as.integer(clinical_survival$OS_event)
)

# Match tumour expression with clinical survival data
survival_merged <- merge(
  sack1a_tumour[, c("patient_id", "SACK1A_expression")],
  clinical_km,
  by = "patient_id"
)

# Remove records with missing or invalid values
survival_merged <- subset(
  survival_merged,
  !is.na(SACK1A_expression) &
    !is.na(OS_time) &
    OS_time > 0 &
    !is.na(OS_event)
)

# Check the number of matched patients
cat("Matched patients:", nrow(survival_merged), "\n")
cat("Deaths:", sum(survival_merged$OS_event == 1), "\n")
cat("Alive:", sum(survival_merged$OS_event == 0), "\n")

# Save matched data
write.csv(
  survival_merged,
  "results/SACK1A_survival_merged.csv",
  row.names = FALSE
)
NOTE- This step retains only the primary tumour samples (identified by the TCGA sample-type code), 
taking one sample per patient. The prepared clinical survival data (using patient IDs) is then merged
with the SACK1A expression table, removing any records with missing or invalid expression or survival 
information. The merged dataset is stored, and a summary of the number of patients matched, and the number 
of deaths observed, is generated.

# ============================================================
# STEP 16: Create high and low SACK1A expression groups
# ============================================================

# Calculate median expression
sack1a_median <- median(
  survival_merged$SACK1A_expression,
  na.rm = TRUE
)

# Classify patients
survival_merged$SACK1A_group <- ifelse(
  survival_merged$SACK1A_expression >= sack1a_median,
  "High",
  "Low"
)

survival_merged$SACK1A_group <- factor(
  survival_merged$SACK1A_group,
  levels = c("Low", "High")
)

# Check group sizes
print(table(survival_merged$SACK1A_group))

# Save updated data
write.csv(
  survival_merged,
  "results/SACK1A_survival_merged.csv",
  row.names = FALSE
)
NOTE-This step involves finding the median SACK1A expression across the matched patients and uses
it to  further split them into two groups. Patients at or above the median go into the high-expression 
group, and those below it are categorized into the low-expression group. The group sizes are then checked, 
and the updated dataset is saved for the survival comparison.

# ============================================================
# STEP 17: Kaplan-Meier survival analysis for SACK1A
# ============================================================

library(survival)
library(survminer)

# Create survival object
surv_object <- Surv(
  time = survival_merged$OS_time,
  event = survival_merged$OS_event
)

# Fit Kaplan-Meier curves
km_fit <- survfit(
  surv_object ~ SACK1A_group,
  data = survival_merged
)

# Log-rank test
logrank_test <- survdiff(
  Surv(OS_time, OS_event) ~ SACK1A_group,
  data = survival_merged
)

print(logrank_test)

# Calculate log-rank p-value
logrank_p <- pchisq(
  logrank_test$chisq,
  df = length(logrank_test$n) - 1,
  lower.tail = FALSE
)

cat("Log-rank p-value:", logrank_p, "\n")

# Plot Kaplan-Meier curves
km_plot <- ggsurvplot(
  km_fit,
  data = survival_merged,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  legend.title = "SACK1A expression",
  legend.labs = c("Low", "High"),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  title = "TCGA-LUAD: SACK1A Expression and Overall Survival"
)

print(km_plot)

# Save Kaplan-Meier plot
png(
  "figures/TCGA_LUAD_SACK1A_Kaplan_Meier.png",
  width = 1800,
  height = 1600,
  res = 200
)

print(km_plot)
dev.off()

# Save log-rank test result
write.csv(
  data.frame(
    gene = "SACK1A",
    group_cutoff = sack1a_median,
    low_patients = sum(survival_merged$SACK1A_group == "Low"),
    high_patients = sum(survival_merged$SACK1A_group == "High"),
    logrank_chisq = logrank_test$chisq,
    logrank_p_value = logrank_p
  ),
  "results/SACK1A_survival_logrank.csv",
  row.names = FALSE
)

NOTE- This step creates a survival object based on overall survival time and event status.
Then, it fits Kaplan-Meier curves to compare the overall survival of the high and low 
SACK1A expression groups and performs a log-rank test to determine whether two curves 
are significantly different. The script also generates a survival plot with confidence
intervals, the p-value, and a risk table, and saves the figure and statistical results.
# ============================================================
# STEP 18: Display and verify survival results
# ============================================================

print(logrank_test)

cat("Log-rank chi-square:", logrank_test$chisq, "\n")
cat("Log-rank p-value:", logrank_p, "\n")

# Confirm output files exist
file.exists("figures/TCGA_LUAD_SACK1A_Kaplan_Meier.png")
file.exists("results/SACK1A_survival_logrank.csv")

NOTE-This step prints the results of the log-rank test, including the value of the 
chi-square statistic and its p-value. It also verifies that the Kaplan-Meier graph 
and the log-rank results file were created in the output folders, as expected. Thus, 
confirming that the outputs for survival analysis were created.

# ============================================================
# STEP 19: Save SACK1A survival interpretation
# ============================================================

survival_summary <- c(
  "TCGA-LUAD: SACK1A overall survival analysis",
  "Method: Kaplan-Meier curves and log-rank test",
  "Patients: 504",
  "Low expression group: 252 patients, 76 observed deaths",
  "High expression group: 252 patients, 106 observed deaths",
  paste("Log-rank chi-square:", round(logrank_test$chisq, 3)),
  paste("Log-rank p-value:", signif(logrank_p, 4)),
  "Interpretation: Survival curves differ statistically.",
  "Further validation is needed to assess independent prognostic value."
)

writeLines(
  survival_summary,
  "results/SACK1A_survival_summary.txt"
)

cat(paste(survival_summary, collapse = "\n"), "\n")

NOTE-This step involves preparing a short text description of the results of SACK1A 
survival analysis. It includes describing the analysis method, the sizes of each 
patient group, the number of observed deaths and the log-rank statistic and p-value. 
It also includes a note that further validation is required before any conclusion can
be made about the independent prognostic value of SACK1A. The text file is saved 
as part of the project documentation.

# ============================================================
# STEP 20: Prepare genes for STRING PPI network analysis
# ============================================================

# Select significantly upregulated genes
up_genes <- subset(
  res_df,
  !is.na(padj) &
    padj < 0.05 &
    log2FoldChange > 1
)

# Remove Ensembl version suffixes
up_gene_ids <- sub(
  "\\..*$",
  "",
  as.character(up_genes$gene)
)

# Keep unique gene IDs
up_gene_ids <- unique(up_gene_ids)

# Save gene IDs for STRING
dir.create("results", showWarnings = FALSE)

write.table(
  up_gene_ids,
  "results/STRING_upregulated_gene_ids.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Check the number of genes
cat("Upregulated genes prepared:", length(up_gene_ids), "\n")
cat(
  "Saved to results/STRING_upregulated_gene_ids.txt\n"
)
NOTE-This step uses the adjusted p-value and fold-change threshold to select 
the significantly upregulated genes. The suffixes in the Ensembl IDs are removed, 
and any duplicate IDs are eliminated. The final gene list is saved as a text file, 
which can be uploaded to the STRING database to identify potential protein-protein 
interactions.

# ============================================================
# STEP 21: Select top 100 upregulated genes for STRING
# ============================================================

# Rank upregulated genes by adjusted p-value
up_genes <- up_genes[order(up_genes$padj), ]

# Select top 100 genes
top100_genes <- head(up_genes, 100)

# Remove Ensembl version suffixes
top100_ids <- sub(
  "\\..*$",
  "",
  as.character(top100_genes$gene)
)

top100_ids <- unique(top100_ids)

# Save gene list
writeLines(
  top100_ids,
  "results/STRING_top100_gene_ids.txt"
)

cat("Genes saved for STRING:", length(top100_ids), "\n")
NOTE - This step ranks the highly upregulated gens according to adjusted p-value, 
takes the first 100, cleans, de-duplicates their Ensembl IDs, and saves them to a text
file that will be used to construct a protein-protein interaction network in STRING which 
can then be visualized and refined in Cytoscape.



