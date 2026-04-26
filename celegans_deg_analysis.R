# =============================================================================
# Differential Gene Expression Analysis: C. elegans Transcription Factor Mutants
# Author: Rhea Charles | University of South Florida
#
# Description:
#   DESeq2-based differential expression analysis comparing wildtype C. elegans
#   against five transcription factor mutant conditions:
#     - pqm-1(rhd90), pqm-1(ok485)       — PQM-1 TF loss-of-function alleles
#     - ceh-60                             — CEH-60 TF knockout
#     - ceh-60;pqm-1(rhd90/ok485)         — double mutants
#
#   Includes count filtering, DESeq2 normalization, significance testing,
#   Pearson correlation heatmap, and hierarchical clustering.
#
# Input:
#   all_samples_exons.txt     — raw exon-level read counts (tab-separated)
#   aa_samples_exons.fpkm.txt — FPKM-normalized expression values
#
# Output:
#   Histogram of count distributions (pre/post filtering)
#   Boxplot of significantly DEGs across conditions
#   Pearson correlation heatmap
#   Hierarchical clustering dendrogram
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

# Install Bioconductor packages if needed
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("DESeq2", "GenomeInfoDb", "GO.db")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}
for (pkg in c("gplots", "RColorBrewer")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(DESeq2)
library(ggplot2)
library(reshape2)

# Update these paths to point to your data directory
counts_path <- "data/all_samples_exons.txt"
fpkm_path   <- "data/aa_samples_exons.fpkm.txt"

# Sample labels — wildtype vs. TF mutant conditions
SAMPLE_LABELS <- c("Wildtype", "pqm(rhd90)", "pqm(ok485)", "ceh60", "ceh60_pqm(rhd90)", "ceh60_pqm(ok485)")


# -----------------------------------------------------------------------------
# 1. Load and inspect raw counts
# -----------------------------------------------------------------------------

df <- read.table(counts_path, header = TRUE, sep = "\t", row.names = 1)

# Select the 6 sample columns (columns 6–11)
df <- df[, 6:11]

cat("Raw count distribution per sample:\n")
print(apply(df, 2, function(x) quantile(as.numeric(x))))

par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))
apply(df, 2, function(x) hist(log2(x + 1), breaks = 100, main = ""))


# -----------------------------------------------------------------------------
# 2. Filter low-expression genes
# -----------------------------------------------------------------------------
# Retain genes with counts > 20 in at least one sample.
# Removes noise from lowly expressed or unexpressed genes.

expressed_mask <- apply(df, 1, function(x) any(x > 20))
df_expressed   <- df[expressed_mask, ]

cat("\nGenes retained after expression filter:", nrow(df_expressed), "/", nrow(df), "\n")

cat("\nPost-filter count distribution:\n")
print(apply(df_expressed, 2, function(x) quantile(as.numeric(x))))

par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))
apply(df_expressed, 2, function(x) hist(log2(x + 1), breaks = 100, main = ""))


# -----------------------------------------------------------------------------
# 3. DESeq2 — normalization and differential expression
# -----------------------------------------------------------------------------

sample_metadata <- data.frame(
  name      = colnames(df_expressed),
  condition = c("control", rep("treatment", 5))
)

dds <- DESeqDataSetFromMatrix(
  countData = df_expressed,
  colData   = sample_metadata,
  design    = ~ condition
)
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE)
results_all <- results(dds)

cat("\nDESeq2 results summary:\n")
summary(results_all)


# -----------------------------------------------------------------------------
# 4. Extract significant DEGs (padj < 0.01)
# -----------------------------------------------------------------------------

results_sig <- results_all[which(results_all$padj < 0.01), ]
cat("\nSignificant DEGs (padj < 0.01):", nrow(results_sig), "\n")


# -----------------------------------------------------------------------------
# 5. Visualize significant DEGs — boxplot across conditions
# -----------------------------------------------------------------------------

sig_counts           <- norm_counts[rownames(results_sig), ]
colnames(sig_counts) <- SAMPLE_LABELS

boxplot(
  log2(sig_counts + 1),
  outline = FALSE,
  col     = c("purple", "orange", "red", "navy", "blue", "green"),
  border  = c("purple", "orange", "red", "navy", "blue", "green"),
  medcol  = rep("white", 6),
  ylab    = "log2(normalized counts + 1)",
  main    = "Significant DEGs by Condition"
)


# -----------------------------------------------------------------------------
# 6. Pearson correlation heatmap
# -----------------------------------------------------------------------------

fpkm           <- read.table(fpkm_path, header = TRUE, sep = "\t", row.names = 1)
cor_matrix     <- cor(fpkm)
colnames(cor_matrix) <- SAMPLE_LABELS
rownames(cor_matrix) <- SAMPLE_LABELS

# Keep upper triangle only to avoid redundancy
get_upper_tri <- function(mat) {
  mat[lower.tri(mat)] <- NA
  mat
}

melted_cor <- melt(get_upper_tri(cor_matrix), na.rm = TRUE)

ggplot(melted_cor, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low      = "navy",
    high     = "red",
    mid      = "white",
    midpoint = 0.9,
    limits   = c(0.8, 1),
    name     = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12)) +
  coord_fixed() +
  labs(title = "Sample-level Pearson Correlation", x = "", y = "")


# -----------------------------------------------------------------------------
# 7. Hierarchical clustering on Pearson correlation distance
# -----------------------------------------------------------------------------
# Uses 1 - correlation as the distance metric (Ward D2 linkage).
# Samples that cluster together have similar expression profiles.

fpkm_matrix <- as.matrix(fpkm)
cor_dist    <- as.dist(1 - cor(fpkm_matrix))

set.seed(42)
hc <- hclust(cor_dist, method = "ward.D2")

plot(
  hc,
  labels = SAMPLE_LABELS,
  main   = "Hierarchical Clustering — C. elegans TF Mutants",
  xlab   = "",
  sub    = "Distance: 1 - Pearson correlation | Linkage: Ward D2"
)

# =============================================================================
# Differential Gene Expression Analysis: C. elegans TF Mutants — Pairwise
# Author: Rhea Charles | University of South Florida
#
# Description:
#   Pairwise DESeq2 comparisons between wildtype C. elegans and five
#   transcription factor mutant conditions using GEO dataset GSE112978.
#   Each comparison isolates the effect of a single mutant allele against
#   wildtype to identify condition-specific differentially expressed genes.
#
#   Comparisons (all vs. wildtype):
#     1. All mutants combined
#     2. ceh-60(ok1485)
#     3. pqm-1(ok485)
#     4. pqm-1(rhd90)
#     5. ceh-60(ok1485);pqm-1(ok485)
#     6. ceh-60(ok1485);pqm-1(rhd90)
#
# Data:
#   GEO accession: GSE112978
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112978
#
# Input:
#   GSE112978_gene_read_counts.txt      — raw gene-level read counts
#   GSE112978_gene_read_counts_pqm1.txt — pqm-1 allele counts
#
# Output:
#   DESeq2 results per comparison (log2FoldChange, padj)
#   res_combined — merged results table across all six conditions
# =============================================================================

setwd('C:\\Users/riocx/Documents/Masters new/Spring 2023/Applied computational genomics/Applied comp project')

# Update path to your local data directory
input <- "data/GSE112978_gene_read_counts.txt"
dfgene <- read.table(input, header = T, sep = '\t', row.names = 1)

dfgene_exp1     <- dfgene[, c(6,7,8,12,13,14)]
dfgene_wildType <- dfgene[, c(6,7,8)]

library(DESeq2)

# Comparison 1: All mutants vs. wildtype
table2 <- data.frame(name = colnames(dfgene_exp1), condition = c('control','control','control','treatment','treatment','treatment'))
dds_1  <- DESeqDataSetFromMatrix(dfgene_exp1, colData = table2, design = ~ condition)
dds_1  <- DESeq(dds_1)
norCounts_1 <- counts(dds_1, normalized = TRUE)
res_1 <- results(dds_1, alpha = 0.01, contrast = c("condition", "treatment", "control"))

# Comparison 2: Wildtype vs. ceh-60(ok1485)
dfgene_ok1485 <- dfgene[, c(12,13,14)]
dfgene_exp2   <- cbind.data.frame(dfgene_wildType, dfgene_ok1485)

table3 <- data.frame(name = colnames(dfgene_exp2), condition = c('control','control','control','treatment','treatment','treatment'))
dds_2  <- DESeqDataSetFromMatrix(dfgene_exp2, colData = table3, design = ~ condition)
dds_2  <- DESeq(dds_2)
norCounts_2 <- counts(dds_2, normalized = TRUE)
res_2 <- results(dds_2, alpha = 0.01, contrast = c("condition", "treatment", "control"))

# Comparison 3: Wildtype vs. pqm-1(ok485)
input_pqm <- "data/GSE112978_gene_read_counts_pqm1.txt"
dfgene2   <- read.table(input_pqm, header = T, sep = '\t', row.names = 1)

dfgene_ok485 <- dfgene[, c(6,7,8)]
dfgene_exp3  <- cbind.data.frame(dfgene_wildType, dfgene_ok485)

table4 <- data.frame(name = colnames(dfgene_exp3), condition = c('control','control','control','treatment','treatment','treatment'))
dds_3  <- DESeqDataSetFromMatrix(dfgene_exp3, colData = table4, design = ~ condition)
dds_3  <- DESeq(dds_3)
norCounts_3 <- counts(dds_3, normalized = TRUE)
res_3 <- results(dds_3, alpha = 0.01, contrast = c("condition", "treatment", "control"))
head(res_3)

# Comparison 4: Wildtype vs. pqm-1(rhd90)
dfgene_rhd90 <- dfgene[, c(9,10,11)]
dfgene_exp4  <- cbind.data.frame(dfgene_wildType, dfgene_rhd90)

table5 <- data.frame(name = colnames(dfgene_exp4), condition = c('control','control','control','treatment','treatment','treatment'))
dds_4  <- DESeqDataSetFromMatrix(dfgene_exp4, colData = table5, design = ~ condition)
dds_4  <- DESeq(dds_4)
norCounts_4 <- counts(dds_4, normalized = TRUE)
res_4 <- results(dds_4, alpha = 0.01, contrast = c("condition", "treatment", "control"))
head(res_4)

# Comparison 5: Wildtype vs. ceh-60(ok1485);pqm-1(ok485)
dfgene_ok485_ok1485 <- dfgene[, c(12,13,14)]
dfgene_exp5         <- cbind.data.frame(dfgene_wildType, dfgene_ok485_ok1485)

table6 <- data.frame(name = colnames(dfgene_exp5), condition = c('control','control','control','treatment','treatment','treatment'))
dds_5  <- DESeqDataSetFromMatrix(dfgene_exp5, colData = table6, design = ~ condition)
dds_5  <- DESeq(dds_5)
norCounts_5 <- counts(dds_5, normalized = TRUE)
res_5 <- results(dds_5, alpha = 0.01, contrast = c("condition", "treatment", "control"))
head(res_5)

# Comparison 6: Wildtype vs. ceh-60(ok1485);pqm-1(rhd90)
dfgene_rhd90_ok1485 <- dfgene[, c(15,16,17)]
dfgene_exp6         <- cbind.data.frame(dfgene_wildType, dfgene_rhd90_ok1485)

table7 <- data.frame(name = colnames(dfgene_exp6), condition = c('control','control','control','treatment','treatment','treatment'))
dds_6  <- DESeqDataSetFromMatrix(dfgene_exp6, colData = table7, design = ~ condition)
dds_6  <- DESeq(dds_6)
norCounts_6 <- counts(dds_6, normalized = TRUE)
res_6 <- results(dds_6, alpha = 0.01, contrast = c("condition", "treatment", "control"))
head(res_6)

# Combine results across all comparisons
res_combined <- Reduce(function(a, b) merge(a, b, by = "row.names"),
                       list(as.data.frame(res_1), as.data.frame(res_2),
                            as.data.frame(res_3), as.data.frame(res_4),
                            as.data.frame(res_5), as.data.frame(res_6)))

# Extract log2FoldChange and padj from each comparison
res1_df <- as.data.frame(res_1)[, c("log2FoldChange", "padj")]
res2_df <- as.data.frame(res_2)[, c("log2FoldChange", "padj")]
res3_df <- as.data.frame(res_3)[, c("log2FoldChange", "padj")]
res4_df <- as.data.frame(res_4)[, c("log2FoldChange", "padj")]
res5_df <- as.data.frame(res_5)[, c("log2FoldChange", "padj")]
res6_df <- as.data.frame(res_6)[, c("log2FoldChange", "padj")]
