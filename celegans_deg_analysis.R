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
