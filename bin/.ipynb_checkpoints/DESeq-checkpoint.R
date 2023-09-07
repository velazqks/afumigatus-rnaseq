# http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual#de

# Install required packages
#install.packages("BiocManager")
#library(BiocManager)
#install.packages("devtools")
#devtools::install_github('alyssafrazee/RSkittleBrewer')
packages_to_install <- c("ballgown","genefilter","dplyr","DESeq2", "ggplot2")
for (package in packages_to_install) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
  } else {
    message(paste(package, "is already installed."))
  }
}

# Import libraries
library(DESeq2)
library(ggplot2)

# Load data and set index column
countData <- as.matrix(read.csv("../data/transcript_count_matrix.csv", row.names="transcript_id"))
colData <- read.csv("../data/phenodata.csv", row.names="id")

# Check that all IDs in colData are in CountData in the same order
all(rownames(colData) %in% colnames(countData))

# Keep only the samples names matching the sample names from colData.
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

# Create a DESeqDataSet from count matrix and labels
# The ~ 1 formula specifies the overall mean (intercept) in all samples without independent variables. 
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ 1)

# Run the default analysis for DESeq and create a table
dds <- DESeq(dds)
res <- results(dds)

# Sort by adjusted p-value and display
(resOrdered <- res[order(res$padj), ])

# Convert DESeqResults object to a data frame
res_df <- data.frame(
  transcript = row.names(resOrdered),
  log2fc = resOrdered$log2FoldChange,
  pval = resOrdered$pvalue,
  padj = resOrdered$padj
)

# DFR adjustment
res_df$adjusted_pval <- p.adjust(res_df$pval, method = "BH")

# Specify if each transcript is UP- or DOWN- regulated
res_df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
res_df$diffexpressed[res_df$log2fc > 0.6 & res_df$pval < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_df$diffexpressed[res_df$log2fc < -0.6 & res_df$pval < 0.05] <- "DOWN"
# Explore a bit
head(res_df[order(res_df$padj) & res_df$diffexpressed == 'DOWN', ])


# Volcano plot
ggplot(data = res_df, aes(x = log2fc, y = -log10(adjusted_pval), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("NO" = "grey", "UP" = "#FFDB6D", "DOWN" = "#00AFBB"),
                     labels = c("NO", "UP", "DOWN"))



# -------------------------------------------------------------------



