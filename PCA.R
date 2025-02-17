library(tidyverse)
library(tximport)
library(rtracklayer)
library(DESeq2)
library(GenomicFeatures)
library(ggplot2)

samples <- c(
  "SRR892995", "SRR892996", "SRR892997", "SRR892998", "SRR892999", "SRR893000", 
  "SRR893001", "SRR893002")

names(samples) <- c('young', 'young', 'young', 'young', 'old', 'old', 'old', 'old')

sampleTable <- data.frame(
  sampleName = samples,
  fileName = quant_filepaths,
  condition = rep(c("young", "old"), each = 4) 
)

sampleTable$condition <- factor(sampleTable$condition, levels = c("young", "old"))

# load transcript-gene csv
tx2gene1 <- read.csv("/cellfile/datapublic/ahyseni/cryption/m.musculus/data/reference_omes/tx2gene.csv")

tx2gene1$TXNAME <- sub("\\..*", "", tx2gene1$TXNAME)

# import Salmon data and summarize to genes
txi <- tximport(quant_filepaths, type = "salmon", tx2gene = tx2gene1, ignoreTxVersion = TRUE)

# Create DESeq2 dataset
dds <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~condition)

# Remove genes with very low counts
dds <- dds[rowSums(counts(dds)) > 5, ]

# Run DESeq normalization
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)


pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

