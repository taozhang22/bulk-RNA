rm(list = ls()); gc()
library(data.table)
library(tidyverse)
library(DESeq2)
setwd("D:\\research\\bioinformatics\\bulk_RNA_analysis\\colorectal_cancer")

exp_data <- fread("D:\\research\\bioinformatics\\database\\TCGA\\COADREAD\\TCGA_COADREAD_RNAseq_count_expr_matrix.csv") |>
    column_to_rownames(var = "V1")
exp_data <- exp_data[rowSums(exp_data) > 0, ]

colData <- data.frame(condition = word(colnames(exp_data), -1, sep = "-"), row.names = colnames(exp_data))
unique(colData$condition)
colData$condition <- factor(ifelse(colData$condition == "11A", "Normal", "Tumor"), levels = c("Normal", "Tumor"))
dds <- DESeqDataSetFromMatrix(countData = exp_data,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05) |>
    as.data.frame()
filter(res, rownames(res) == "GPX4")
