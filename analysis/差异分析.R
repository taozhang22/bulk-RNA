############################################################################################################# 
# DESeq2包差异分析（使用count表达矩阵）
#############################################################################################################
rm(list = ls()); gc()
library(data.table)
library(tidyverse)
library(DESeq2)
setwd("D:\\research\\bioinformatics\\bulk_RNA_analysis\\colorectal_cancer")

exp_data <- fread("D:\\research\\bioinformatics\\database\\TCGA\\COADREAD\\TCGA_COADREAD_RNAseq_count_expr_matrix.csv") |> # 读取count矩阵
    column_to_rownames(var = "V1")
nrow(exp_data)
exp_data <- exp_data[rowSums(exp_data) > 0, ] # 去除在所有的样本中表达值都为0的基因
nrow(exp_data)

colData <- data.frame(condition = word(colnames(exp_data), -1, sep = "-"), row.names = colnames(exp_data)) # 制作分类表
unique(colData$condition)
colData$condition <- factor(ifelse(colData$condition == "11A", "Normal", "Tumor"), levels = c("Normal", "Tumor"))
colData <- arrange(colData, condition)
exp_data <- exp_data[, rownames(colData)]  # 让表达矩阵中Normal在前，Tumor在后

dds <- DESeqDataSetFromMatrix(countData = exp_data,
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05) |> as.data.frame()
