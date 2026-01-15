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

############################################################################################################# 
# 绘制箱线图
#############################################################################################################
gene = "GPX4"
data <- merge(x = as.data.frame(t(exp_data[gene, ])),
              y = colData,
              by = 0) |>
    group_by(condition) |>
    mutate(q1 = quantile(GPX4, 0.25, na.rm = TRUE),
           q3 = quantile(GPX4, 0.75, na.rm = TRUE),
           IQR = q3 - q1,
           lower = q1 - 1.5 * IQR,
           upper = q3 + 1.5 * IQR) |>
    ungroup()
# 获取差异分析的结果
p <- res[gene, "padj"]
star <- case_when(
  is.na(p)    ~ "NA",
  p < 1e-4    ~ "****",
  p < 1e-3    ~ "***",
  p < 1e-2    ~ "**",
  p < 5e-2    ~ "*",
  TRUE        ~ "ns"
)


p <- ggplot(data = data, aes(x = condition, y = GPX4)) +
    geom_boxplot(aes(fill = condition), width = 0.5, outlier.shape = NA) +
    geom_jitter(data = subset(data, GPX4 >= lower & GPX4 <= upper), width = 0.2, size = 0.5, alpha = 0.5) +
    annotate("segment", x = 1, xend = 2, y = 15000, yend = 15000, linewidth = 0.6) +
    annotate("text", x = 1.5, y = 15500, label = star, size = 6) +
    labs(x = "", y = "GPX4") +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", color = NA),
          axis.line = element_line(color = "black")) +
    scale_fill_manual(values = c("Normal" = "#1f77b4", "Tumor" = "#ff7f0e")) +
    coord_cartesian(ylim = c(0, 15000))

    
  


ggplot(data = data, aes(x = condition, y = GPX4)) +
  geom_boxplot(aes(fill = condition), width = 0.5, outlier.shape = NA) +
  geom_jitter(data = data_no_outliers, width = 0.2, size = 0.5, alpha = 0.5) +
  labs(x = "Condition", y = "GPX4") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA) # 背景为白色
        )) +
  scale_fill_manual(values = c("Normal" = "#1f77b4", "Tumor" = "#ff7f0e"))


