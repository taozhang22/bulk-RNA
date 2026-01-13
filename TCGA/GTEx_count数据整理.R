rm(list=ls()); gc()
library(data.table)
library(tidyverse)
setwd("D:\\research\\bioinformatics\\database\\GTEx")

count <- fread("GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_reads.gct") |>
    distinct(Description, .keep_all = TRUE) |>
    column_to_rownames("Description") |>
    dplyr::select(-Name) |>
    setNames(word(colnames(.), start = 1, end = 3, sep = "-"))
