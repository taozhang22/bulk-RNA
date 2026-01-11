# 将TCGA的RNA seq表达数据进行标准化
acquire_tcga_rnaseq_expr_matrix_data <- function(dtype = "unstranded", outfile) {
  sheet <- fread(list.files("expr_data", pattern = "sheet", recursive = TRUE, full.names = TRUE))
  filenames <- list.files(pattern = "augmented_star_gene_counts", recursive = TRUE, full.names = TRUE)
  data_list <- list()
  for (filename in filenames) {
    data <- fread(filename) |>
      filter(!gene_id %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")) |>
      distinct(gene_name, .keep_all = TRUE) |>
      column_to_rownames(var = "gene_name") |>
      dplyr::select(all_of(dtype))
    colnames(data) <- sheet$`Sample ID`[match(basename(filename), sheet$`File Name`)]
    data_list[[filename]] <- data
  }
  data <- do.call(cbind, data_list)
  fwrite(data, file = outfile, row.names = TRUE)
}
