# 将TCGA的临床数据进行标准化
acquire_tcga_clinical_data <- function(outfile) {
  filenames <- list.files(pattern = "*xml$", recursive = TRUE)
  data_list <- list()
  for (filename in filenames) {
    data <- xmlParse(filename)
    data <- xmlRoot(data)
    data <- xmlToDataFrame(data[2])
    data_list[[filename]] <- data
  }
  data <- do.call(rbind, data_list) |>
    mutate(survival_time = ifelse(vital_status == "Alive", days_to_last_followup, days_to_death),
    stage_event = str_squish(stage_event)) |>
    extract(stage_event,
            into = c("stage", "T", "N", "M"),
            regex = "^(?:.*Stage)\\s+(IV|III|II|I)(?:[ABC]?)T([0-9]|X)(?:[a-z]?)N([0-9]|X)(?:[a-z]?)M([0-9]|X)(?:[a-z]?)$",
            remove = FALSE)
  fwrite(data, file = outfile, row.names = FALSE)
}
