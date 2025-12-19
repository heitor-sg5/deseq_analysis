library(here)

source("scripts/parser.R")
source("scripts/run_deseq2.R")
source("scripts/summaries.R")
source("scripts/volcano_plot.R")
source("scripts/pca_plot.R")
source("scripts/ma_plot.R")

input_file <- here("data", "example_data.xlsx")
fig_dir <- "output/figures"
tab_dir <- "output/tables"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

data <- read_rnaseq_xlsx(input_file)

count_matrix <- data$count_matrix
sample_metadata <- data$sample_metadata

deseq <- run_deseq2(count_matrix, sample_metadata)

dds <- deseq$dds
results <- deseq$results

write.csv(
  as.data.frame(results),
  file = file.path(tab_dir, "deseq2_results.csv")
)

run_summaries(results, count_matrix)

plot_volcano(results, file.path(fig_dir, "volcano_plot.pdf"))

plot_pca(count_matrix, sample_metadata, file.path(fig_dir, "pca_plot.pdf"))

plot_ma(results,file.path(fig_dir, "ma_plot.pdf"))