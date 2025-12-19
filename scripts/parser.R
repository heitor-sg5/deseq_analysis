library(openxlsx)

read_rnaseq_xlsx <- function(path) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path)
  }
  
  sheets <- getSheetNames(path)
  required_sheets <- c("counts", "metadata")
  
  missing_sheets <- setdiff(required_sheets, sheets)
  if (length(missing_sheets) > 0) {
    stop("Missing required sheet(s): ", paste(missing_sheets, collapse = ", "))
  }
  
  counts <- as.data.frame(
    read.xlsx(path, sheet = "counts", colNames = TRUE)
  )
  
  meta <- as.data.frame(
    read.xlsx(path, sheet = "metadata", colNames = TRUE)
  )
  
  if (!"gene" %in% colnames(counts)) {
    stop("Counts sheet must contain a 'gene' column")
  }
  
  if (ncol(counts) < 2) {
    stop("Counts sheet must contain at least one sample column")
  }
  
  if (anyDuplicated(counts$gene)) {
    stop("Duplicate gene identifiers found in counts sheet")
  }
  
  rownames(counts) <- counts$gene
  counts$gene <- NULL
  
  count_matrix <- data.matrix(counts)
  
  if (anyNA(count_matrix)) {
    stop("Counts matrix contains non-numeric or missing values")
  }
  
  if (any(count_matrix < 0)) {
    stop("Counts matrix contains negative values")
  }
  
  required_meta_cols <- c("sample", "condition")
  missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))
  
  if (length(missing_meta_cols) > 0) {
    stop(
      "Metadata sheet missing column(s): ",
      paste(missing_meta_cols, collapse = ", ")
    )
  }
  
  if (anyDuplicated(meta$sample)) {
    stop("Duplicate sample names found in metadata")
  }
  
  if (!setequal(colnames(count_matrix), meta$sample)) {
    stop(
      "Sample names in counts and metadata do not match\n",
      "Counts: ", paste(colnames(count_matrix), collapse = ", "), "\n",
      "Metadata: ", paste(meta$sample, collapse = ", ")
    )
  }
  
  meta <- meta[match(colnames(count_matrix), meta$sample), ]
  
  list(count_matrix = count_matrix, sample_metadata = meta)
}