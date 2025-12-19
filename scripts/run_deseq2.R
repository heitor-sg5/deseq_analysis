library(DESeq2)

run_deseq2 <- function(count_matrix, sample_metadata) {
  rownames(sample_metadata) <- sample_metadata$sample
  stopifnot(all(colnames(count_matrix) == sample_metadata$sample))
  
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_metadata,
    design = ~ condition
  )
  
  dds <- DESeq(dds)
  
  results <- results(
    dds,
    contrast = c("condition", "Mutant", "WT")
  )
  
  results$padj <- results$padj + 1e-10
  
  list(dds = dds, results = results)
}