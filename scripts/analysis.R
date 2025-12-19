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

run_summaries <- function(results, count_matrix) {
  cat("\n===== SUMMARY =====\n")
  
  sig <- results[!is.na(results$padj) & results$padj < 0.05, ]
  up <- sig[sig$log2FoldChange > 0, ]
  down <- sig[sig$log2FoldChange < 0, ]
  
  cat("Significant genes (padj < 0.05):", nrow(sig), "\n")
  cat("Upregulated:", nrow(up), "\n")
  cat("Downregulated:", nrow(down), "\n")
  cat("Up/Down ratio:", round(nrow(up) / nrow(down), 2), "\n\n")
  
  strongest_up <- rownames(up)[which.max(up$log2FoldChange)]
  strongest_down <- rownames(down)[which.min(down$log2FoldChange)]
  
  cat("Strongest upregulated gene:", strongest_up, "\n")
  cat("Strongest downregulated gene:", strongest_down, "\n")
  
  counts_up <- count_matrix[strongest_up, ]
  cat(
    "Sample with the highest expression:",
    names(which.max(counts_up)), "\n"
  )
  
  cat("\nMean expression per sample:\n")
  print(colMeans(count_matrix))
  
  vars <- apply(count_matrix, 2, var)
  cat("\nSample with highest variance:", names(which.max(vars)), "\n")
  
  gene_means <- rowMeans(count_matrix)
  
  cat(
    "\nHighest expressed gene overall:",
    names(which.max(gene_means)), "\n"
  )
  cat(
    "Lowest expressed gene overall:",
    names(which.min(gene_means)), "\n"
  )
  
  filtered <- down[down$baseMean > 100, ]
  if (nrow(filtered) > 0) {
    cat(
      "\nLowest log2FC (baseMean > 100):",
      round(min(filtered$log2FoldChange), 2), "\n"
    )
  }
}