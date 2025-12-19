library(ggplot2)

plot_pca <- function(count_matrix, sample_metadata, outfile) {
  gene_var <- apply(count_matrix, 1, var)
  count_matrix <- count_matrix[gene_var > 0, ]
  
  mat <- t(log2(count_matrix + 1))
  pca <- prcomp(mat, scale. = TRUE)
  
  var_exp <- summary(pca)$importance[2, ] * 100
  
  df <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    sample = rownames(pca$x)
  )
  
  df$condition <- sample_metadata$condition[
    match(df$sample, sample_metadata$sample)
  ]
  
  p <- ggplot(df, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2])
    ) +
    theme_minimal()
  
  ggsave(outfile, p, width = 8, height = 6)
}

plot_ma <- function(results, outfile) {
  df <- as.data.frame(results)
  
  df$A <- log2(df$baseMean)
  df$M <- df$log2FoldChange
  df$significant <- !is.na(df$padj) & df$padj < 0.05
  
  p <- ggplot(df, aes(A, M, color = significant)) +
    geom_point(alpha = 0.5, size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(
      x = "mean expression (log2)",
      y = "log2 fold change"
    )
  
  ggsave(outfile, p, width = 8, height = 6)
}

plot_volcano <- function(results, outfile) {
  
  df <- as.data.frame(results)
  df <- df[!is.na(df$padj), ]
  
  df$threshold <- ifelse(
    df$padj < 0.05 & df$log2FoldChange > 0, "Upregulated",
    ifelse(
      df$padj < 0.05 & df$log2FoldChange < 0,
      "Downregulated",
      "Not Significant"
    )
  )
  p <- ggplot(df, aes(log2FoldChange, -log10(padj), color = threshold)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(
      values = c(
        "Upregulated" = "darkred",
        "Downregulated" = "darkblue",
        "Not Significant" = "darkgrey"
      )
    ) +
    theme_minimal() +
    labs(
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = "Status"
    )
  ggsave(outfile, p, width = 7, height = 6)
}