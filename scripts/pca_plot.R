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