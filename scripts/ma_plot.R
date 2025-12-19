library(ggplot2)

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