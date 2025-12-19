library(ggplot2)

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