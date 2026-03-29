
############################# Volcano Plot and heatmap #########################


library(ggplot2)
library(ggrepel)
library(pheatmap)
library(scales)
library(readr)
library(dplyr)
library(grid)

plot_from_hgnc_csv <- function(csv_path, mat1, mat2, group1, group2, output_dir) {
  
  result_df <- read_csv(csv_path, show_col_types = FALSE)
  result_df <- as.data.frame(result_df)
  rownames(result_df) <- make.unique(result_df$hgnc_symbol)
  
  result_df$diffexpressed <- "Not significant"
  result_df$diffexpressed[result_df$padj < 0.05 & result_df$log2FoldChange > 1.5]  <- "Upregulated"
  result_df$diffexpressed[result_df$padj < 0.05 & result_df$log2FoldChange < -1.5] <- "Downregulated"
  
  top_up   <- result_df %>% filter(diffexpressed == "Upregulated")   %>% arrange(desc(log2FoldChange)) %>% head(30)
  top_down <- result_df %>% filter(diffexpressed == "Downregulated") %>% arrange(log2FoldChange)       %>% head(15)
  top_genes <- rbind(top_up, top_down)
  
  volcano_plot <- ggplot(result_df, aes(log2FoldChange, -log10(padj), fill = diffexpressed)) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_point(shape = 21, size = 2, alpha = 0.6, color = "black") +
    scale_fill_manual(values = c("Upregulated" = "#bb0c00", "Downregulated" = "#00AFBB", "Not significant" = "grey")) +
    coord_cartesian(xlim = c(-4, 4), ylim = c(0, 5)) +
    labs(title = paste(group1, "vs", group2),
         x = expression(log[2]~FC),
         y = expression(-log[10]~p)) +
    theme_classic(base_size = 20) +
    geom_label_repel(data = top_genes, aes(label = hgnc_symbol), size = 3)
  
  ## TIFF
  tiff(file.path(output_dir, "volcano_plot_top_genes.tiff"),
       width = 10, height = 8, units = "in", res = 300)
  print(volcano_plot)
  dev.off()
  
  ## SVG
  svg(file.path(output_dir, "volcano_plot_top_genes.svg"),
      width = 10, height = 8)
  print(volcano_plot)
  dev.off()
  
  
  ################################ Heatmap #######################################
  
  
  expr_matrix <- cbind(mat1, mat2)
  rownames(expr_matrix) <- make.unique(rownames(expr_matrix))
  
  sig_res <- result_df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1.5)
  
  sig_res <- sig_res[rownames(sig_res) %in% rownames(expr_matrix), ]
  
  if (nrow(sig_res) < 2) {
    message("Not enough significant genes for heatmap: ", group1, " vs ", group2)
    return(NULL)
  }
  
  sig_genes <- rownames(sig_res)
  expr_sig <- expr_matrix[sig_genes, , drop = FALSE]
  expr_sig[is.na(expr_sig)] <- 0
  
  scaled_expr <- t(scale(t(expr_sig)))
  scaled_expr <- scaled_expr[apply(scaled_expr, 1, function(x) all(is.finite(x))), ]
  
  annotation_col <- data.frame(
    Group = c(rep(group1, ncol(mat1)), rep(group2, ncol(mat2)))
  )
  rownames(annotation_col) <- colnames(expr_sig)
  
  ph <- pheatmap(
    scaled_expr,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    fontsize_row = 6,
    main = paste("Heatmap:", group1, "vs", group2),
    silent = TRUE
  )
  
  ## TIFF
  tiff(file.path(output_dir, "heatmap_fc1.5_padj0.05.tiff"),
       width = 8, height = 12, units = "in", res = 300)
  grid.newpage(); grid.draw(ph$gtable)
  dev.off()
  
  ## SVG
  svg(file.path(output_dir, "heatmap_fc1.5_padj0.05.svg"),
      width = 8, height = 12)
  grid.newpage(); grid.draw(ph$gtable)a
  dev.off()
}



######################## Apply Plots to Contrasts ##############################


contrast_plot_list <- list(
  over_weight_vs_normal_weight = c("over_weight", "normal_weight"),
  obese_vs_normal_weight       = c("obese", "normal_weight"),
  obese_vs_over_weight         = c("obese", "over_weight")
)

for (nm in names(contrast_plot_list)) {
  
  group1 <- contrast_plot_list[[nm]][1]
  group2 <- contrast_plot_list[[nm]][2]
  
  output_dir <- file.path("results", paste0(group1, "_vs_", group2))
  csv_path <- file.path(output_dir, "deseq2_results_with_HGNC.csv")
  
  if (!file.exists(csv_path)) next
  
  plot_from_hgnc_csv(
    csv_path,
    group_matrices[[group1]],
    group_matrices[[group2]],
    group1,
    group2,
    output_dir
  )
}
