
###################### GSEA_Custon immport signture ###############
#Gene signature were downloaded from https://www.dev.immport.org/shared/genelists accessed on 14-04-2025, 12:41 AM, IST
library(GSEABase)
library(GSVA)
library(DESeq2)
library(pheatmap)
library(matrixStats)
library(org.Hs.eg.db)
library(clusterProfiler)
# Load custom ImmPort gene set from GMT
gmt_file <- "immport_signature/immport_signature.gmt"
gene_sets <- getGmt(gmt_file)

# Convert to named list
gene_set_list <- geneIds(gene_sets)
names(gene_set_list) <- names(gene_sets)


run_custom_gmt_ssgsea <- function(dds_object, gmt_path, prefix = "ImmPort", out_dir = "results/ssGSEA", top_n = 20) {

  
  # Output path
  out_path <- file.path(out_dir, prefix)
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  
  # VST normalization
  vsd_mat <- assay(vst(dds_object, blind = FALSE))
  
  # Ensembl to HGNC conversion
  ensembl_ids <- rownames(vsd_mat)
  mapping <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
  common <- intersect(mapping$ENSEMBL, rownames(vsd_mat))
  vsd_mat <- vsd_mat[common, ]
  rownames(vsd_mat) <- mapping$SYMBOL[match(common, mapping$ENSEMBL)]
  vsd_mat <- vsd_mat[!duplicated(rownames(vsd_mat)), ]
  
  # Load GMT file
  gene_sets <- getGmt(gmt_path)
  gene_set_list <- geneIds(gene_sets)
  names(gene_set_list) <- names(gene_sets)
  
  # Run ssGSEA
  message("Running ssGSEA with ImmPort GMT...")
  ssgsea_param <- ssgseaParam(expr = vsd_mat, geneSets = gene_set_list)
  ssgsea_scores <- gsva(param = ssgsea_param)
  
  # Save scores
  write.csv(ssgsea_scores, file = file.path(out_path, paste0("ssGSEA_scores_", prefix, ".csv")))
  
  # Top variable pathways
  top_paths <- head(order(rowVars(ssgsea_scores), decreasing = TRUE), top_n)
  top_matrix <- ssgsea_scores[top_paths, , drop = FALSE]
  
  # Heatmap
  tiff(file.path(out_path, paste0("ssGSEA_top", top_n, "_heatmap_", prefix, ".tiff")),
       width = 10, height = 10, units = "in", res = 300)
  pheatmap(
    top_matrix,
    scale = "row",
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 7,
    main = paste("ssGSEA Top", top_n, "Pathways (ImmPort)")
  )
  dev.off()
  
  message("ssGSEA with custom GMT completed: ", prefix)
}



run_custom_gmt_ssgsea(
  dds_object = dds_all,
  gmt_path = "immport_signature/immport_signature.gmt",
  prefix = "ImmPort_custom",
  out_dir = "results/ssGSEA_custom/3group_bmi"
)


#################### Kruskal-Wallis or ANOVA (plus post-hoc tests) #############
library(dplyr)
library(ggpubr)
library(FSA)      # For dunnTest
library(rstatix)  # For kruskal_test

run_kruskal_posthoc <- function(ssgsea_matrix, bmi_groups, output_csv) {
  # Ensure sample names match
  names(bmi_groups) <- gsub("[-;]", ".", names(bmi_groups))
  colnames(ssgsea_matrix) <- gsub("[-;]", ".", colnames(ssgsea_matrix))
  
  # Filter to common samples
  common_samples <- intersect(colnames(ssgsea_matrix), names(bmi_groups))
  bmi_groups <- bmi_groups[common_samples]
  ssgsea_matrix <- ssgsea_matrix[, common_samples]
  
  # Store results
  stat_results <- data.frame()
  
  for (pathway in rownames(ssgsea_matrix)) {
    scores <- ssgsea_matrix[pathway, ]
    group_df <- data.frame(score = as.numeric(scores), group = as.factor(bmi_groups[names(scores)]))
    
    # Run Kruskal-Wallis if multiple groups exist
    if (nlevels(group_df$group) >= 2 && all(table(group_df$group) >= 2)) {
      kruskal <- kruskal_test(score ~ group, data = group_df)
      
      if (kruskal$p < 0.05) {
        # Run posthoc Dunn's test
        posthoc <- tryCatch({
          dunnTest(score ~ group, data = group_df, method = "bh")
        }, error = function(e) {
          message("Dunn's test failed for pathway: ", pathway)
          return(NULL)
        })
        
        # Process and store results
        if (!is.null(posthoc) && "res" %in% names(posthoc)) {
          posthoc_result <- posthoc$res %>%
            filter(P.adj < 0.05) %>%
            mutate(Signature = pathway)
          
          stat_results <- bind_rows(stat_results, posthoc_result)
        }
      }
    }
  }
  
  # Save result
  write.csv(stat_results, output_csv, row.names = FALSE)
  message("Post-hoc comparison results saved to: ", output_csv)
}


ssgsea_matrix_3group <- read.csv("results/ssGSEA_custom/3group_bmi/ImmPort_custom/ssGSEA_scores_ImmPort_custom.csv",
                                 row.names = 1, check.names = FALSE)

bmi_groups_3group <- colData(dds_all)$BMI_group
names(bmi_groups_3group) <- colnames(dds_all)


# Run for 3-group comparison
run_kruskal_posthoc(
  ssgsea_matrix = ssgsea_matrix_3group,
  bmi_groups = bmi_groups_3group,
  output_csv = "results/ssGSEA_custom/3group_bmi/statistical_comparisons_posthoc.csv"
)


########################## violin plot combined ################################
library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)

# Load posthoc statistical results (e.g., Dunn test output)
posthoc_df <- read.csv("results/ssGSEA_custom/3group_bmi/statistical_comparisons_posthoc.csv")

# Load ssGSEA matrix and group labels
ssgsea_matrix_3group <- read.csv("results/ssGSEA_custom/3group_bmi/ImmPort_custom/ssGSEA_scores_ImmPort_custom.csv", row.names = 1, check.names = FALSE)
bmi_groups_3group <- colData(dds_all)$BMI_group
names(bmi_groups_3group) <- colnames(dds_all)

# Normalize IDs
normalize_ids <- function(x) gsub("(C3[LN])[.-](\\d{5}).*", "\\1-\\2", x)
colnames(ssgsea_matrix_3group) <- normalize_ids(colnames(ssgsea_matrix_3group))
names(bmi_groups_3group) <- normalize_ids(names(bmi_groups_3group))

# Match samples
common_ids <- intersect(colnames(ssgsea_matrix_3group), names(bmi_groups_3group))
ssgsea_filtered <- ssgsea_matrix_3group[, common_ids]
bmi_groups_filtered <- bmi_groups_3group[common_ids]

# Output directory
plot_dir <- "results/ssGSEA_custom/3group_bmi/final_violin_plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

plot_list <- list()

for (i in 1:nrow(posthoc_df)) {
  sig <- posthoc_df$Signature[i]
  comparison <- posthoc_df$Comparison[i]
  pval <- round(posthoc_df$P.adj[i], 4)
  
  # Split group comparison
  groups <- strsplit(comparison, " - ")[[1]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  group_samples <- names(bmi_groups_filtered)[bmi_groups_filtered %in% c(group1, group2)]
  scores <- as.numeric(ssgsea_filtered[sig, group_samples])
  
  plot_df <- data.frame(
    Score = scores,
    Group = factor(bmi_groups_filtered[group_samples], levels = c(group1, group2))
  )
  
  # Format long signature names
  display_title <- stringr::str_wrap(sig, width = 35)
  
  # Plot
  p <- ggplot(plot_df, aes(x = Group, y = Score, fill = Group)) +
    geom_violin(trim = FALSE, color = "black", alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.7) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c(group1, group2)),
      label = "p.signif",
      label.y.npc = "top",
      size = 12,
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
        symbols = c("****", "***", "**", "*", "ns")
      )
    ) +
    labs(
      title = stringr::str_wrap(display_title, width = 35),
      subtitle = paste(group1, "vs", group2, "| adj.p =", pval),
      x = NULL, y = "ssGSEA Score"
    ) +
    scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
    theme_pubr(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5, size = 18),
      legend.position = "none",
      panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      panel.grid.major.y = element_line(linetype = "dashed", color = "#b4aea9"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "grey40"),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 20, face = "bold"),
      axis.text = element_text(size = 20)
    )
  
  # Save individual TIFF
  filename <- paste0(plot_dir, "/", gsub("[/: ]", "_", sig), "_", group1, "_vs_", group2, ".tiff")
  ggsave(filename, plot = p, width = 7, height = 6.5, dpi = 300)
  
  # Store for combined
  plot_list[[i]] <- p
}


# Save combined plot
tiff("results/ssGSEA_custom/3group_bmi/final_violin_plots/combined_ssGSEA_violin_grid.tiff", 
     width = 18, height = 12, units = "in", res = 300)
print(plot_grid(plotlist = plot_list, ncol = 3, align = "hv"))
dev.off()
