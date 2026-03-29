
####################### GSEA Analysis Using DESeq2 Stat Ranking ################

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

run_gsea_and_save <- function(res_df, out_dir, comparison_label) {
  
  # Prepare ranked gene list
  res_df <- res_df[!is.na(res_df$stat), ]
  res_df <- res_df[!duplicated(res_df$hgnc_symbol), ]
  rownames(res_df) <- res_df$hgnc_symbol
  
  gene_list <- res_df$stat
  names(gene_list) <- rownames(res_df)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Run GSEA
  gsea_obj <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.1,
    by = "fgsea",
    verbose = FALSE
  )
  
  # Save outputs
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(as.data.frame(gsea_obj), file = file.path(out_dir, "GO_GSEA_results.csv"), row.names = FALSE)
  saveRDS(gene_list, file = file.path(out_dir, "gene_list.rds"))
  saveRDS(gsea_obj, file = file.path(out_dir, "GO_GSEA_object.rds"))
  
  # Dotplot
  tiff(file.path(out_dir, "GO_dotplot.tiff"), width = 7, height = 10, units = "in", res = 300)
  print(dotplot(gsea_obj, showCategory = 20) +
          ggtitle(paste("GSEA-GO:", comparison_label)) +
          theme(axis.text.y = element_text(size = 10, hjust = 1)))
  dev.off()
  
  svg(file.path(out_dir, "GO_dotplot.svg"), width = 7, height = 10)
  print(dotplot(gsea_obj, showCategory = 20) +
          ggtitle(paste("GSEA-GO:", comparison_label)) +
          theme(axis.text.y = element_text(size = 10, hjust = 1)))
  dev.off()
}


# Load DESeq2 results with HGNC
res_ow_nw <- read.csv("/results/over_weight_vs_normal_weight/deseq2_results_with_HGNC.csv")
res_ob_nw <- read.csv("/results/obese_vs_normal_weight/deseq2_results_with_HGNC.csv")
res_ob_ow <- read.csv("/results/obese_vs_over_weight/deseq2_results_with_HGNC.csv")

# Run and save GSEA for each comparison
run_gsea_and_save(res_bmi, "results/GSEA_sets/high_bmi_vs_normal_bmi/GO", "high_bmi_vs_normal_bmi")
run_gsea_and_save(res_ow_nw, "results/GSEA_sets/over_weight_vs_normal_weight/GO", "over_weight_vs_normal_weight")
run_gsea_and_save(res_ob_nw, "results/GSEA_sets/obese_vs_normal_weight/GO", "obese_vs_normal_weight")
run_gsea_and_save(res_ob_ow, "results/GSEA_sets/obese_vs_over_weight/GO", "obese_vs_over_weight")






####################### GSEA (Entrez id input) + Plot function #################

# ======================== Entrez id mapping ================================= #
prepare_entrez_from_ensembl <- function(res_df) {
  library(org.Hs.eg.db)
  library(clusterProfiler)
  
  res_df <- res_df[!is.na(res_df$stat), ]
  res_df <- res_df[!duplicated(res_df$ensembl_id), ]
  
  entrez_df <- bitr(res_df$ensembl_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  merged_df <- merge(res_df, entrez_df, by.x = "ensembl_id", by.y = "ENSEMBL")
  merged_df <- merged_df[!duplicated(merged_df$ENTREZID), ]
  
  gene_list <- merged_df$stat
  names(gene_list) <- merged_df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

# ======================== GSEA Entrez id set enrichment ===================== #
run_gsea_entrez_sets <- function(gene_list, comparison_name,
                                 output_base = "/results/GSEA_sets") {
  library(clusterProfiler)
  library(ReactomePA)
  library(msigdbr)
  library(enrichplot)
  
  make_plot <- function(gsea_res, name, out_dir) {
    if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
      write.csv(as.data.frame(gsea_res),
                file = file.path(out_dir, paste0(name, "_GSEA_results.csv")))
      
      tiff(file.path(out_dir, paste0(name, "_dotplot.tiff")),
           width = 7, height = 10, units = "in", res = 300)
      print(dotplot(gsea_res, showCategory = 20) + ggtitle(name))
      dev.off()
      
      svg(file.path(out_dir, paste0(name, "_dotplot.svg")),
          width = 7, height = 10)
      print(dotplot(gsea_res, showCategory = 20) + ggtitle(name))
      dev.off()
    }
  }
  
  categories <- list(
    KEGG = list(fun = gseKEGG, args = list(organism = "hsa", keyType = "kegg")),
    Reactome = list(fun = gsePathway, args = list(organism = "human")),
    Hallmark = list(msig_category = "H"),
    C4 = list(msig_category = "C4"),
    C5 = list(msig_category = "C5"),
    C6 = list(msig_category = "C6"),
    C7 = list(msig_category = "C7"),
    WikiPathways = list(msig_category = "C2", subcategory = "CP:WIKIPATHWAYS")
  )
  
  for (cat in names(categories)) {
    out_dir <- file.path(output_base, comparison_name, cat)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (!is.null(categories[[cat]]$fun)) {
      res <- tryCatch({
        do.call(categories[[cat]]$fun,
                c(list(geneList = gene_list,
                       pvalueCutoff = 0.05,
                       minGSSize = 10,
                       maxGSSize = 500),
                  categories[[cat]]$args))
      }, error = function(e) NULL)
      
      make_plot(res, cat, out_dir)
      
    } else {
      msig <- msigdbr(
        species = "Homo sapiens",
        category = categories[[cat]]$msig_category,
        subcategory = categories[[cat]]$subcategory
      )
      
      gene_sets <- split(msig$entrez_gene, msig$gs_name)
      subset_list <- gene_list[names(gene_list) %in% msig$entrez_gene]
      
      res <- tryCatch({
        GSEA(subset_list, TERM2GENE = gene_sets, pvalueCutoff = 0.05)
      }, error = function(e) NULL)
      
      make_plot(res, cat, out_dir)
    }
  }
}

# Run GSEA for each comparison
# 1. Prepare gene list from res_3a and run GSEA for over_weight_vs_normal_weight
gene_list_ov_nw <- prepare_entrez_from_ensembl(res_3a)
run_gsea_entrez_sets(gene_list_ov_nw, "over_weight_vs_normal_weight")

# 2. Prepare gene list from res_3b and run GSEA for obese_vs_normal_weight
gene_list_ob_nw <- prepare_entrez_from_ensembl(res_3b)
run_gsea_entrez_sets(gene_list_ob_nw, "obese_vs_normal_weight")

# 3. Prepare gene list from res_3c and run GSEA for obese_vs_over_weight
gene_list_ob_ow <- prepare_entrez_from_ensembl(res_3c)
run_gsea_entrez_sets(gene_list_ob_ow, "obese_vs_over_weight")


####################### GSEA (Entrez id input) + Plot function #################
# ======================== entrez id mapping ================================= #
prepare_entrez_from_ensembl <- function(res_df) {
  library(org.Hs.eg.db)
  library(clusterProfiler)
  
  res_df <- res_df[!is.na(res_df$stat), ]
  res_df <- res_df[!duplicated(res_df$ensembl_id), ]
  
  entrez_df <- bitr(res_df$ensembl_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  merged_df <- merge(res_df, entrez_df, by.x = "ensembl_id", by.y = "ENSEMBL")
  merged_df <- merged_df[!duplicated(merged_df$ENTREZID), ]
  
  gene_list <- merged_df$stat
  names(gene_list) <- merged_df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}

# ====================== Bubble Plot Function (Facet + Mirror) =============== #
plot_gsea_bubble_facet_mirror <- function(comparison_name, gsea_type = c("GO", "Reactome", "Hallmark", "WikiPathways"), base_path = "C:/Users/arunv/OneDrive/Desktop/CPTAC analysis 06-04-2025/results/GSEA_sets") {
  library(ggplot2)
  library(dplyr)
  library(stringr)
  
  for (etype in gsea_type) {
    file_path <- file.path(base_path, comparison_name, etype, paste0(etype, "_GSEA_results.csv"))
    if (!file.exists(file_path)) {
      message("Skipping: ", etype, " (file not found)")
      next
    }
    
    gsea_df <- read.csv(file_path)
    if (!"NES" %in% colnames(gsea_df)) {
      message("Skipping: ", etype, " (missing NES column)")
      next
    }
    
    gsea_df$direction <- ifelse(gsea_df$NES > 0, "Activated", "Suppressed")
    gsea_df$GeneCount <- sapply(strsplit(as.character(gsea_df$core_enrichment), "/"), length)
    gsea_df$GeneRatio <- gsea_df$GeneCount / gsea_df$setSize
    gsea_df$GeneRatio_mirror <- ifelse(gsea_df$direction == "Suppressed", -gsea_df$GeneRatio, gsea_df$GeneRatio)
    
    top_activated <- head(gsea_df[gsea_df$direction == "Activated", ], 15)
    top_suppressed <- head(gsea_df[gsea_df$direction == "Suppressed", ], 15)
    top_combined <- rbind(top_activated, top_suppressed)
    top_combined$Description <- stringr::str_wrap(top_combined$Description, width = 30)
    
    out_dir <- file.path(base_path, comparison_name, "BubblePlots")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Facet Plot
    p_facet <- ggplot(top_combined, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
      geom_point(aes(size = GeneCount, fill = -log10(p.adjust)), color = "black", shape = 21, stroke = 0.5) +
      scale_fill_gradient(low = "#4575b4", high = "#d73027", name = expression(-log[10]~"(adj. p-value)")) +
      scale_size(range = c(3, 10), name = "Gene Count") +
      facet_wrap(~direction, scales = "free_x", nrow = 1, strip.position = "top") +
      theme_minimal(base_size = 14) +
      labs(title = paste("GSEA:", etype, "\n", comparison_name), x = "GeneRatio", y = NULL) +
      theme(
        strip.background = element_rect(fill = "grey90", color = "grey20"),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.spacing.x = unit(1.5, "lines"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)
      )
    
    ggsave(file.path(out_dir, paste0("Bubble_Facetcomp_", etype, ".tiff")), plot = p_facet, width = 10, height = 20, dpi = 300)
    
    # Mirror Plot
    p_mirror <- ggplot(top_combined, aes(x = GeneRatio_mirror, y = reorder(Description, GeneRatio))) +
      geom_point(aes(size = GeneCount, fill = -log10(p.adjust)), color = "black", shape = 21, stroke = 0.5) +
      scale_fill_gradient(low = "#4575b4", high = "#d73027", name = expression(-log[10]~"(adj. p-value)")) +
      scale_size(range = c(3, 10), name = "Gene Count") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
      theme_minimal(base_size = 14) +
      labs(title = paste("GSEA:", etype, "\n", comparison_name), 
           x = "\u2190 Suppressed    GeneRatio    Activated \u2192", y = NULL) +
      theme(
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.title = element_text(face = "bold")
      )
    
    ggsave(file.path(out_dir, paste0("Bubble_Mirrorcomp_", etype, ".tiff")), plot = p_mirror, width = 12, height = 20, dpi = 300)
  }
}

# ================ run bubble plot function for all comparisons ============== #
comparisons <- c("over_weight_vs_normal_weight", "obese_vs_normal_weight", 
                 "obese_vs_over_weight")

for (comp in comparisons) {
  plot_gsea_bubble_facet_mirror(comp)
}

# ================= Bubble Plot Function (Facet + Mirror) ==================== #
plot_gsea_bubble_facet_mirror <- function(comparison_name, 
                                          gsea_type = c("GO", "Reactome", "Hallmark", "WikiPathways"), 
                                          base_path = "results/GSEA_sets") {
  library(ggplot2)
  library(dplyr)
  library(stringr)
  
  for (etype in gsea_type) {
    file_path <- file.path(base_path, comparison_name, etype, paste0(etype, "_GSEA_results.csv"))
    if (!file.exists(file_path)) {
      message("Skipping: ", etype, " (file not found)")
      next
    }
    
    gsea_df <- read.csv(file_path)
    if (!"NES" %in% colnames(gsea_df)) {
      message("Skipping: ", etype, " (missing NES column)")
      next
    }
    
    gsea_df$direction <- ifelse(gsea_df$NES > 0, "Activated", "Suppressed")
    gsea_df$GeneCount <- sapply(strsplit(as.character(gsea_df$core_enrichment), "/"), length)
    gsea_df$GeneRatio <- gsea_df$GeneCount / gsea_df$setSize
    gsea_df$GeneRatio_mirror <- ifelse(gsea_df$direction == "Suppressed", -gsea_df$GeneRatio, gsea_df$GeneRatio)
    
    top_activated <- head(gsea_df %>% filter(direction == "Activated") %>% arrange(p.adjust), 15)
    top_suppressed <- head(gsea_df %>% filter(direction == "Suppressed") %>% arrange(p.adjust), 15)
    top_combined <- rbind(top_activated, top_suppressed)
    top_combined$Description <- stringr::str_wrap(top_combined$Description, width = 30)
    
    out_dir <- file.path(base_path, comparison_name, "BubblePlots")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Facet Plot
    p_facet <- ggplot(top_combined, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
      geom_point(aes(size = GeneCount, fill = -log10(p.adjust)), color = "black", shape = 21, stroke = 0.5) +
      scale_fill_gradient(low = "#4575b4", high = "#d73027", name = expression(-log[10]~"(adj. p-value)")) +
      scale_size(range = c(3, 10), name = "Gene Count") +
      facet_wrap(~direction, scales = "free_x", nrow = 1, strip.position = "top") +
      theme_minimal(base_size = 14) +
      labs(title = paste("GSEA:", etype, "\n", comparison_name), x = "GeneRatio", y = NULL) +
      theme(
        strip.background = element_rect(fill = "grey90", color = "grey20"),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.spacing.x = unit(1, "lines"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75)
      )
    
    ggsave(file.path(out_dir, paste0("Bubble_Facetcomp_top15_", etype, ".tiff")),
           plot = p_facet, width = 10, height = 16, dpi = 300)
    
    # Mirror Plot
    p_mirror <- ggplot(top_combined, aes(x = GeneRatio_mirror, y = reorder(Description, GeneRatio))) +
      geom_point(aes(size = GeneCount, fill = -log10(p.adjust)), color = "black", shape = 21, stroke = 0.5) +
      scale_fill_gradient(low = "#4575b4", high = "#d73027", name = expression(-log[10]~"(adj. p-value)")) +
      scale_size(range = c(3, 10), name = "Gene Count") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
      theme_minimal(base_size = 14) +
      labs(title = paste("GSEA:", etype, "\n", comparison_name), 
           x = "\u2190 Suppressed    GeneRatio    Activated \u2192", y = NULL) +
      theme(
        axis.text.y = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")
      )
    
    ggsave(file.path(out_dir, paste0("Bubble_Mirrorcomp_top15_", etype, ".tiff")),
           plot = p_mirror, width = 10, height = 16, dpi = 300)
  }
}


# ================== run bubble plot function for all comparisons ============ #
comparisons <- c("over_weight_vs_normal_weight", "obese_vs_normal_weight", 
                 "obese_vs_over_weight")

for (comp in comparisons) {
  plot_gsea_bubble_facet_mirror(comp)
}
# ======================== Bubble Plot Function (Facet + Mirror) ============= # 
plot_gsea_bubble_facet_mirror <- function(
    comparison_name,
    gsea_type = c("GO", "Reactome", "Hallmark", "WikiPathways", "KEGG"),
    base_path = "results/GSEA_sets"
) {
  library(ggplot2)
  library(dplyr)
  library(stringr)
  
  for (etype in gsea_type) {
    file_path <- file.path(base_path, comparison_name, etype, paste0(etype, "_GSEA_results.csv"))
    if (!file.exists(file_path)) {
      message("Skipping: ", etype, " (file not found)")
      next
    }
    
    gsea_df <- read.csv(file_path)
    if (!"NES" %in% colnames(gsea_df)) {
      message("Skipping: ", etype, " (missing NES column)")
      next
    }
    
    gsea_df$direction <- ifelse(gsea_df$NES > 0, "Activated", "Suppressed")
    gsea_df$GeneCount <- sapply(strsplit(as.character(gsea_df$core_enrichment), "/"), length)
    gsea_df$GeneRatio <- gsea_df$GeneCount / gsea_df$setSize
    gsea_df$GeneRatio_mirror <- ifelse(
      gsea_df$direction == "Suppressed",
      -gsea_df$GeneRatio,
      gsea_df$GeneRatio
    )
    
    top_activated <- head(gsea_df[gsea_df$direction == "Activated", ] %>% arrange(p.adjust), 10)
    top_suppressed <- head(gsea_df[gsea_df$direction == "Suppressed", ] %>% arrange(p.adjust), 10)
    top_combined <- rbind(top_activated, top_suppressed)
    top_combined$Description <- stringr::str_wrap(top_combined$Description, width = 40)
    
    out_dir <- file.path(base_path, comparison_name, "BubblePlots")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Facet plot
    p_facet <- ggplot(top_combined, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
      geom_point(aes(size = GeneCount, fill = -log10(p.adjust)),
                 color = "black", shape = 21, stroke = 0.5) +
      scale_fill_gradient(low = "#4575b4", high = "#d73027",
                          name = expression(-log[10]~"(adj. p-value)")) +
      scale_size(range = c(3, 10), name = "Gene Count") +
      facet_wrap(~direction, scales = "free_x", nrow = 1, strip.position = "top") +
      theme_minimal(base_size = 14) +
      labs(title = paste("GSEA:", etype, "\n", comparison_name),
           x = "GeneRatio", y = NULL) +
      theme(
        strip.background = element_rect(fill = "grey90", color = "grey20"),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.spacing.x = unit(1.5, "lines"),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        legend.position = "right"
      )
    
    ggsave(
      file.path(out_dir, paste0("Bubble_Facetcomp_top20_", etype, ".tiff")),
      plot = p_facet, width = 12, height = 10, dpi = 300
    )
    
    ggsave(
      file.path(out_dir, paste0("Bubble_Facetcomp_top20_", etype, ".svg")),
      plot = p_facet, width = 12, height = 10
    )
    
    # Mirror plot
    p_mirror <- ggplot(top_combined,
                       aes(x = GeneRatio_mirror, y = reorder(Description, GeneRatio))) +
      geom_point(aes(size = GeneCount, fill = -log10(p.adjust)),
                 color = "black", shape = 21, stroke = 0.5) +
      scale_fill_gradient(low = "#4575b4", high = "#d73027",
                          name = expression(-log[10]~"(adj. p-value)")) +
      scale_size(range = c(3, 10), name = "Gene Count") +
      geom_vline(xintercept = 0, linetype = "dashed",
                 color = "gray40", size = 0.6) +
      theme_minimal(base_size = 14) +
      labs(
        title = paste("GSEA:", etype, "\n", comparison_name),
        x = "\u2190 Suppressed    GeneRatio    Activated \u2192",
        y = NULL
      ) +
      theme(
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.title = element_text(face = "bold"),
        legend.position = "right"
      )
    
    ggsave(
      file.path(out_dir, paste0("Bubble_Mirrorcomp_top20_", etype, ".tiff")),
      plot = p_mirror, width = 12, height = 16, dpi = 300
    )
    
    ggsave(
      file.path(out_dir, paste0("Bubble_Mirrorcomp_top20_", etype, ".svg")),
      plot = p_mirror, width = 12, height = 16
    )
  }
}


# =============== run bubble plot function for all comparisons =============== #
comparisons <- c("over_weight_vs_normal_weight", "obese_vs_normal_weight", 
                 "obese_vs_over_weight")

for (comp in comparisons) {
  plot_gsea_bubble_facet_mirror(comp)
}
