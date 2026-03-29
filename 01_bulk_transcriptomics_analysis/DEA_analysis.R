#################################################### Load libraries #################################################
library(DESeq2)
library(dplyr)
library(readr)

#################################################### Load and preprocess data #################################################

# Load the starcount data
normal_weight <- read_csv("results/expression_matrix/pdac_normal_weight.csv")
over_weight   <- read_csv("results/expression_matrix/pdac_over_weight.csv")
obese         <- read_csv("results/expression_matrix/pdac_Obese.csv")

# Remove Ensembl version
remove_ensembl_version <- function(df) {
  df %>%
    mutate(Ensembl_ID = sub("\\..*", "", df[[1]])) %>%
    dplyr::select(Ensembl_ID, everything(), -1)
}

normal_weight <- remove_ensembl_version(normal_weight)
over_weight   <- remove_ensembl_version(over_weight)
obese         <- remove_ensembl_version(obese)

################################## Convert to matrices #########################

normal_weight_matrix <- as.matrix(normal_weight[, -1])
rownames(normal_weight_matrix) <- normal_weight$Ensembl_ID

over_weight_matrix <- as.matrix(over_weight[, -1])
rownames(over_weight_matrix) <- over_weight$Ensembl_ID

obese_matrix <- as.matrix(obese[, -1])
rownames(obese_matrix) <- obese$Ensembl_ID

################################# Prepare group info ###########################

group_matrices <- list(
  normal_weight = normal_weight_matrix,
  over_weight   = over_weight_matrix,
  obese         = obese_matrix
)

################################### FIXED CONTRASTS (AS REQUESTED) #############

contrast_list <- list(
  over_weight_vs_normal_weight = c("over_weight", "normal_weight"),
  obese_vs_normal_weight       = c("obese", "normal_weight"),
  obese_vs_over_weight         = c("obese", "over_weight")
)

################################### DESeq2 runner function #####################

run_deseq2 <- function(group1, group2, group_matrices, output_dir = "results") {
  
  mat1 <- group_matrices[[group1]]
  mat2 <- group_matrices[[group2]]
  
  combined <- cbind(mat1, mat2)
  
  condition <- factor(
    c(rep(group1, ncol(mat1)), rep(group2, ncol(mat2))),
    levels = c(group1, group2)
  )
  
  coldata <- data.frame(condition = condition)
  rownames(coldata) <- colnames(combined)
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(combined),
    colData   = coldata,
    design    = ~ condition
  )
  
  dds <- DESeq(dds)
  
  # IMPORTANT:
  # Positive log2FC = higher in FIRST group (group1)
  # Negative log2FC = higher in SECOND group (group2)
  res <- results(dds, contrast = c("condition", group1, group2))
  
  # Embed contrast metadata
  res_df <- as.data.frame(res)
  res_df$contrast_group1 <- group1
  res_df$contrast_group2 <- group2
  
  out_path <- file.path(output_dir, paste0(group1, "_vs_", group2))
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(
    res_df,
    file = file.path(out_path,
                     paste0("deseq2_results_", group1, "_vs_", group2, ".csv"))
  )
  
  return(res_df)
}

########################### Loop over FIXED contrasts ##########################

all_results <- list()

for (nm in names(contrast_list)) {
  
  group1 <- contrast_list[[nm]][1]
  group2 <- contrast_list[[nm]][2]
  
  message("Running DESeq2 for: ", group1, " vs ", group2)
  
  res <- run_deseq2(group1, group2, group_matrices)
  all_results[[paste0(group1, "_vs_", group2)]] <- res
}

######################### HGNC conversion ######################################

library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

convert_ensembl_to_hgnc <- function(res_df) {
  
  res_df$ensembl_id <- rownames(res_df)
  
  conversion <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = res_df$ensembl_id,
    mart       = ensembl
  )
  
  res_df <- merge(
    res_df,
    conversion,
    by.x = "ensembl_id",
    by.y = "ensembl_gene_id",
    all.x = TRUE
  )
  
  res_df$hgnc_symbol <- ifelse(
    is.na(res_df$hgnc_symbol) | res_df$hgnc_symbol == "",
    res_df$ensembl_id,
    res_df$hgnc_symbol
  )
  
  rownames(res_df) <- make.unique(res_df$hgnc_symbol)
  return(res_df)
}

######################### HGNC conversion for all contrasts ####################

for (nm in names(contrast_list)) {
  
  group1 <- contrast_list[[nm]][1]
  group2 <- contrast_list[[nm]][2]
  
  result_name <- paste0(group1, "_vs_", group2)
  res_df <- all_results[[result_name]]
  
  if (is.null(res_df)) next
  
  res_hgnc <- convert_ensembl_to_hgnc(res_df)
  
  out_dir <- file.path("results", result_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(
    res_hgnc,
    file = file.path(out_dir, "deseq2_results_with_HGNC.csv")
  )
}

