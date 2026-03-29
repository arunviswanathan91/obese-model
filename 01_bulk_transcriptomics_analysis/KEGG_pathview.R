############### Pathway view KEGG ###############
#🔴 Red = Upregulated in normal weight (positive log2FC).
#🟢 Green = Downregulated in normal weight (i.e., up in overweight).

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(AnnotationDbi)
library(pathview)
library(ggplot2)

# Use DESeq2 result
res_df <- res_ow_nw  # or any other results
res_df <- na.omit(res_df)

# Map Ensembl IDs to Entrez IDs
conversion <- bitr(res_df$ensembl_id,
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# Merge converted IDs back to the DESeq2 result
res_mapped <- merge(res_df, conversion, by.x = "ensembl_id", by.y = "ENSEMBL")
res_mapped <- res_mapped[!duplicated(res_mapped$ENTREZID), ]

# Prepare named vector of log2FC with Entrez IDs
gene_list <- res_mapped$log2FoldChange
names(gene_list) <- res_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# Create output directory 
output_dir <- "/results/pathview_overweight_vs_normalweight"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

# Run Pathview for T cell receptor pathway
pathview(
  gene.data   = gene_list,
  pathway.id  = "hsa04660",
  species     = "hsa",
  out.suffix  = "TCR_normalweight_vs_overweight",
  kegg.native = TRUE,      
  xml.file    = TRUE,       
  map.null    = TRUE,      
  both.dirs   = TRUE      
)


# Run Pathview for ECM-receptor interaction
pathview(
  gene.data   = gene_list,
  pathway.id  = "hsa04512",
  species     = "hsa",
  out.suffix  = "ECM_normalweight_vs_overweight",
  kegg.native = TRUE
)

# Run Pathview for Focal adhesion pathway
pathview(
  gene.data   = gene_list,
  pathway.id  = "hsa04510",
  species     = "hsa",
  out.suffix  = "FocalAdhesion_normalweight_vs_overweight",
  kegg.native = TRUE
)

# Run Pathview for Cholesterol metabolism
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04979",        
  species     = "hsa",             
  out.suffix  = "Cholesterol metabolism_normalweight_vs_overweight",     
  kegg.native = TRUE
)
# Run Pathview for RNA transport
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa03013",        
  species     = "hsa",             
  out.suffix  = "RNA transport_normalweight_vs_overweight",     
  kegg.native = TRUE
)

# Run Pathview for Spliceosome
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa03040",        
  species     = "hsa",             
  out.suffix  = "Spliceosome_normalweight_vs_overweight",     
  kegg.native = TRUE
)

# Run Pathview for Antigen binding
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04612",        
  species     = "hsa",             
  out.suffix  = "Antigen processing and presentation_normalweight_vs_overweight",     
  kegg.native = TRUE
)

# Run Pathview for mRNA surveillance pathway
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa03015",        
  species     = "hsa",             
  out.suffix  = "mRNA surveillance pathway_normalweight_vs_overweight",     
  kegg.native = TRUE
)


# Run Pathview for Primary immunodeficiency
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa05340",        
  species     = "hsa",             
  out.suffix  = "Primary immunodeficiency",     
  kegg.native = TRUE
)

# Run Pathview for Citrate cycle (TCA cycle)

pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa00020",        
  species     = "hsa",             
  out.suffix  = "Citrate cycle (TCA cycle)",     
  kegg.native = TRUE
)


# Run Pathview for Neutrophil extracellular trap formation
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04613",        
  species     = "hsa",             
  out.suffix  = "Neutrophil extracellular trap formation",     
  kegg.native = TRUE
)


# Run Pathview for Oxidative phosphorylation

pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa00190",        
  species     = "hsa",             
  out.suffix  = "Oxidative phosphorylation",     
  kegg.native = TRUE
)

# Run Pathview for Toll-like receptor signaling pathway
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04620",        
  species     = "hsa",             
  out.suffix  = "Toll-like receptor signaling pathway",     
  kegg.native = TRUE
)

# Run Pathview for Pancreatic secretion
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04972",        
  species     = "hsa",             
  out.suffix  = "Pancreatic secretion",     
  kegg.native = TRUE
)

# Run Pathview for Pancreatic secretion
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04060",        
  species     = "hsa",             
  out.suffix  = "Cytokine-cytokine receptor interaction",     
  kegg.native = TRUE
)

# Run Pathview for Pancreatic secretion
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04064",        
  species     = "hsa",             
  out.suffix  = "NF-kappa B signaling pathway",     
  kegg.native = TRUE
)
############### Pathway view KEGG ###############

# Use DESeq2 result
res_df <- res_ob_nw  
res_df <- na.omit(res_df)

# Map Ensembl IDs (from column) to Entrez IDs
conversion <- bitr(res_df$ensembl_id,
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# Merge converted IDs back to the DESeq2 result
res_mapped <- merge(res_df, conversion, by.x = "ensembl_id", by.y = "ENSEMBL")
res_mapped <- res_mapped[!duplicated(res_mapped$ENTREZID), ]

# Prepare named vector of log2FC with Entrez IDs
gene_list <- res_mapped$log2FoldChange
names(gene_list) <- res_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# Create output directory 
output_dir <- "/results/pathview_obese_vs_normalweight"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

# Run Pathview for T cell receptor pathway
pathview(
  gene.data   = gene_list,
  pathway.id  = "hsa04660",
  species     = "hsa",
  out.suffix  = "TCR_normalweight_vs_obese",
  kegg.native = TRUE
)

# Run Pathview for ECM-receptor interaction
pathview(
  gene.data   = gene_list,
  pathway.id  = "hsa04512",
  species     = "hsa",
  out.suffix  = "ECM_normalweight_vs_obese",
  kegg.native = TRUE
)

# Run Pathview for Focal adhesion pathway
pathview(
  gene.data   = gene_list,
  pathway.id  = "hsa04510",
  species     = "hsa",
  out.suffix  = "FocalAdhesion_normalweight_vs_obese",
  kegg.native = TRUE
)

# Run Pathview for Cholesterol metabolism
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04979",        
  species     = "hsa",             
  out.suffix  = "Cholesterol metabolism_normalweight_vs_obese",     
  kegg.native = TRUE
)

# Run Pathview for RNA transport
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa03013",        
  species     = "hsa",             
  out.suffix  = "RNA transport_normalweight_vs_obese",     
  kegg.native = TRUE
)

# Run Pathview for Spliceosome
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa03040",        
  species     = "hsa",             
  out.suffix  = "Spliceosome_normalweight_vs_obese",     
  kegg.native = TRUE
)

# Run Pathview for Antigen binding
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04612",        
  species     = "hsa",             
  out.suffix  = "Antigen processing and presentation_normalweight_vs_obese",     
  kegg.native = TRUE
)

# Run Pathview for mRNA surveillance pathway
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa03015",        
  species     = "hsa",             
  out.suffix  = "mRNA surveillance pathway_normalweight_vs_obese",     
  kegg.native = TRUE
)

# Run Pathview for Primary immunodeficiency
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa05340",        
  species     = "hsa",             
  out.suffix  = "Primary immunodeficiency",     
  kegg.native = TRUE
)

# Run Pathview for Citrate cycle (TCA cycle)

pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa00020",        
  species     = "hsa",             
  out.suffix  = "Citrate cycle (TCA cycle)",     
  kegg.native = TRUE
)


# Run Pathview for Neutrophil extracellular trap formation
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04613",        
  species     = "hsa",             
  out.suffix  = "Neutrophil extracellular trap formation",     
  kegg.native = TRUE
)


# Run Pathview for Oxidative phosphorylation

pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa00190",        
  species     = "hsa",             
  out.suffix  = "Oxidative phosphorylation",     
  kegg.native = TRUE
)

# Run Pathview for Toll-like receptor signaling pathway
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04620",        
  species     = "hsa",             
  out.suffix  = "Toll-like receptor signaling pathway",     
  kegg.native = TRUE
)

# Run Pathview for Pancreatic secretion
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04972",        
  species     = "hsa",             
  out.suffix  = "Pancreatic secretion",     
  kegg.native = TRUE
)

# Run Pathview for Pancreatic secretion
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04060",        
  species     = "hsa",             
  out.suffix  = "Cytokine-cytokine receptor interaction",     
  kegg.native = TRUE
)

# Run Pathview for Pancreatic secretion
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04064",        
  species     = "hsa",             
  out.suffix  = "NF-kappa B signaling pathway",     
  kegg.native = TRUE
)

############### Pathway view KEGG ###############

# Use DESeq2 result
res_df <- res_ob_ow  
res_df <- na.omit(res_df)

# Map Ensembl IDs to Entrez IDs
conversion <- bitr(res_df$ensembl_id,
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# Merge converted IDs back to the DESeq2 result
res_mapped <- merge(res_df, conversion, by.x = "ensembl_id", by.y = "ENSEMBL")
res_mapped <- res_mapped[!duplicated(res_mapped$ENTREZID), ]

# Prepare named vector of log2FC with Entrez IDs
gene_list <- res_mapped$log2FoldChange
names(gene_list) <- res_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# Create output directory 
output_dir <- "/results/pathview_obese_vs_overweight"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

# Run Pathview for T cell receptor pathway
pathview(
  gene.data   = gene_list,
  pathway.id  = "hsa04660",
  species     = "hsa",
  out.suffix  = "TCR_overweight_vs_obese",
  kegg.native = TRUE
)

# Run Pathview for ECM-receptor interaction
pathview(
  gene.data   = gene_list,
  pathway.id  = "hsa04512",
  species     = "hsa",
  out.suffix  = "ECM_overweight_vs_obese",
  kegg.native = TRUE
)

# Run Pathview for Focal adhesion pathway
pathview(
  gene.data   = gene_list,
  pathway.id  = "hsa04510",
  species     = "hsa",
  out.suffix  = "FocalAdhesion_overweight_vs_obese",
  kegg.native = TRUE
)

# Run Pathview for Cholesterol metabolism
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04979",        
  species     = "hsa",             
  out.suffix  = "Cholesterol metabolism_overweight_vs_obese",     
  kegg.native = TRUE
)

# Run Pathview for RNA transport
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa03013",        
  species     = "hsa",             
  out.suffix  = "RNA transport_overweight_vs_obese",     
  kegg.native = TRUE
)

# Run Pathview for Spliceosome
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa03040",        
  species     = "hsa",             
  out.suffix  = "Spliceosome_overweight_vs_obese",     
  kegg.native = TRUE
)

# Run Pathview for Antigen binding
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04612",        
  species     = "hsa",             
  out.suffix  = "Antigen processing and presentation_overweight_vs_obese",     
  kegg.native = TRUE
)

# Run Pathview for mRNA surveillance pathway
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa03015",        
  species     = "hsa",             
  out.suffix  = "mRNA surveillance pathway_overweight_vs_obese",     
  kegg.native = TRUE
)

# Run Pathview for Primary immunodeficiency
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa05340",        
  species     = "hsa",             
  out.suffix  = "Primary immunodeficiency",     
  kegg.native = TRUE
)

# Run Pathview for Citrate cycle (TCA cycle)

pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa00020",        
  species     = "hsa",             
  out.suffix  = "Citrate cycle (TCA cycle)",     
  kegg.native = TRUE
)


# Run Pathview for Neutrophil extracellular trap formation
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04613",        
  species     = "hsa",             
  out.suffix  = "Neutrophil extracellular trap formation",     
  kegg.native = TRUE
)


# Run Pathview for Oxidative phosphorylation

pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa00190",        
  species     = "hsa",             
  out.suffix  = "Oxidative phosphorylation",     
  kegg.native = TRUE
)

# Run Pathview for Toll-like receptor signaling pathway
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04620",        
  species     = "hsa",             
  out.suffix  = "Toll-like receptor signaling pathway",     
  kegg.native = TRUE
)

# Run Pathview for Pancreatic secretion
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04972",        
  species     = "hsa",             
  out.suffix  = "Pancreatic secretion",     
  kegg.native = TRUE
)

# Run Pathview for Pancreatic secretion
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04060",        
  species     = "hsa",             
  out.suffix  = "Cytokine-cytokine receptor interaction",     
  kegg.native = TRUE
)

# Run Pathview for Pancreatic secretion
pathview(
  gene.data   = gene_list,         
  pathway.id  = "hsa04064",        
  species     = "hsa",             
  out.suffix  = "NF-kappa B signaling pathway",     
  kegg.native = TRUE
)
