# Bulk RNA-seq Analysis

This directory contains R scripts for differential expression analysis, pathway enrichment, and immune signature scoring of CPTAC-PDAC bulk RNA-seq data stratified by BMI.

## Scripts

| Script                                        | Purpose                                                                        |
| --------------------------------------------- | ------------------------------------------------------------------------------ |
| `DEA_analysis.R`                              | DESeq2 differential expression and HGNC gene symbol conversion                 |
| `GSEA.R`                                      | Gene set enrichment analysis across GO, KEGG, Reactome, and MSigDB collections |
| `immport_custom_GSEA.R`                       | ssGSEA scoring using ImmPort gene sets + Kruskal-Wallis testing                |
| `ssGSEA using ImmPort genes & Violin plots.R` | ssGSEA scoring and violin plot generation for significant signatures           |
| `KEGG_pathview.R`                             | KEGG pathway maps overlaid with log2 fold changes                              |

---

## 1. Data and Patient Stratification

Raw STAR count RNA-seq data were obtained from CPTAC-PDAC primary tumors via the Genomic Data Commons (GDC) using `TCGAbiolinks`. Patients with available BMI data (n = 140) were grouped by WHO classification:

| BMI Group     | BMI Range       | n   |
| ------------- | --------------- | --- |
| Normal weight | 18.5–24.9 kg/m² | 51  |
| Overweight    | 25.0–29.9 kg/m² | 58  |
| Obese         | ≥ 30.0 kg/m²    | 18  |

---

## 2. Differential Expression Analysis (`DEA_analysis.R`)

Raw counts were processed with DESeq2. Three pairwise contrasts were run:

| Contrast                     | Interpretation                         |
| ---------------------------- | -------------------------------------- |
| Overweight vs. Normal weight | Positive log2FC = higher in overweight |
| Obese vs. Normal weight      | Positive log2FC = higher in obese      |
| Obese vs. Overweight         | Positive log2FC = higher in obese      |

Ensembl version suffixes were stripped before analysis. Ensembl IDs were mapped to HGNC symbols using BioMart (`biomaRt`). Genes with no symbol were retained under their Ensembl ID.

---

## 3. Gene Set Enrichment Analysis (`GSEA.R`)

Genes were ranked by DESeq2 Wald test statistic (`stat` column) in descending order. GSEA was run per contrast using `clusterProfiler`.

**Gene set collections:**

| Collection          | Tool / Source             | ID type |
| ------------------- | ------------------------- | ------- |
| GO (ALL ontologies) | `gseGO` (org.Hs.eg.db)    | SYMBOL  |
| KEGG                | `gseKEGG`                 | Entrez  |
| Reactome            | `gsePathway` (ReactomePA) | Entrez  |

**GSEA parameters:**

| Parameter        | GO    | KEGG / Reactome |
| ---------------- | ----- | --------------- |
| `pvalueCutoff`   | 0.05  | 0.05            |
| `minGSSize`      | 10    | 10              |
| `maxGSSize`      | 500   | 500             |
| `by` (algorithm) | fgsea | default         |

---

## 4. ssGSEA with ImmPort Gene Sets (`immport_custom_GSEA.R`, `ssGSEA using ImmPort genes & Violin plots.R`)

Single-sample GSEA (ssGSEA) was run across all 140 samples using immune gene sets from ImmPort.

**Gene sets:** Downloaded from [immport.org](https://www.immport.org/home) on 14-04-2025 (`immport_signature.gmt`)

**ssGSEA parameters:**

| Parameter        | Value                                      |
| ---------------- | ------------------------------------------ |
| Input expression | VST-normalised (`vst(blind = FALSE)`)      |
| ID mapping       | Ensembl → HGNC symbol (bitr, org.Hs.eg.db) |
| Method           | ssGSEA via `GSVA::ssgseaParam` + `gsva`    |

**Statistical testing (3-group):**

| Step              | Method                      | Threshold    |
| ----------------- | --------------------------- | ------------ |
| Overall test      | Kruskal-Wallis              | p < 0.05     |
| Post-hoc pairwise | Dunn's test (BH correction) | adj.p < 0.05 |

---

## 5. KEGG Pathway Maps (`KEGG_pathview.R`)

Selected KEGG pathways were visualised using the `pathview` package, with gene-level log2 fold changes overlaid. Maps were generated for all three contrasts.

**Pathways visualised (per contrast):**

| KEGG ID  | Pathway                                 |
| -------- | --------------------------------------- |
| hsa04660 | T cell receptor signaling               |
| hsa04512 | ECM-receptor interaction                |
| hsa04510 | Focal adhesion                          |
| hsa04979 | Cholesterol metabolism                  |
| hsa03013 | RNA transport                           |
| hsa03040 | Spliceosome                             |
| hsa04612 | Antigen processing and presentation     |
| hsa03015 | mRNA surveillance pathway               |
| hsa05340 | Primary immunodeficiency                |
| hsa00020 | Citrate cycle (TCA cycle)               |
| hsa04613 | Neutrophil extracellular trap formation |
| hsa00190 | Oxidative phosphorylation               |
| hsa04620 | Toll-like receptor signaling            |
| hsa04972 | Pancreatic secretion                    |
| hsa04060 | Cytokine-cytokine receptor interaction  |
| hsa04064 | NF-kappa B signaling pathway            |
