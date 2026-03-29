# 00_data_acquisition

## Overview

This module downloads bulk RNA-seq data for pancreatic ductal adenocarcinoma (PDAC) patients from the Clinical Proteomic Tumor Analysis Consortium (CPTAC-3) via the Genomic Data Commons (GDC). The data are retrieved using the `TCGAbiolinks` R package and exported as gene-level count matrices for downstream analyses.

---

## Data Source

- **Project:** CPTAC-3
- **Portal:** Genomic Data Commons (GDC)
- **Data Type:** Gene Expression Quantification
- **Workflow:** STAR - Counts

---

## Cohort Definition

Samples are stratified based on Body Mass Index (BMI):

- **Normal weight:** BMI 18.5–24.9
- **Overweight:** BMI 25.0–29.9
- **Obese:** BMI ≥ 30.0

Predefined sample barcodes corresponding to each group are included in the script.

---

## Script

- `tcga_cptac_data_pull.R`  
  Performs:
  - Querying GDC
  - Downloading RNA-seq data
  - Preparing expression matrices
  - Exporting count data

---

## Output

Expression matrices are saved to: "results/expression_matrix/pdac_Obese.csv"
