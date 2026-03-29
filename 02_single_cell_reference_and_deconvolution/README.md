# Single-Cell Reference Construction and BayesPrism Deconvolution

This directory contains notebooks for building scRNA-seq reference matrices and running BayesPrism deconvolution on CPTAC-PDAC bulk RNA-seq data, stratified by BMI group (normal-weight, overweight, obese).

## Notebooks

| Notebook                      | Reference Type                            | Source Dataset                                                |
| ----------------------------- | ----------------------------------------- | ------------------------------------------------------------- |
| `CD45+BayesPrism.ipynb`       | Immune (CD45+) cell reference             | GSM7502530 (GSE235452), n=3, 6,553 cells post-QC              |
| `non-immune_BayesPrism.ipynb` | Whole-tumour stromal/epithelial reference | GSE242230 (GSM7755911–GSM7755935), n=25, ~35,000 cells pre-QC |

## Analysis Overview

### 1. Quality Control

Adaptive, percentile-based thresholds were used in each dataset independently.

**CD45+ immune reference:**

- Feature counts: `max(200, 1st percentile)` to `min(6,000, 99th percentile)`
- Total UMI counts: below 98th percentile
- Mitochondrial content: below `min(20%, 95th percentile)`

**Whole-tumour reference:**

- Feature counts: `max(200, 0.5th percentile)` to `min(8,000, 99.5th percentile)`
- Total UMI counts: below 99th percentile
- Mitochondrial content: below `min(25%, 95th percentile)`

### 2. Preprocessing

Both datasets were processed with the standard Seurat workflow:

- Log normalisation (`NormalizeData`)
- 2,000 highly variable features (VST method)
- Scaling with regression of `percent.mt` and `nCount_RNA`
- PCA → UMAP (dims 1:30), Louvain clustering (resolution 0.4–0.5)

### 3. Cell Type Annotation

**CD45+ immune reference (`CD45+BayesPrism.ipynb`):**

- Immune gating on PTPRC expression, with epithelial and stromal cells excluded by marker scoring
- Automated annotation using SingleR against four references: Monaco, DICE, Blueprint, Novershtern
- Macrophages vs. monocytes distinguished by a macrophage-to-monocyte delta score (z-score ≥ 1.0–1.8 adaptive threshold) plus expression of ≥ 3 TAM hallmark genes
- Final reference: 21 immune cell types with ≥ 20 cells each (e.g. classical monocytes, TAMs, CD8+ exhausted T cells, NK cells)

**Whole-tumour reference (`non-immune_BayesPrism.ipynb`):**

- Epithelial cells identified by EPCAM/KRT marker scoring; malignant vs. normal epithelial cells separated by CopyKAT copy number inference
- Malignant cells subtyped as basal-like or classical PDAC using published gene signatures
- Stromal subtypes (myCAF, iCAF, apCAF, endothelial, pericyte, acinar) assigned by module scoring with lineage-specific markers
- Immune contamination excluded before finalising the stromal reference

### 4. BayesPrism Deconvolution

BayesPrism was run in three phases on CPTAC-PDAC bulk RNA-seq data (normal-weight, overweight, obese cohorts):

- **Phase 1 — Immune coarse:** CD45+ reference with broad immune cell types (Monaco main-level)
- **Phase 2 — Immune fine:** CD45+ reference with fine-grained immune subtypes (21 cell types)
- **Phase 3 — Non-immune:** Whole-tumour stromal/epithelial reference

All three phases used identical Gibbs sampling parameters:

| Parameter          | Value |
| ------------------ | ----- |
| `chain.length`     | 2,000 |
| `burn.in`          | 500   |
| `thinning`         | 2     |
| `outlier.cut`      | 0.01  |
| `outlier.fraction` | 0.05  |
| `input.type`       | count |

Outputs include per-sample cell type proportion estimates (theta), coefficient of variation (theta.cv), and reference expression profiles (phi).

---

## Reference Tables

### Table 1 — Quality Control Thresholds

| Parameter               | CD45+ Immune Reference      | Whole-Tumour Reference        |
| ----------------------- | --------------------------- | ----------------------------- |
| Feature count (lower)   | max(200, 1st percentile)    | max(200, 0.5th percentile)    |
| Feature count (upper)   | min(6,000, 99th percentile) | min(8,000, 99.5th percentile) |
| Total UMI count (upper) | < 98th percentile           | < 99th percentile             |
| Mitochondrial % (upper) | < min(20%, 95th percentile) | < min(25%, 95th percentile)   |

---

### Table 2 — Epithelial and Malignancy Classification Markers

| Cell Type                 | Markers                                                                               | Threshold                              |
| ------------------------- | ------------------------------------------------------------------------------------- | -------------------------------------- |
| Epithelial (CD45+ gating) | EPCAM, KRT8, KRT18, KRT19, KRT7, TACSTD2                                              | Score sum > 70th percentile per sample |
| Stromal (CD45+ gating)    | COL1A1, COL1A2, COL3A1, DCN, BGN, SPARC                                               | Score sum > threshold                  |
| Immune gate               | PTPRC, CD3E, CD19, CD14, FCGR3A, CD68, MS4A1                                          | PTPRC > 0 or ≥ 25% markers expressed   |
| Basal-like PDAC           | S100A2, KRT6A, KRT17, SERPINB4, LY6D, HMGA2, TP63, KRT5, DSC3, PKP1                   | Score > 0.20 and > classical score     |
| Classical PDAC            | TFF1, LGALS4, TSPAN8, REG4, ST6GALNAC1, CTSE, GATA6, HNF1A, HNF4A, FOXA2, FOXA3, PDX1 | Score > 0.20                           |
| Normal ductal epithelial  | CFTR, BICC1, SLC4A4, GLIS3, SCTR, AMBP, CA2                                           | Module score > 0.40                    |

CopyKAT parameters for malignancy inference: `ngene.chr=5`, `win.size=25`, `KS.cut=0.1`

---

### Table 3 — Stromal Cell Type Markers and Score Thresholds

| Cell Type              | Markers                                                                                  | Score Threshold |
| ---------------------- | ---------------------------------------------------------------------------------------- | --------------- |
| Fibroblast (quiescent) | DPT, DCN, LUM, COL1A1, COL1A2, COL3A1, BGN                                               | > 0.15          |
| myCAF                  | ACTA2, TAGLN, MYL9, POSTN, MMP11, HOPX, MYH11, CNN1, MYLK                                | > 0.30          |
| iCAF                   | IL6, CXCL12, CXCL14, CCL2, PDGFRA, CFD, PLA2G2A, CD34                                    | > 0.30          |
| apCAF                  | HLA-DRA, HLA-DPB1, CD74, CIITA, SLPI, HLA-DPA1                                           | > 0.30          |
| Blood endothelial      | PECAM1, VWF, CDH5, KDR, FLT1, FLI1, ERG                                                  | > 0.35          |
| Lymphatic endothelial  | PROX1, LYVE1, FLT4, PDPN, CCL21                                                          | > 0.35          |
| Pericyte / SMC         | RGS5, PDGFRB, CSPG4, MCAM, NOTCH3, ABCC9, ACTA2, MYH11, TAGLN, CNN1, MYLK, CASQ2, KCNAB1 | > 0.35          |
| Pancreatic stellate    | RGS5, PDGFRB, CSPG4, MCAM, NOTCH3, ABCC9, GFAP, SPARC, VIM, ACTA2, TIMP1, COL1A1, DES    | > 0.30          |
| Acinar                 | PRSS1, CPA1, CTRB1, CELA3A, AMY2A, PNLIP, CTRB2, CTRC                                    | > 0.30          |

---

### Table 4 — Immune Cell Annotation (CD45+ Reference)

**SingleR parameters:**

| Setting       | Main level                           | Fine level                           |
| ------------- | ------------------------------------ | ------------------------------------ |
| `de.method`   | wilcox                               | wilcox                               |
| `de.n`        | 50                                   | 20                                   |
| `fine.tune`   | TRUE                                 | TRUE                                 |
| `tune.thresh` | 0.05                                 | 0.10                                 |
| `sd.thresh`   | 1.0                                  | 1.5                                  |
| References    | Monaco, DICE, Blueprint, Novershtern | Monaco, DICE, Blueprint, Novershtern |

**TAM identification criteria:**

| Criterion                               | Value                                                                            |
| --------------------------------------- | -------------------------------------------------------------------------------- |
| Delta_MM z-score                        | ≥ 1.0–1.8 (adaptive, reduced in 0.1 steps until ≥ 30 TAMs found)                 |
| TAM hallmark genes (≥ 3 above median)   | C1QA, C1QB, C1QC, APOE, TREM2, SPP1, SIGLEC1, MARCO, MRC1, CD163                 |
| Monocyte signature                      | LYZ, S100A8, S100A9, FCN1, MS4A7, CTSS, VCAN, TYROBP, CCR2, LST1, FCGR3A, CX3CR1 |
| Exclusion: DC markers (top 10%)         | FLT3, IRF8, CLEC9A, XCR1, CD1C, FCER1A, LAMP3, CCR7                              |
| Exclusion: Neutrophil markers (top 10%) | CXCR2, FCGR3B, MMP8, CSF3R, ELANE, MPO                                           |

**CD8+ exhausted T cell criteria:**

| Criterion          | Value                                                                       |
| ------------------ | --------------------------------------------------------------------------- |
| Cell type          | CD8+ T cell subset                                                          |
| Exhaustion markers | PDCD1, HAVCR2, LAG3, TIGIT, TOX                                             |
| Score threshold    | Mean expression > 80th percentile among CD8+ T cells AND module score > 0.5 |
