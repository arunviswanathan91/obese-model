# TimiGP Analysis

This directory contains code and outputs supporting the TimiGP-based immune cell interaction analysis described in the supplementary methods of the associated manuscript.

## Overview

The TimiGP framework was applied to characterise immune cell gene pair interactions and infer cell-cell communication networks from bulk RNA-seq data. Analysis was restricted to normal-weight and overweight patient cohorts; the obese group was excluded due to insufficient sample size for robust statistical inference (n = 18).

## Immune Signature Sets

Three curated immune signature collections available in TimiGP were used:

| Signature Set      | Cell Types                 | Source                             |
| ------------------ | -------------------------- | ---------------------------------- |
| Bindea et al. 2013 | 28 immune cell types       | Colorectal cancer immune landscape |
| Newman LM22        | 22 immune cell types       | CIBERSORT reference matrix         |
| Zheng et al. 2021  | Pan-cancer T cell subtypes | scRNA-seq across 21 cancer types   |

## Analysis Pipeline

### 1. Cox Regression and Gene Pair Selection

- Gene pairs derived from signature marker genes
- Cox proportional hazards regression performed against survival endpoints
- Gene pairs filtered at raw p < 0.01; top 5% ranked by p-value retained for downstream analysis

### 3. Cell-Cell Interaction Enrichment

- `TimiEnrich` applied with Benjamini-Hochberg FDR (adjusted p < 0.05)
- Background gene pair distributions computed via `TimiBG`
- Cell pair definitions generated with `TimiCellPair` (multi-core)
- Permutation-based FDR estimated using `TimiPermFDR` (100 iterations)

### 4. Network Construction and Favourability Scoring

- Interaction networks built with `TimiCellNetwork` using signature-specific annotations
- Prognostic favourability scores computed via `TimiFS`
- Network topology exported to Cytoscape for visualisation
- Outputs include chord diagrams and dot plots stratified by BMI group
