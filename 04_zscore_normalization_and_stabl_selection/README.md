# Z-Score Calculation and STABL Feature Selection

This notebook computes immune signature scores from deconvoluted bulk RNA-seq data and applies the STABL machine learning framework to identify signatures associated with BMI.

**Notebook:** `zscore_calcualtion_and_stabl.ipynb`

---

## 1. Signature Score Calculation

Z-scores were computed per gene across samples, then averaged within each signature to produce a single score per sample per signature.

**Key steps:**

- Gene expression values were mean-centred and scaled across samples
- Z-scores were capped at ±3 (winsorization) to limit the influence of extreme values
- Missing gene values were imputed with 0 (neutral contribution)
- Signatures with fewer than 4 detectable genes were excluded

**Parameters:**

| Parameter                              | Value |
| -------------------------------------- | ----- |
| Winsorization cap                      | ±3    |
| Minimum detectable genes per signature | 4     |
| Missing gene imputation                | 0     |

---

## 2. STABL Feature Selection

STABL (v1.0.0) was used to identify signatures robustly associated with BMI from high-dimensional, small-sample data. It combines stability selection with Knockoff-based FDR control.

**Pre-processing (applied equally to train and test):**

- Samples with missing BMI metadata were removed before splitting
- Features with > 50% missing values were dropped
- Remaining missing values were filled by median imputation (`SimpleImputer`)
- Near-zero variance features were removed (`VarianceThreshold = 1e-4`)

**Train/test split:**

- 80% train / 20% test (`TEST_SIZE = 0.20`)
- Stratified by BMI group to preserve class proportions
- Analysis repeated across 8 random seeds: [11, 22, 33, 44, 55, 66, 77, 88]

**STABL modes:**

| Mode        | Base estimator         | Regularisation grid        | Purpose                            |
| ----------- | ---------------------- | -------------------------- | ---------------------------------- |
| Categorical | L1 Logistic Regression | C: logspace(−3, 3, 30)     | Distinguish BMI groups             |
| Continuous  | Lasso                  | alpha: logspace(−7, 1, 50) | Track signatures changing with BMI |

**Validation estimators (applied to held-out test set):**

| Mode        | Estimator                                        |
| ----------- | ------------------------------------------------ |
| Categorical | L2 Logistic Regression (balanced, max_iter=1000) |
| Continuous  | Ridge Regression (alpha=1.0)                     |

**STABL parameters:**

| Parameter             | Value                    |
| --------------------- | ------------------------ |
| `n_bootstraps`        | 500                      |
| `sample_fraction`     | 0.5                      |
| `artificial_type`     | knockoff                 |
| `fdr_threshold_range` | 0.10 to 0.98 (step 0.02) |
| Parallel workers      | 8                        |

**Output:**

Features selected in categorical or continuous mode across seeds were merged into a single consensus feature set for downstream analysis.
