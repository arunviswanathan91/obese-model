# Continuous Bayesian Hierarchical Model

This directory contains notebooks for fitting a Bayesian hierarchical model in PyMC that quantifies how immune and stromal signature scores change continuously with BMI.

## Notebooks

| Notebook                           | Purpose                                |
| ---------------------------------- | -------------------------------------- |
| `Continues_Bayesian_Model.ipynb`   | Model fitting and posterior inference  |
| `Heatmap_Continous.ipynb`          | Heatmap visualisation of BMI slopes    |
| `Ridge_Plot_Continoues.ipynb`      | Ridge plots of posterior distributions |
| `Stacked_Bar_Plot_Continues.ipynb` | Stacked bar plots of credible effects  |
| `Diagnostic_Plots_Continous.ipynb` | Convergence and diagnostic plots       |

---

## Model Structure

The continuous model estimates a per-feature BMI slope using partial pooling across signatures within each cell type:

```
Z_ij = baseline_j + slope_j × BMI_standardized_i
```

BMI was z-score standardised before fitting. Slopes reflect the expected change in signature score per standard deviation of BMI.

The hierarchical structure has two levels:

- **Cell-type level:** shared slope distribution per lineage (`celltype_bmi_slope`)
- **Feature level:** individual feature slopes as deviations from the cell-type mean (`feature_bmi_slope`)

Non-centered parameterisation was used at the feature level to avoid sampler divergences.

---

## Priors

Priors are compartment-specific to reflect differences in expected signal magnitude and noise:

| Parameter                            | Non-immune | Immune fine | Immune coarse |
| ------------------------------------ | ---------- | ----------- | ------------- |
| Cell-type slope σ (`celltype_sigma`) | 0.20       | 0.25        | 0.18          |
| Feature slope σ (`feature_sigma`)    | 0.30       | 0.40        | 0.28          |
| Baseline σ (`baseline_sigma`)        | 1.5        | 1.5         | 1.5           |
| Observation noise σ (`obs_sigma`)    | 1.0        | 1.0         | 1.0           |

---

## Sampling (NUTS)

| Parameter         | Value                    |
| ----------------- | ------------------------ |
| Sampler           | NUTS (No-U-Turn Sampler) |
| Chains            | 4                        |
| Tuning steps      | 2,000                    |
| Sampling draws    | 2,000                    |
| Target acceptance | 0.99                     |

**Convergence criteria:** R-hat < 1.01 for all parameters.

---

## Effect Assessment

Statistical credibility and practical significance were assessed using HDI and ROPE:

**Credibility:** 95% HDI excludes zero (Tier 3)

**ROPE analysis** — five thresholds tested on the standardised slope scale:

| Threshold              | Values                            |
| ---------------------- | --------------------------------- |
| ROPE thresholds tested | ±0.05, ±0.10, ±0.15, ±0.20, ±0.30 |

**Effect tiers:**

| Tier            | Criteria                                                      |
| --------------- | ------------------------------------------------------------- |
| Tier 1 (large)  | Credible AND \|slope\| > 0.2 with > 95% posterior probability |
| Tier 2 (medium) | Credible AND \|slope\| > 0.1 with > 95% posterior probability |
| Tier 3 (any)    | 95% HDI excludes zero                                         |
