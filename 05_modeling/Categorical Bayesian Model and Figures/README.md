# Categorical Bayesian Hierarchical Model

This directory contains notebooks for fitting a Bayesian hierarchical model in PyMC that estimates the effect of BMI group (normal-weight, overweight, obese) on immune and stromal signature scores.

## Notebooks

| Notebook                             | Purpose                                |
| ------------------------------------ | -------------------------------------- |
| `Categorical_Bayesian_Model.ipynb`   | Model fitting and posterior inference  |
| `Heatmap_Categorical.ipynb`          | Heatmap of group effects               |
| `Adaptive_Heatmap_Categorical.ipynb` | Adaptive heatmap visualisation         |
| `Ridge_plots_categorical.ipynb`      | Ridge plots of posterior distributions |
| `Stacked_barplots_categorical.ipynb` | Stacked bar plots of credible effects  |
| `Diagnostic_Plots_Categorical.ipynb` | Convergence and diagnostic plots       |

---

## Model Structure

The categorical model estimates separate effects for overweight and obese groups relative to normal-weight as the reference:

```
Z_ij = baseline_j + β_overweight[k[j]] × I(BMI = overweight)
                  + β_obese[k[j]]      × I(BMI = obese)
```

This also allows a derived contrast: obese vs. overweight (`feature_effect_obese_vs_overweight`), computed directly from the posterior.

The hierarchical structure has two levels:

- **Cell-type level:** shared group effects per lineage
- **Feature level:** individual deviations from the cell-type mean

Non-centered parameterisation was used at the feature level to avoid sampler divergences.

---

## Priors

Priors are compartment-specific:

| Parameter                             | Non-immune | Immune fine | Immune coarse |
| ------------------------------------- | ---------- | ----------- | ------------- |
| Cell-type effect σ (`celltype_sigma`) | 0.20       | 0.25        | 0.18          |
| Feature effect σ (`feature_sigma`)    | 0.30       | 0.40        | 0.28          |
| Baseline σ (`baseline_sigma`)         | 1.5        | 1.5         | 1.5           |
| Observation noise σ (`obs_sigma`)     | 1.0        | 1.0         | 1.0           |

---

## Sampling (NUTS)

| Parameter         | Value                    |
| ----------------- | ------------------------ |
| Sampler           | NUTS (No-U-Turn Sampler) |
| Chains            | 4                        |
| Tuning steps      | 2,000                    |
| Sampling draws    | 2,000                    |
| Target acceptance | 0.99                     |

**Convergence criterion:** R-hat < 1.01

---

## Effect Assessment

Credibility is assessed using the 95% HDI for three contrasts per feature:

| Contrast              | Credible if           |
| --------------------- | --------------------- |
| Overweight vs. Normal | 95% HDI excludes zero |
| Obese vs. Normal      | 95% HDI excludes zero |
| Obese vs. Overweight  | 95% HDI excludes zero |
