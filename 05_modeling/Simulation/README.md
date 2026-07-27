# Cohort-Adequacy Simulation

This directory contains a simulation study that checks whether the Bayesian hierarchical model can pick up a real BMI effect at the actual CPTAC-3 cohort sizes, and whether it avoids inventing an effect when none exists. It is an in-silico adequacy check, not external replication.

## Notebook

| Notebook | Purpose |
| --- | --- |
| `BHM_simulation.ipynb` | Power (recovery) and calibration (false-positive) simulation |

---

## What It Tests

Two arms, run over `N_REPLICATES = 4` independently simulated datasets each:

- **POWER arm** — a true overweight-vs-normal effect (`TRUE_EFFECT_SIZE = 0.5` SD) is injected into 50% of features (`FRAC_ACTIVE`). Measures the fraction of truly-active features whose 95% HDI correctly excludes 0 with the right sign (detection rate).
- **NULL arm** — no effect is injected (flat across BMI). Measures the fraction of truly-null features whose 95% HDI wrongly excludes 0 (false-positive rate); should sit at or below the nominal ~5%.

## Simulation Setup

The simulation re-uses the same model architecture, priors, and sampler settings as the manuscript model:

- Three-level hierarchy: cell type > feature (signature) > observation, plus a non-centered patient-level random intercept for pseudoreplication control
- Real cohort sizes preserved: N = 51 normal, OW = 58 overweight, OB = 18 obese
- Structure approximates the immune-coarse compartment: 9 cell types x 6 signatures/cell type = 54 features
- Priors copied from `PRIOR_CONFIGS["immune_coarse"]`

| Parameter | Value |
| --- | --- |
| Cell-type σ (`celltype_sigma`) | 0.15 |
| Feature σ (`feature_sigma`) | 0.40 |
| Patient intercept σ (`patient_sigma`) | 0.50 |
| Baseline σ (`baseline_sigma`) | 1.5 |
| Observation noise σ (`obs_sigma`) | 1.0 |

Sampler (NUTS): 4 chains, 2,000 tuning steps, 2,000 draws, target acceptance 0.99 — identical to the manuscript model.

**Credibility rule:** 95% HDI excludes 0, R-hat < 1.01, bulk-ESS > 400.

---

## Results

| Metric | Value |
| --- | --- |
| Convergence rate (both arms) | 100% |
| Detection rate, power arm (HDI excludes 0) | 59.3% ± 3.7% |
| Detection rate, power arm (correct sign) | 59.3% ± 3.7% |
| False-positive rate, null arm | 0.0% ± 0.0% (0/216 tests) |

At the real cohort sizes, the model has ~59% power to recover a true overweight-vs-normal effect of 0.5 SD, and does not invent effects from noise — the false-positive rate stays well below the nominal 5%.
