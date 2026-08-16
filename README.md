# Simulation Study for a Bayesian Collective Reserving Framework

This repository contains the supplementary simulation study for the paper

**A Unified Bayesian Collective Reserving Framework with Hierarchical Temporal Smoothing**.

The purpose of the simulation is to provide a controlled assessment of the
proposed reserving architecture beyond the empirical insurance portfolio used
in the paper.

The study is not intended to establish that one fixed Bayesian specification is
universally optimal. Instead, it evaluates whether the proposed collective
Bayesian framework can recover and predict outstanding liabilities when the
underlying claims process contains several simultaneous sources of temporal,
developmental, and claim-type heterogeneity.

## Methods Compared

Three reserving approaches are evaluated:

1. **Full Bayesian model**  
   A joint count-composition-payment model with claim-type decomposition,
   hierarchical payment effects, evolving reporting patterns, type-specific
   settlement development, and calendar effects.

2. **Baseline Bayesian model**  
   Uses the same total-claim-count structure but models payments at the aggregate
   level without claim-type decomposition or the richer calendar and
   type-specific settlement structure.

3. **Paid Chain Ladder**  
   An aggregate paid-loss Chain Ladder fitted on the accident-period by total
   development-period triangle.

## Simulation Design

The principal `rich` scenario is deliberately not generated from an exact copy
of either fitted Bayesian model.

The synthetic claims process includes:

- a nonlinear accident-period claim-frequency pattern;
- overdispersed reported claim counts;
- a reporting-delay distribution that changes over accident time;
- unequal and time-varying claim-type composition;
- reporting-delay selection in claim type;
- a hurdle payment process with zero and positive payments;
- right-skewed Log-normal positive payments;
- heterogeneous claim-type payment levels;
- type-specific settlement trajectories;
- gradual calendar inflation; and
- a temporary calendar shock.

The simulated grid contains:

- 36 accident periods;
- 8 reporting-delay periods;
- 10 settlement-delay periods; and
- 4 claim types.

One claim type is deliberately uncommon, accounting for approximately 3.5% of
claims at baseline.

An optional `simple` scenario removes the evolving reporting-profile,
claim-composition, claim-type payment, type-specific settlement, and calendar
effects. It is included as an additional implementation and robustness check.

## Evaluation Design

For each independently generated synthetic portfolio, only information
observable at the valuation date is supplied to the fitted reserving models.

Because the data-generating process is known, the simulation permits two
distinct forms of evaluation.

### Expected-liability recovery

The primary target is the **oracle expected outstanding liability** under the
known data-generating mechanism.

The RBNS expectation is calculated exactly conditional on the reported claim
counts. The IBNR expectation integrates over the unknown future claim counts and
claim-type composition using Monte Carlo integration.

For this target, we report:

- relative bias;
- mean absolute percentage error (MAPE);
- root mean squared error (RMSE); and
- normalized RMSE.

### Posterior predictive performance

A separate realized future runoff is simulated in each replication. For the
Bayesian models we additionally evaluate:

- predictive-median error;
- empirical coverage of 95% posterior predictive intervals; and
- relative predictive interval width.

## Main Results: Rich Scenario

The principal experiment consists of **50 independently simulated portfolios**.

### Total Outstanding Liability

| Method | Relative Bias (%) | MAPE vs. Oracle (%) | Normalized RMSE (%) | MAPE vs. Realized Runoff (%) |
|---|---:|---:|---:|---:|
| Full Bayesian | -0.21 | 6.65 | 9.03 | 8.26 |
| Baseline Bayesian | 5.91 | 7.65 | 9.77 | 10.35 |
| Paid Chain Ladder | -1.26 | 15.38 | 22.62 | 14.95 |

The full Bayesian specification is essentially unbiased on average and reduces
MAPE relative to both the baseline Bayesian specification and Paid Chain
Ladder. Its oracle-expected-liability MAPE is **6.65%**, compared with **7.65%**
for the baseline Bayesian model and **15.38%** for Paid Chain Ladder.

The improvement from the baseline to the full Bayesian specification is also
visible in dispersion: normalized RMSE decreases from **9.77%** to **9.03%**,
while the median relative predictive interval width decreases from **70.0%** to
**52.1%**.

Both Bayesian specifications achieve 100% empirical coverage of the realized
total runoff by their 95% posterior predictive intervals across the 50
replications. Given the finite number of replications, this should not be
interpreted as evidence of exact nominal calibration.

### Component-Level Recovery

| Model | Quantity | Relative Bias (%) | MAPE vs. Oracle (%) | Predictive-Median MAPE vs. Realized (%) | 95% Predictive Coverage (%) |
|---|---|---:|---:|---:|---:|
| Full Bayesian | RBNS | -0.05 | 6.37 | 8.67 | 98 |
| Full Bayesian | IBNR | -0.42 | 10.51 | 12.81 | 98 |
| Full Bayesian | IBNR count | 1.63 | 5.69 | 7.45 | 98 |
| Full Bayesian | Total | -0.21 | 6.65 | 8.19 | 100 |
| Baseline Bayesian | RBNS | 5.91 | 7.79 | 9.34 | 100 |
| Baseline Bayesian | IBNR | 6.29 | 10.20 | 15.39 | 98 |
| Baseline Bayesian | IBNR count | 1.67 | 5.66 | 7.42 | 98 |
| Baseline Bayesian | Total | 5.91 | 7.65 | 9.41 | 100 |

The count component performs similarly under the two Bayesian specifications,
as expected because they use essentially the same aggregate claim-count
structure. The main gains from the richer specification occur in the payment
component: the full model substantially reduces RBNS bias and improves total
reserve accuracy while retaining strong IBNR and count recovery.

These results should not be interpreted as proving universal superiority of the
full specification. Simulation conclusions are necessarily conditional on the
chosen data-generating mechanism. The experiment instead demonstrates that the
proposed Bayesian collective architecture can recover reliable reserve
predictions when several forms of temporal evolution, development structure,
claim-type heterogeneity, and calendar effects operate simultaneously.

## Convergence

Across the 50 rich-scenario replications:

- the Full Bayesian model produced no divergent transitions;
- the Baseline Bayesian model produced no divergent transitions;
- the maximum observed R-hat was approximately 1.024 for the Full Bayesian
  model and 1.015 for the Baseline Bayesian model.

## Repository Structure

```text
.
├── README.md
├── run_simulation.sh
├── R/
│   ├── chain_ladder.R
│   ├── run_replication.R
│   ├── simulate_dataset.R
│   └── summarize_results.R
├── stan/
│   ├── baseline_model.stan
│   └── full_model.stan
└── results/
    ├── rich/
    │   └── summary/
    │       ├── component_recovery.csv
    │       └── total_reserve_comparison.csv
    └── simple/
        └── summary/
            ├── component_recovery.csv
            └── total_reserve_comparison.csv
```

## Reproducing the Simulation

The analysis requires R, `cmdstanr`, and a working CmdStan installation.

To run the principal rich scenario with 50 replications:

```bash
bash run_simulation.sh 50 rich 2 final
```

The final fitting profile uses four chains per Bayesian model, with 500 warm-up
and 500 retained sampling iterations per chain.

The optional simpler scenario can be run using:

```bash
bash run_simulation.sh 30 simple 2 final
```

The summary files are written to:

```text
results/rich/summary/total_reserve_comparison.csv
results/rich/summary/component_recovery.csv
```

with corresponding outputs for the simple scenario.

## Reproducibility

Replication-specific random seeds are deterministic. The repository contains
the code required to generate the synthetic portfolios, fit the Full and
Baseline Bayesian specifications, calculate the Paid Chain Ladder benchmark,
and reproduce the reported summary statistics.
