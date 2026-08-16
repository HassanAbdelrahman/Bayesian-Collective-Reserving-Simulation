# Targeted reserving simulation study

This study compares three methods on synthetic portfolios:

1. **Full Bayesian model**: joint count/composition + hierarchical hurdle-lognormal payment model.
2. **Baseline Bayesian model**: same total-count process, but aggregate hurdle-lognormal severity without claim-type structure or calendar enrichment.
3. **Paid Chain Ladder**: aggregate paid-loss chain ladder on accident month × total development month.

The synthetic DGP is deliberately not identical to either fitted Bayesian model. The `rich` scenario contains changing reporting patterns, changing claim-type composition, type-specific severity/settlement effects, and calendar effects. The optional `simple` scenario removes those heterogeneity/calendar features.

## Primary target

The primary simulation truth is the **oracle expected outstanding liability** under the known DGP. RBNS expectation is calculated exactly conditional on reported counts; IBNR expectation is integrated over unknown counts and type composition by Monte Carlo. A separately simulated realized runoff is retained for predictive-interval coverage and predictive-median diagnostics.

## Pilot

```bash
cd reserving_simulation_study
bash run_simulation.sh 3 rich 1 pilot
```

Pilot profile: 2 chains, 300 warmup, 300 sampling iterations per Bayesian model.

## Final

```bash
bash run_simulation.sh 50 rich 2 final
```

Final profile: 4 chains, 500 warmup, 500 sampling iterations. Each replication fits **two** Bayesian models, so choose `MAX_PARALLEL` conservatively.

Optional robustness scenario:

```bash
bash run_simulation.sh 30 simple 2 final
```

## Main output

`results/rich/summary/total_reserve_comparison.csv`

reports, for each method:

- relative bias against oracle expected liability;
- MAPE against oracle expected liability;
- RMSE and normalized RMSE;
- error against one realized runoff;
- 95% predictive coverage and interval width for Bayesian models;
- fitting time and basic convergence diagnostics.

`component_recovery.csv` additionally reports RBNS, IBNR, total-reserve and IBNR-count recovery for the two Bayesian models.
