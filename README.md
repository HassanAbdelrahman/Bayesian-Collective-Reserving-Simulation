# Simulation Study for a Bayesian Collective Reserving Framework

This repository contains the supplementary simulation study for the paper

**A Unified Bayesian Collective Reserving Framework with Hierarchical Temporal Smoothing**.

The purpose of the simulation is to provide a controlled assessment of the
proposed reserving architecture beyond the empirical insurance portfolio used
in the paper.

The study is not intended to establish that one fixed Bayesian specification is
universally optimal. Instead, it evaluates whether the proposed collective
Bayesian framework can recover expected outstanding liabilities and generate
well-calibrated predictive distributions when the underlying claims process
contains several simultaneous sources of temporal, development, and claim-type
heterogeneity.

## Methods Compared

Three reserving approaches are evaluated.

1. **Full Bayesian model**  
   A joint count-composition-payment model with claim-type decomposition,
   evolving reporting patterns, hierarchical claim-type payment effects,
   type-specific settlement development, and calendar effects.

2. **Baseline Bayesian model**  
   Uses the same total-count model, but aggregates payments across claim types
   and omits claim-type payment effects, type-specific settlement development,
   and calendar effects.

3. **Paid Chain Ladder**  
   An aggregate paid-loss Chain Ladder fitted to an accident-period by total
   development-period triangle.

The purpose of the Baseline Bayesian model is to provide a reduced Bayesian
benchmark with the same general count machinery but substantially less payment
structure. The Paid Chain Ladder provides a conventional aggregate reserving
benchmark.

## Simulation Design

### Indexing and valuation structure

Each synthetic portfolio contains

- 36 accident periods;
- 8 reporting-delay periods;
- 10 settlement-delay periods; and
- 4 claim types.

Let \(t=1,\ldots,36\) denote accident period, \(d=0,\ldots,7\) reporting
delay, \(s=0,\ldots,9\) settlement delay, and \(j=1,\ldots,4\) claim type.
The valuation date is the end of accident period 36.

A claim-count cell is observed when

\[
t+d \leq 36,
\]

and is IBNR when \(t+d>36\). A payment cell is observed when

\[
t+d+s \leq 36.
\]

Future payments from already reported claims form the RBNS reserve, whereas
future payments associated with cells satisfying \(t+d>36\) form the IBNR
reserve. Because the entire future is simulated, the realized RBNS, IBNR, and
total reserves are known exactly in each replication.

### 1. Exposure and reported claim counts

The exposure level varies smoothly over accident time:

\[
e_t
=
\exp\left[
0.08\sin\left\{\frac{2\pi(t-1)}{12}\right\}
+0.04u_t
\right],
\]

where \(u_t\) is a linearly scaled accident-time index ranging from \(-1\) to
\(1\).

Reported claim counts follow a Negative Binomial distribution,

\[
N_{t,d}
\sim
\operatorname{NegBin}_2(\lambda_{t,d},\phi_N),
\]

with

\[
\log \lambda_{t,d}
=
\log e_t
+\mu_N
+\alpha_t^{(N)}
+\beta_d^{(N)}
+\kappa_d u_t.
\]

The dispersion parameter is \(\phi_N=35\), and the intercept is chosen to
produce approximately 300 claims per accident period before temporal
variation.

The accident-period effect is deliberately nonlinear:

\[
\alpha_t^{(N)}
=
\operatorname{center}\left[
0.16\sin\left\{\frac{2\pi(t-1)}{18}\right\}
+
0.13\tanh\left\{\frac{t-17}{4}\right\}
\right].
\]

Thus, the simulated frequency path contains both cyclical movement and a
smooth level transition rather than being generated from a random-walk prior.

The baseline reporting-delay probabilities are proportional to

\[
(0.58,\;0.20,\;0.09,\;0.05,\;0.035,\;0.020,\;0.015,\;0.010),
\]

and \(\beta_d^{(N)}\) is obtained by centering their logarithms.

In the rich scenario, the reporting pattern also changes over accident time.
The delay-specific slopes are the centered version of

\[
(0.24,\;0.11,\;0.03,\;-0.04,\;-0.10,\;-0.15,\;-0.20,\;-0.25).
\]

Consequently, early and late reporting delays evolve differently over time.

### 2. Claim-type composition

Conditional on the total count,

\[
(N_{t,d,1},\ldots,N_{t,d,4})
\mid N_{t,d}
\sim
\operatorname{Multinomial}
\left(N_{t,d},\boldsymbol{\pi}_{t,d}\right).
\]

The baseline type shares are

\[
(0.49,\;0.28,\;0.195,\;0.035),
\]

so the fourth claim type is deliberately sparse.

Type 1 is the reference category. For types \(j=2,3,4\),

\[
\log\frac{\pi_{t,d,j}}{\pi_{t,d,1}}
=
\gamma_j
+a_j u_t
+c_{d,j},
\]

where the baseline log-odds \(\gamma_j\) reproduce the shares above. The
accident-time slopes are

\[
(a_2,a_3,a_4)=(-0.18,\;0.30,\;0.12).
\]

The reporting-delay tilts increase linearly with normalized reporting delay.
Their terminal magnitudes are \(0.22\), \(-0.18\), and \(0.50\) for types
2, 3, and 4, respectively, with zero tilt at reporting delay \(d=0\).

The rich scenario therefore contains both changing portfolio composition over
accident time and selection of claim type by reporting delay.

### 3. Payment occurrence

Payments are generated on the
\((t,d,s,j)\) grid. For a cell containing \(n\) claims,

\[
Z_{t,d,s,j}
\sim
\operatorname{Bernoulli}(p_{t,d,s,j}),
\]

where \(Z=1\) indicates a positive aggregate payment and

\[
\operatorname{logit}(p_{t,d,s,j})
=
\mu_p
+\alpha_t^{(p)}
+\beta_d^{(p)}
+\zeta_s^{(p)}
+\tau_j^{(p)}
+0.48\log n.
\]

The global intercept is \(\mu_p=-2.25\).

The accident-period payment-incidence effect is again deterministic and
nonlinear:

\[
\alpha_t^{(p)}
=
\operatorname{center}\left[
0.22\sin\left\{\frac{2\pi(t-1)}{15}\right\}
+
0.17\exp\left\{
-\frac{1}{2}\left(\frac{t-20}{4}\right)^2
\right\}
\right].
\]

The reporting-delay main effect is the centered version of

\[
(0.22,\;0.14,\;0.05,\;-0.05,\;-0.16,\;-0.27,\;-0.38,\;-0.49),
\]

and the settlement-delay main effect is the centered version of

\[
(1.35,\;1.00,\;0.62,\;0.25,\;-0.10,\;-0.48,\;-0.86,\;-1.24,\;-1.62,\;-2.00).
\]

In the rich scenario, claim types also differ in payment incidence through the
centered type effects

\[
(0.25,\;0.02,\;-0.12,\;-0.38).
\]

### 4. Positive payment amounts

Conditional on a positive payment,

\[
X_{t,d,s,j}\mid Z_{t,d,s,j}=1
\sim
\operatorname{Lognormal}
\left(
m_{t,d,s,j}-\frac{\sigma^2}{2},
\sigma
\right),
\]

with \(\sigma=0.55\). This parameterization implies

\[
E[X_{t,d,s,j}\mid Z_{t,d,s,j}=1]
=
\exp(m_{t,d,s,j}).
\]

The log conditional mean is

\[
m_{t,d,s,j}
=
\log(900)
+\alpha_t^{(X)}
+\beta_d^{(X)}
+\zeta_s^{(X)}
+\tau_j^{(X)}
+h_{j,s}
+\delta_1 C^{(1)}_{t+d+s}
+\delta_2 C^{(2)}_{t+d+s}
+\log n.
\]

The accident-period severity effect is

\[
\alpha_t^{(X)}
=
\operatorname{center}\left[
0.12\sin\left\{\frac{2\pi(t-1)}{20}\right\}
+0.15u_t
\right].
\]

The reporting-delay effect is the centered version of

\[
(-0.10,\;-0.04,\;0.02,\;0.07,\;0.12,\;0.17,\;0.21,\;0.25),
\]

while the settlement-delay effect is generated from

\[
-0.13s
-0.18\log(1+s)
+
0.24\exp\left\{
-\frac{1}{2}\left(\frac{s-2}{1.25}\right)^2
\right\},
\]

and is centered across settlement periods.

The rich scenario also contains claim-type severity effects given by the
centered values

\[
(-0.18,\;0.04,\;0.24,\;0.58).
\]

The term \(h_{j,s}\) introduces type-specific settlement trajectories. These
are generated from four distinct deterministic shapes and then double-centered
across claim type and settlement delay so that they represent interaction
deviations rather than shifts in the corresponding main effects.

Finally, two calendar covariates enter the positive-payment mean:

- a centered linear calendar-time trend representing gradual inflation; and
- a centered temporary-shock indicator active over four consecutive calendar
  periods.

Their coefficients are

\[
(\delta_1,\delta_2)=(0.22,\;0.30).
\]

### 5. Rich and simple scenarios

The **rich scenario** uses all of the structures described above.

The optional **simple scenario** retains the nonlinear accident-time effects,
reporting-delay main effects, settlement-delay main effects, overdispersed
counts, hurdle payments, and unequal baseline claim-type shares, but removes

- evolution of the reporting-delay profile over accident time;
- time variation in claim-type composition;
- reporting-delay selection in claim type;
- claim-type differences in payment incidence;
- claim-type differences in positive-payment levels;
- type-specific settlement trajectories; and
- calendar inflation and temporary-shock effects.

The simple scenario is therefore a reduced-heterogeneity robustness experiment.
It is not intended to make the Baseline Bayesian model the exact
data-generating model.

### 6. Relationship between the DGP and fitted Bayesian models

The simulation is deliberately not a prior-predictive self-recovery exercise.
The deterministic temporal curves above are fixed functions chosen by the
simulation designer; they are not generated by drawing random-walk effects from
the fitted Bayesian models' priors.

The fitted Bayesian models must infer these temporal patterns from the
valuation-date data. In particular, the Bayesian models represent accident-time
effects using regularized first-order random walks rather than being supplied
with the true sinusoidal, hyperbolic-tangent, Gaussian-bump, or linear
functional forms used by the simulation engine.

The Full Bayesian model additionally estimates claim-type composition,
claim-type payment effects, type-specific settlement deviations, and calendar
effects. The Baseline Bayesian model aggregates payments across claim types and
omits these richer payment structures.

This design makes the simulation a test of whether the proposed architecture can
adapt to a structured synthetic claims environment, rather than simply recover
parameters drawn from its own prior distribution.

## Evaluation Design

For each independently generated synthetic portfolio, only information
observable at the valuation date is supplied to the fitted reserving models.

Because the data-generating process is known, the simulation permits two
different validation targets.

### Expected-liability recovery

The primary target is the **oracle expected outstanding liability** under the
known DGP.

The RBNS expectation is evaluated exactly conditional on the reported counts.
The IBNR expectation integrates over future unreported claim counts and
claim-type composition using Monte Carlo integration under the known DGP.

For this target, we report

- relative bias;
- mean absolute percentage error (MAPE);
- RMSE; and
- normalized RMSE.

### Posterior predictive performance

A separate future runoff is generated in every replication. For the Bayesian
models we additionally report

- predictive-median MAPE against realized runoff;
- empirical coverage of 95% posterior predictive intervals; and
- median relative predictive interval width.

The distinction is important: oracle-based criteria assess recovery of the
conditional expected liability, whereas realized-runoff criteria additionally
reflect future process variation.

## Main Results: Rich Scenario

The principal experiment consists of **50 independently simulated portfolios**.

### Total Outstanding Liability

| Method | Relative Bias (%) | MAPE vs. Oracle (%) | Normalized RMSE (%) | MAPE vs. Realized Runoff (%) |
|---|---:|---:|---:|---:|
| Full Bayesian | -0.21 | 6.65 | 9.03 | 8.26 |
| Baseline Bayesian | 5.91 | 7.65 | 9.77 | 10.35 |
| Paid Chain Ladder | -1.26 | 15.38 | 22.62 | 14.95 |

The Full Bayesian model is essentially unbiased on average for the oracle
expected total liability. Its MAPE is **6.65%**, compared with **7.65%** for
the Baseline Bayesian model and **15.38%** for Paid Chain Ladder.

The Full Bayesian model also reduces normalized RMSE from **9.77%** under the
Baseline Bayesian model to **9.03%**. Its median relative 95% predictive
interval width is **52.1%**, compared with **70.0%** for the Baseline Bayesian
model.

Both Bayesian specifications contain the realized total runoff in their 95%
posterior predictive intervals in all 50 rich-scenario replications. With only
50 replications, this should not be interpreted as evidence of exact nominal
coverage.

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

The IBNR-count results are nearly identical across the two Bayesian models,
which is expected because they use the same general aggregate count structure.
The additional gains from the Full Bayesian model arise primarily in the
payment layer. In particular, the Full model removes most of the positive RBNS
and total-reserve bias observed under the reduced payment specification while
retaining accurate IBNR-count prediction.

## Additional Results: Simple Scenario

The reduced-heterogeneity scenario was also evaluated over **50 independently
simulated portfolios**.

### Total Outstanding Liability

| Method | Relative Bias (%) | MAPE vs. Oracle (%) | Normalized RMSE (%) | MAPE vs. Realized Runoff (%) |
|---|---:|---:|---:|---:|
| Full Bayesian | -0.13 | 6.02 | 7.05 | 8.35 |
| Baseline Bayesian | 12.54 | 13.24 | 15.54 | 15.89 |
| Paid Chain Ladder | -0.87 | 15.61 | 20.23 | 17.40 |

The Full Bayesian model again recovers the expected total liability with
negligible average bias and the lowest MAPE and normalized RMSE. The purpose of
this scenario is not to make the Baseline Bayesian model correctly specified;
rather, it removes the principal time-varying claim-type and calendar
heterogeneity while retaining a collective hurdle-payment process generated at
the claim-type cell level.

The 95% posterior predictive interval covers realized total runoff in all 50
replications for the Full Bayesian model and in 98% of replications for the
Baseline Bayesian model.

### Component-Level Recovery

| Model | Quantity | Relative Bias (%) | MAPE vs. Oracle (%) | Predictive-Median MAPE vs. Realized (%) | 95% Predictive Coverage (%) |
|---|---|---:|---:|---:|---:|
| Full Bayesian | RBNS | 0.35 | 6.16 | 9.46 | 96 |
| Full Bayesian | IBNR | -1.29 | 7.34 | 13.33 | 100 |
| Full Bayesian | IBNR count | -0.50 | 4.68 | 8.81 | 98 |
| Full Bayesian | Total | -0.13 | 6.02 | 8.47 | 100 |
| Baseline Bayesian | RBNS | 13.18 | 14.00 | 15.75 | 98 |
| Baseline Bayesian | IBNR | 11.03 | 13.14 | 17.20 | 100 |
| Baseline Bayesian | IBNR count | -0.53 | 4.69 | 8.80 | 98 |
| Baseline Bayesian | Total | 12.54 | 13.24 | 14.18 | 98 |

The simple-scenario count results are again almost identical across the two
Bayesian specifications. Differences in reserve performance therefore arise
from the payment model rather than from materially different IBNR-count
predictions.

## Interpretation

The simulation results are not intended to establish universal superiority of
the Full Bayesian specification. Any simulation conclusion is necessarily
conditional on the chosen DGP.

The principal conclusions are narrower:

1. the proposed framework can recover expected RBNS, IBNR, and total
   outstanding liabilities with low average bias under a claims process
   containing several simultaneous forms of temporal and claim-type
   heterogeneity;
2. the IBNR-count component is recovered accurately under both Bayesian
   specifications;
3. richer payment structure materially improves reserve recovery when the
   synthetic process contains claim-type and calendar structure; and
4. the framework continues to perform well in the reduced-heterogeneity
   scenario.

The simulation therefore complements, rather than replaces, the rolling
out-of-sample empirical validation reported in the paper.

## Convergence and Computation

For the 50 rich-scenario replications:

- the Full Bayesian model produced zero divergent transitions in total;
- the Baseline Bayesian model produced zero divergent transitions in total;
- the maximum recorded R-hat was approximately 1.024 for the Full Bayesian
  model and 1.015 for the Baseline Bayesian model.

For the 50 simple-scenario replications:

- both Bayesian models again produced zero divergent transitions in total;
- the maximum recorded R-hat was approximately 1.018 for the Full Bayesian
  model and 1.014 for the Baseline Bayesian model.

The final fitting profile uses four chains per Bayesian model, with 500 warm-up
and 500 retained sampling iterations per chain.

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

To run the reduced-heterogeneity scenario:

```bash
bash run_simulation.sh 50 simple 2 final
```

The summary files are written to

```text
results/rich/summary/total_reserve_comparison.csv
results/rich/summary/component_recovery.csv
```

with corresponding outputs under `results/simple/summary/`.

## Reproducibility

Replication-specific data-generating and MCMC seeds are deterministic functions
of the replication index. The repository contains the code required to

- generate both synthetic scenarios;
- construct the observed, RBNS, and IBNR regions at the valuation date;
- calculate the oracle expected liability under the known DGP;
- fit the Full and Baseline Bayesian models;
- calculate the Paid Chain Ladder benchmark; and
- reproduce the summary statistics reported above.
