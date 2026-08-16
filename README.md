# Simulation Study for a Bayesian Collective Reserving Framework

This repository contains the supplementary simulation study for the paper

**A Unified Bayesian Collective Reserving Framework with Hierarchical Temporal Smoothing**.

The purpose of the simulation is to provide a controlled assessment of the
proposed reserving architecture beyond the empirical insurance portfolio used
in the paper.

The study is not intended to establish that one fixed Bayesian specification is
universally optimal. Instead, it evaluates whether the proposed collective
Bayesian framework can recover expected outstanding liabilities and generate
reliable predictive distributions when the underlying claims process contains
several simultaneous sources of temporal, development, and claim-type
heterogeneity.

## Methods Compared

Three reserving approaches are evaluated.

1. **Full Bayesian simulation model**  
   A simulation analogue of the paper's proposed Bayesian collective reserving
   architecture. It jointly models total claim counts, claim-type composition,
   payment occurrence, and positive payment amounts, with temporal
   regularization, hierarchical claim-type effects, type-specific settlement
   development, and calendar effects.

2. **Baseline Bayesian simulation model**  
   Uses the same total-count model, but aggregates payments across claim types
   and omits claim-type payment effects, type-specific settlement development,
   and calendar effects.

3. **Paid Chain Ladder**  
   An aggregate paid-loss Chain Ladder fitted to an accident-period by total
   development-period triangle.

The Full Bayesian simulation model is intentionally **not an exact copy of the
final empirical specification used in the paper**. It preserves the core
count--composition--payment architecture and posterior reserve propagation, but
uses related, deliberately simplified structural choices. The Baseline Bayesian
model provides a reduced Bayesian benchmark with the same general count
machinery but substantially less payment structure. Paid Chain Ladder provides
a conventional aggregate reserving benchmark.

## Relationship to the Empirical Paper Specification

The purpose of the simulation is to assess the proposed **modeling
architecture**, rather than to generate data from the exact application-specific
model and then refit that same specification.

The Full Bayesian simulation model and the final empirical model share the
following central features:

- a Negative Binomial model for reported claim counts with an exposure offset
  and regularized accident-period dynamics;
- a multinomial decomposition of predicted claims across claim types;
- separate Bernoulli payment-incidence and Lognormal positive-payment
  components;
- hierarchical or partially pooled claim-type effects;
- explicit modeling of reporting and settlement development;
- posterior propagation from claim counts through claim type and payments; and
- separate posterior predictions of RBNS, IBNR, and total outstanding
  liabilities.

The exact structural implementations differ in several respects.

### Claim-count layer

In the empirical paper, the claim-count model uses a first-order random walk for
the accident-period effect, a hybrid reporting-delay representation with a
spline-regularized tail, calendar-related covariates acting on the reporting
profile, and an additional tail-only accident-period effect for long reporting
delays.

In the simulation model, the accident-period effect is also estimated through a
first-order random walk, but the reporting-delay structure is represented more
simply through centered delay effects together with delay-specific linear
changes over scaled accident time.

### Claim-type composition layer

In the final empirical model, baseline claim-type log-odds are decomposed into
persistent type-specific levels and first-order random-walk temporal deviations.
Reporting-delay selection is then represented through a regularized
body-and-tail delay-tilt structure.

In the simulation model, claim-type composition uses type-specific baseline
log-odds, a linear accident-time change in relative type prevalence, and
delay-specific tilts relative to reporting delay zero. Thus, the simulation
preserves changing claim-type composition and delay selection without using the
exact random-walk and body-tail parameterization of the empirical model.

### Payment layer

Both models use a hurdle Lognormal payment model. The empirical specification
is more application-specific: it allows hierarchical claim-type deviations in
the payment predictors, uses a hybrid fixed--spline settlement-delay structure,
includes payment-month seasonality, and contains a smooth binned
reporting-delay by settlement-delay interaction in positive-payment severity.

The simulation model instead uses regularized global accident-period effects,
centered claim-type payment effects, a type-by-settlement interaction, and two
known calendar covariates representing gradual inflation and a temporary
shock. It does not reproduce the paper's exact payment-month seasonal term,
hybrid settlement spline, or binned reporting--payment interaction.

### Dimension and purpose

The synthetic experiment uses 36 accident periods, 8 reporting delays,
10 settlement delays, and 4 claim types so that repeated full posterior fits are
computationally feasible. The empirical application is larger and contains six
claim types.

These differences are deliberate. Several components of the final empirical
specification---such as the hybrid settlement spline, the payment-month seasonal
effect, and the particular reporting-delay interaction---were selected to
represent features of the observed insurance portfolio. Reproducing all of
those exact choices in the simulation DGP would make the exercise closer to a
self-recovery experiment. Instead, the simulation retains the principal
Bayesian collective reserving architecture while changing the functional form
of several important components. This provides a controlled assessment of
whether the architecture can adapt to synthetic claims processes that are
related to, but not identical to, the empirical specification.

# Simulation Design

## Indexing and valuation structure

Each synthetic portfolio contains:

- 36 accident periods;
- 8 reporting-delay periods;
- 10 settlement-delay periods; and
- 4 claim types.

Let $t=1,\ldots,36$ denote accident period, $d=0,\ldots,7$ reporting delay,
$s=0,\ldots,9$ settlement delay, and $j=1,\ldots,4$ claim type. The valuation
date is the end of accident period 36.

A claim-count cell is observed when

```math
t+d \le 36
```

and is IBNR when

```math
t+d > 36.
```

A payment cell is observed when

```math
t+d+s \le 36.
```

Future payments from already reported claims form the RBNS reserve. Future
payments associated with claims not yet reported at valuation form the IBNR
reserve. Because the complete future runoff is simulated, realized RBNS, IBNR,
and total outstanding liabilities are known in every replication.

## 1. Exposure and reported claim counts

### Exposure

Define a scaled accident-time index

```math
u_t = -1 + \frac{2(t-1)}{35}.
```

Exposure varies smoothly over accident time according to

```math
e_t =
\exp\left[
0.08\sin\left(\frac{2\pi(t-1)}{12}\right)
+0.04u_t
\right].
```

Thus exposure contains a modest seasonal component together with a gradual
time trend.

### Total reported claim counts

Reported claim counts are generated from an overdispersed count distribution:

```math
N_{t,d} \sim \mathrm{NB2}(\lambda_{t,d},\phi_N),
\qquad
\phi_N = 35.
```

The conditional mean is

```math
\log(\lambda_{t,d})
=
\log(e_t)
+\mu_N
+\alpha_t^{(N)}
+\beta_d^{(N)}
+\kappa_d u_t.
```

The intercept $\mu_N$ is calibrated so that the portfolio contains
approximately 300 claims per accident period before temporal variation.

### Nonlinear accident-period frequency pattern

The raw accident-period frequency curve is

```math
f_t^{(N)}
=
0.16\sin\left(\frac{2\pi(t-1)}{18}\right)
+
0.13\tanh\left(\frac{t-17}{4}\right).
```

It is centered before entering the count predictor:

```math
\alpha_t^{(N)}
=
f_t^{(N)}
-
\frac{1}{36}
\sum_{r=1}^{36} f_r^{(N)}.
```

This produces a frequency path containing both cyclical movement and a smooth
level transition. Importantly, the path is a fixed deterministic function in
the simulation engine; it is not generated from the random-walk prior used by
the fitted Bayesian model.

### Reporting-delay structure

The baseline reporting-delay probability vector is proportional to

```math
(0.58,\;0.20,\;0.09,\;0.05,\;0.035,\;0.020,\;0.015,\;0.010).
```

If these normalized probabilities are denoted by $q_d$, the reporting-delay
main effects are

```math
\beta_d^{(N)}
=
\log(q_d)
-
\frac{1}{8}
\sum_{r=0}^{7}\log(q_r).
```

In the rich scenario, the reporting-delay profile also changes with accident
time. Start from

```math
\kappa^{*}
=
(0.24,\;0.11,\;0.03,\;-0.04,\;-0.10,\;-0.15,\;-0.20,\;-0.25),
```

and center the vector:

```math
\kappa_d
=
\kappa_d^{*}
-
\frac{1}{8}
\sum_{r=0}^{7}\kappa_r^{*}.
```

The interaction term $\kappa_d u_t$ therefore makes early and late reporting
delays evolve differently over accident time.

## 2. Claim-type composition

Conditional on the total number of claims in an accident-reporting cell, claim
types are generated from

```math
(N_{t,d,1},\ldots,N_{t,d,4})
\mid N_{t,d}
\sim
\mathrm{Multinomial}
\left(
N_{t,d},
\pi_{t,d,1},\ldots,\pi_{t,d,4}
\right).
```

The baseline type-share vector is

```math
(0.49,\;0.28,\;0.195,\;0.035).
```

Thus the fourth claim type is deliberately sparse.

Type 1 is the reference category. For $j=2,3,4$,

```math
\log\left(
\frac{\pi_{t,d,j}}{\pi_{t,d,1}}
\right)
=
\gamma_j
+
a_j u_t
+
c_{d,j}.
```

The baseline log-odds $\gamma_j$ are chosen to reproduce the baseline shares
above. The accident-time slopes are

```math
(a_2,a_3,a_4)
=
(-0.18,\;0.30,\;0.12).
```

Hence the relative prevalence of the claim types changes over accident time.

Reporting delay also affects claim-type composition. The delay tilt for each
non-reference type increases linearly from zero at $d=0$ to its terminal value
at $d=7$. The terminal values are

```math
(c_{7,2},c_{7,3},c_{7,4})
=
(0.22,\;-0.18,\;0.50).
```

Thus the rich scenario contains both time-varying claim-type composition and
reporting-delay selection in claim type.

## 3. Payment occurrence

Payments are generated on the $(t,d,s,j)$ grid. For a cell containing $n$
claims, let $Z_{t,d,s,j}$ indicate whether the aggregate payment in that cell is
positive.

```math
Z_{t,d,s,j}
\sim
\mathrm{Bernoulli}(p_{t,d,s,j}).
```

The payment-incidence probability satisfies

```math
\mathrm{logit}(p_{t,d,s,j})
=
\mu_p
+
\alpha_t^{(p)}
+
\beta_d^{(p)}
+
\zeta_s^{(p)}
+
\tau_j^{(p)}
+
0.48\log(n),
```

with

```math
\mu_p=-2.25.
```

### Nonlinear accident-period payment-incidence effect

Define

```math
f_t^{(p)}
=
0.22\sin\left(\frac{2\pi(t-1)}{15}\right)
+
0.17\exp\left[
-\frac{1}{2}
\left(\frac{t-20}{4}\right)^2
\right].
```

The centered effect used in the predictor is

```math
\alpha_t^{(p)}
=
f_t^{(p)}
-
\frac{1}{36}
\sum_{r=1}^{36}f_r^{(p)}.
```

This combines a cyclical pattern with a localized temporary increase in payment
incidence.

### Reporting- and settlement-delay effects

The reporting-delay effect is obtained by centering

```math
(0.22,\;0.14,\;0.05,\;-0.05,\;-0.16,\;-0.27,\;-0.38,\;-0.49).
```

The settlement-delay effect is obtained by centering

```math
(1.35,\;1.00,\;0.62,\;0.25,\;-0.10,\;-0.48,\;-0.86,\;-1.24,\;-1.62,\;-2.00).
```

Consequently, the probability of a positive payment generally decreases with
settlement delay.

In the rich scenario, claim types also differ in payment incidence. The
claim-type effect is the centered version of

```math
(0.25,\;0.02,\;-0.12,\;-0.38).
```

## 4. Positive payment amounts

Conditional on a positive aggregate payment,

```math
X_{t,d,s,j}
\mid
Z_{t,d,s,j}=1
\sim
\mathrm{Lognormal}
\left(
m_{t,d,s,j}-\frac{\sigma^2}{2},
\sigma
\right),
```

where

```math
\sigma=0.55.
```

Under this parameterization,

```math
E
\left[
X_{t,d,s,j}
\mid
Z_{t,d,s,j}=1
\right]
=
\exp(m_{t,d,s,j}).
```

The log conditional mean is

```math
m_{t,d,s,j}
=
\log(900)
+
\alpha_t^{(X)}
+
\beta_d^{(X)}
+
\zeta_s^{(X)}
+
\tau_j^{(X)}
+
h_{j,s}
+
\delta_1 C^{(1)}_{t+d+s}
+
\delta_2 C^{(2)}_{t+d+s}
+
\log(n).
```

The final $\log(n)$ term makes the conditional mean aggregate payment scale
with the number of claims in the cell.

### Accident-period severity effect

Define

```math
f_t^{(X)}
=
0.12\sin\left(\frac{2\pi(t-1)}{20}\right)
+
0.15u_t.
```

The centered effect is

```math
\alpha_t^{(X)}
=
f_t^{(X)}
-
\frac{1}{36}
\sum_{r=1}^{36}f_r^{(X)}.
```

### Reporting-delay severity effect

The reporting-delay severity effect is the centered version of

```math
(-0.10,\;-0.04,\;0.02,\;0.07,\;0.12,\;0.17,\;0.21,\;0.25).
```

### Settlement-delay severity effect

For settlement delay $s$, define

```math
g_s
=
-0.13s
-
0.18\log(1+s)
+
0.24\exp\left[
-\frac{1}{2}
\left(\frac{s-2}{1.25}\right)^2
\right].
```

The effect entering the predictor is

```math
\zeta_s^{(X)}
=
g_s
-
\frac{1}{10}
\sum_{r=0}^{9}g_r.
```

This produces an overall declining settlement profile together with a modest
early-settlement bump.

### Claim-type effects and type-specific settlement development

In the rich scenario, the claim-type positive-payment effects are obtained by
centering

```math
(-0.18,\;0.04,\;0.24,\;0.58).
```

The interaction $h_{j,s}$ introduces four different claim-type settlement
shapes. The raw deterministic type-specific settlement curves are
double-centered across claim type and settlement delay before entering the
positive-payment predictor. This makes $h_{j,s}$ an interaction deviation rather
than another claim-type or settlement-delay main effect.

### Calendar effects

Let $k=t+d+s$ denote calendar period.

Two calendar variables affect positive payment size:

1. a centered linear calendar-time index, representing gradual inflation; and
2. a centered indicator for a temporary four-period calendar shock.

Their coefficients are

```math
(\delta_1,\delta_2)
=
(0.22,\;0.30).
```

Thus positive payments are affected both by gradual calendar inflation and by a
temporary external disturbance.

## 5. Rich and simple scenarios

The **rich scenario** uses all structures described above.

The **simple scenario** retains:

- nonlinear accident-time count effects;
- nonlinear accident-time payment effects;
- reporting-delay main effects;
- settlement-delay main effects;
- overdispersed claim counts;
- the hurdle payment mechanism; and
- unequal baseline claim-type shares.

It removes:

- evolution of the reporting-delay profile over accident time;
- time variation in claim-type composition;
- reporting-delay selection in claim type;
- claim-type differences in payment incidence;
- claim-type differences in positive-payment levels;
- type-specific settlement trajectories; and
- calendar inflation and temporary-shock effects.

The simple scenario is therefore a reduced-heterogeneity robustness experiment.
It is not constructed to make the Baseline Bayesian model exactly equal to the
data-generating mechanism.

## 6. Why the DGP is not the fitted simulation model

The simulation is deliberately not a prior-predictive self-recovery exercise.

The temporal curves in the simulation engine are fixed deterministic functions.
For example, the accident-period count process combines a sine curve and a
hyperbolic-tangent transition, while the payment-incidence process combines a
sine curve and a localized Gaussian-shaped bump.

The fitted Bayesian simulation models are not given these functional forms.
Instead, they must learn the temporal structure from the observed
valuation-date data using regularized first-order random-walk effects.

Similarly, the Full Bayesian simulation model must estimate claim-type
composition, claim-type payment effects, type-specific settlement deviations,
and calendar effects from the observed data rather than being supplied with
their generating values.

The simulation therefore evaluates whether a simulation implementation of the
proposed architecture can adapt to a structured synthetic claims environment,
rather than merely recover parameters generated from its own priors or from an
identical fitted specification.

# Evaluation Design

For each independently generated synthetic portfolio, only information
observable at the valuation date is supplied to the fitted reserving models.

Because the DGP is known, two distinct validation targets are available.

## Expected-liability recovery

The primary target is the **oracle expected outstanding liability** under the
known DGP.

The RBNS expectation is calculated conditional on the observed reported claim
counts. The IBNR expectation integrates over future unreported claim counts and
claim-type composition using Monte Carlo integration under the known DGP.

For oracle expected liability, we report:

- relative bias;
- mean absolute percentage error (MAPE);
- RMSE; and
- normalized RMSE.

## Posterior predictive performance

A separate future runoff is simulated in each replication. For the Bayesian
models we also report:

- predictive-median MAPE against realized runoff;
- empirical coverage of 95% posterior predictive intervals; and
- median relative predictive interval width.

Oracle-based criteria assess recovery of the conditional expected liability.
Realized-runoff criteria additionally contain future process variation.

# Main Results: Rich Scenario

The principal experiment consists of **50 independently simulated portfolios**.

## Total Outstanding Liability

| Method | Relative Bias (%) | MAPE vs. Oracle (%) | Normalized RMSE (%) | Expected-Reserve MAPE vs. Realized Runoff (%) |
|---|---:|---:|---:|---:|
| Full Bayesian | -0.21 | 6.65 | 9.03 | 8.26 |
| Baseline Bayesian | 5.91 | 7.65 | 9.77 | 10.35 |
| Paid Chain Ladder | -1.26 | 15.38 | 22.62 | 14.95 |

The Full Bayesian simulation model is essentially unbiased on average for the oracle
expected total liability. Its oracle MAPE is **6.65%**, compared with **7.65%**
for the Baseline Bayesian model and **15.38%** for Paid Chain Ladder.

The Full Bayesian simulation model also reduces normalized RMSE from **9.77%** under the
Baseline Bayesian model to **9.03%**. Its median relative 95% predictive
interval width is **52.1%**, compared with **70.0%** for the Baseline Bayesian
model.

Both Bayesian specifications contain the realized total runoff in their 95%
posterior predictive intervals in all 50 rich-scenario replications. With 50
replications, this high empirical coverage should not be interpreted as
evidence of exact nominal calibration.

## Component-Level Recovery

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

The IBNR-count results are nearly identical across the two Bayesian
specifications, as expected because they use the same general aggregate count
structure.

The richer Full Bayesian payment model primarily improves the monetary reserve
components. It reduces RBNS oracle MAPE from **7.79%** to **6.37%** and removes
most of the positive RBNS bias. For IBNR, oracle MAPE is similar under the two
models (**10.51%** versus **10.20%**), but the Full model reduces relative bias
from **6.29%** to **-0.42%**, lowers RMSE, and improves predictive-median MAPE
against realized runoff from **15.39%** to **12.81%**.

Thus the simulation does not show uniform dominance on every component metric,
but the Full specification provides the strongest overall recovery of the total
outstanding liability.

# Additional Results: Simple Scenario

The reduced-heterogeneity scenario was also evaluated using **50 independently
simulated portfolios**.

## Total Outstanding Liability

| Method | Relative Bias (%) | MAPE vs. Oracle (%) | Normalized RMSE (%) | Expected-Reserve MAPE vs. Realized Runoff (%) |
|---|---:|---:|---:|---:|
| Full Bayesian | -0.13 | 6.02 | 7.05 | 8.35 |
| Baseline Bayesian | 12.54 | 13.24 | 15.54 | 15.89 |
| Paid Chain Ladder | -0.87 | 15.61 | 20.23 | 17.40 |

The Full Bayesian simulation model again recovers the expected total liability with
negligible average bias and the lowest oracle MAPE and normalized RMSE.

The purpose of this scenario is not to make the Baseline Bayesian model
correctly specified. Instead, it removes the principal time-varying claim-type
and calendar heterogeneity while retaining a collective hurdle-payment process
generated on the claim-type cell grid.

The 95% posterior predictive interval contains the realized total runoff in all
50 replications for the Full Bayesian model and in 98% of replications for the
Baseline Bayesian model.

## Component-Level Recovery

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
from the payment layer rather than from materially different IBNR-count
predictions.

# Interpretation

The simulation is not intended to establish universal superiority of the Full
Bayesian specification. Any simulation conclusion is necessarily conditional on
the chosen DGP.

The principal conclusions are narrower:

1. the proposed framework can recover expected RBNS, IBNR, and total
   outstanding liabilities with low average bias under a synthetic claims
   process containing simultaneous temporal, development, claim-type, and
   calendar structure;
2. IBNR claim counts are recovered accurately under both Bayesian
   specifications;
3. the richer payment architecture improves overall reserve recovery,
   particularly for RBNS and total outstanding liabilities; and
4. the framework continues to perform well in the reduced-heterogeneity
   scenario.

The simulation therefore complements, rather than replaces, the rolling
out-of-sample empirical validation reported in the paper.

# Convergence and Computation

For the 50 rich-scenario replications:

- the Full Bayesian model produced zero divergent transitions in total;
- the Baseline Bayesian model produced zero divergent transitions in total;
- the maximum recorded R-hat was approximately 1.024 for the Full Bayesian
  model and 1.015 for the Baseline Bayesian model.

For the 50 simple-scenario replications:

- both Bayesian models produced zero divergent transitions in total;
- the maximum recorded R-hat was approximately 1.018 for the Full Bayesian
  model and 1.014 for the Baseline Bayesian model.

The final fitting profile uses four chains per Bayesian model, with 500 warm-up
and 500 retained sampling iterations per chain.

# Repository Structure

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

# Reproducing the Simulation

The analysis requires R, `cmdstanr`, `posterior`, `dplyr`, and a working CmdStan installation.

To run the rich scenario with 50 replications:

```bash
bash run_simulation.sh 50 rich 2 final
```

To run the simple scenario with 50 replications:

```bash
bash run_simulation.sh 50 simple 2 final
```

The summary files are written to:

```text
results/rich/summary/total_reserve_comparison.csv
results/rich/summary/component_recovery.csv
```

with corresponding outputs under `results/simple/summary/`.

# Reproducibility

Replication-specific data-generating and MCMC seeds are deterministic functions
of the replication index.

The repository contains the code required to:

- generate both synthetic scenarios;
- construct the observed, RBNS, and IBNR regions at the valuation date;
- calculate the oracle expected liability under the known DGP;
- fit the Full and Baseline Bayesian simulation models;
- calculate the Paid Chain Ladder benchmark; and
- reproduce the summary statistics reported above.
