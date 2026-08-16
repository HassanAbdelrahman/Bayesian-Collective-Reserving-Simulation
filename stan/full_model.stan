// ============================================================================
// Illustrative simulation model for a rich collective reserving process
//
// Components:
//   1. Negative-binomial reported claim counts on an (accident, reporting delay)
//      grid, including a smooth accident-time effect and a time-varying reporting
//      delay profile.
//   2. Multinomial claim-type composition with changing portfolio mix and
//      reporting-delay selection.
//   3. Hierarchical hurdle-lognormal aggregate payments on
//      (accident, reporting delay, settlement delay, claim type) cells.
//   4. Type-specific settlement deviations and known calendar covariates.
//   5. Posterior-predictive RBNS, IBNR, and total reserve simulation.
//
// The severity predictor is parameterized in terms of the conditional mean of a
// positive payment. The Lognormal location used in the likelihood is therefore
// log_mean_positive - sigma^2 / 2.
// ============================================================================

data {
  int<lower=1> T;
  int<lower=2> D;
  int<lower=2> S;
  int<lower=2> J;
  int<lower=1> C;
  int<lower=1> P_cal;

  vector<lower=0>[T] exposure;
  vector[T] time_scaled;
  matrix[C, P_cal] calendar_cov;

  // Observed count and composition cells: t + d0 <= valuation
  int<lower=1> N_count_obs;
  array[N_count_obs] int<lower=1, upper=T> t_count_obs;
  array[N_count_obs] int<lower=1, upper=D> d_count_obs;
  array[N_count_obs] int<lower=0> y_total_obs;
  array[N_count_obs, J] int<lower=0> y_type_obs;

  // Unobserved count cells used for IBNR prediction
  int<lower=0> N_count_pred;
  array[N_count_pred] int<lower=1, upper=T> t_count_pred;
  array[N_count_pred] int<lower=1, upper=D> d_count_pred;

  // Observed payment cells
  int<lower=1> N_sev_obs;
  array[N_sev_obs] int<lower=1, upper=T> t_sev_obs;
  array[N_sev_obs] int<lower=1, upper=D> d_sev_obs;
  array[N_sev_obs] int<lower=1, upper=S> s_sev_obs;
  array[N_sev_obs] int<lower=1, upper=J> j_sev_obs;
  vector<lower=0>[N_sev_obs] payment_obs;
  array[N_sev_obs] int<lower=0, upper=1> z_obs;
  vector[N_sev_obs] log_n_obs;

  // Future payment cells for already reported claims
  int<lower=0> N_rbns;
  array[N_rbns] int<lower=1, upper=T> t_rbns;
  array[N_rbns] int<lower=1, upper=D> d_rbns;
  array[N_rbns] int<lower=1, upper=S> s_rbns;
  array[N_rbns] int<lower=1, upper=J> j_rbns;
  vector[N_rbns] log_n_rbns;
}

parameters {
  // Count model
  real mu_count;
  vector[T] alpha_count_raw;
  real<lower=0> sigma_alpha_count;
  vector[D] beta_count_raw;
  vector[D] report_slope_raw;
  real<lower=0> phi_count;

  // Claim-type composition; type 1 and delay 1 are references
  vector[J - 1] comp_intercept;
  vector[J - 1] comp_time_slope;
  matrix[D - 1, J - 1] comp_delay_raw;
  real<lower=0> sigma_comp_delay;

  // Hurdle probability component
  real mu_prob;
  vector[T] alpha_prob_raw;
  real<lower=0> sigma_alpha_prob;
  vector[D] beta_prob_raw;
  vector[S] zeta_prob_raw;
  vector[J] type_prob_raw;
  real<lower=0> sigma_type_prob;
  real theta_log_n_prob;

  // Positive-payment mean component
  real mu_sev;
  vector[T] alpha_sev_raw;
  real<lower=0> sigma_alpha_sev;
  vector[D] beta_sev_raw;
  vector[S] zeta_sev_raw;
  vector[J] type_sev_raw;
  real<lower=0> sigma_type_sev;

  // Hierarchical type-by-settlement deviations
  matrix[J, S] settle_type_raw;
  real<lower=0> sigma_settle_type;

  vector[P_cal] beta_calendar;
  real<lower=0> sigma_lognormal;
}

transformed parameters {
  vector[T] alpha_count;
  vector[T] alpha_prob;
  vector[T] alpha_sev;

  vector[D] beta_count =
    beta_count_raw - rep_vector(mean(beta_count_raw), D);
  vector[D] report_slope =
    report_slope_raw - rep_vector(mean(report_slope_raw), D);

  vector[D] beta_prob =
    beta_prob_raw - rep_vector(mean(beta_prob_raw), D);
  vector[D] beta_sev =
    beta_sev_raw - rep_vector(mean(beta_sev_raw), D);

  vector[S] zeta_prob =
    zeta_prob_raw - rep_vector(mean(zeta_prob_raw), S);
  vector[S] zeta_sev =
    zeta_sev_raw - rep_vector(mean(zeta_sev_raw), S);

  vector[J] type_prob;
  vector[J] type_sev;

  matrix[D, J - 1] comp_delay_tilt;
  matrix[J, S] settle_uncentered =
    sigma_settle_type * settle_type_raw;
  matrix[J, S] settle_type;

  vector[J] settle_row_mean;
  vector[S] settle_col_mean;
  real settle_grand_mean;

  alpha_count[1] =
    sigma_alpha_count * alpha_count_raw[1];
  alpha_prob[1] =
    sigma_alpha_prob * alpha_prob_raw[1];
  alpha_sev[1] =
    sigma_alpha_sev * alpha_sev_raw[1];

  for (t in 2:T) {
    alpha_count[t] =
      alpha_count[t - 1]
      + sigma_alpha_count * alpha_count_raw[t];
    alpha_prob[t] =
      alpha_prob[t - 1]
      + sigma_alpha_prob * alpha_prob_raw[t];
    alpha_sev[t] =
      alpha_sev[t - 1]
      + sigma_alpha_sev * alpha_sev_raw[t];
  }

  alpha_count -= rep_vector(mean(alpha_count), T);
  alpha_prob -= rep_vector(mean(alpha_prob), T);
  alpha_sev -= rep_vector(mean(alpha_sev), T);

  type_prob =
    sigma_type_prob
    * (type_prob_raw - rep_vector(mean(type_prob_raw), J));
  type_sev =
    sigma_type_sev
    * (type_sev_raw - rep_vector(mean(type_sev_raw), J));

  comp_delay_tilt[1] = rep_row_vector(0, J - 1);
  for (d in 2:D) {
    comp_delay_tilt[d] =
      sigma_comp_delay * comp_delay_raw[d - 1];
  }

  for (j in 1:J) {
    real acc = 0;
    for (s in 1:S) {
      acc += settle_uncentered[j, s];
    }
    settle_row_mean[j] = acc / S;
  }

  for (s in 1:S) {
    real acc = 0;
    for (j in 1:J) {
      acc += settle_uncentered[j, s];
    }
    settle_col_mean[s] = acc / J;
  }

  settle_grand_mean = mean(to_vector(settle_uncentered));

  for (j in 1:J) {
    for (s in 1:S) {
      settle_type[j, s] =
        settle_uncentered[j, s]
        - settle_row_mean[j]
        - settle_col_mean[s]
        + settle_grand_mean;
    }
  }
}

model {
  // Count priors
  mu_count ~ normal(log(40), 1.5);
  alpha_count_raw ~ std_normal();
  sigma_alpha_count ~ normal(0, 0.35);
  beta_count_raw ~ normal(0, 1.25);
  report_slope_raw ~ normal(0, 0.35);
  phi_count ~ gamma(2, 0.05);

  // Composition priors
  comp_intercept ~ normal(0, 1.5);
  comp_time_slope ~ normal(0, 0.5);
  to_vector(comp_delay_raw) ~ std_normal();
  sigma_comp_delay ~ normal(0, 0.5);

  // Payment-incidence priors
  mu_prob ~ normal(-2, 1.5);
  alpha_prob_raw ~ std_normal();
  sigma_alpha_prob ~ normal(0, 0.35);
  beta_prob_raw ~ normal(0, 0.75);
  zeta_prob_raw ~ normal(0, 1);
  type_prob_raw ~ std_normal();
  sigma_type_prob ~ normal(0, 0.75);
  theta_log_n_prob ~ normal(0.4, 0.35);

  // Positive-payment priors
  mu_sev ~ normal(log(1000), 1.5);
  alpha_sev_raw ~ std_normal();
  sigma_alpha_sev ~ normal(0, 0.35);
  beta_sev_raw ~ normal(0, 0.75);
  zeta_sev_raw ~ normal(0, 1);
  type_sev_raw ~ std_normal();
  sigma_type_sev ~ normal(0, 0.75);

  to_vector(settle_type_raw) ~ std_normal();
  sigma_settle_type ~ normal(0, 0.5);

  beta_calendar ~ normal(0, 0.5);
  sigma_lognormal ~ normal(0, 0.75);

  // Count and composition likelihood
  for (n in 1:N_count_obs) {
    int tt = t_count_obs[n];
    int dd = d_count_obs[n];
    real log_mu_count =
      log(exposure[tt])
      + mu_count
      + alpha_count[tt]
      + beta_count[dd]
      + report_slope[dd] * time_scaled[tt];
    vector[J] composition_logits;
    vector[J] composition_prob;

    composition_logits[1] = 0;
    for (j in 2:J) {
      composition_logits[j] =
        comp_intercept[j - 1]
        + comp_time_slope[j - 1] * time_scaled[tt]
        + comp_delay_tilt[dd, j - 1];
    }
    composition_prob = softmax(composition_logits);

    target += neg_binomial_2_lpmf(
      y_total_obs[n] | exp(log_mu_count), phi_count
    );
    target += multinomial_lpmf(
      y_type_obs[n] | composition_prob
    );
  }

  // Hurdle-lognormal likelihood
  for (n in 1:N_sev_obs) {
    int tt = t_sev_obs[n];
    int dd = d_sev_obs[n];
    int ss = s_sev_obs[n];
    int jj = j_sev_obs[n];
    int cc = tt + dd + ss - 2;

    real eta_prob =
      mu_prob
      + alpha_prob[tt]
      + beta_prob[dd]
      + zeta_prob[ss]
      + type_prob[jj]
      + theta_log_n_prob * log_n_obs[n];

    target += bernoulli_logit_lpmf(
      z_obs[n] | eta_prob
    );

    if (z_obs[n] == 1) {
      real log_mean_positive =
        mu_sev
        + alpha_sev[tt]
        + beta_sev[dd]
        + zeta_sev[ss]
        + type_sev[jj]
        + settle_type[jj, ss]
        + calendar_cov[cc] * beta_calendar
        + log_n_obs[n];

      real lognormal_location =
        log_mean_positive
        - 0.5 * square(sigma_lognormal);

      target += lognormal_lpdf(
        payment_obs[n]
        | lognormal_location, sigma_lognormal
      );
    }
  }
}

generated quantities {
  real rbns_reserve_draw = 0;
  real ibnr_reserve_draw = 0;
  real total_reserve_draw = 0;

  real rbns_expected_draw = 0;
  real ibnr_expected_given_counts_draw = 0;

  int ibnr_count_draw = 0;

  // RBNS: future payments associated with observed reported counts
  for (n in 1:N_rbns) {
    int tt = t_rbns[n];
    int dd = d_rbns[n];
    int ss = s_rbns[n];
    int jj = j_rbns[n];
    int cc = tt + dd + ss - 2;

    real eta_prob =
      mu_prob
      + alpha_prob[tt]
      + beta_prob[dd]
      + zeta_prob[ss]
      + type_prob[jj]
      + theta_log_n_prob * log_n_rbns[n];

    real log_mean_positive =
      mu_sev
      + alpha_sev[tt]
      + beta_sev[dd]
      + zeta_sev[ss]
      + type_sev[jj]
      + settle_type[jj, ss]
      + calendar_cov[cc] * beta_calendar
      + log_n_rbns[n];

    real p_payment = inv_logit(eta_prob);

    rbns_expected_draw +=
      p_payment * exp(log_mean_positive);

    if (bernoulli_rng(p_payment) == 1) {
      rbns_reserve_draw +=
        lognormal_rng(
          log_mean_positive
          - 0.5 * square(sigma_lognormal),
          sigma_lognormal
        );
    }
  }

  // IBNR: first simulate unreported counts and type composition, then payments
  for (n in 1:N_count_pred) {
    int tt = t_count_pred[n];
    int dd = d_count_pred[n];

    real log_mu_count =
      log(exposure[tt])
      + mu_count
      + alpha_count[tt]
      + beta_count[dd]
      + report_slope[dd] * time_scaled[tt];

    vector[J] composition_logits;
    vector[J] composition_prob;
    int n_total;
    array[J] int n_type;

    composition_logits[1] = 0;
    for (j in 2:J) {
      composition_logits[j] =
        comp_intercept[j - 1]
        + comp_time_slope[j - 1] * time_scaled[tt]
        + comp_delay_tilt[dd, j - 1];
    }
    composition_prob = softmax(composition_logits);

    n_total =
      neg_binomial_2_rng(
        exp(log_mu_count), phi_count
      );
    n_type = multinomial_rng(composition_prob, n_total);
    ibnr_count_draw += n_total;

    for (jj in 1:J) {
      if (n_type[jj] > 0) {
        real log_n = log(n_type[jj]);

        for (ss in 1:S) {
          int cc = tt + dd + ss - 2;

          real eta_prob =
            mu_prob
            + alpha_prob[tt]
            + beta_prob[dd]
            + zeta_prob[ss]
            + type_prob[jj]
            + theta_log_n_prob * log_n;

          real log_mean_positive =
            mu_sev
            + alpha_sev[tt]
            + beta_sev[dd]
            + zeta_sev[ss]
            + type_sev[jj]
            + settle_type[jj, ss]
            + calendar_cov[cc] * beta_calendar
            + log_n;

          real p_payment = inv_logit(eta_prob);

          ibnr_expected_given_counts_draw +=
            p_payment * exp(log_mean_positive);

          if (bernoulli_rng(p_payment) == 1) {
            ibnr_reserve_draw +=
              lognormal_rng(
                log_mean_positive
                - 0.5 * square(sigma_lognormal),
                sigma_lognormal
              );
          }
        }
      }
    }
  }

  total_reserve_draw =
    rbns_reserve_draw + ibnr_reserve_draw;
}
