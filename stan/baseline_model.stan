data {
  int<lower=1> T;
  int<lower=2> D;
  int<lower=2> S;
  vector<lower=0>[T] exposure;
  vector[T] time_scaled;

  int<lower=1> N_count_obs;
  array[N_count_obs] int<lower=1, upper=T> t_count_obs;
  array[N_count_obs] int<lower=1, upper=D> d_count_obs;
  array[N_count_obs] int<lower=0> y_total_obs;

  int<lower=0> N_count_pred;
  array[N_count_pred] int<lower=1, upper=T> t_count_pred;
  array[N_count_pred] int<lower=1, upper=D> d_count_pred;

  int<lower=1> N_sev_obs;
  array[N_sev_obs] int<lower=1, upper=T> t_sev_obs;
  array[N_sev_obs] int<lower=1, upper=D> d_sev_obs;
  array[N_sev_obs] int<lower=1, upper=S> s_sev_obs;
  vector<lower=0>[N_sev_obs] payment_obs;
  array[N_sev_obs] int<lower=0, upper=1> z_obs;
  vector[N_sev_obs] log_n_obs;

  int<lower=0> N_rbns;
  array[N_rbns] int<lower=1, upper=T> t_rbns;
  array[N_rbns] int<lower=1, upper=D> d_rbns;
  array[N_rbns] int<lower=1, upper=S> s_rbns;
  vector[N_rbns] log_n_rbns;
}
parameters {
  real mu_count;
  vector[T] alpha_count_raw;
  real<lower=0> sigma_alpha_count;
  vector[D] beta_count_raw;
  vector[D] report_slope_raw;
  real<lower=0> phi_count;

  real mu_prob;
  vector[T] alpha_prob_raw;
  real<lower=0> sigma_alpha_prob;
  vector[D] beta_prob_raw;
  vector[S] zeta_prob_raw;
  real theta_log_n_prob;

  real mu_sev;
  vector[T] alpha_sev_raw;
  real<lower=0> sigma_alpha_sev;
  vector[D] beta_sev_raw;
  vector[S] zeta_sev_raw;
  real<lower=0> sigma_lognormal;
}
transformed parameters {
  vector[T] alpha_count;
  vector[T] alpha_prob;
  vector[T] alpha_sev;
  vector[D] beta_count = beta_count_raw - rep_vector(mean(beta_count_raw), D);
  vector[D] report_slope = report_slope_raw - rep_vector(mean(report_slope_raw), D);
  vector[D] beta_prob = beta_prob_raw - rep_vector(mean(beta_prob_raw), D);
  vector[D] beta_sev = beta_sev_raw - rep_vector(mean(beta_sev_raw), D);
  vector[S] zeta_prob = zeta_prob_raw - rep_vector(mean(zeta_prob_raw), S);
  vector[S] zeta_sev = zeta_sev_raw - rep_vector(mean(zeta_sev_raw), S);

  alpha_count[1] = sigma_alpha_count * alpha_count_raw[1];
  alpha_prob[1] = sigma_alpha_prob * alpha_prob_raw[1];
  alpha_sev[1] = sigma_alpha_sev * alpha_sev_raw[1];
  for (t in 2:T) {
    alpha_count[t] = alpha_count[t - 1] + sigma_alpha_count * alpha_count_raw[t];
    alpha_prob[t] = alpha_prob[t - 1] + sigma_alpha_prob * alpha_prob_raw[t];
    alpha_sev[t] = alpha_sev[t - 1] + sigma_alpha_sev * alpha_sev_raw[t];
  }
  alpha_count -= rep_vector(mean(alpha_count), T);
  alpha_prob -= rep_vector(mean(alpha_prob), T);
  alpha_sev -= rep_vector(mean(alpha_sev), T);
}
model {
  mu_count ~ normal(log(40), 1.5);
  alpha_count_raw ~ std_normal();
  sigma_alpha_count ~ normal(0, 0.35);
  beta_count_raw ~ normal(0, 1.25);
  report_slope_raw ~ normal(0, 0.35);
  phi_count ~ gamma(2, 0.05);

  mu_prob ~ normal(-2, 1.5);
  alpha_prob_raw ~ std_normal();
  sigma_alpha_prob ~ normal(0, 0.35);
  beta_prob_raw ~ normal(0, 0.75);
  zeta_prob_raw ~ normal(0, 1);
  theta_log_n_prob ~ normal(0.4, 0.35);

  mu_sev ~ normal(log(1000), 1.5);
  alpha_sev_raw ~ std_normal();
  sigma_alpha_sev ~ normal(0, 0.35);
  beta_sev_raw ~ normal(0, 0.75);
  zeta_sev_raw ~ normal(0, 1);
  sigma_lognormal ~ normal(0, 0.75);

  for (n in 1:N_count_obs) {
    int tt = t_count_obs[n];
    int dd = d_count_obs[n];
    real log_mu = log(exposure[tt]) + mu_count + alpha_count[tt] +
      beta_count[dd] + report_slope[dd] * time_scaled[tt];
    y_total_obs[n] ~ neg_binomial_2(exp(log_mu), phi_count);
  }

  for (n in 1:N_sev_obs) {
    int tt = t_sev_obs[n];
    int dd = d_sev_obs[n];
    int ss = s_sev_obs[n];
    real eta_prob = mu_prob + alpha_prob[tt] + beta_prob[dd] + zeta_prob[ss] +
      theta_log_n_prob * log_n_obs[n];
    z_obs[n] ~ bernoulli_logit(eta_prob);
    if (z_obs[n] == 1) {
      real log_mean_positive = mu_sev + alpha_sev[tt] + beta_sev[dd] + zeta_sev[ss] +
        log_n_obs[n];
      payment_obs[n] ~ lognormal(log_mean_positive - 0.5 * square(sigma_lognormal), sigma_lognormal);
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

  for (n in 1:N_rbns) {
    int tt = t_rbns[n];
    int dd = d_rbns[n];
    int ss = s_rbns[n];
    real eta_prob = mu_prob + alpha_prob[tt] + beta_prob[dd] + zeta_prob[ss] +
      theta_log_n_prob * log_n_rbns[n];
    real log_mean_positive = mu_sev + alpha_sev[tt] + beta_sev[dd] + zeta_sev[ss] + log_n_rbns[n];
    real p = inv_logit(eta_prob);
    rbns_expected_draw += p * exp(log_mean_positive);
    if (bernoulli_rng(p))
      rbns_reserve_draw += lognormal_rng(log_mean_positive - 0.5 * square(sigma_lognormal), sigma_lognormal);
  }

  for (n in 1:N_count_pred) {
    int tt = t_count_pred[n];
    int dd = d_count_pred[n];
    real log_mu = log(exposure[tt]) + mu_count + alpha_count[tt] +
      beta_count[dd] + report_slope[dd] * time_scaled[tt];
    int n_total = neg_binomial_2_rng(exp(log_mu), phi_count);
    ibnr_count_draw += n_total;
    if (n_total > 0) {
      real log_n = log(n_total);
      for (ss in 1:S) {
        real eta_prob = mu_prob + alpha_prob[tt] + beta_prob[dd] + zeta_prob[ss] +
          theta_log_n_prob * log_n;
        real log_mean_positive = mu_sev + alpha_sev[tt] + beta_sev[dd] + zeta_sev[ss] + log_n;
        real p = inv_logit(eta_prob);
        ibnr_expected_given_counts_draw += p * exp(log_mean_positive);
        if (bernoulli_rng(p))
          ibnr_reserve_draw += lognormal_rng(log_mean_positive - 0.5 * square(sigma_lognormal), sigma_lognormal);
      }
    }
  }
  total_reserve_draw = rbns_reserve_draw + ibnr_reserve_draw;
}
