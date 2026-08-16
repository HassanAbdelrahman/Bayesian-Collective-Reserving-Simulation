center_vector <- function(x) x - mean(x)

softmax_r <- function(x) {
  z <- exp(x - max(x))
  z / sum(z)
}

double_center_matrix <- function(x) {
  x - rowMeans(x) - rep(colMeans(x), each = nrow(x)) + mean(x)
}

simulation_dimensions <- function() {
  list(T = 36L, D = 8L, S = 10L, J = 4L, valuation = 36L)
}

make_dgp_parameters <- function(scenario = c("rich", "simple")) {
  scenario <- match.arg(scenario)
  dims <- simulation_dimensions()
  T <- dims$T; D <- dims$D; S <- dims$S; J <- dims$J
  t_index <- seq_len(T)
  d0 <- 0:(D - 1L)
  s0 <- 0:(S - 1L)
  time_scaled <- seq(-1, 1, length.out = T)

  exposure <- exp(
    0.08 * sin(2 * pi * (t_index - 1) / 12) +
      0.04 * time_scaled
  )

  alpha_count <- center_vector(
    0.16 * sin(2 * pi * (t_index - 1) / 18) +
      0.13 * tanh((t_index - 17) / 4)
  )
  reporting_prob <- c(0.58, 0.20, 0.09, 0.05, 0.035, 0.020, 0.015, 0.010)
  reporting_prob <- reporting_prob / sum(reporting_prob)
  beta_count <- center_vector(log(reporting_prob))
  report_slope <- center_vector(c(0.24, 0.11, 0.03, -0.04, -0.10, -0.15, -0.20, -0.25))
  if (scenario == "simple") report_slope[] <- 0

  target_monthly_claims <- 300
  mu_count <- log(target_monthly_claims / sum(exp(beta_count)))
  phi_count <- 35

  base_share <- c(0.49, 0.28, 0.195, 0.035)
  comp_intercept <- log(base_share[-1] / base_share[[1]])
  comp_time_slope <- c(-0.18, 0.30, 0.12)
  delay_fraction <- d0 / max(d0)
  comp_delay_tilt <- matrix(0, D, J - 1L)
  comp_delay_tilt[, 1] <- 0.22 * delay_fraction
  comp_delay_tilt[, 2] <- -0.18 * delay_fraction
  comp_delay_tilt[, 3] <- 0.50 * delay_fraction
  comp_delay_tilt[1, ] <- 0
  if (scenario == "simple") {
    comp_time_slope[] <- 0
    comp_delay_tilt[,] <- 0
  }

  alpha_prob <- center_vector(
    0.22 * sin(2 * pi * (t_index - 1) / 15) +
      0.17 * exp(-0.5 * ((t_index - 20) / 4)^2)
  )
  beta_prob <- center_vector(c(0.22, 0.14, 0.05, -0.05, -0.16, -0.27, -0.38, -0.49))
  zeta_prob <- center_vector(c(1.35, 1.00, 0.62, 0.25, -0.10, -0.48, -0.86, -1.24, -1.62, -2.00))
  type_prob <- center_vector(c(0.25, 0.02, -0.12, -0.38))
  mu_prob <- -2.25
  theta_log_n_prob <- 0.48

  alpha_sev <- center_vector(
    0.12 * sin(2 * pi * (t_index - 1) / 20) + 0.15 * time_scaled
  )
  beta_sev <- center_vector(c(-0.10, -0.04, 0.02, 0.07, 0.12, 0.17, 0.21, 0.25))
  zeta_sev <- center_vector(
    -0.13 * s0 - 0.18 * log1p(s0) +
      0.24 * exp(-0.5 * ((s0 - 2) / 1.25)^2)
  )
  type_sev <- center_vector(c(-0.18, 0.04, 0.24, 0.58))

  raw_settle_type <- rbind(
    seq(0.28, -0.28, length.out = S),
    0.06 * sin(2 * pi * s0 / max(1, S - 1)),
    seq(-0.22, 0.25, length.out = S),
    0.38 * exp(-0.5 * ((s0 - 4) / 1.5)^2) - 0.08
  )
  settle_type <- double_center_matrix(raw_settle_type)

  C <- T + D + S - 2L
  calendar_index <- seq_len(C)
  calendar_linear <- center_vector(seq(-1, 1, length.out = C))
  calendar_shock <- center_vector(as.numeric(calendar_index >= 16 & calendar_index <= 19))
  calendar_cov <- cbind(
    linear_inflation = calendar_linear,
    temporary_shock = calendar_shock
  )
  beta_calendar <- c(0.22, 0.30)

  if (scenario == "simple") {
    type_prob[] <- 0
    type_sev[] <- 0
    settle_type[,] <- 0
    beta_calendar[] <- 0
  }

  list(
    dims = dims,
    exposure = exposure,
    time_scaled = time_scaled,
    mu_count = mu_count,
    phi_count = phi_count,
    alpha_count = alpha_count,
    beta_count = beta_count,
    report_slope = report_slope,
    comp_intercept = comp_intercept,
    comp_time_slope = comp_time_slope,
    comp_delay_tilt = comp_delay_tilt,
    mu_prob = mu_prob,
    alpha_prob = alpha_prob,
    beta_prob = beta_prob,
    zeta_prob = zeta_prob,
    type_prob = type_prob,
    theta_log_n_prob = theta_log_n_prob,
    mu_sev = log(900),
    alpha_sev = alpha_sev,
    beta_sev = beta_sev,
    zeta_sev = zeta_sev,
    type_sev = type_sev,
    settle_type = settle_type,
    calendar_cov = calendar_cov,
    beta_calendar = beta_calendar,
    sigma_lognormal = 0.55
  )
}

count_mean_true <- function(par, t, d) {
  par$exposure[[t]] * exp(
    par$mu_count + par$alpha_count[[t]] + par$beta_count[[d]] +
      par$report_slope[[d]] * par$time_scaled[[t]]
  )
}

composition_prob_true <- function(par, t, d) {
  logits <- c(
    0,
    par$comp_intercept +
      par$comp_time_slope * par$time_scaled[[t]] +
      par$comp_delay_tilt[d, ]
  )
  softmax_r(logits)
}

payment_moments_true <- function(par, t, d, s, j, n) {
  if (n <= 0) return(c(prob = 0, mean_positive = 0, mean = 0))
  cc <- t + d + s - 2L
  log_n <- log(n)
  eta_prob <- par$mu_prob + par$alpha_prob[[t]] + par$beta_prob[[d]] +
    par$zeta_prob[[s]] + par$type_prob[[j]] +
    par$theta_log_n_prob * log_n
  log_mean_positive <- par$mu_sev + par$alpha_sev[[t]] + par$beta_sev[[d]] +
    par$zeta_sev[[s]] + par$type_sev[[j]] + par$settle_type[j, s] +
    sum(par$calendar_cov[cc, ] * par$beta_calendar) + log_n
  p <- plogis(eta_prob)
  mp <- exp(log_mean_positive)
  c(prob = p, mean_positive = mp, mean = p * mp)
}

simulate_payment_true <- function(par, t, d, s, j, n) {
  mom <- payment_moments_true(par, t, d, s, j, n)
  if (runif(1) > mom[["prob"]]) return(0)
  rlnorm(
    1,
    meanlog = log(mom[["mean_positive"]]) - 0.5 * par$sigma_lognormal^2,
    sdlog = par$sigma_lognormal
  )
}

oracle_expected_reserve <- function(sim, mc_draws = 5000L, seed = 99173L) {
  par <- sim$truth
  dims <- par$dims
  valuation <- dims$valuation
  S <- dims$S
  J <- dims$J

  # Exact conditional mean for RBNS: reported counts are observed at valuation.
  rbns_expected <- 0
  if (nrow(sim$rbns_grid) > 0) {
    for (i in seq_len(nrow(sim$rbns_grid))) {
      r <- sim$rbns_grid[i, ]
      rbns_expected <- rbns_expected + payment_moments_true(
        par, r$t, r$d, r$s, r$j, r$N
      )[["mean"]]
    }
  }

  # IBNR expected liability integrates over the unknown future counts and type mix.
  set.seed(seed + sim$replication_id)
  ibnr_mc <- numeric(mc_draws)
  count_mc <- numeric(mc_draws)
  pred <- sim$count_pred

  for (m in seq_len(mc_draws)) {
    acc_loss <- 0
    acc_count <- 0
    if (nrow(pred) > 0) {
      for (i in seq_len(nrow(pred))) {
        tt <- pred$t[[i]]
        dd <- pred$d[[i]]
        n_total <- rnbinom(1, mu = count_mean_true(par, tt, dd), size = par$phi_count)
        acc_count <- acc_count + n_total
        if (n_total > 0) {
          n_type <- as.integer(rmultinom(1, n_total, composition_prob_true(par, tt, dd))[, 1])
          for (jj in seq_len(J)) {
            if (n_type[[jj]] > 0) {
              for (ss in seq_len(S)) {
                acc_loss <- acc_loss + payment_moments_true(
                  par, tt, dd, ss, jj, n_type[[jj]]
                )[["mean"]]
              }
            }
          }
        }
      }
    }
    ibnr_mc[[m]] <- acc_loss
    count_mc[[m]] <- acc_count
  }

  c(
    RBNS = rbns_expected,
    IBNR = mean(ibnr_mc),
    Total = rbns_expected + mean(ibnr_mc),
    IBNR_count = mean(count_mc)
  )
}

make_full_stan_data <- function(sim) {
  par <- sim$truth
  dims <- par$dims
  list(
    T = dims$T, D = dims$D, S = dims$S, J = dims$J,
    C = nrow(par$calendar_cov), P_cal = ncol(par$calendar_cov),
    exposure = as.numeric(par$exposure),
    time_scaled = as.numeric(par$time_scaled),
    calendar_cov = unname(par$calendar_cov),
    N_count_obs = nrow(sim$count_obs),
    t_count_obs = as.integer(sim$count_obs$t),
    d_count_obs = as.integer(sim$count_obs$d),
    y_total_obs = as.integer(sim$count_obs$N_total),
    y_type_obs = unname(as.matrix(sim$count_obs[, paste0("N_type_", seq_len(dims$J)), drop = FALSE])),
    N_count_pred = nrow(sim$count_pred),
    t_count_pred = as.integer(sim$count_pred$t),
    d_count_pred = as.integer(sim$count_pred$d),
    N_sev_obs = nrow(sim$sev_obs),
    t_sev_obs = as.integer(sim$sev_obs$t),
    d_sev_obs = as.integer(sim$sev_obs$d),
    s_sev_obs = as.integer(sim$sev_obs$s),
    j_sev_obs = as.integer(sim$sev_obs$j),
    payment_obs = as.numeric(sim$sev_obs$X),
    z_obs = as.integer(sim$sev_obs$X > 0),
    log_n_obs = as.numeric(log(sim$sev_obs$N)),
    N_rbns = nrow(sim$rbns_grid),
    t_rbns = as.integer(sim$rbns_grid$t),
    d_rbns = as.integer(sim$rbns_grid$d),
    s_rbns = as.integer(sim$rbns_grid$s),
    j_rbns = as.integer(sim$rbns_grid$j),
    log_n_rbns = as.numeric(log(sim$rbns_grid$N))
  )
}

make_baseline_stan_data <- function(sim) {
  dims <- sim$truth$dims
  agg_obs <- aggregate(
    X ~ t + d + s + d0 + s0 + report_month + payment_month,
    data = sim$payment_grid[sim$payment_grid$payment_month <= dims$valuation, ],
    FUN = sum
  )
  counts_lookup <- sim$count_grid[, c("t", "d", "N_total")]
  agg_obs <- merge(agg_obs, counts_lookup, by = c("t", "d"), all.x = TRUE, sort = FALSE)

  rbns_agg <- aggregate(
    X ~ t + d + s + d0 + s0 + report_month + payment_month,
    data = sim$payment_grid[
      sim$payment_grid$report_month <= dims$valuation &
        sim$payment_grid$payment_month > dims$valuation,
    ],
    FUN = sum
  )
  rbns_agg <- merge(rbns_agg, counts_lookup, by = c("t", "d"), all.x = TRUE, sort = FALSE)

  list(
    T = dims$T, D = dims$D, S = dims$S,
    exposure = as.numeric(sim$truth$exposure),
    time_scaled = as.numeric(sim$truth$time_scaled),
    N_count_obs = nrow(sim$count_obs),
    t_count_obs = as.integer(sim$count_obs$t),
    d_count_obs = as.integer(sim$count_obs$d),
    y_total_obs = as.integer(sim$count_obs$N_total),
    N_count_pred = nrow(sim$count_pred),
    t_count_pred = as.integer(sim$count_pred$t),
    d_count_pred = as.integer(sim$count_pred$d),
    N_sev_obs = nrow(agg_obs),
    t_sev_obs = as.integer(agg_obs$t),
    d_sev_obs = as.integer(agg_obs$d),
    s_sev_obs = as.integer(agg_obs$s),
    payment_obs = as.numeric(agg_obs$X),
    z_obs = as.integer(agg_obs$X > 0),
    log_n_obs = as.numeric(log(pmax(agg_obs$N_total, 1))),
    N_rbns = nrow(rbns_agg),
    t_rbns = as.integer(rbns_agg$t),
    d_rbns = as.integer(rbns_agg$d),
    s_rbns = as.integer(rbns_agg$s),
    log_n_rbns = as.numeric(log(pmax(rbns_agg$N_total, 1)))
  )
}

simulate_reserving_data <- function(replication_id, scenario = c("rich", "simple"), seed_base = 812347L,
                                    oracle_mc_draws = 5000L) {
  scenario <- match.arg(scenario)
  par <- make_dgp_parameters(scenario)
  dims <- par$dims
  T <- dims$T; D <- dims$D; S <- dims$S; J <- dims$J; valuation <- dims$valuation
  seed <- as.integer(seed_base + 1009L * replication_id + ifelse(scenario == "simple", 100000L, 0L))
  set.seed(seed)

  count_grid <- expand.grid(t = seq_len(T), d0 = 0:(D - 1L), KEEP.OUT.ATTRS = FALSE)
  count_grid <- count_grid[order(count_grid$t, count_grid$d0), ]
  count_grid$d <- count_grid$d0 + 1L
  count_grid$report_month <- count_grid$t + count_grid$d0
  count_grid$mu_count <- mapply(count_mean_true, MoreArgs = list(par = par),
                                t = count_grid$t, d = count_grid$d)
  count_grid$N_total <- rnbinom(nrow(count_grid), mu = count_grid$mu_count, size = par$phi_count)

  y_type <- matrix(0L, nrow(count_grid), J)
  for (i in seq_len(nrow(count_grid))) {
    y_type[i, ] <- as.integer(rmultinom(
      1, count_grid$N_total[[i]], composition_prob_true(par, count_grid$t[[i]], count_grid$d[[i]])
    )[, 1])
  }
  colnames(y_type) <- paste0("N_type_", seq_len(J))
  count_grid <- cbind(count_grid, y_type)

  type_grid <- count_grid[rep(seq_len(nrow(count_grid)), each = J), c("t", "d0", "d", "report_month")]
  type_grid$j <- rep(seq_len(J), times = nrow(count_grid))
  type_grid$N <- as.integer(as.vector(t(y_type)))
  positive_type_grid <- type_grid[type_grid$N > 0, ]

  payment_grid <- positive_type_grid[rep(seq_len(nrow(positive_type_grid)), each = S), ]
  payment_grid$s0 <- rep(0:(S - 1L), times = nrow(positive_type_grid))
  payment_grid$s <- payment_grid$s0 + 1L
  payment_grid$payment_month <- payment_grid$t + payment_grid$d0 + payment_grid$s0
  payment_grid$X <- 0
  payment_grid$true_mean <- 0
  for (i in seq_len(nrow(payment_grid))) {
    r <- payment_grid[i, ]
    payment_grid$true_mean[[i]] <- payment_moments_true(par, r$t, r$d, r$s, r$j, r$N)[["mean"]]
    payment_grid$X[[i]] <- simulate_payment_true(par, r$t, r$d, r$s, r$j, r$N)
  }

  count_obs <- count_grid[count_grid$report_month <= valuation, ]
  count_pred <- count_grid[count_grid$report_month > valuation, ]
  sev_obs <- payment_grid[payment_grid$payment_month <= valuation, ]
  rbns_grid <- payment_grid[payment_grid$report_month <= valuation & payment_grid$payment_month > valuation, ]
  ibnr_grid <- payment_grid[payment_grid$report_month > valuation, ]

  sim <- list(
    replication_id = replication_id,
    scenario = scenario,
    seed = seed,
    truth = par,
    count_grid = count_grid,
    type_grid = type_grid,
    payment_grid = payment_grid,
    count_obs = count_obs,
    count_pred = count_pred,
    sev_obs = sev_obs,
    rbns_grid = rbns_grid,
    ibnr_grid = ibnr_grid,
    actual_reserve = c(
      RBNS = sum(rbns_grid$X),
      IBNR = sum(ibnr_grid$X),
      Total = sum(rbns_grid$X) + sum(ibnr_grid$X)
    ),
    actual_ibnr_count = sum(count_pred$N_total)
  )

  sim$oracle_expected <- oracle_expected_reserve(
    sim,
    mc_draws = oracle_mc_draws,
    seed = seed_base + 700000L
  )
  sim$full_stan_data <- make_full_stan_data(sim)
  sim$baseline_stan_data <- make_baseline_stan_data(sim)
  sim
}
