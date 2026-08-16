suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
rep_id <- if (length(args) >= 1) as.integer(args[[1]]) else 1L
scenario <- if (length(args) >= 2) args[[2]] else "rich"
if (!scenario %in% c("rich", "simple")) stop("scenario must be rich or simple")

root <- normalizePath(Sys.getenv("SIM_ROOT", unset = "."), mustWork = TRUE)
source(file.path(root, "R", "simulate_dataset.R"))
source(file.path(root, "R", "chain_ladder.R"))

profile <- Sys.getenv("SIM_PROFILE", unset = "pilot")
if (profile == "pilot") {
  chains <- 2L; parallel_chains <- 2L; warmup <- 300L; sampling <- 300L
  adapt_delta <- 0.90; max_treedepth <- 12L; oracle_mc <- 2000L
} else {
  chains <- 4L; parallel_chains <- 4L; warmup <- 500L; sampling <- 500L
  adapt_delta <- 0.95; max_treedepth <- 12L; oracle_mc <- 5000L
}

seed_base <- 812347L
sim <- simulate_reserving_data(
  replication_id = rep_id,
  scenario = scenario,
  seed_base = seed_base,
  oracle_mc_draws = oracle_mc
)

out_dir <- file.path(root, "results", scenario, sprintf("rep_%03d", rep_id))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(sim, file.path(out_dir, "simulated_data.rds"), compress = TRUE)

full_model <- cmdstan_model(file.path(root, "stan", "full_model.stan"), force_recompile = FALSE)
baseline_model <- cmdstan_model(file.path(root, "stan", "baseline_model.stan"), force_recompile = FALSE)

fit_one <- function(model, data, model_name, seed) {
  csv_dir <- file.path(out_dir, paste0("csv_", model_name))
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  fit <- model$sample(
    data = data,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = warmup,
    iter_sampling = sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = 0,
    output_dir = csv_dir
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  vars <- c(
    "rbns_reserve_draw", "ibnr_reserve_draw", "total_reserve_draw",
    "rbns_expected_draw", "ibnr_expected_given_counts_draw", "ibnr_count_draw"
  )
  draws <- as.data.frame(as_draws_matrix(fit$draws(variables = vars)))

  diag_sum <- fit$summary()
  max_rhat <- suppressWarnings(max(diag_sum$rhat, na.rm = TRUE))
  min_ess_bulk <- suppressWarnings(min(diag_sum$ess_bulk, na.rm = TRUE))
  ds <- tryCatch(fit$diagnostic_summary(), error = function(e) NULL)
  divergences <- if (!is.null(ds) && "num_divergent" %in% names(ds)) sum(ds$num_divergent) else NA_integer_

  list(draws = draws, elapsed_seconds = elapsed, max_rhat = max_rhat,
       min_ess_bulk = min_ess_bulk, divergences = divergences)
}

message("Fitting full Bayesian model for replication ", rep_id)
full <- fit_one(full_model, sim$full_stan_data, "full", seed_base + 500000L + rep_id)
message("Fitting baseline Bayesian model for replication ", rep_id)
baseline <- fit_one(baseline_model, sim$baseline_stan_data, "baseline", seed_base + 600000L + rep_id)

cl <- paid_chain_ladder(
  sim$payment_grid,
  valuation = sim$truth$dims$valuation,
  max_dev = sim$truth$dims$D + sim$truth$dims$S - 2L
)

summarize_bayes <- function(obj, model_name) {
  d <- obj$draws
  rbns_exp_draw <- d$rbns_expected_draw
  ibnr_exp_draw <- d$ibnr_expected_given_counts_draw
  total_exp_draw <- rbns_exp_draw + ibnr_exp_draw

  qfun <- function(x) unname(quantile(x, c(.025, .5, .975), na.rm = TRUE))
  make_row <- function(quantity, expected_draw, predictive_draw, oracle, realized) {
    q <- qfun(predictive_draw)
    data.frame(
      replication = rep_id,
      scenario = scenario,
      model = model_name,
      quantity = quantity,
      oracle_expected = oracle,
      realized = realized,
      expected_point = mean(expected_draw, na.rm = TRUE),
      predictive_median = q[[2]],
      pred_lower_95 = q[[1]],
      pred_upper_95 = q[[3]],
      predictive_covered_realized = realized >= q[[1]] && realized <= q[[3]],
      rel_error_expected_vs_oracle = (mean(expected_draw, na.rm = TRUE) - oracle) / oracle,
      rel_error_median_vs_realized = (q[[2]] - realized) / realized,
      predictive_width_relative = (q[[3]] - q[[1]]) / oracle,
      elapsed_seconds = obj$elapsed_seconds,
      max_rhat = obj$max_rhat,
      min_ess_bulk = obj$min_ess_bulk,
      divergences = obj$divergences
    )
  }
  rbind(
    make_row("RBNS", rbns_exp_draw, d$rbns_reserve_draw,
             sim$oracle_expected[["RBNS"]], sim$actual_reserve[["RBNS"]]),
    make_row("IBNR", ibnr_exp_draw, d$ibnr_reserve_draw,
             sim$oracle_expected[["IBNR"]], sim$actual_reserve[["IBNR"]]),
    make_row("Total", total_exp_draw, d$total_reserve_draw,
             sim$oracle_expected[["Total"]], sim$actual_reserve[["Total"]]),
    data.frame(
      replication = rep_id, scenario = scenario, model = model_name, quantity = "IBNR count",
      oracle_expected = sim$oracle_expected[["IBNR_count"]], realized = sim$actual_ibnr_count,
      expected_point = mean(d$ibnr_count_draw), predictive_median = median(d$ibnr_count_draw),
      pred_lower_95 = quantile(d$ibnr_count_draw, .025), pred_upper_95 = quantile(d$ibnr_count_draw, .975),
      predictive_covered_realized = sim$actual_ibnr_count >= quantile(d$ibnr_count_draw, .025) &&
        sim$actual_ibnr_count <= quantile(d$ibnr_count_draw, .975),
      rel_error_expected_vs_oracle = (mean(d$ibnr_count_draw) - sim$oracle_expected[["IBNR_count"]]) /
        sim$oracle_expected[["IBNR_count"]],
      rel_error_median_vs_realized = (median(d$ibnr_count_draw) - sim$actual_ibnr_count) /
        sim$actual_ibnr_count,
      predictive_width_relative = (quantile(d$ibnr_count_draw, .975) - quantile(d$ibnr_count_draw, .025)) /
        sim$oracle_expected[["IBNR_count"]],
      elapsed_seconds = obj$elapsed_seconds, max_rhat = obj$max_rhat,
      min_ess_bulk = obj$min_ess_bulk, divergences = obj$divergences
    )
  )
}

res <- rbind(
  summarize_bayes(full, "Full Bayesian"),
  summarize_bayes(baseline, "Baseline Bayesian"),
  data.frame(
    replication = rep_id,
    scenario = scenario,
    model = "Paid Chain Ladder",
    quantity = "Total",
    oracle_expected = sim$oracle_expected[["Total"]],
    realized = sim$actual_reserve[["Total"]],
    expected_point = cl$reserve,
    predictive_median = NA_real_,
    pred_lower_95 = NA_real_, pred_upper_95 = NA_real_,
    predictive_covered_realized = NA,
    rel_error_expected_vs_oracle = (cl$reserve - sim$oracle_expected[["Total"]]) / sim$oracle_expected[["Total"]],
    rel_error_median_vs_realized = NA_real_,
    predictive_width_relative = NA_real_,
    elapsed_seconds = NA_real_, max_rhat = NA_real_, min_ess_bulk = NA_real_, divergences = NA_integer_
  )
)

write.csv(res, file.path(out_dir, "performance.csv"), row.names = FALSE)
saveRDS(list(sim = sim, full = full, baseline = baseline, chain_ladder = cl, performance = res),
        file.path(out_dir, "replication_result.rds"), compress = TRUE)
message("Saved replication ", rep_id, " to ", out_dir)
