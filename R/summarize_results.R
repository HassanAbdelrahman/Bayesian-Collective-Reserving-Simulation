suppressPackageStartupMessages(library(dplyr))

root <- normalizePath(Sys.getenv("SIM_ROOT", unset = "."), mustWork = TRUE)
scenario <- Sys.getenv("SIM_SCENARIO", unset = "rich")
files <- list.files(file.path(root, "results", scenario), pattern = "performance\\.csv$",
                    recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No performance.csv files found")
all <- bind_rows(lapply(files, read.csv))

summary_expected <- all %>%
  filter(quantity == "Total") %>%
  group_by(model) %>%
  summarise(
    n_replications = n(),
    relative_bias_percent = 100 * mean(rel_error_expected_vs_oracle, na.rm = TRUE),
    MAPE_oracle_percent = 100 * mean(abs(rel_error_expected_vs_oracle), na.rm = TRUE),
    RMSE_oracle = sqrt(mean((expected_point - oracle_expected)^2, na.rm = TRUE)),
    normalized_RMSE_percent = 100 * RMSE_oracle / mean(oracle_expected, na.rm = TRUE),
    MAPE_realized_percent = 100 * mean(abs(expected_point - realized) / pmax(abs(realized), 1), na.rm = TRUE),
    predictive_coverage_95_percent = 100 * mean(predictive_covered_realized, na.rm = TRUE),
    median_relative_predictive_width_percent = 100 * median(predictive_width_relative, na.rm = TRUE),
    mean_elapsed_seconds = mean(elapsed_seconds, na.rm = TRUE),
    max_rhat = suppressWarnings(max(max_rhat, na.rm = TRUE)),
    total_divergences = sum(divergences, na.rm = TRUE),
    .groups = "drop"
  )

summary_components <- all %>%
  filter(model %in% c("Full Bayesian", "Baseline Bayesian")) %>%
  group_by(model, quantity) %>%
  summarise(
    n_replications = n(),
    relative_bias_percent = 100 * mean(rel_error_expected_vs_oracle, na.rm = TRUE),
    MAPE_oracle_percent = 100 * mean(abs(rel_error_expected_vs_oracle), na.rm = TRUE),
    RMSE_oracle = sqrt(mean((expected_point - oracle_expected)^2, na.rm = TRUE)),
    predictive_median_MAPE_realized_percent = 100 * mean(abs(rel_error_median_vs_realized), na.rm = TRUE),
    predictive_coverage_95_percent = 100 * mean(predictive_covered_realized, na.rm = TRUE),
    median_relative_predictive_width_percent = 100 * median(predictive_width_relative, na.rm = TRUE),
    .groups = "drop"
  )

out <- file.path(root, "results", scenario, "summary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(all, file.path(out, "all_replications.csv"), row.names = FALSE)
write.csv(summary_expected, file.path(out, "total_reserve_comparison.csv"), row.names = FALSE)
write.csv(summary_components, file.path(out, "component_recovery.csv"), row.names = FALSE)

cat("\n=== TOTAL RESERVE: ORACLE EXPECTED LIABILITY ===\n")
print(summary_expected)
cat("\n=== BAYESIAN COMPONENT RECOVERY ===\n")
print(summary_components)
