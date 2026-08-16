paid_chain_ladder <- function(payment_grid, valuation, max_dev = NULL) {
  x <- payment_grid
  x$dev <- x$d0 + x$s0
  if (is.null(max_dev)) max_dev <- max(x$dev)

  inc <- matrix(NA_real_, nrow = valuation, ncol = max_dev + 1L)
  rownames(inc) <- seq_len(valuation)
  colnames(inc) <- 0:max_dev

  obs <- aggregate(X ~ t + dev, data = x[x$payment_month <= valuation, ], FUN = sum)
  for (tt in seq_len(valuation)) {
    max_obs_dev <- min(max_dev, valuation - tt)
    if (max_obs_dev >= 0) inc[tt, seq_len(max_obs_dev + 1L)] <- 0
  }
  for (i in seq_len(nrow(obs))) {
    if (obs$t[[i]] <= valuation && obs$dev[[i]] <= max_dev) {
      inc[obs$t[[i]], obs$dev[[i]] + 1L] <- obs$X[[i]]
    }
  }

  cum <- inc
  for (tt in seq_len(nrow(cum))) {
    seen <- which(!is.na(cum[tt, ]))
    if (length(seen)) cum[tt, seen] <- cumsum(cum[tt, seen])
  }

  factors <- rep(1, max_dev)
  if (max_dev >= 1) {
    for (k in seq_len(max_dev)) {
      eligible <- which(!is.na(cum[, k]) & !is.na(cum[, k + 1L]))
      den <- sum(cum[eligible, k], na.rm = TRUE)
      num <- sum(cum[eligible, k + 1L], na.rm = TRUE)
      factors[[k]] <- if (is.finite(den) && den > 0 && is.finite(num)) num / den else 1
      if (!is.finite(factors[[k]]) || factors[[k]] < 1) factors[[k]] <- 1
    }
  }

  latest <- numeric(valuation)
  ultimate <- numeric(valuation)
  latest_dev <- integer(valuation)
  for (tt in seq_len(valuation)) {
    seen <- which(!is.na(cum[tt, ]))
    if (!length(seen)) next
    jj <- max(seen)
    latest[[tt]] <- cum[tt, jj]
    latest_dev[[tt]] <- jj - 1L
    tail_factor <- if (jj <= max_dev) prod(factors[jj:max_dev]) else 1
    ultimate[[tt]] <- latest[[tt]] * tail_factor
  }

  reserve <- sum(ultimate - latest)
  list(
    reserve = reserve,
    factors = factors,
    incremental = inc,
    cumulative = cum,
    latest = latest,
    ultimate = ultimate,
    latest_dev = latest_dev
  )
}
