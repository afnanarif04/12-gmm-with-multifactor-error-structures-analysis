###########################################################################
# Simple Monte Carlo simulation (single self-contained file)
#
# A teaching and replication version. It contains everything in one file:
# the data generating process, the benchmark estimator, the factor-aware
# estimator, and a short Monte Carlo loop that compares them. It is kept
# compact so that any researcher can read it top to bottom and reproduce
# the central simulation finding in a minute or two.
#
# Central finding: when the data contain an unobserved common factor, the
# standard first-difference estimator is biased upward in the persistence
# parameter, and the factor-aware estimator removes most of that bias.
#
# Run with:  Rscript simple_01_monte_carlo.R
# Dependencies: base R only.
###########################################################################

set.seed(2024)

# =================================================================
# 1. DATA GENERATING PROCESS
# =================================================================
# Dynamic panel:
#   y_it = gamma * y_{i,t-1} + beta * x_it + alpha_i + eta_i * f_t + e_it
# The term eta_i * f_t is the common factor. It is included only when m = 1
# and is the single source of cross-sectional dependence.

simulate_panel <- function(N, T_obs, gamma, beta, m = 0) {
  alpha <- rnorm(N)                                  # fixed effects
  x     <- matrix(rnorm(N * (T_obs + 1)), N, T_obs + 1)
  f     <- as.numeric(arima.sim(list(ar = 0.5), n = T_obs + 1))  # common factor
  eta   <- rnorm(N, mean = 1)                        # factor loadings
  y     <- matrix(0, N, T_obs + 1)
  y[, 1] <- alpha + rnorm(N)                         # initial condition
  for (t in 2:(T_obs + 1)) {
    y[, t] <- gamma * y[, t - 1] + beta * x[, t] + alpha + rnorm(N)
    if (m == 1) y[, t] <- y[, t] + eta * f[t]
  }
  list(y = y, x = x, N = N, T_obs = T_obs)
}

# =================================================================
# 2. FIRST-DIFFERENCE INSTRUMENTAL VARIABLE ESTIMATOR
# =================================================================
# Differences out the fixed effect, then stacks all periods into one system.
# Instruments: the twice-lagged level of y for the differenced lagged
# dependent variable, and the differenced regressor instruments itself.
#
# Setting factor_aware = TRUE adds the cross-sectional averages of y and x at
# both the current and the previous period as extra regressors and
# instruments. Including both periods means the differenced common factor is
# spanned by the cross-sectional averages, which proxies the factor and
# restores consistency.

estimate_fd_iv <- function(d, factor_aware = FALSE) {
  y <- d$y; x <- d$x; N <- d$N; T_obs <- d$T_obs
  ybar <- colMeans(y); xbar <- colMeans(x)

  dep <- reg_lag <- reg_x <- inst_lag <- numeric(0)
  csa_y0 <- csa_y1 <- csa_x0 <- csa_x1 <- numeric(0)
  for (t in 3:(T_obs + 1)) {
    dep      <- c(dep,      y[, t]     - y[, t - 1])   # dy_t
    reg_lag  <- c(reg_lag,  y[, t - 1] - y[, t - 2])   # dy_{t-1}
    reg_x    <- c(reg_x,    x[, t]     - x[, t - 1])   # dx_t
    inst_lag <- c(inst_lag, y[, t - 2])                # instrument for dy_{t-1}
    csa_y0   <- c(csa_y0,   rep(ybar[t],     N))       # CSA of y, current period
    csa_y1   <- c(csa_y1,   rep(ybar[t - 1], N))       # CSA of y, previous period
    csa_x0   <- c(csa_x0,   rep(xbar[t],     N))       # CSA of x, current period
    csa_x1   <- c(csa_x1,   rep(xbar[t - 1], N))       # CSA of x, previous period
  }

  if (factor_aware) {
    X <- cbind(reg_lag, reg_x, csa_y0, csa_y1, csa_x0, csa_x1)
    Z <- cbind(inst_lag, reg_x, csa_y0, csa_y1, csa_x0, csa_x1)
  } else {
    X <- cbind(reg_lag, reg_x)
    Z <- cbind(inst_lag, reg_x)
  }

  # Two-stage least squares with a small ridge term for numerical stability.
  p   <- ncol(X)
  ZZi <- solve(crossprod(Z) + 1e-6 * diag(ncol(Z)))
  XZ  <- crossprod(X, Z)
  coef <- solve(XZ %*% ZZi %*% t(XZ) + 1e-8 * diag(p)) %*% XZ %*% ZZi %*% crossprod(Z, dep)
  c(gamma = coef[1], beta = coef[2])
}

estimate_benchmark    <- function(d) estimate_fd_iv(d, factor_aware = FALSE)
estimate_factor_aware <- function(d) estimate_fd_iv(d, factor_aware = TRUE)

# =================================================================
# 3. MONTE CARLO LOOP
# =================================================================
# The single twice-lagged instrument is weak in a small share of draws, which
# can produce an extreme estimate. We therefore report the median estimate,
# the median bias, and the median absolute error, all of which are robust to
# such draws and convey the central comparison cleanly.

run_simulation <- function(R = 500, N = 200, T_obs = 6,
                            gamma = 0.5, beta = 1.0, m = 1) {
  bench <- matrix(NA, R, 2); facaw <- matrix(NA, R, 2)
  for (r in 1:R) {
    d <- simulate_panel(N, T_obs, gamma, beta, m = m)
    bench[r, ] <- tryCatch(estimate_benchmark(d),    error = function(e) c(NA, NA))
    facaw[r, ] <- tryCatch(estimate_factor_aware(d), error = function(e) c(NA, NA))
  }
  summarise <- function(est, truth) {
    est <- est[is.finite(est)]
    c(median = median(est), bias = median(est) - truth, mae = median(abs(est - truth)))
  }
  list(bench_gamma = summarise(bench[, 1], gamma),
       facaw_gamma = summarise(facaw[, 1], gamma),
       bench_beta  = summarise(bench[, 2], beta),
       facaw_beta  = summarise(facaw[, 2], beta))
}

# =================================================================
# 4. RUN AND REPORT
# =================================================================
cat("Monte Carlo: dynamic panel with one unobserved common factor\n")
cat("True gamma = 0.5, true beta = 1.0, N = 200, T = 6, replications = 500\n")
cat("Reported: median estimate, median bias, median absolute error\n")
cat(strrep("-", 64), "\n")

res <- run_simulation(R = 500, N = 200, T_obs = 6, gamma = 0.5, beta = 1.0, m = 1)

cat("\nPersistence parameter gamma (true value 0.5):\n")
cat(sprintf("  Benchmark (ignores factor):  median=%.3f  bias=%+.3f  mae=%.3f\n",
            res$bench_gamma["median"], res$bench_gamma["bias"], res$bench_gamma["mae"]))
cat(sprintf("  Factor-aware estimator:      median=%.3f  bias=%+.3f  mae=%.3f\n",
            res$facaw_gamma["median"], res$facaw_gamma["bias"], res$facaw_gamma["mae"]))

cat("\nSlope parameter beta (true value 1.0):\n")
cat(sprintf("  Benchmark (ignores factor):  median=%.3f  bias=%+.3f  mae=%.3f\n",
            res$bench_beta["median"], res$bench_beta["bias"], res$bench_beta["mae"]))
cat(sprintf("  Factor-aware estimator:      median=%.3f  bias=%+.3f  mae=%.3f\n",
            res$facaw_beta["median"], res$facaw_beta["bias"], res$facaw_beta["mae"]))

out <- data.frame(
  estimator  = c("Benchmark", "Factor-aware"),
  gamma_bias = c(res$bench_gamma["bias"], res$facaw_gamma["bias"]),
  gamma_mae  = c(res$bench_gamma["mae"],  res$facaw_gamma["mae"]),
  beta_bias  = c(res$bench_beta["bias"],  res$facaw_beta["bias"]),
  beta_mae   = c(res$bench_beta["mae"],   res$facaw_beta["mae"])
)
write.csv(out, "simple_monte_carlo_results.csv", row.names = FALSE)
cat("\nResults saved to simple_monte_carlo_results.csv\n")
cat("Done.\n")
