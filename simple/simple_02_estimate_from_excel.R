###########################################################################
# Simple estimation on your own panel data (single self-contained file)
#
# This script shows any researcher how to apply the benchmark and the
# factor-aware estimator to their own panel stored in an Excel file. It reads
# the example file shipped with this folder, but you can point it at your own
# Excel file as long as it has the same column layout (see below).
#
# Required Excel columns (one row per country-period):
#   iso3c        country identifier (text)
#   period       time index, 1, 2, 3, ... (integer, equally spaced)
#   <dependent>  the dependent variable (for example log_gdp)
#   <regressor>  one or more regressors (for example log_energy, log_invest)
#
# Run with:  Rscript simple_02_estimate_from_excel.R
# Dependencies: readxl (install.packages("readxl"))
###########################################################################

library(readxl)

# =================================================================
# 1. SETTINGS — edit these three lines to use your own data
# =================================================================
excel_file <- "example_panel_data.xlsx"   # path to your Excel file
dep_var    <- "log_gdp"                    # name of the dependent variable
regressors <- c("log_energy", "log_invest")# names of the regressors

# =================================================================
# 2. READ AND RESHAPE THE PANEL
# =================================================================
df <- as.data.frame(read_excel(excel_file))
cat(sprintf("Loaded %d rows, %d countries, periods %d to %d\n",
            nrow(df), length(unique(df$iso3c)),
            min(df$period), max(df$period)))

# Reshape the long data frame into the (y, x) array format the estimator uses.
# y is N x (T+1); x is N x (T+1) x k. Only countries observed in every period
# are kept, so the panel is balanced.
to_arrays <- function(df, dep, regs) {
  ids   <- sort(unique(df$iso3c))
  times <- sort(unique(df$period))
  # keep only balanced units
  ok <- names(which(table(df$iso3c) == length(times)))
  ids <- ids[ids %in% ok]
  N <- length(ids); Tp1 <- length(times); k <- length(regs)
  y <- matrix(NA_real_, N, Tp1); x <- array(NA_real_, c(N, Tp1, k))
  for (i in seq_along(ids)) {
    sub <- df[df$iso3c == ids[i], ]
    sub <- sub[order(sub$period), ]
    y[i, ] <- sub[[dep]]
    for (j in seq_along(regs)) x[i, , j] <- sub[[regs[j]]]
  }
  list(y = y, x = x, N = N, T_obs = Tp1 - 1, k = k, ids = ids)
}

d <- to_arrays(df, dep_var, regressors)
cat(sprintf("Balanced panel: %d countries x %d periods\n", d$N, d$T_obs + 1))

# =================================================================
# 3. THE TWO ESTIMATORS (same logic as the simulation file)
# =================================================================
# Benchmark: first-difference instrumental variables.
# Factor-aware: adds current and previous cross-sectional averages of the
# dependent variable and the regressors as extra controls, proxying any
# unobserved common factor (such as a global business cycle).

estimate_fd_iv <- function(d, factor_aware = FALSE) {
  y <- d$y; x <- d$x; N <- d$N; T_obs <- d$T_obs; k <- d$k
  ybar <- colMeans(y, na.rm = TRUE)
  xbar <- apply(x, c(2, 3), mean, na.rm = TRUE)   # (T+1) x k

  dep <- reg_lag <- inst_lag <- numeric(0)
  reg_x <- matrix(0, 0, k)
  csa <- matrix(0, 0, 2 * (k + 1))
  for (t in 3:(T_obs + 1)) {
    dep      <- c(dep,      y[, t]     - y[, t - 1])
    reg_lag  <- c(reg_lag,  y[, t - 1] - y[, t - 2])
    inst_lag <- c(inst_lag, y[, t - 2])
    dxt <- x[, t, ] - x[, t - 1, ]
    reg_x <- rbind(reg_x, matrix(dxt, N, k))
    # current and previous CSAs of y and of each regressor
    block <- cbind(rep(ybar[t], N), rep(ybar[t - 1], N))
    for (j in 1:k) block <- cbind(block, rep(xbar[t, j], N), rep(xbar[t - 1, j], N))
    csa <- rbind(csa, block)
  }

  if (factor_aware) {
    X <- cbind(reg_lag, reg_x, csa)
    Z <- cbind(inst_lag, reg_x, csa)
  } else {
    X <- cbind(reg_lag, reg_x)
    Z <- cbind(inst_lag, reg_x)
  }

  p   <- ncol(X)
  ZZi <- solve(crossprod(Z) + 1e-6 * diag(ncol(Z)))
  XZ  <- crossprod(X, Z)
  coef <- solve(XZ %*% ZZi %*% t(XZ) + 1e-8 * diag(p)) %*% XZ %*% ZZi %*% crossprod(Z, dep)
  names_out <- c("gamma (persistence)", regressors)
  setNames(coef[1:(1 + k)], names_out)
}

# =================================================================
# 4. RUN AND REPORT
# =================================================================
cat("\nDependent variable:", dep_var, "\n")
cat("Regressors:", paste(regressors, collapse = ", "), "\n")
cat(strrep("-", 60), "\n")

bench <- estimate_fd_iv(d, factor_aware = FALSE)
facaw <- estimate_fd_iv(d, factor_aware = TRUE)

res <- data.frame(
  parameter    = names(bench),
  benchmark    = round(as.numeric(bench), 4),
  factor_aware = round(as.numeric(facaw), 4)
)
print(res, row.names = FALSE)

write.csv(res, "simple_estimation_results.csv", row.names = FALSE)
cat("\nResults saved to simple_estimation_results.csv\n")
cat("\nInterpretation: if the persistence estimate (gamma) falls noticeably\n")
cat("when moving from the benchmark column to the factor-aware column, your\n")
cat("data likely contain an unobserved common factor that the benchmark\n")
cat("estimator does not handle.\n")
cat("Done.\n")
