###########################################################################
# Estimator 3 of 5: GMM with common correlated effects and instrumental variables
#
# Requires 02_estimator_ab_gmm.R in the same directory.
#
# Run from the command line with:  Rscript <this_file>.R
# Or open in an editor, set the working directory to the code/ folder,
# and run interactively.
#
# Dependencies: base R and the MASS package only (install.packages("MASS")).
###########################################################################

e <- new.env()
sys.source("02_estimator_ab_gmm.R", envir = e)
generate_dgp     <- e$generate_dgp
first_difference  <- e$first_difference
arellano_bond_gmm <- e$arellano_bond_gmm

library(MASS)

# =================================================================
# 1. GMM-CCE-IV
# =================================================================
gmm_cce_iv <- function(data, step = 2) {
  y <- data$y; x <- data$x
  N <- data$N; k <- data$k; T_obs <- data$T_obs
  fd <- first_difference(data)
  dy <- fd$dy; dx <- fd$dx
  n_eq <- T_obs - 1; p_s <- 1 + k

  dep  <- if (n_eq == 1) matrix(dy[, 2], ncol=1) else dy[, 2:T_obs]
  dlag <- if (n_eq == 1) matrix(dy[, 1], ncol=1) else dy[, 1:(T_obs-1)]
  dxs  <- array(0, c(N, n_eq, k))
  for (j in 1:k) {
    sl <- if (n_eq == 1) matrix(dx[, 2, j], ncol=1) else dx[, 2:T_obs, j]
    dxs[, , j] <- sl
  }

  # Level cross-sectional averages
  ybar <- colMeans(y)                     # (T+1) vector
  xbar <- apply(x, c(2, 3), mean)        # (T+1) x k matrix

  # CSA control matrix: [ybar_t, xbar_t, ybar_{t-1}, xbar_{t-1}]
  k_w <- 1 + k; n_csa <- k_w * 2; p_aug <- p_s + n_csa
  W_ctrl <- matrix(0, n_eq, n_csa)
  for (t_idx in 1:n_eq) {
    t_act <- t_idx + 2  # 1-indexed in y
    W_ctrl[t_idx, 1] <- ybar[t_act]
    for (j in 1:k) W_ctrl[t_idx, 1 + j] <- xbar[t_act, j]
    W_ctrl[t_idx, k_w + 1] <- ybar[t_act - 1]
    for (j in 1:k) W_ctrl[t_idx, k_w + 1 + j] <- xbar[t_act - 1, j]
  }

  # Regressors: [dlag, dxs, W_ctrl]
  reg <- array(0, c(N, n_eq, p_aug))
  reg[, , 1] <- dlag
  for (j in 1:k) reg[, , 1 + j] <- dxs[, , j]
  for (i in 1:N) reg[i, , (p_s+1):p_aug] <- W_ctrl

  # Instruments: [y_{i,t-2}, dxs, W_ctrl]
  n_iv_pe <- 1 + k + n_csa; n_iv <- n_eq * n_iv_pe
  Z <- array(0, c(N, n_eq, n_iv))
  for (t_idx in 1:n_eq) {
    t_act <- t_idx + 2
    cb <- (t_idx - 1) * n_iv_pe
    Z[, t_idx, cb + 1] <- y[, max(t_act - 2, 1)]
    for (j in 1:k) Z[, t_idx, cb + 1 + j] <- dxs[, t_idx, j]
    for (i in 1:N) Z[i, t_idx, (cb + 1 + k + 1):(cb + 1 + k + n_csa)] <- W_ctrl[t_idx, ]
  }

  # GMM
  XZ <- matrix(0, p_aug, n_iv); Zy <- rep(0, n_iv)
  for (i in 1:N) {
    ri <- matrix(reg[i, , ], n_eq, p_aug)
    zi <- matrix(Z[i, , ], n_eq, n_iv)
    XZ <- XZ + t(ri) %*% zi
    Zy <- Zy + t(zi) %*% dep[i, ]
  }

  H <- 2 * diag(n_eq)
  if (n_eq > 1) for (ii in 1:(n_eq-1)) { H[ii, ii+1] <- -1; H[ii+1, ii] <- -1 }
  ZHZ <- matrix(0, n_iv, n_iv)
  for (i in 1:N) {
    zi <- matrix(Z[i, , ], n_eq, n_iv)
    ZHZ <- ZHZ + t(zi) %*% H %*% zi
  }
  ZHZ <- ZHZ / N
  W1 <- tryCatch(solve(ZHZ + 1e-6 * diag(n_iv)), error = function(e) ginv(ZHZ))
  th1 <- tryCatch(solve(XZ %*% W1 %*% t(XZ) + 1e-10 * diag(p_aug), XZ %*% W1 %*% Zy),
                   error = function(e) ginv(XZ %*% W1 %*% t(XZ)) %*% (XZ %*% W1 %*% Zy))

  if (step >= 2) {
    ZeZ <- matrix(0, n_iv, n_iv)
    for (i in 1:N) {
      ri <- matrix(reg[i, , ], n_eq, p_aug)
      zi <- matrix(Z[i, , ], n_eq, n_iv)
      ei <- dep[i, ] - ri %*% th1
      ze <- t(zi) %*% ei; ZeZ <- ZeZ + ze %*% t(ze)
    }
    ZeZ <- ZeZ / N
    W2 <- tryCatch(solve(ZeZ + 1e-6 * diag(n_iv)), error = function(e) ginv(ZeZ))
    th <- tryCatch(solve(XZ %*% W2 %*% t(XZ) + 1e-10 * diag(p_aug), XZ %*% W2 %*% Zy),
                   error = function(e) ginv(XZ %*% W2 %*% t(XZ)) %*% (XZ %*% W2 %*% Zy))
    Wf <- W2
  } else { th <- th1; Wf <- W1 }

  # Sandwich variance (structural params only)
  ef <- matrix(0, N, n_eq)
  for (i in 1:N) {
    ri <- matrix(reg[i, , ], n_eq, p_aug)
    ef[i, ] <- dep[i, ] - ri %*% th
  }
  An <- XZ %*% Wf %*% t(XZ) / N
  Om <- matrix(0, n_iv, n_iv)
  for (i in 1:N) {
    zi <- matrix(Z[i, , ], n_eq, n_iv)
    ze <- t(zi) %*% ef[i, ]; Om <- Om + ze %*% t(ze)
  }
  Om <- Om / N
  Ai <- tryCatch(solve(An), error = function(e) ginv(An))
  V <- Ai %*% (XZ %*% Wf %*% Om %*% Wf %*% t(XZ) / N^2) %*% Ai
  se <- sqrt(pmax(diag(V), 0))

  list(theta = as.numeric(th[1:p_s]),
       gamma_hat = as.numeric(th[1]),
       beta_hat = as.numeric(th[2:p_s]),
       se = as.numeric(se[1:p_s]),
       method = "GMM-CCE-IV", N = N, T_obs = T_obs)
}

# =================================================================
# 2. STANDALONE TEST
# =================================================================
cat("=============================================================\n")
cat("GMM-3: GMM-CCE-IV — STANDALONE TEST\n")
cat("=============================================================\n")

cat("\n--- Test 1: m=0 -> CCE should be close to AB ---\n")
d0 <- generate_dgp(500, 6, 0.5, list(c(0.8, 1.2)), m = 0, seed = 42)
ab <- arellano_bond_gmm(d0)
cce <- gmm_cce_iv(d0)
cat(sprintf("  AB:  gamma=%.4f, beta=[%.4f, %.4f]\n",
            ab$gamma_hat, ab$beta_hat[1], ab$beta_hat[2]))
cat(sprintf("  CCE: gamma=%.4f, beta=[%.4f, %.4f]\n",
            cce$gamma_hat, cce$beta_hat[1], cce$beta_hat[2]))

cat("\n--- Test 2: m=1 -> CCE should reduce RMSE vs AB ---\n")
d1 <- generate_dgp(500, 6, 0.5, list(c(0.8, 1.2)), m = 1, sigma_eta = 0.5, seed = 42)
ab1 <- arellano_bond_gmm(d1)
cce1 <- gmm_cce_iv(d1)
cat(sprintf("  AB:  gamma=%.4f (bias=%+.4f)\n", ab1$gamma_hat, ab1$gamma_hat - 0.5))
cat(sprintf("  CCE: gamma=%.4f (bias=%+.4f)\n", cce1$gamma_hat, cce1$gamma_hat - 0.5))

cat("\n--- Test 3: m=2 -> CCE adapts (no m needed) ---\n")
d2 <- generate_dgp(500, 8, 0.5, list(c(0.8, 1.2)), m = 2, sigma_eta = 0.5, seed = 42)
ab2 <- arellano_bond_gmm(d2)
cce2 <- gmm_cce_iv(d2)
cat(sprintf("  AB:  gamma=%.4f (bias=%+.4f)\n", ab2$gamma_hat, ab2$gamma_hat - 0.5))
cat(sprintf("  CCE: gamma=%.4f (bias=%+.4f)\n", cce2$gamma_hat, cce2$gamma_hat - 0.5))

cat("\n--- Test 4: Mini MC (R=100, N=200, T=6, m=1) ---\n")
Rmc <- 100; g_ab <- g_cce <- numeric(Rmc)
for (r in 1:Rmc) {
  d <- generate_dgp(200, 6, 0.5, list(c(0.8, 1.2)), m = 1, sigma_eta = 0.5, seed = r)
  g_ab[r]  <- tryCatch(arellano_bond_gmm(d)$gamma_hat, error = function(e) NA)
  g_cce[r] <- tryCatch(gmm_cce_iv(d)$gamma_hat, error = function(e) NA)
}
g_ab <- na.omit(g_ab); g_cce <- na.omit(g_cce)
cat(sprintf("  AB:  Bias=%+.4f, RMSE=%.4f (n=%d)\n",
            mean(g_ab) - 0.5, sqrt(mean((g_ab - 0.5)^2)), length(g_ab)))
cat(sprintf("  CCE: Bias=%+.4f, RMSE=%.4f (n=%d)\n",
            mean(g_cce) - 0.5, sqrt(mean((g_cce - 0.5)^2)), length(g_cce)))
cat(sprintf("  RMSE reduction: %.0f%%\n",
            100 * (1 - sqrt(mean((g_cce - 0.5)^2)) / sqrt(mean((g_ab - 0.5)^2)))))

cat("\n--- Test 5: Strong persistence (gamma=0.8, m=1) ---\n")
g_ab8 <- g_cce8 <- numeric(Rmc)
for (r in 1:Rmc) {
  d <- generate_dgp(200, 6, 0.8, list(c(0.8, 1.2)), m = 1, sigma_eta = 0.5, seed = r)
  g_ab8[r]  <- tryCatch(arellano_bond_gmm(d)$gamma_hat, error = function(e) NA)
  g_cce8[r] <- tryCatch(gmm_cce_iv(d)$gamma_hat, error = function(e) NA)
}
g_ab8 <- na.omit(g_ab8); g_cce8 <- na.omit(g_cce8)
cat(sprintf("  AB:  Bias=%+.4f, RMSE=%.4f\n",
            mean(g_ab8) - 0.8, sqrt(mean((g_ab8 - 0.8)^2))))
cat(sprintf("  CCE: Bias=%+.4f, RMSE=%.4f\n",
            mean(g_cce8) - 0.8, sqrt(mean((g_cce8 - 0.8)^2))))

cat("\nGMM-3 COMPLETE.\n")
