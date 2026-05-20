###########################################################################
# Estimator 1 of 5: Arellano-Bond first-difference GMM (benchmark)
#
# Also defines the data generating process and first-difference helper used by all other scripts.
#
# Run from the command line with:  Rscript <this_file>.R
# Or open in an editor, set the working directory to the code/ folder,
# and run interactively.
#
# Dependencies: base R and the MASS package only (install.packages("MASS")).
###########################################################################

library(MASS)

# =================================================================
# 1. DGP GENERATOR
# =================================================================
generate_dgp <- function(N, T_obs, gamma, beta_groups, m = 0,
                         sigma_u = 1.0, sigma_eta = 0.5,
                         corr_eta_x = FALSE, heterosked = FALSE,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  G <- length(beta_groups)
  k <- length(beta_groups[[1]])
  groups <- rep(0:(G-1), length.out = N)
  groups <- sample(groups)
  alpha <- rnorm(N, 0, 1)
  delta <- rnorm(T_obs + 1, 0, 0.3)
  f_mat <- matrix(0, T_obs + 1, max(m, 1))
  if (m > 0) {
    for (j in 1:m) {
      f_mat[1, j] <- rnorm(1)
      for (t in 2:(T_obs + 1)) f_mat[t, j] <- 0.5 * f_mat[t-1, j] + rnorm(1)
    }
  }
  eta <- matrix(0, N, max(m, 1))
  if (m > 0) eta <- matrix(rnorm(N * m, 1.0, sigma_eta), N, m)
  rho_x <- 0.3
  x_array <- array(0, c(N, T_obs + 1, k))
  for (i in 1:N) {
    for (j in 1:k) {
      x_array[i, 1, j] <- rnorm(1)
      for (t in 2:(T_obs + 1)) {
        x_array[i, t, j] <- rho_x * x_array[i, t-1, j] + rnorm(1)
        if (corr_eta_x && m > 0) x_array[i, t, j] <- x_array[i, t, j] + 0.3 * eta[i, 1]
      }
    }
  }
  u <- matrix(rnorm(N * (T_obs + 1), 0, sigma_u), N, T_obs + 1)
  if (heterosked) {
    sig_i <- 0.5 + runif(N)
    for (i in 1:N) u[i, ] <- rnorm(T_obs + 1, 0, sig_i[i])
  }
  y <- matrix(0, N, T_obs + 1)
  T_burn <- 50
  for (i in 1:N) {
    g <- groups[i] + 1
    beta_i <- beta_groups[[g]]
    yp <- rnorm(1)
    for (tb in 1:T_burn) yp <- gamma * yp + alpha[i] + rnorm(1, 0, sigma_u)
    y[i, 1] <- yp
    for (t in 2:(T_obs + 1)) {
      y[i, t] <- alpha[i] + delta[t] + gamma * y[i, t-1] + sum(x_array[i, t, ] * beta_i)
      if (m > 0) y[i, t] <- y[i, t] + sum(eta[i, 1:m] * f_mat[t, 1:m])
      y[i, t] <- y[i, t] + u[i, t]
    }
  }
  list(y = y, x = x_array, alpha = alpha, delta = delta,
       eta = eta, f = f_mat, groups = groups,
       beta_true = beta_groups, gamma = gamma, u = u,
       N = N, T_obs = T_obs, k = k, m = m, G = G)
}

# =================================================================
# 2. FIRST DIFFERENCING
# =================================================================
first_difference <- function(data) {
  y <- data$y; x <- data$x
  N <- data$N; T_obs <- data$T_obs; k <- data$k
  dy <- y[, 2:(T_obs+1)] - y[, 1:T_obs]
  dx <- array(0, c(N, T_obs, k))
  for (j in 1:k) dx[, , j] <- x[, 2:(T_obs+1), j] - x[, 1:T_obs, j]
  list(dy = dy, dx = dx)
}

# =================================================================
# 3. ARELLANO-BOND GMM (2-step)
# =================================================================
arellano_bond_gmm <- function(data, step = 2) {
  y <- data$y; N <- data$N; T_obs <- data$T_obs; k <- data$k
  fd <- first_difference(data)
  dy <- fd$dy; dx <- fd$dx
  n_eq <- T_obs - 1; p <- 1 + k

  dep  <- if (n_eq == 1) matrix(dy[, 2], ncol=1) else dy[, 2:T_obs]
  dlag <- if (n_eq == 1) matrix(dy[, 1], ncol=1) else dy[, 1:(T_obs-1)]
  dxs  <- array(0, c(N, n_eq, k))
  for (j in 1:k) {
    sl <- if (n_eq == 1) matrix(dx[, 2, j], ncol=1) else dx[, 2:T_obs, j]
    dxs[, , j] <- sl
  }

  n_iv_pe <- 1 + k; n_iv <- n_eq * n_iv_pe
  reg <- array(0, c(N, n_eq, p))
  Z   <- array(0, c(N, n_eq, n_iv))
  for (t_idx in 1:n_eq) {
    t_act <- t_idx + 2
    reg[, t_idx, 1] <- dlag[, t_idx]
    for (j in 1:k) reg[, t_idx, 1+j] <- dxs[, t_idx, j]
    cb <- (t_idx - 1) * n_iv_pe
    Z[, t_idx, cb + 1] <- y[, max(t_act - 2, 1)]
    for (j in 1:k) Z[, t_idx, cb + 1 + j] <- dxs[, t_idx, j]
  }

  XZ <- matrix(0, p, n_iv); Zy <- rep(0, n_iv)
  for (i in 1:N) {
    ri <- matrix(reg[i, , ], n_eq, p)
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
  W1 <- tryCatch(solve(ZHZ + 1e-8 * diag(n_iv)), error = function(e) ginv(ZHZ))
  th1 <- tryCatch(solve(XZ %*% W1 %*% t(XZ) + 1e-10 * diag(p), XZ %*% W1 %*% Zy),
                   error = function(e) ginv(XZ %*% W1 %*% t(XZ)) %*% (XZ %*% W1 %*% Zy))

  if (step >= 2) {
    ZeZ <- matrix(0, n_iv, n_iv)
    for (i in 1:N) {
      ri <- matrix(reg[i, , ], n_eq, p)
      zi <- matrix(Z[i, , ], n_eq, n_iv)
      ei <- dep[i, ] - ri %*% th1
      ze <- t(zi) %*% ei; ZeZ <- ZeZ + ze %*% t(ze)
    }
    ZeZ <- ZeZ / N
    W2 <- tryCatch(solve(ZeZ + 1e-8 * diag(n_iv)), error = function(e) ginv(ZeZ))
    th <- tryCatch(solve(XZ %*% W2 %*% t(XZ) + 1e-10 * diag(p), XZ %*% W2 %*% Zy),
                   error = function(e) ginv(XZ %*% W2 %*% t(XZ)) %*% (XZ %*% W2 %*% Zy))
    Wf <- W2
  } else { th <- th1; Wf <- W1 }

  # Sandwich variance
  ef <- matrix(0, N, n_eq)
  for (i in 1:N) {
    ri <- matrix(reg[i, , ], n_eq, p)
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

  # J-test
  g_bar <- rep(0, n_iv)
  for (i in 1:N) {
    zi <- matrix(Z[i, , ], n_eq, n_iv)
    g_bar <- g_bar + t(zi) %*% ef[i, ]
  }
  g_bar <- g_bar / N
  J <- as.numeric(N * t(g_bar) %*% Wf %*% g_bar)

  list(theta = as.numeric(th), gamma_hat = as.numeric(th[1]),
       beta_hat = as.numeric(th[2:p]), se = as.numeric(se),
       J_stat = J, df_J = n_iv - p, N = N, T_obs = T_obs)
}

# =================================================================
# 4. STANDALONE TEST
# =================================================================
cat("=============================================================\n")
cat("GMM-1: Arellano-Bond GMM — STANDALONE TEST\n")
cat("=============================================================\n")

cat("\n--- Test 1: m=0, G=1 (no factors, no groups) ---\n")
cat("True: gamma=0.5, beta=[0.8, 1.2]\n")
d0 <- generate_dgp(500, 6, 0.5, list(c(0.8, 1.2)), m = 0, seed = 42)
r0 <- arellano_bond_gmm(d0)
cat(sprintf("  gamma = %.4f (se=%.4f)\n", r0$gamma_hat, r0$se[1]))
cat(sprintf("  beta  = [%.4f, %.4f]\n", r0$beta_hat[1], r0$beta_hat[2]))
cat(sprintf("  J-stat = %.3f (df=%d)\n", r0$J_stat, r0$df_J))

cat("\n--- Test 2: m=1 (AB should be biased by factors) ---\n")
d1 <- generate_dgp(500, 6, 0.5, list(c(0.8, 1.2)), m = 1, sigma_eta = 0.5, seed = 42)
r1 <- arellano_bond_gmm(d1)
cat(sprintf("  gamma = %.4f (bias=%+.4f)\n", r1$gamma_hat, r1$gamma_hat - 0.5))

cat("\n--- Test 3: Mini MC (R=100, N=200, T=6, m=0) ---\n")
R_mc <- 100; g_est <- numeric(R_mc)
for (r in 1:R_mc) {
  d <- generate_dgp(200, 6, 0.5, list(c(0.8, 1.2)), m = 0, seed = r)
  g_est[r] <- arellano_bond_gmm(d)$gamma_hat
}
cat(sprintf("  gamma: Mean=%.4f, Bias=%+.4f, RMSE=%.4f\n",
            mean(g_est), mean(g_est) - 0.5, sqrt(mean((g_est - 0.5)^2))))

cat("\n--- Test 4: Mini MC (R=100, N=200, T=6, m=1) ---\n")
g_est2 <- numeric(R_mc)
for (r in 1:R_mc) {
  d <- generate_dgp(200, 6, 0.5, list(c(0.8, 1.2)), m = 1, sigma_eta = 0.5, seed = r)
  g_est2[r] <- tryCatch(arellano_bond_gmm(d)$gamma_hat, error = function(e) NA)
}
g_est2 <- na.omit(g_est2)
cat(sprintf("  gamma: Mean=%.4f, Bias=%+.4f, RMSE=%.4f (n=%d)\n",
            mean(g_est2), mean(g_est2) - 0.5, sqrt(mean((g_est2 - 0.5)^2)), length(g_est2)))

cat("\nGMM-1 COMPLETE.\n")
