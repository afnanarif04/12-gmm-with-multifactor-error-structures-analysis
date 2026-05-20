###########################################################################
# Estimator 2 of 5: GMM with interactive fixed effects (Mundlak-Chamberlain projection)
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
# If that fails because GMM-1 runs tests at top level, just copy functions.
# Alternatively: manually source only the function definitions.
# For safety, re-attach the key functions:
generate_dgp       <- e$generate_dgp
first_difference    <- e$first_difference
arellano_bond_gmm   <- e$arellano_bond_gmm

library(MASS)

# =================================================================
# 1. GMM-IFE PROJECTION (Mundlak-Chamberlain augmentation)
# =================================================================
gmm_ife_proj <- function(data, m_assumed = 1, step = 2) {
  if (m_assumed == 0) return(arellano_bond_gmm(data, step))

  y <- data$y; x <- data$x
  N <- data$N; k <- data$k; T_obs <- data$T_obs
  fd <- first_difference(data)
  dy <- fd$dy; dx <- fd$dx
  n_eq <- T_obs - 1
  p_s <- 1 + k; p_aug <- p_s + k

  dep  <- if (n_eq == 1) matrix(dy[, 2], ncol=1) else dy[, 2:T_obs]
  dlag <- if (n_eq == 1) matrix(dy[, 1], ncol=1) else dy[, 1:(T_obs-1)]
  dxs  <- array(0, c(N, n_eq, k))
  for (j in 1:k) {
    sl <- if (n_eq == 1) matrix(dx[, 2, j], ncol=1) else dx[, 2:T_obs, j]
    dxs[, , j] <- sl
  }

  # Cross-section demean
  dep_dm  <- sweep(dep, 2, colMeans(dep))
  dlag_dm <- sweep(dlag, 2, colMeans(dlag))
  dxs_dm  <- array(0, c(N, n_eq, k))
  for (j in 1:k) {
    dxs_dm[, , j] <- sweep(matrix(dxs[, , j], N, n_eq), 2,
                            colMeans(matrix(dxs[, , j], N, n_eq)))
  }

  # Mundlak device: individual time-means of x
  xbar <- matrix(0, N, k)
  for (j in 1:k) xbar[, j] <- rowMeans(x[, 2:(T_obs+1), j])
  xbar_dm <- sweep(xbar, 2, colMeans(xbar))

  # Augmented regressors: [dlag_dm, dxs_dm, xbar_dm]
  reg <- array(0, c(N, n_eq, p_aug))
  reg[, , 1] <- dlag_dm
  for (j in 1:k) reg[, , 1 + j] <- dxs_dm[, , j]
  for (t_idx in 1:n_eq) reg[, t_idx, (p_s+1):p_aug] <- xbar_dm

  # Instruments: [y_{i,t-2} demeaned, dxs_dm, xbar_dm]
  n_iv_pe <- 1 + k + k; n_iv <- n_eq * n_iv_pe
  Z <- array(0, c(N, n_eq, n_iv))
  for (t_idx in 1:n_eq) {
    t_act <- t_idx + 2
    cb <- (t_idx - 1) * n_iv_pe
    Z[, t_idx, cb + 1] <- y[, max(t_act - 2, 1)] - mean(y[, max(t_act - 2, 1)])
    for (j in 1:k) Z[, t_idx, cb + 1 + j] <- dxs_dm[, t_idx, j]
    Z[, t_idx, (cb + 1 + k + 1):(cb + 1 + k + k)] <- xbar_dm
  }

  # GMM estimation
  XZ <- matrix(0, p_aug, n_iv); Zy <- rep(0, n_iv)
  for (i in 1:N) {
    ri <- matrix(reg[i, , ], n_eq, p_aug)
    zi <- matrix(Z[i, , ], n_eq, n_iv)
    XZ <- XZ + t(ri) %*% zi
    Zy <- Zy + t(zi) %*% dep_dm[i, ]
  }

  # Step 1 weight (H matrix)
  H <- 2 * diag(n_eq)
  if (n_eq > 1) for (ii in 1:(n_eq-1)) { H[ii, ii+1] <- -1; H[ii+1, ii] <- -1 }
  ZHZ <- matrix(0, n_iv, n_iv)
  for (i in 1:N) {
    zi <- matrix(Z[i, , ], n_eq, n_iv)
    ZHZ <- ZHZ + t(zi) %*% H %*% zi
  }
  ZHZ <- ZHZ / N
  W1 <- tryCatch(solve(ZHZ + 1e-8 * diag(n_iv)), error = function(e) ginv(ZHZ))
  th1 <- tryCatch(solve(XZ %*% W1 %*% t(XZ) + 1e-10 * diag(p_aug),
                         XZ %*% W1 %*% Zy),
                   error = function(e) ginv(XZ %*% W1 %*% t(XZ)) %*% (XZ %*% W1 %*% Zy))

  if (step >= 2) {
    ZeZ <- matrix(0, n_iv, n_iv)
    for (i in 1:N) {
      ri <- matrix(reg[i, , ], n_eq, p_aug)
      zi <- matrix(Z[i, , ], n_eq, n_iv)
      ei <- dep_dm[i, ] - ri %*% th1
      ze <- t(zi) %*% ei; ZeZ <- ZeZ + ze %*% t(ze)
    }
    ZeZ <- ZeZ / N
    W2 <- tryCatch(solve(ZeZ + 1e-8 * diag(n_iv)), error = function(e) ginv(ZeZ))
    th <- tryCatch(solve(XZ %*% W2 %*% t(XZ) + 1e-10 * diag(p_aug),
                         XZ %*% W2 %*% Zy),
                   error = function(e) ginv(XZ %*% W2 %*% t(XZ)) %*% (XZ %*% W2 %*% Zy))
    Wf <- W2
  } else { th <- th1; Wf <- W1 }

  # Sandwich variance (structural parameters only)
  ef <- matrix(0, N, n_eq)
  for (i in 1:N) {
    ri <- matrix(reg[i, , ], n_eq, p_aug)
    ef[i, ] <- dep_dm[i, ] - ri %*% th
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
       method = "GMM-IFE-Proj", N = N, T_obs = T_obs)
}

# Wrapper
gmm_ife <- function(data, m_assumed = 1, method = "proj", step = 2) {
  if (method == "proj") return(gmm_ife_proj(data, m_assumed, step))
  stop(paste("Unknown method:", method))
}

# =================================================================
# 2. STANDALONE TEST
# =================================================================
cat("=============================================================\n")
cat("GMM-2: GMM-IFE (Projection) — STANDALONE TEST\n")
cat("=============================================================\n")

cat("\n--- Test 1: m=0 -> should match AB ---\n")
d0 <- generate_dgp(500, 6, 0.5, list(c(0.8, 1.2)), m = 0, seed = 42)
ab <- arellano_bond_gmm(d0)
ife <- gmm_ife(d0, m_assumed = 0)
cat(sprintf("  AB:  gamma=%.4f\n", ab$gamma_hat))
cat(sprintf("  IFE: gamma=%.4f (should match AB)\n", ife$gamma_hat))

cat("\n--- Test 2: m=1 -> IFE should reduce bias vs AB ---\n")
d1 <- generate_dgp(500, 6, 0.5, list(c(0.8, 1.2)), m = 1, sigma_eta = 0.5, seed = 42)
ab1 <- arellano_bond_gmm(d1)
ife1 <- gmm_ife(d1, m_assumed = 1)
cat(sprintf("  AB:  gamma=%.4f (bias=%+.4f)\n", ab1$gamma_hat, ab1$gamma_hat - 0.5))
cat(sprintf("  IFE: gamma=%.4f (bias=%+.4f)\n", ife1$gamma_hat, ife1$gamma_hat - 0.5))
cat(sprintf("  IFE: beta =[%.4f, %.4f]\n", ife1$beta_hat[1], ife1$beta_hat[2]))

cat("\n--- Test 3: Mini MC (R=100, N=200, T=6, m=1) ---\n")
Rmc <- 100
g_ab <- g_ife <- numeric(Rmc)
for (r in 1:Rmc) {
  d <- generate_dgp(200, 6, 0.5, list(c(0.8, 1.2)), m = 1, sigma_eta = 0.5, seed = r)
  g_ab[r]  <- tryCatch(arellano_bond_gmm(d)$gamma_hat, error = function(e) NA)
  g_ife[r] <- tryCatch(gmm_ife(d, m_assumed = 1)$gamma_hat, error = function(e) NA)
}
g_ab  <- na.omit(g_ab); g_ife <- na.omit(g_ife)
cat(sprintf("  AB:  Bias=%+.4f, RMSE=%.4f (n=%d)\n",
            mean(g_ab) - 0.5, sqrt(mean((g_ab - 0.5)^2)), length(g_ab)))
cat(sprintf("  IFE: Bias=%+.4f, RMSE=%.4f (n=%d)\n",
            mean(g_ife) - 0.5, sqrt(mean((g_ife - 0.5)^2)), length(g_ife)))
cat(sprintf("  RMSE reduction: %.0f%%\n", 100 * (1 - sqrt(mean((g_ife - 0.5)^2)) / sqrt(mean((g_ab - 0.5)^2)))))

cat("\nGMM-2 COMPLETE.\n")
