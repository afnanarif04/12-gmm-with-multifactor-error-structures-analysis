###########################################################################
# Estimator 5 of 5: GMM with interactive fixed effects and latent group structure
#
# Requires 02_estimator_ab_gmm.R, 03_estimator_gmm_ife.R, 05_estimator_gmm_lgs.R.
#
# Run from the command line with:  Rscript <this_file>.R
# Or open in an editor, set the working directory to the code/ folder,
# and run interactively.
#
# Dependencies: base R and the MASS package only (install.packages("MASS")).
###########################################################################

e1 <- new.env(); sys.source("02_estimator_ab_gmm.R", envir = e1)
generate_dgp       <- e1$generate_dgp
first_difference    <- e1$first_difference
arellano_bond_gmm   <- e1$arellano_bond_gmm

e2 <- new.env(); sys.source("03_estimator_gmm_ife.R", envir = e2)
gmm_ife <- e2$gmm_ife

e4 <- new.env(); sys.source("05_estimator_gmm_lgs.R", envir = e4)
individual_gmm          <- e4$individual_gmm
kmeans_classify         <- e4$kmeans_classify
group_gmm               <- e4$group_gmm
classification_accuracy <- e4$classification_accuracy
gmm_lgs                 <- e4$gmm_lgs

library(MASS)

# =================================================================
# 1. FACTOR-PURGED INDIVIDUAL ESTIMATES (PC method)
# =================================================================
individual_gmm_purged_pc <- function(data, m_assumed = 1) {
  if (m_assumed == 0) return(individual_gmm(data))

  y <- data$y; x <- data$x
  N <- data$N; k <- data$k; T_obs <- data$T_obs
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

  # Initial pooled AB estimate
  ab <- tryCatch(arellano_bond_gmm(data, step = 1),
                  error = function(e) list(gamma_hat = 0, beta_hat = rep(0, k)))

  # Compute residuals
  resid <- dep - ab$gamma_hat * dlag
  for (j in 1:k) resid <- resid - ab$beta_hat[j] * dxs[, , j]

  # Cross-section demean
  resid_dm <- sweep(resid, 2, colMeans(resid))

  # Extract m factors via SVD
  m_use <- min(m_assumed, min(dim(resid_dm)) - 1)
  sv <- svd(resid_dm, nu = m_use, nv = m_use)
  eta_h <- sv$u[, 1:m_use, drop = FALSE] * rep(sv$d[1:m_use], each = N)
  df_h  <- t(sv$v[, 1:m_use, drop = FALSE])   # m x n_eq
  fc    <- eta_h %*% df_h                      # N x n_eq (factor component in dep)

  # Factor component in lagged dep
  dlag_dm <- sweep(dlag, 2, colMeans(dlag))
  ete <- t(eta_h) %*% eta_h + 1e-8 * diag(m_use)
  df_lag_h <- solve(ete, t(eta_h) %*% dlag_dm)
  fc_lag <- eta_h %*% df_lag_h

  # Purged data
  dep_p  <- dep - fc
  dlag_p <- dlag - fc_lag

  # Re-estimate individually on purged data
  theta_i <- matrix(NA, N, p)
  valid   <- rep(FALSE, N)

  for (i in 1:N) {
    Xi <- matrix(0, n_eq, p)
    Xi[, 1] <- dlag_p[i, ]
    for (j in 1:k) Xi[, 1 + j] <- dxs[i, , j]

    Zi <- matrix(0, n_eq, p)
    for (t_idx in 1:n_eq) {
      t_act <- t_idx + 2
      Zi[t_idx, 1] <- y[i, max(t_act - 2, 1)]
      for (j in 1:k) Zi[t_idx, 1 + j] <- dxs[i, t_idx, j]
    }

    if (n_eq < p) next

    res <- tryCatch({
      ZX <- t(Zi) %*% Xi
      if (qr(ZX)$rank < p) NULL
      else {
        th <- solve(ZX + 1e-8 * diag(p), t(Zi) %*% dep_p[i, ])
        if (abs(th[1]) < 3.0 && all(abs(th[2:p]) < 20)) {
          theta_i[i, ] <- as.numeric(th)
          valid[i] <- TRUE
        }
      }
    }, error = function(e) NULL)
  }

  list(theta_i = theta_i, valid = valid)
}

# =================================================================
# 2. POST-CLASSIFICATION GMM-IFE PER GROUP
# =================================================================
group_gmm_ife <- function(data, labels, G, m_assumed = 1) {
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

  gg <- rep(NA_real_, G); gb <- matrix(NA_real_, G, k); gn <- integer(G)

  for (g in 1:G) {
    mg <- which(labels == (g - 1)); Ng <- length(mg); gn[g] <- Ng
    if (Ng < max(5, p_s + 2)) next

    dep_g  <- dep[mg, , drop = FALSE]
    dlag_g <- dlag[mg, , drop = FALSE]
    dxs_g  <- dxs[mg, , , drop = FALSE]
    y_g    <- y[mg, , drop = FALSE]

    # Mundlak augmentation within group
    p_aug <- p_s + k
    xbar_g <- matrix(0, Ng, k)
    for (j in 1:k) xbar_g[, j] <- rowMeans(x[mg, 2:(T_obs+1), j])
    xbar_g_dm <- sweep(xbar_g, 2, colMeans(xbar_g))

    dep_dm  <- sweep(dep_g, 2, colMeans(dep_g))
    dlag_dm <- sweep(dlag_g, 2, colMeans(dlag_g))
    dxs_dm  <- array(0, c(Ng, n_eq, k))
    for (j in 1:k) {
      dxs_dm[, , j] <- sweep(matrix(dxs_g[, , j], Ng, n_eq), 2,
                               colMeans(matrix(dxs_g[, , j], Ng, n_eq)))
    }

    reg <- array(0, c(Ng, n_eq, p_aug))
    reg[, , 1] <- dlag_dm
    for (j in 1:k) reg[, , 1 + j] <- dxs_dm[, , j]
    for (ti in 1:n_eq) reg[, ti, (p_s+1):p_aug] <- xbar_g_dm

    n_iv_pe <- 1 + k + k; n_iv <- n_eq * n_iv_pe
    Z <- array(0, c(Ng, n_eq, n_iv))
    for (ti in 1:n_eq) {
      ta <- ti + 2; cb <- (ti - 1) * n_iv_pe
      Z[, ti, cb + 1] <- y_g[, max(ta - 2, 1)] - mean(y_g[, max(ta - 2, 1)])
      for (j in 1:k) Z[, ti, cb + 1 + j] <- dxs_dm[, ti, j]
      Z[, ti, (cb + 1 + k + 1):(cb + 1 + k + k)] <- xbar_g_dm
    }
    dep_use <- dep_dm

    # 2-step GMM
    XZ <- matrix(0, p_aug, n_iv); Zy <- rep(0, n_iv)
    for (i in 1:Ng) {
      ri <- matrix(reg[i, , ], n_eq, p_aug)
      zi <- matrix(Z[i, , ], n_eq, n_iv)
      XZ <- XZ + t(ri) %*% zi; Zy <- Zy + t(zi) %*% dep_use[i, ]
    }

    H <- 2 * diag(n_eq)
    if (n_eq > 1) for (ii in 1:(n_eq-1)) { H[ii, ii+1] <- -1; H[ii+1, ii] <- -1 }
    ZHZ <- matrix(0, n_iv, n_iv)
    for (i in 1:Ng) {
      zi <- matrix(Z[i, , ], n_eq, n_iv)
      ZHZ <- ZHZ + t(zi) %*% H %*% zi
    }
    ZHZ <- ZHZ / Ng
    W1 <- tryCatch(solve(ZHZ + 1e-6 * diag(n_iv)), error = function(e) ginv(ZHZ))
    th1 <- tryCatch(solve(XZ %*% W1 %*% t(XZ) + 1e-10 * diag(p_aug), XZ %*% W1 %*% Zy),
                     error = function(e) ginv(XZ %*% W1 %*% t(XZ)) %*% (XZ %*% W1 %*% Zy))

    ZeZ <- matrix(0, n_iv, n_iv)
    for (i in 1:Ng) {
      ri <- matrix(reg[i, , ], n_eq, p_aug)
      zi <- matrix(Z[i, , ], n_eq, n_iv)
      ei <- dep_use[i, ] - ri %*% th1; ze <- t(zi) %*% ei
      ZeZ <- ZeZ + ze %*% t(ze)
    }
    ZeZ <- ZeZ / Ng
    W2 <- tryCatch(solve(ZeZ + 1e-6 * diag(n_iv)), error = function(e) ginv(ZeZ))
    th <- tryCatch(solve(XZ %*% W2 %*% t(XZ) + 1e-10 * diag(p_aug), XZ %*% W2 %*% Zy),
                   error = function(e) ginv(XZ %*% W2 %*% t(XZ)) %*% (XZ %*% W2 %*% Zy))

    gg[g] <- as.numeric(th[1])
    gb[g, ] <- as.numeric(th[2:p_s])
  }

  list(group_gamma = gg, group_beta = gb, group_N = gn, labels = labels, G = G)
}

# =================================================================
# 3. FULL GMM-IFE-LGS WRAPPER
# =================================================================
gmm_ife_lgs <- function(data, m_assumed = 1, G = 2) {
  N <- data$N; k <- data$k

  # Stage 1: factor-purged individual estimates
  ind <- individual_gmm_purged_pc(data, m_assumed)
  theta_i <- ind$theta_i; valid <- ind$valid
  nv <- sum(valid)

  cat(sprintf("    Purged valid: %d / %d\n", nv, N))

  if (nv < 10) {
    cat("    Falling back to pooled IFE\n")
    res <- gmm_ife(data, m_assumed, method = "proj")
    return(list(group_gamma = res$gamma_hat,
                group_beta = matrix(res$beta_hat, 1, k),
                group_N = N, labels = rep(0L, N), G = 1))
  }

  # Stage 2a: classify purged estimates
  cl <- kmeans_classify(theta_i, valid, G)
  labels <- cl$labels
  labels[labels == -1] <- 0L

  # Stage 2b: post-classification GMM-IFE
  res <- group_gmm_ife(data, labels, G, m_assumed)
  res$n_valid <- nv; res$centers <- cl$centers
  res
}

# =================================================================
# 4. STANDALONE TEST
# =================================================================
cat("=============================================================\n")
cat("GMM-5: GMM-IFE-LGS (Combined) — STANDALONE TEST\n")
cat("=============================================================\n")

bg2 <- list(c(0.5, 1.5), c(1.5, 0.5))

cat("\n--- Test 1: m=1, G=2, N=500, T=6 (full model) ---\n")
d1 <- generate_dgp(500, 6, 0.5, bg2, m = 1, sigma_eta = 0.5, seed = 42)
ab <- arellano_bond_gmm(d1)
cat(sprintf("  AB (pooled): gamma=%.4f, beta=[%.4f, %.4f]\n",
            ab$gamma_hat, ab$beta_hat[1], ab$beta_hat[2]))
cat("  Running IFE-LGS...\n")
res1 <- gmm_ife_lgs(d1, m_assumed = 1, G = 2)
for (g in 1:res1$G) {
  if (!is.na(res1$group_gamma[g]))
    cat(sprintf("  G%d (n=%d): gamma=%.4f, beta=[%.4f, %.4f]\n",
                g, res1$group_N[g], res1$group_gamma[g],
                res1$group_beta[g, 1], res1$group_beta[g, 2]))
}
acc1 <- classification_accuracy(res1$labels, d1$groups, 2)
cat(sprintf("  Accuracy=%.3f, ARI=%.3f\n", acc1$accuracy, acc1$adj_rand))

cat("\n--- Test 2: Compare LGS vs IFE-LGS (m=1, G=2) ---\n")
cat("  LGS (no factor purging)...\n")
r_lgs <- gmm_lgs(d1, G = 2)
for (g in 1:2) {
  if (!is.na(r_lgs$group_gamma[g]))
    cat(sprintf("    LGS G%d (n=%d): gamma=%.4f\n",
                g, r_lgs$group_N[g], r_lgs$group_gamma[g]))
}
acc_lgs <- classification_accuracy(r_lgs$labels, d1$groups, 2)
cat(sprintf("    LGS Accuracy=%.3f\n", acc_lgs$accuracy))

cat("  IFE-LGS (with factor purging)...\n")
cat(sprintf("    IFE-LGS Accuracy=%.3f\n", acc1$accuracy))
cat(sprintf("    -> Factor purging %s classification\n",
            ifelse(acc1$accuracy > acc_lgs$accuracy, "IMPROVES", "does not improve")))

cat("\n--- Test 3: m=0, G=2 (no factors, groups only) ---\n")
d3 <- generate_dgp(500, 6, 0.5, bg2, m = 0, seed = 42)
cat("  Running IFE-LGS with m=0...\n")
res3 <- gmm_ife_lgs(d3, m_assumed = 0, G = 2)
for (g in 1:res3$G) {
  if (!is.na(res3$group_gamma[g]))
    cat(sprintf("  G%d (n=%d): gamma=%.4f, beta=[%.4f, %.4f]\n",
                g, res3$group_N[g], res3$group_gamma[g],
                res3$group_beta[g, 1], res3$group_beta[g, 2]))
}

cat("\n--- Test 4: Strong persistence (gamma=0.8, m=1, G=2) ---\n")
d4 <- generate_dgp(300, 6, 0.8, bg2, m = 1, sigma_eta = 0.5, seed = 42)
ab4 <- arellano_bond_gmm(d4)
cat(sprintf("  AB: gamma=%.4f (bias=%+.4f, should be ~-0.2)\n",
            ab4$gamma_hat, ab4$gamma_hat - 0.8))
cat("  Running IFE-LGS...\n")
res4 <- gmm_ife_lgs(d4, m_assumed = 1, G = 2)
w4 <- res4$group_N / sum(res4$group_N)
g4_mg <- sum(w4 * res4$group_gamma, na.rm = TRUE)
cat(sprintf("  IFE-LGS MG: gamma=%.4f (bias=%+.4f)\n", g4_mg, g4_mg - 0.8))

cat("\n--- Test 5: Mini MC (R=30, N=300, T=6, m=1, G=2) ---\n")
Rmc <- 30
g_ab <- g_lgs <- g_ife_lgs <- acc_lgs_v <- acc_ife_v <- numeric(Rmc)
for (r in 1:Rmc) {
  d <- generate_dgp(300, 6, 0.5, bg2, m = 1, sigma_eta = 0.5, seed = r)

  g_ab[r] <- tryCatch(arellano_bond_gmm(d)$gamma_hat, error = function(e) NA)

  tryCatch({
    rl <- gmm_lgs(d, G = 2)
    wl <- rl$group_N / sum(rl$group_N)
    g_lgs[r] <- sum(wl * rl$group_gamma, na.rm = TRUE)
    acc_lgs_v[r] <- classification_accuracy(rl$labels, d$groups, 2)$accuracy
  }, error = function(e) { g_lgs[r] <<- NA; acc_lgs_v[r] <<- NA })

  tryCatch({
    ri <- gmm_ife_lgs(d, m_assumed = 1, G = 2)
    wi <- ri$group_N / sum(ri$group_N)
    g_ife_lgs[r] <- sum(wi * ri$group_gamma, na.rm = TRUE)
    acc_ife_v[r] <- classification_accuracy(ri$labels, d$groups, 2)$accuracy
  }, error = function(e) { g_ife_lgs[r] <<- NA; acc_ife_v[r] <<- NA })

  if (r %% 10 == 0) cat(sprintf("    r=%d/%d\n", r, Rmc))
}
g_ab <- na.omit(g_ab); g_lgs <- na.omit(g_lgs); g_ife_lgs <- na.omit(g_ife_lgs)

cat(sprintf("\n  AB:       Bias=%+.4f, RMSE=%.4f\n",
            mean(g_ab) - 0.5, sqrt(mean((g_ab - 0.5)^2))))
cat(sprintf("  LGS:      Bias=%+.4f, RMSE=%.4f, Acc=%.3f\n",
            mean(g_lgs) - 0.5, sqrt(mean((g_lgs - 0.5)^2)),
            mean(acc_lgs_v, na.rm = TRUE)))
cat(sprintf("  IFE-LGS:  Bias=%+.4f, RMSE=%.4f, Acc=%.3f\n",
            mean(g_ife_lgs) - 0.5, sqrt(mean((g_ife_lgs - 0.5)^2)),
            mean(acc_ife_v, na.rm = TRUE)))

rmse_ab <- sqrt(mean((g_ab - 0.5)^2))
rmse_il <- sqrt(mean((g_ife_lgs - 0.5)^2))
cat(sprintf("\n  RMSE reduction (IFE-LGS vs AB): %.0f%%\n", 100 * (1 - rmse_il / rmse_ab)))

cat("\nGMM-5 COMPLETE.\n")
