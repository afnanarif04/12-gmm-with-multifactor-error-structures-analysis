###########################################################################
# Estimator 4 of 5: GMM with latent group structure
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
generate_dgp       <- e$generate_dgp
first_difference    <- e$first_difference
arellano_bond_gmm   <- e$arellano_bond_gmm

library(MASS)

# =================================================================
# 1. INDIVIDUAL GMM (unit-by-unit IV)
# =================================================================
individual_gmm <- function(data) {
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

  theta_i <- matrix(NA, N, p)
  valid   <- rep(FALSE, N)

  for (i in 1:N) {
    Xi <- matrix(0, n_eq, p)
    Xi[, 1] <- dlag[i, ]
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
        th <- solve(ZX + 1e-8 * diag(p), t(Zi) %*% dep[i, ])
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
# 2. K-MEANS CLASSIFICATION (on slopes only)
# =================================================================
kmeans_classify <- function(theta_i, valid, G, n_init = 50) {
  N <- nrow(theta_i); p <- ncol(theta_i); k <- p - 1
  beta_valid <- theta_i[valid, 2:p, drop = FALSE]
  nv <- sum(valid)

  if (nv < G) {
    labels <- rep(0L, N); labels[!valid] <- -1L
    return(list(labels = labels, centers = matrix(0, G, k), inertia = Inf))
  }

  if (G == 1) {
    ctr <- matrix(colMeans(beta_valid), 1, k)
    wss <- sum(sweep(beta_valid, 2, ctr)^2)
    labels <- rep(0L, N); labels[!valid] <- -1L
    return(list(labels = labels, centers = ctr, inertia = wss))
  }

  set.seed(42)
  best_wss <- Inf; best_lab <- NULL; best_ctr <- NULL

  for (init in 1:n_init) {
    idx <- sample(nv, G)
    ctr <- beta_valid[idx, , drop = FALSE]

    empty <- FALSE
    for (iter in 1:200) {
      dists <- matrix(0, nv, G)
      for (g in 1:G) dists[, g] <- rowSums(sweep(beta_valid, 2, ctr[g, ])^2)
      lab <- apply(dists, 1, which.min) - 1L

      new_ctr <- matrix(0, G, k)
      for (g in 1:G) {
        mg <- lab == (g - 1)
        if (sum(mg) == 0) { empty <- TRUE; break }
        new_ctr[g, ] <- colMeans(beta_valid[mg, , drop = FALSE])
      }
      if (empty) break
      if (max(abs(new_ctr - ctr)) < 1e-8) { ctr <- new_ctr; break }
      ctr <- new_ctr
    }
    if (empty) next

    wss <- 0
    for (g in 1:G) {
      mg <- lab == (g - 1)
      if (sum(mg) > 0) wss <- wss + sum(sweep(beta_valid[mg, , drop = FALSE], 2, ctr[g, ])^2)
    }
    if (wss < best_wss) { best_wss <- wss; best_lab <- lab; best_ctr <- ctr }
  }

  if (is.null(best_lab)) {
    labels <- rep(0L, N); labels[!valid] <- -1L
    return(list(labels = labels, centers = matrix(colMeans(beta_valid), 1, k), inertia = Inf))
  }

  labels <- rep(-1L, N)
  vi <- which(valid)
  for (j in seq_along(vi)) labels[vi[j]] <- best_lab[j]

  list(labels = labels, centers = best_ctr, inertia = best_wss)
}

# =================================================================
# 3. POST-CLASSIFICATION GROUP GMM
# =================================================================
group_gmm <- function(data, labels, G) {
  y <- data$y; N <- data$N; k <- data$k; T_obs <- data$T_obs
  fd <- first_difference(data); dy <- fd$dy; dx <- fd$dx
  n_eq <- T_obs - 1; p <- 1 + k

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
    if (Ng < max(5, p + 1)) next

    dep_g  <- dep[mg, , drop = FALSE]
    dlag_g <- dlag[mg, , drop = FALSE]
    dxs_g  <- dxs[mg, , , drop = FALSE]
    y_g    <- y[mg, , drop = FALSE]

    n_iv_pe <- 1 + k; n_iv <- n_eq * n_iv_pe
    reg_g <- array(0, c(Ng, n_eq, p))
    Z_g   <- array(0, c(Ng, n_eq, n_iv))

    for (t_idx in 1:n_eq) {
      t_act <- t_idx + 2
      cb <- (t_idx - 1) * n_iv_pe
      reg_g[, t_idx, 1] <- dlag_g[, t_idx]
      for (j in 1:k) reg_g[, t_idx, 1 + j] <- dxs_g[, t_idx, j]
      Z_g[, t_idx, cb + 1] <- y_g[, max(t_act - 2, 1)]
      for (j in 1:k) Z_g[, t_idx, cb + 1 + j] <- dxs_g[, t_idx, j]
    }

    XZ <- matrix(0, p, n_iv); Zy <- rep(0, n_iv)
    for (i in 1:Ng) {
      ri <- matrix(reg_g[i, , ], n_eq, p)
      zi <- matrix(Z_g[i, , ], n_eq, n_iv)
      XZ <- XZ + t(ri) %*% zi; Zy <- Zy + t(zi) %*% dep_g[i, ]
    }

    H <- 2 * diag(n_eq)
    if (n_eq > 1) for (ii in 1:(n_eq-1)) { H[ii, ii+1] <- -1; H[ii+1, ii] <- -1 }
    ZHZ <- matrix(0, n_iv, n_iv)
    for (i in 1:Ng) {
      zi <- matrix(Z_g[i, , ], n_eq, n_iv)
      ZHZ <- ZHZ + t(zi) %*% H %*% zi
    }
    ZHZ <- ZHZ / Ng
    W1 <- tryCatch(solve(ZHZ + 1e-6 * diag(n_iv)), error = function(e) ginv(ZHZ))
    th1 <- tryCatch(solve(XZ %*% W1 %*% t(XZ) + 1e-10 * diag(p), XZ %*% W1 %*% Zy),
                     error = function(e) ginv(XZ %*% W1 %*% t(XZ)) %*% (XZ %*% W1 %*% Zy))

    # Step 2
    ZeZ <- matrix(0, n_iv, n_iv)
    for (i in 1:Ng) {
      ri <- matrix(reg_g[i, , ], n_eq, p)
      zi <- matrix(Z_g[i, , ], n_eq, n_iv)
      ei <- dep_g[i, ] - ri %*% th1
      ze <- t(zi) %*% ei; ZeZ <- ZeZ + ze %*% t(ze)
    }
    ZeZ <- ZeZ / Ng
    W2 <- tryCatch(solve(ZeZ + 1e-6 * diag(n_iv)), error = function(e) ginv(ZeZ))
    th <- tryCatch(solve(XZ %*% W2 %*% t(XZ) + 1e-10 * diag(p), XZ %*% W2 %*% Zy),
                   error = function(e) ginv(XZ %*% W2 %*% t(XZ)) %*% (XZ %*% W2 %*% Zy))

    gg[g] <- as.numeric(th[1])
    gb[g, ] <- as.numeric(th[2:p])
  }

  list(group_gamma = gg, group_beta = gb, group_N = gn, labels = labels, G = G)
}

# =================================================================
# 4. FULL GMM-LGS WRAPPER
# =================================================================
gmm_lgs <- function(data, G = 2) {
  ind <- individual_gmm(data)
  theta_i <- ind$theta_i; valid <- ind$valid
  cl <- kmeans_classify(theta_i, valid, G)
  labels <- cl$labels
  labels[labels == -1] <- 0L
  res <- group_gmm(data, labels, G)
  res$theta_i <- theta_i; res$valid <- valid; res$centers <- cl$centers
  res
}

# =================================================================
# 5. CLASSIFICATION ACCURACY (best permutation)
# =================================================================
classification_accuracy <- function(labels_est, labels_true, G) {
  N <- length(labels_est)
  if (G == 1) return(list(accuracy = 1.0, adj_rand = 1.0))

  # Generate all permutations of 1:G
  perm_list <- function(n) {
    if (n == 1) return(list(1))
    out <- list()
    for (i in 1:n) {
      sub <- perm_list(n - 1)
      for (s in sub) {
        p <- integer(n)
        p[1] <- i
        idx <- 2
        for (v in s) {
          if (v >= i) p[idx] <- v + 1 else p[idx] <- v
          idx <- idx + 1
        }
        # Fix: simpler approach
      }
    }
    # Use a simple recursive permutation
    if (n <= 6) {
      perms <- combinat_permn(n)
    }
    perms
  }

  # Simple permutation generator (no external package)
  gen_perms <- function(n) {
    if (n == 1) return(list(c(1)))
    if (n == 2) return(list(c(1, 2), c(2, 1)))
    result <- list()
    for (i in 1:n) {
      rest <- setdiff(1:n, i)
      sub_perms <- gen_perms(n - 1)
      for (sp in sub_perms) {
        mapped <- rest[sp]
        result <- c(result, list(c(i, mapped)))
      }
    }
    result
  }

  best_acc <- 0
  for (perm in gen_perms(G)) {
    relab <- rep(-1L, N)
    for (ge in 1:G) relab[labels_est == (ge - 1)] <- perm[ge] - 1
    acc <- mean(relab == labels_true, na.rm = TRUE)
    best_acc <- max(best_acc, acc)
  }

  # ARI
  ct <- sort(unique(labels_true))
  cp <- sort(unique(labels_est[labels_est >= 0]))
  nij <- matrix(0, length(ct), length(cp))
  for (i in seq_along(ct))
    for (j in seq_along(cp))
      nij[i, j] <- sum(labels_true == ct[i] & labels_est == cp[j])
  ai <- rowSums(nij); bj <- colSums(nij)
  c2 <- function(n) n * (n - 1) / 2
  s_nij <- sum(sapply(nij, c2))
  s_ai <- sum(sapply(ai, c2)); s_bj <- sum(sapply(bj, c2))
  cn <- c2(N)
  if (cn == 0) { ari <- 0 } else {
    ex <- s_ai * s_bj / cn; mx <- 0.5 * (s_ai + s_bj)
    ari <- if (mx == ex) { if (s_nij == ex) 1.0 else 0.0 } else (s_nij - ex) / (mx - ex)
  }

  list(accuracy = best_acc, adj_rand = ari)
}

# =================================================================
# 6. STANDALONE TEST
# =================================================================
cat("=============================================================\n")
cat("GMM-4: GMM-LGS — STANDALONE TEST\n")
cat("=============================================================\n")

cat("\n--- Test 1: Individual estimates (m=0, G=1, N=500, T=6) ---\n")
d1 <- generate_dgp(500, 6, 0.5, list(c(0.8, 1.2)), m = 0, seed = 42)
ind <- individual_gmm(d1)
nv <- sum(ind$valid)
cat(sprintf("  Valid: %d / %d\n", nv, 500))
if (nv > 10) {
  g_ind <- ind$theta_i[ind$valid, 1]
  b_ind <- ind$theta_i[ind$valid, 2]
  cat(sprintf("  gamma_i: mean=%.4f, sd=%.4f\n", mean(g_ind), sd(g_ind)))
  cat(sprintf("  beta1_i: mean=%.4f, sd=%.4f\n", mean(b_ind), sd(b_ind)))
}

cat("\n--- Test 2: G=2 groups (well-separated, m=0) ---\n")
bg2 <- list(c(0.5, 1.5), c(1.5, 0.5))
d2 <- generate_dgp(500, 6, 0.5, bg2, m = 0, seed = 42)
r2 <- gmm_lgs(d2, G = 2)
for (g in 1:2) {
  cat(sprintf("  G%d (n=%d): gamma=%.4f, beta=[%.4f, %.4f]\n",
              g, r2$group_N[g], r2$group_gamma[g], r2$group_beta[g, 1], r2$group_beta[g, 2]))
}
acc2 <- classification_accuracy(r2$labels, d2$groups, 2)
cat(sprintf("  Accuracy=%.3f, ARI=%.3f\n", acc2$accuracy, acc2$adj_rand))

cat("\n--- Test 3: G=2 with factors (m=1) -> LGS should struggle ---\n")
d3 <- generate_dgp(500, 6, 0.5, bg2, m = 1, sigma_eta = 0.5, seed = 42)
r3 <- gmm_lgs(d3, G = 2)
for (g in 1:2) {
  if (!is.na(r3$group_gamma[g]))
    cat(sprintf("  G%d (n=%d): gamma=%.4f, beta=[%.4f, %.4f]\n",
                g, r3$group_N[g], r3$group_gamma[g], r3$group_beta[g, 1], r3$group_beta[g, 2]))
}
acc3 <- classification_accuracy(r3$labels, d3$groups, 2)
cat(sprintf("  Accuracy=%.3f (degraded by factors)\n", acc3$accuracy))

cat("\n--- Test 4: G=1 (should recover pooled AB) ---\n")
d4 <- generate_dgp(500, 6, 0.5, list(c(0.8, 1.2)), m = 0, seed = 42)
r4 <- gmm_lgs(d4, G = 1)
ab4 <- arellano_bond_gmm(d4)
cat(sprintf("  LGS(G=1): gamma=%.4f, beta=[%.4f, %.4f]\n",
            r4$group_gamma[1], r4$group_beta[1, 1], r4$group_beta[1, 2]))
cat(sprintf("  AB:       gamma=%.4f, beta=[%.4f, %.4f]\n",
            ab4$gamma_hat, ab4$beta_hat[1], ab4$beta_hat[2]))

cat("\n--- Test 5: Mini MC (R=50, N=300, T=6, G=2, m=0) ---\n")
Rmc <- 50
acc_vec <- numeric(Rmc); g_rmse <- numeric(Rmc)
for (r in 1:Rmc) {
  d <- generate_dgp(300, 6, 0.5, bg2, m = 0, seed = r)
  res <- tryCatch({
    rr <- gmm_lgs(d, G = 2)
    w <- rr$group_N / sum(rr$group_N)
    gm <- sum(w * rr$group_gamma, na.rm = TRUE)
    ac <- classification_accuracy(rr$labels, d$groups, 2)$accuracy
    list(gm = gm, ac = ac)
  }, error = function(e) list(gm = NA, ac = NA))
  g_rmse[r] <- res$gm; acc_vec[r] <- res$ac
}
g_rmse <- na.omit(g_rmse); acc_vec <- na.omit(acc_vec)
cat(sprintf("  gamma MG: Bias=%+.4f, RMSE=%.4f\n",
            mean(g_rmse) - 0.5, sqrt(mean((g_rmse - 0.5)^2))))
cat(sprintf("  Accuracy: mean=%.3f, sd=%.3f\n", mean(acc_vec), sd(acc_vec)))

cat("\nGMM-4 COMPLETE.\n")
