###########################################################################
# Monte Carlo engine: runs all five estimators across the simulation designs
#
# Requires 02 through 06 in the same directory. Produces the simulation results table.
#
# Run from the command line with:  Rscript <this_file>.R
# Or open in an editor, set the working directory to the code/ folder,
# and run interactively.
#
# Dependencies: base R and the MASS package only (install.packages("MASS")).
###########################################################################

cat("Loading estimators...\n")
e1 <- new.env(); sys.source("02_estimator_ab_gmm.R", envir = e1)
generate_dgp       <- e1$generate_dgp
first_difference    <- e1$first_difference
arellano_bond_gmm   <- e1$arellano_bond_gmm

e2 <- new.env(); sys.source("03_estimator_gmm_ife.R", envir = e2)
gmm_ife <- e2$gmm_ife

e3 <- new.env(); sys.source("04_estimator_gmm_cce.R", envir = e3)
gmm_cce_iv <- e3$gmm_cce_iv

e4 <- new.env(); sys.source("05_estimator_gmm_lgs.R", envir = e4)
individual_gmm          <- e4$individual_gmm
kmeans_classify         <- e4$kmeans_classify
group_gmm               <- e4$group_gmm
gmm_lgs                 <- e4$gmm_lgs
classification_accuracy <- e4$classification_accuracy

e5 <- new.env(); sys.source("06_estimator_gmm_ife_lgs.R", envir = e5)
individual_gmm_purged_pc <- e5$individual_gmm_purged_pc
group_gmm_ife            <- e5$group_gmm_ife
gmm_ife_lgs              <- e5$gmm_ife_lgs

library(MASS)
cat("All estimators loaded.\n")

# =================================================================
# 1. DGP CONFIGURATIONS
# =================================================================
make_dgp_configs <- function() {
  list(
    DGP1 = list(name = "DGP1: Baseline (m=0, G=1)",
                m = 0, G = 1, gamma = 0.5,
                betas = list(c(0.8, 1.2)),
                sigma_eta = 0),
    DGP2 = list(name = "DGP2: Factors only (m=1, G=1)",
                m = 1, G = 1, gamma = 0.5,
                betas = list(c(0.8, 1.2)),
                sigma_eta = 0.5),
    DGP3 = list(name = "DGP3: Groups only (m=0, G=2)",
                m = 0, G = 2, gamma = 0.5,
                betas = list(c(0.5, 1.5), c(1.5, 0.5)),
                sigma_eta = 0),
    DGP4 = list(name = "DGP4: Full model (m=1, G=2, gamma=0.5)",
                m = 1, G = 2, gamma = 0.5,
                betas = list(c(0.5, 1.5), c(1.5, 0.5)),
                sigma_eta = 0.5),
    DGP5 = list(name = "DGP5: Full model (m=1, G=2, gamma=0.8)",
                m = 1, G = 2, gamma = 0.8,
                betas = list(c(0.5, 1.5), c(1.5, 0.5)),
                sigma_eta = 0.5)
  )
}

# =================================================================
# 2. SINGLE-REPLICATION RUNNER
# =================================================================
run_one_rep <- function(dgp_cfg, N, T_obs, seed) {
  data <- generate_dgp(N, T_obs, dgp_cfg$gamma, dgp_cfg$betas,
                        m = dgp_cfg$m, sigma_eta = dgp_cfg$sigma_eta,
                        seed = seed)
  true_gamma <- dgp_cfg$gamma
  G_true <- dgp_cfg$G
  m_true <- dgp_cfg$m

  results <- list()

  # --- GMM-1: AB ---
  tryCatch({
    r <- arellano_bond_gmm(data)
    results$AB <- list(gamma = r$gamma_hat, beta = r$beta_hat, ok = TRUE)
  }, error = function(e) { results$AB <<- list(gamma = NA, beta = NA, ok = FALSE) })

  # --- GMM-2: IFE-Proj (use m=1 always, tests robustness if m=0) ---
  tryCatch({
    r <- gmm_ife(data, m_assumed = max(m_true, 1), method = "proj")
    results$IFE <- list(gamma = r$gamma_hat, beta = r$beta_hat, ok = TRUE)
  }, error = function(e) { results$IFE <<- list(gamma = NA, beta = NA, ok = FALSE) })

  # --- GMM-3: CCE-IV ---
  tryCatch({
    r <- gmm_cce_iv(data)
    results$CCE <- list(gamma = r$gamma_hat, beta = r$beta_hat, ok = TRUE)
  }, error = function(e) { results$CCE <<- list(gamma = NA, beta = NA, ok = FALSE) })

  # --- GMM-4: LGS (G = G_true or 2) ---
  G_use <- max(G_true, 2)
  tryCatch({
    r <- gmm_lgs(data, G = G_use)
    w <- r$group_N / sum(r$group_N)
    gm <- sum(w * r$group_gamma, na.rm = TRUE)
    acc <- if (G_true > 1) classification_accuracy(r$labels, data$groups, G_use)$accuracy else NA
    results$LGS <- list(gamma = gm, beta = NA, acc = acc, ok = TRUE)
  }, error = function(e) { results$LGS <<- list(gamma = NA, beta = NA, acc = NA, ok = FALSE) })

  # --- GMM-5: IFE-LGS-Proj ---
  tryCatch({
    r <- gmm_ife_lgs(data, m_assumed = max(m_true, 1), G = G_use)
    w <- r$group_N / sum(r$group_N)
    gm <- sum(w * r$group_gamma, na.rm = TRUE)
    acc <- if (G_true > 1) classification_accuracy(r$labels, data$groups, G_use)$accuracy else NA
    results$IFE_LGS <- list(gamma = gm, beta = NA, acc = acc, ok = TRUE)
  }, error = function(e) { results$IFE_LGS <<- list(gamma = NA, beta = NA, acc = NA, ok = FALSE) })

  results
}

# =================================================================
# 3. MONTE CARLO ENGINE
# =================================================================
run_mc <- function(dgp_cfg, R = 1000, N = 300, T_obs = 6) {
  cat(sprintf("\n=== %s === (R=%d, N=%d, T=%d)\n", dgp_cfg$name, R, N, T_obs))
  true_gamma <- dgp_cfg$gamma
  est_names <- c("AB", "IFE", "CCE", "LGS", "IFE_LGS")

  gamma_mat <- matrix(NA, R, length(est_names))
  colnames(gamma_mat) <- est_names
  acc_mat <- matrix(NA, R, 2)
  colnames(acc_mat) <- c("LGS", "IFE_LGS")

  t0 <- proc.time()[3]
  for (r in 1:R) {
    res <- run_one_rep(dgp_cfg, N, T_obs, seed = r)

    for (en in est_names) {
      if (!is.null(res[[en]]) && res[[en]]$ok && !is.na(res[[en]]$gamma)) {
        gamma_mat[r, en] <- res[[en]]$gamma
      }
    }
    if (!is.null(res$LGS) && !is.na(res$LGS$acc)) acc_mat[r, "LGS"] <- res$LGS$acc
    if (!is.null(res$IFE_LGS) && !is.na(res$IFE_LGS$acc)) acc_mat[r, "IFE_LGS"] <- res$IFE_LGS$acc

    if (r %% 25 == 0) {
      elapsed <- proc.time()[3] - t0
      cat(sprintf("  r=%d/%d  (%.1fs elapsed, ~%.0fs remaining)\n",
                  r, R, elapsed, elapsed / r * (R - r)))
    }
  }

  # Summarise
  summary_df <- data.frame(Estimator = est_names, stringsAsFactors = FALSE)
  summary_df$N_ok  <- apply(gamma_mat, 2, function(x) sum(!is.na(x)))
  summary_df$Bias  <- apply(gamma_mat, 2, function(x) mean(x, na.rm = TRUE) - true_gamma)
  summary_df$SD    <- apply(gamma_mat, 2, function(x) sd(x, na.rm = TRUE))
  summary_df$RMSE  <- apply(gamma_mat, 2, function(x) sqrt(mean((x - true_gamma)^2, na.rm = TRUE)))
  summary_df$Acc   <- c(NA, NA, NA,
                         mean(acc_mat[, "LGS"], na.rm = TRUE),
                         mean(acc_mat[, "IFE_LGS"], na.rm = TRUE))

  cat(sprintf("\n  %-12s %5s %8s %8s %8s %8s\n",
              "Estimator", "N_ok", "Bias", "SD", "RMSE", "Acc"))
  cat(sprintf("  %s\n", paste(rep("-", 55), collapse = "")))
  for (i in 1:nrow(summary_df)) {
    acc_str <- ifelse(is.na(summary_df$Acc[i]), "   ---", sprintf("%8.3f", summary_df$Acc[i]))
    cat(sprintf("  %-12s %5d %+8.4f %8.4f %8.4f %s\n",
                summary_df$Estimator[i], summary_df$N_ok[i],
                summary_df$Bias[i], summary_df$SD[i], summary_df$RMSE[i], acc_str))
  }
  elapsed <- proc.time()[3] - t0
  cat(sprintf("  Time: %.1fs\n", elapsed))

  list(summary = summary_df, gamma_mat = gamma_mat, acc_mat = acc_mat,
       dgp = dgp_cfg$name, R = R, N = N, T_obs = T_obs)
}

# =================================================================
# 4. MAIN EXECUTION
# =================================================================
cat("\n")
cat("###############################################################\n")
cat("#  GMM-6: FULL MONTE CARLO ENGINE                             #\n")
cat("###############################################################\n")

configs <- make_dgp_configs()

# Control parameters — reduce R for testing, increase for publication
R_rep <- 1000   # Use R=500 or 1000 for final tables
N_sim <- 300
T_sim <- 6

all_results <- list()
for (dgp_name in names(configs)) {
  all_results[[dgp_name]] <- run_mc(configs[[dgp_name]], R = R_rep, N = N_sim, T_obs = T_sim)
}

# =================================================================
# 5. COMBINED TABLE
# =================================================================
cat("\n\n")
cat("###############################################################\n")
cat("#  COMBINED RESULTS TABLE                                     #\n")
cat("###############################################################\n")

combined <- data.frame()
for (dgp_name in names(all_results)) {
  s <- all_results[[dgp_name]]$summary
  s$DGP <- dgp_name
  combined <- rbind(combined, s)
}

cat(sprintf("\n  %-6s %-12s %5s %+8s %8s %8s %8s\n",
            "DGP", "Estimator", "N_ok", "Bias", "SD", "RMSE", "Acc"))
cat(sprintf("  %s\n", paste(rep("-", 60), collapse = "")))
for (i in 1:nrow(combined)) {
  acc_str <- ifelse(is.na(combined$Acc[i]), "   ---", sprintf("%8.3f", combined$Acc[i]))
  cat(sprintf("  %-6s %-12s %5d %+8.4f %8.4f %8.4f %s\n",
              combined$DGP[i], combined$Estimator[i], combined$N_ok[i],
              combined$Bias[i], combined$SD[i], combined$RMSE[i], acc_str))
}

# =================================================================
# 6. SAVE TO CSV
# =================================================================
csv_path <- "MC_Tables_R.csv"
write.csv(combined, csv_path, row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", csv_path))

# Try xlsx if openxlsx is available
tryCatch({
  library(openxlsx)
  xlsx_path <- "MC_Tables_R.xlsx"
  wb <- createWorkbook()

  # Sheet 1: Combined
  addWorksheet(wb, "Combined")
  writeData(wb, "Combined", combined)

  # Sheets per DGP
  for (dgp_name in names(all_results)) {
    addWorksheet(wb, dgp_name)
    writeData(wb, dgp_name, all_results[[dgp_name]]$summary)
  }

  # Sheet: Summary comparison (Table 5 format)
  summary5 <- data.frame(
    DGP = character(), gamma_true = numeric(),
    AB_RMSE = numeric(), IFE_RMSE = numeric(), CCE_RMSE = numeric(),
    LGS_RMSE = numeric(), IFE_LGS_RMSE = numeric(),
    LGS_Acc = numeric(), IFE_LGS_Acc = numeric(),
    stringsAsFactors = FALSE
  )
  for (dgp_name in names(configs)) {
    s <- all_results[[dgp_name]]$summary
    row <- data.frame(
      DGP = dgp_name,
      gamma_true = configs[[dgp_name]]$gamma,
      AB_RMSE      = s$RMSE[s$Estimator == "AB"],
      IFE_RMSE     = s$RMSE[s$Estimator == "IFE"],
      CCE_RMSE     = s$RMSE[s$Estimator == "CCE"],
      LGS_RMSE     = s$RMSE[s$Estimator == "LGS"],
      IFE_LGS_RMSE = s$RMSE[s$Estimator == "IFE_LGS"],
      LGS_Acc      = s$Acc[s$Estimator == "LGS"],
      IFE_LGS_Acc  = s$Acc[s$Estimator == "IFE_LGS"],
      stringsAsFactors = FALSE
    )
    summary5 <- rbind(summary5, row)
  }
  addWorksheet(wb, "Table5_Summary")
  writeData(wb, "Table5_Summary", summary5)

  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  cat(sprintf("Excel saved to: %s\n", xlsx_path))
}, error = function(e) {
  cat("(openxlsx not available — skipping Excel output)\n")
  cat("Install with: install.packages('openxlsx')\n")
})

# =================================================================
# 7. FINAL SUMMARY
# =================================================================
cat("\n")
cat("###############################################################\n")
cat("#  KEY FINDINGS                                               #\n")
cat("###############################################################\n")

for (dgp_name in names(all_results)) {
  s <- all_results[[dgp_name]]$summary
  best_idx <- which.min(s$RMSE)
  cat(sprintf("  %s: Best = %s (RMSE=%.4f)\n",
              dgp_name, s$Estimator[best_idx], s$RMSE[best_idx]))
}

# DGP4 comparison
if ("DGP4" %in% names(all_results)) {
  s4 <- all_results[["DGP4"]]$summary
  rmse_ab <- s4$RMSE[s4$Estimator == "AB"]
  rmse_il <- s4$RMSE[s4$Estimator == "IFE_LGS"]
  cat(sprintf("\n  DGP4 RMSE reduction (IFE-LGS vs AB): %.0f%%\n",
              1000 * (1 - rmse_il / rmse_ab)))

  rmse_lgs <- s4$RMSE[s4$Estimator == "LGS"]
  cat(sprintf("  DGP4 LGS without purging: RMSE=%.4f (%s than AB=%.4f)\n",
              rmse_lgs,
              ifelse(rmse_lgs > rmse_ab, "WORSE", "better"),
              rmse_ab))
}

cat("\nGMM-6 COMPLETE.\n")
