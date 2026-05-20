###########################################################################
# Empirical applications: energy-growth nexus and the environmental
# Kuznets curve for carbon dioxide emissions
#
# Loads the seven World Bank CSV files from the data/ folder, builds two
# balanced five-year-average panels, applies all five estimators, and writes
# the coefficient estimates to the results/ folder.
#
# Requires 02 through 06 in the same directory and the CSV files in ../data.
#
# Run with:  Rscript 08_empirical_applications.R
# Dependencies: base R and MASS only.
###########################################################################

# --- Load the five estimators -------------------------------------------
# Each estimator file runs a short self-test when sourced; the messages can
# be ignored. We only need the functions they define.
e1 <- new.env(); sys.source("02_estimator_ab_gmm.R",      envir = e1)
e2 <- new.env(); sys.source("03_estimator_gmm_ife.R",     envir = e2)
e3 <- new.env(); sys.source("04_estimator_gmm_cce.R",     envir = e3)
e4 <- new.env(); sys.source("05_estimator_gmm_lgs.R",     envir = e4)
e5 <- new.env(); sys.source("06_estimator_gmm_ife_lgs.R", envir = e5)

arellano_bond_gmm <- e1$arellano_bond_gmm
gmm_ife           <- e2$gmm_ife
gmm_cce_iv        <- e3$gmm_cce_iv
gmm_lgs           <- e4$gmm_lgs
gmm_ife_lgs       <- e5$gmm_ife_lgs

DATA_DIR    <- file.path("..", "data")
RESULTS_DIR <- file.path("..", "results")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# =================================================================
# 1. LOAD AND RESHAPE THE WORLD BANK CSV FILES
# =================================================================
# The World Bank CSV files have four metadata rows at the top, then one row
# per country with one column per year. We read each file, drop the metadata,
# and reshape to long format (country, year, value).

read_wdi_csv <- function(path, value_name) {
  raw <- read.csv(path, skip = 4, check.names = FALSE, stringsAsFactors = FALSE)
  year_cols <- as.character(1990:2019)
  year_cols <- year_cols[year_cols %in% names(raw)]
  keep <- c("Country Name", "Country Code", year_cols)
  raw <- raw[, keep]
  long <- reshape(raw,
                  varying   = year_cols,
                  v.names   = value_name,
                  timevar   = "year",
                  times     = as.integer(year_cols),
                  idvar     = c("Country Name", "Country Code"),
                  direction = "long")
  names(long)[names(long) == "Country Code"] <- "iso3c"
  rownames(long) <- NULL
  long[, c("iso3c", "year", value_name)]
}

gdp   <- read_wdi_csv(file.path(DATA_DIR, "wdi_gdp_per_capita.csv"),     "gdp_pc")
energy<- read_wdi_csv(file.path(DATA_DIR, "wdi_energy_per_capita.csv"),  "energy_pc")
co2   <- read_wdi_csv(file.path(DATA_DIR, "wdi_co2_per_capita.csv"),     "co2_pc")
invest<- read_wdi_csv(file.path(DATA_DIR, "wdi_investment.csv"),         "invest_pct")
trade <- read_wdi_csv(file.path(DATA_DIR, "wdi_trade_openness.csv"),     "trade_pct")
labor <- read_wdi_csv(file.path(DATA_DIR, "wdi_labour_participation.csv"),"labor_pct")
renew <- read_wdi_csv(file.path(DATA_DIR, "wdi_renewable_share.csv"),    "renew_pct")

merge_all <- function(...) Reduce(function(a, b) merge(a, b, by = c("iso3c", "year")), list(...))
panel <- merge_all(gdp, energy, co2, invest, trade, labor, renew)

# Five-year periods: 1990-94, 1995-99, ..., 2015-19
panel$period <- ((panel$year - 1990) %/% 5) + 1

# =================================================================
# 2. BUILD A BALANCED FIVE-YEAR-AVERAGE PANEL
# =================================================================
# Average each variable within each five-year window, then keep only the
# countries that have a complete record for all six windows.

build_balanced_panel <- function(df, vars) {
  agg <- aggregate(df[, vars], by = list(iso3c = df$iso3c, period = df$period),
                   FUN = function(z) mean(z, na.rm = TRUE))
  agg <- agg[complete.cases(agg), ]
  agg <- agg[is.finite(rowSums(agg[, vars, drop = FALSE])), ]
  counts <- table(agg$iso3c)
  keep   <- names(counts[counts == 6])
  agg    <- agg[agg$iso3c %in% keep, ]
  agg[order(agg$iso3c, agg$period), ]
}

# Reshape a long balanced panel into the (y, x) array format expected by the
# estimators: y is N x (T+1); x is N x (T+1) x k.
to_estimator_format <- function(df, dv, regs) {
  ids   <- sort(unique(df$iso3c))
  times <- sort(unique(df$period))
  N <- length(ids); Tp1 <- length(times); k <- length(regs)
  y <- matrix(NA_real_, N, Tp1)
  x <- array(NA_real_, c(N, Tp1, k))
  for (i in seq_along(ids)) {
    sub <- df[df$iso3c == ids[i], ]
    sub <- sub[order(sub$period), ]
    y[i, ] <- sub[[dv]]
    for (j in seq_along(regs)) x[i, , j] <- sub[[regs[j]]]
  }
  list(y = y, x = x, N = N, T_obs = Tp1 - 1, k = k, m = 0, G = 1,
       ids = ids)
}

# --- Application 1: energy and economic growth --------------------------
p1 <- panel
p1$log_gdp    <- log(p1$gdp_pc)
p1$log_energy <- log(p1$energy_pc)
p1$log_invest <- log(p1$invest_pct)
p1$log_trade  <- log(p1$trade_pct)
p1$log_labor  <- log(p1$labor_pct)
vars1 <- c("log_gdp", "log_energy", "log_invest", "log_trade", "log_labor")
p1bal <- build_balanced_panel(p1, vars1)
data1 <- to_estimator_format(p1bal, dv = "log_gdp",
                             regs = c("log_energy", "log_invest", "log_trade", "log_labor"))
cat(sprintf("Application 1: %d countries, %d periods\n", data1$N, data1$T_obs + 1))

# --- Application 2: carbon dioxide and the Kuznets curve ----------------
p2 <- panel
p2$log_co2    <- log(pmax(p2$co2_pc, 1e-3))
p2$log_gdp    <- log(p2$gdp_pc)
p2$log_gdp_sq <- p2$log_gdp^2
p2$log_renew  <- log(pmax(p2$renew_pct, 1e-2))
p2$log_trade  <- log(p2$trade_pct)
vars2 <- c("log_co2", "log_gdp", "log_gdp_sq", "log_renew", "log_trade")
p2bal <- build_balanced_panel(p2, vars2)
data2 <- to_estimator_format(p2bal, dv = "log_co2",
                             regs = c("log_gdp", "log_gdp_sq", "log_renew", "log_trade"))
cat(sprintf("Application 2: %d countries, %d periods\n", data2$N, data2$T_obs + 1))

# =================================================================
# 3. RUN ALL FIVE ESTIMATORS ON EACH APPLICATION
# =================================================================
run_all_estimators <- function(data, reg_names) {
  out <- list()
  ab  <- tryCatch(arellano_bond_gmm(data), error = function(e) NULL)
  ife <- tryCatch(gmm_ife(data, m_assumed = 1), error = function(e) NULL)
  cce <- tryCatch(gmm_cce_iv(data), error = function(e) NULL)
  lgs <- tryCatch(gmm_lgs(data, G = 2), error = function(e) NULL)
  il  <- tryCatch(gmm_ife_lgs(data, m_assumed = 1, G = 2), error = function(e) NULL)

  fmt <- function(res) {
    if (is.null(res)) return(rep(NA, 1 + length(reg_names)))
    if (!is.null(res$theta)) return(res$theta)
    c(res$gamma_hat, res$beta_hat)
  }
  tab <- data.frame(parameter = c("gamma (persistence)", reg_names))
  tab$AB_GMM      <- fmt(ab)
  tab$GMM_IFE     <- fmt(ife)
  tab$GMM_CCE_IV  <- fmt(cce)
  # Group estimators report per-group coefficients; report group 1 here.
  g1 <- function(res) {
    if (is.null(res) || is.null(res$group_gamma)) return(rep(NA, 1 + length(reg_names)))
    c(res$group_gamma[1], res$group_beta[1, ])
  }
  tab$GMM_LGS_G1     <- g1(lgs)
  tab$GMM_IFE_LGS_G1 <- g1(il)
  tab
}

cat("\nRunning estimators for Application 1 ...\n")
res1 <- run_all_estimators(data1, c("log_energy", "log_invest", "log_trade", "log_labor"))
print(res1, row.names = FALSE)

cat("\nRunning estimators for Application 2 ...\n")
res2 <- run_all_estimators(data2, c("log_gdp", "log_gdp_sq", "log_renew", "log_trade"))
print(res2, row.names = FALSE)

# =================================================================
# 4. WRITE RESULTS
# =================================================================
write.csv(res1, file.path(RESULTS_DIR, "application_1_energy_growth.csv"), row.names = FALSE)
write.csv(res2, file.path(RESULTS_DIR, "application_2_co2_ekc.csv"),       row.names = FALSE)
cat("\nResults written to the results/ folder.\n")
cat("Done.\n")
