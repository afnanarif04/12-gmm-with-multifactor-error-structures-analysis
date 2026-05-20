###########################################################################
# Data download: World Bank World Development Indicators
#
# Downloads the six macro indicators used in the two empirical applications
# directly from the World Bank API and writes one CSV per indicator into
# the data/ folder. The repository already ships these CSV files, so this
# script is provided only so that a reviewer can reproduce the raw download.
#
# Indicator codes (World Bank WDI):
#   NY.GDP.PCAP.KD        GDP per capita, constant 2015 US dollars
#   EG.USE.PCAP.KG.OE     Energy use per capita, kg of oil equivalent
#   NE.GDI.TOTL.ZS        Gross capital formation, percent of GDP
#   NE.TRD.GNFS.ZS        Trade openness, percent of GDP
#   SL.TLF.CACT.ZS        Labour force participation, percent
#   EG.FEC.RNEW.ZS        Renewable energy share, percent
#   EN.GHG.CO2.PC.CE.AR5  Carbon dioxide emissions per capita (AR5, EDGAR v8.0)
#
# Source pages (for verification):
#   https://data.worldbank.org/indicator/NY.GDP.PCAP.KD
#   https://data.worldbank.org/indicator/EG.USE.PCAP.KG.OE
#   https://data.worldbank.org/indicator/NE.GDI.TOTL.ZS
#   https://data.worldbank.org/indicator/NE.TRD.GNFS.ZS
#   https://data.worldbank.org/indicator/SL.TLF.CACT.ZS
#   https://data.worldbank.org/indicator/EG.FEC.RNEW.ZS
#   https://data.worldbank.org/indicator/EN.GHG.CO2.PC.CE.AR5
#
# Run with:  Rscript 01_download_data.R
# Dependencies: WDI, dplyr, tidyr
###########################################################################

library(WDI)
library(dplyr)
library(tidyr)

cat("=============================================================\n")
cat("Downloading WDI data for Section 5 (FIXED indicators)\n")
cat("=============================================================\n")

# =================================================================
# STEP 1: Download in TWO batches to handle failures gracefully
# =================================================================

# Batch 1: indicators that definitely work
ind_batch1 <- c(
  gdp_pc     = "NY.GDP.PCAP.KD",      # GDP per capita (constant 2015 USD)
  energy_pc  = "EG.USE.PCAP.KG.OE",   # Energy use per capita (kg oil eq.)
  invest_gdp = "NE.GDI.TOTL.ZS",      # Gross capital formation (% GDP)
  trade_gdp  = "NE.TRD.GNFS.ZS",      # Trade openness (% GDP)
  labor_part = "SL.TLF.CACT.ZS",      # Labor force participation (%)
  renew_share = "EG.FEC.RNEW.ZS"      # Renewable energy share (%)
)

cat("Batch 1: core indicators...\n")
raw1 <- WDI(indicator = ind_batch1, country = "all",
            start = 1990, end = 2019, extra = TRUE)
cat(sprintf("  Downloaded: %d rows\n", nrow(raw1)))

# Batch 2: CO2 — try NEW code first, then OLD
cat("Batch 2: CO2 per capita...\n")
raw_co2 <- tryCatch({
  cat("  Trying NEW code: EN.GHG.CO2.PC.CE.AR5\n")
  r <- WDI(indicator = c(co2_pc = "EN.GHG.CO2.PC.CE.AR5"),
           country = "all", start = 1990, end = 2019, extra = TRUE)
  cat(sprintf("  Success: %d rows\n", nrow(r)))
  r
}, error = function(e) {
  cat("  NEW code failed. Trying OLD code: EN.ATM.CO2E.PC\n")
  tryCatch({
    r <- WDI(indicator = c(co2_pc = "EN.ATM.CO2E.PC"),
             country = "all", start = 1990, end = 2019, extra = TRUE)
    cat(sprintf("  Success: %d rows\n", nrow(r)))
    r
  }, error = function(e2) {
    cat("  OLD code also failed.\n")
    cat("  FALLBACK: Will use CO2 intensity (kg per PPP $ GDP) instead.\n")
    r <- WDI(indicator = c(co2_pc = "EN.ATM.CO2E.PP.GD"),
             country = "all", start = 1990, end = 2019, extra = TRUE)
    cat(sprintf("  Downloaded CO2 intensity: %d rows\n", nrow(r)))
    r
  })
})

# Batch 3: Energy intensity — try multiple codes
cat("Batch 3: Energy intensity...\n")
raw_eint <- tryCatch({
  cat("  Trying: EG.EGY.PRIM.PP.KD\n")
  r <- WDI(indicator = c(energy_int = "EG.EGY.PRIM.PP.KD"),
           country = "all", start = 1990, end = 2019, extra = TRUE)
  cat(sprintf("  Success: %d rows\n", nrow(r)))
  r
}, error = function(e) {
  cat("  Failed. Trying: EG.USE.COMM.GD.PP.KD\n")
  tryCatch({
    r <- WDI(indicator = c(energy_int = "EG.USE.COMM.GD.PP.KD"),
             country = "all", start = 1990, end = 2019, extra = TRUE)
    cat(sprintf("  Success: %d rows\n", nrow(r)))
    r
  }, error = function(e2) {
    cat("  Also failed. Will compute energy_int = energy_pc / gdp_pc later.\n")
    NULL
  })
})

# =================================================================
# STEP 2: MERGE ALL BATCHES
# =================================================================

cat("\nMerging batches...\n")
merge_cols <- c("iso2c", "iso3c", "country", "year", "region",
                "capital", "longitude", "latitude", "income", "lending",
                "status")
# Keep only columns that exist in all
common <- intersect(names(raw1), names(raw_co2))

df <- merge(raw1, raw_co2[, c("iso3c", "year", "co2_pc")],
            by = c("iso3c", "year"), all.x = TRUE)

if (!is.null(raw_eint)) {
  df <- merge(df, raw_eint[, c("iso3c", "year", "energy_int")],
              by = c("iso3c", "year"), all.x = TRUE)
} else {
  # Compute energy intensity as ratio
  df$energy_int <- df$energy_pc / df$gdp_pc
}

cat(sprintf("Merged: %d rows\n", nrow(df)))

# =================================================================
# STEP 3: CLEAN — keep countries only
# =================================================================

df <- df %>%
  filter(region != "Aggregates") %>%
  filter(!is.na(iso3c))

cat(sprintf("After removing aggregates: %d rows, %d countries\n",
            nrow(df), length(unique(df$iso3c))))

# =================================================================
# STEP 4: CREATE 5-YEAR AVERAGES
# =================================================================

df <- df %>%
  mutate(period = case_when(
    year >= 1990 & year <= 1994 ~ "1990-94",
    year >= 1995 & year <= 1999 ~ "1995-99",
    year >= 2000 & year <= 2004 ~ "2000-04",
    year >= 2005 & year <= 2009 ~ "2005-09",
    year >= 2010 & year <= 2014 ~ "2010-14",
    year >= 2015 & year <= 2019 ~ "2015-19",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(period))

panel <- df %>%
  group_by(iso3c, country, period, income, region) %>%
  summarise(across(c(gdp_pc, energy_pc, co2_pc, invest_gdp,
                     trade_gdp, labor_part, energy_int, renew_share),
                   ~mean(.x, na.rm = TRUE)),
            n_years = sum(!is.na(gdp_pc)),
            .groups = "drop")

# Replace NaN with NA
panel <- panel %>% mutate(across(where(is.numeric), ~ifelse(is.nan(.x), NA, .x)))

cat(sprintf("Panel: %d obs, %d countries, %d periods\n",
            nrow(panel), length(unique(panel$iso3c)),
            length(unique(panel$period))))

# =================================================================
# STEP 5: BALANCED PANEL
# =================================================================

# For App 1: need gdp_pc, energy_pc, invest_gdp, trade_gdp, labor_part
app1_countries <- panel %>%
  filter(!is.na(gdp_pc) & gdp_pc > 0 &
         !is.na(energy_pc) & energy_pc > 0 &
         !is.na(invest_gdp) &
         !is.na(trade_gdp) &
         !is.na(labor_part)) %>%
  group_by(iso3c) %>%
  filter(n() == 6) %>%
  ungroup() %>%
  pull(iso3c) %>% unique()

# For App 2: need co2_pc, gdp_pc, energy_int (renew_share may have gaps)
app2_countries <- panel %>%
  filter(!is.na(gdp_pc) & gdp_pc > 0 &
         !is.na(co2_pc) & co2_pc > 0) %>%
  group_by(iso3c) %>%
  filter(n() == 6) %>%
  ungroup() %>%
  pull(iso3c) %>% unique()

cat(sprintf("\nBalanced countries: App1=%d, App2=%d\n",
            length(app1_countries), length(app2_countries)))

# =================================================================
# STEP 6: CREATE LOG VARIABLES
# =================================================================

panel <- panel %>%
  mutate(
    log_gdp_pc      = log(pmax(gdp_pc, 1)),
    log_energy_pc   = log(pmax(energy_pc, 1)),
    log_co2_pc      = log(pmax(co2_pc, 0.01)),
    log_invest_gdp  = log(pmax(invest_gdp, 1)),
    log_trade_gdp   = log(pmax(trade_gdp, 1)),
    log_labor_part  = log(pmax(labor_part, 1)),
    log_energy_int  = log(pmax(energy_int, 0.001)),
    log_renew_share = log(pmax(renew_share, 0.1)),
    log_gdp_pc_sq   = log_gdp_pc^2
  )

app1 <- panel %>% filter(iso3c %in% app1_countries) %>% arrange(iso3c, period)
app2 <- panel %>% filter(iso3c %in% app2_countries) %>% arrange(iso3c, period)

N1 <- length(unique(app1$iso3c)); N2 <- length(unique(app2$iso3c))

# =================================================================
# STEP 7: SAVE TO EXCEL
# =================================================================

library(openxlsx)
wb <- createWorkbook()

addWorksheet(wb, "App1_Energy_Growth")
writeData(wb, "App1_Energy_Growth", app1)

addWorksheet(wb, "App2_CO2_EKC")
writeData(wb, "App2_CO2_EKC", app2)

# Descriptive stats
desc_fn <- function(x) {
  x <- x[!is.na(x)]
  c(N=length(x), Mean=mean(x), SD=sd(x), Min=min(x),
    Median=median(x), Max=max(x))
}

stats1 <- data.frame(
  Variable = c("GDP per capita","Energy per capita","Investment/GDP",
               "Trade/GDP","Labor force part."),
  t(sapply(c("gdp_pc","energy_pc","invest_gdp","trade_gdp","labor_part"),
           function(v) desc_fn(app1[[v]])))
)
addWorksheet(wb, "Stats_App1"); writeData(wb, "Stats_App1", stats1)

stats2 <- data.frame(
  Variable = c("CO2 per capita","GDP per capita","Energy intensity",
               "Renewable share","Trade/GDP"),
  t(sapply(c("co2_pc","gdp_pc","energy_int","renew_share","trade_gdp"),
           function(v) desc_fn(app2[[v]])))
)
addWorksheet(wb, "Stats_App2"); writeData(wb, "Stats_App2", stats2)

# Country list
clist <- panel %>%
  select(iso3c, country, income, region) %>%
  distinct() %>%
  filter(iso3c %in% union(app1_countries, app2_countries)) %>%
  arrange(income, country)
addWorksheet(wb, "Country_List"); writeData(wb, "Country_List", clist)

xlsx_path <- "Section5_WDI_RealData.xlsx"
saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# =================================================================
# STEP 8: SUMMARY
# =================================================================

cat("\n=============================================================\n")
cat("DATA DOWNLOAD COMPLETE\n")
cat("=============================================================\n")
cat(sprintf("File: %s\n", xlsx_path))
cat(sprintf("App 1 (Energy-Growth): N=%d countries, T=6\n", N1))
cat(sprintf("App 2 (CO2-EKC):       N=%d countries, T=6\n", N2))
cat(sprintf("Periods: %s\n", paste(sort(unique(panel$period)), collapse=", ")))
cat("\nIncome distribution (App1):\n")
print(table(app1$income[!duplicated(app1$iso3c)]))
cat("\nIncome distribution (App2):\n")
print(table(app2$income[!duplicated(app2$iso3c)]))

cat("\n=============================================================\n")
cat("NEXT: Run R_empirical_estimation.R on this file\n")
cat("=============================================================\n")

# =================================================================
# DIAGNOSTIC: Check which indicators downloaded successfully
# =================================================================

cat("\n--- Indicator coverage check ---\n")
for (v in c("gdp_pc","energy_pc","co2_pc","invest_gdp","trade_gdp",
            "labor_part","energy_int","renew_share")) {
  n_obs <- sum(!is.na(panel[[v]]))
  n_countries <- length(unique(panel$iso3c[!is.na(panel[[v]])]))
  cat(sprintf("  %-15s: %5d obs, %3d countries\n", v, n_obs, n_countries))
}
