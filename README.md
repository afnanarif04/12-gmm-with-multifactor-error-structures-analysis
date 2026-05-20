# Short Dynamic Panel GMM with Common Factors and Latent Groups

Replication code and data for a set of generalised method of moments
estimators for short dynamic panel data models. The estimators handle two
features that often appear together in cross-country panels: unobserved
common factors that make units move together, and discrete differences in
slope coefficients across hidden groups of units.

This repository contains everything needed to reproduce the simulation and
empirical results, plus a separate folder with simplified scripts for anyone
who wants to understand or apply the methods quickly.

## Repository layout

```
.
├── code/      Full replication scripts (run in numbered order)
├── data/      Raw World Bank CSV files used in the empirical applications
├── results/   Output folder (populated when you run the code)
├── simple/    Simplified, self-contained scripts for quick replication
├── README.md
└── LICENSE
```

## Requirements

All scripts use base R (version 4.0 or later). A few packages are needed for
specific scripts:

- `MASS` for the estimator scripts (usually bundled with R).
- `readxl` for the simplified Excel example in `simple/`.
- `WDI`, `dplyr`, `tidyr` only if you want to re-download the raw data with
  `code/01_download_data.R`. The CSV files are already provided in `data/`,
  so this download step is optional.

Install them once with:

```r
install.packages(c("MASS", "readxl", "WDI", "dplyr", "tidyr"))
```

## How to run the full replication

Open R or a terminal, set the working directory to the `code/` folder, and
run the scripts in order. Each script can also be run from the command line
with `Rscript`.

```bash
cd code
Rscript 01_download_data.R          # optional: re-download raw data
Rscript 02_estimator_ab_gmm.R       # benchmark estimator + self-test
Rscript 03_estimator_gmm_ife.R      # interactive fixed effects estimator
Rscript 04_estimator_gmm_cce.R      # common correlated effects estimator
Rscript 05_estimator_gmm_lgs.R      # latent group structure estimator
Rscript 06_estimator_gmm_ife_lgs.R  # combined estimator
Rscript 07_monte_carlo_engine.R     # full simulation across all designs
Rscript 08_empirical_applications.R # two empirical applications
```

Scripts `02` through `06` each define one estimator and run a short
self-test when executed directly, so you can confirm each estimator works on
its own. Script `07` runs the full Monte Carlo study. Script `08` builds the
two empirical panels from the CSV files in `data/` and writes the coefficient
tables to `results/`.

## Data

The empirical applications use seven indicators from the World Bank World
Development Indicators, downloaded for 1990 to 2019. The raw CSV files are in
`data/`. For verification, the original source pages are:

| File | Indicator | Source |
|------|-----------|--------|
| `wdi_gdp_per_capita.csv` | GDP per capita, constant 2015 USD | https://data.worldbank.org/indicator/NY.GDP.PCAP.KD |
| `wdi_energy_per_capita.csv` | Energy use per capita, kg oil eq. | https://data.worldbank.org/indicator/EG.USE.PCAP.KG.OE |
| `wdi_co2_per_capita.csv` | CO2 emissions per capita (AR5) | https://data.worldbank.org/indicator/EN.GHG.CO2.PC.CE.AR5 |
| `wdi_investment.csv` | Gross capital formation, % GDP | https://data.worldbank.org/indicator/NE.GDI.TOTL.ZS |
| `wdi_trade_openness.csv` | Trade openness, % GDP | https://data.worldbank.org/indicator/NE.TRD.GNFS.ZS |
| `wdi_labour_participation.csv` | Labour force participation, % | https://data.worldbank.org/indicator/SL.TLF.CACT.ZS |
| `wdi_renewable_share.csv` | Renewable energy share, % | https://data.worldbank.org/indicator/EG.FEC.RNEW.ZS |

The carbon dioxide series uses the AR5 framework (EDGAR v8.0), which the
World Bank introduced in December 2024 to replace the previous emissions
series. The download script in `code/01_download_data.R` falls back to the
older code automatically if the new one is unavailable.

## Quick start for new users

If you only want to understand the core idea or apply the method to your own
data, start in the `simple/` folder, which has its own short guide. It
contains three small self-contained scripts: a teaching Monte Carlo, an
example of applying the estimator to a panel stored in an Excel file, and a
script that draws a simple comparison figure.

## License

Released under the MIT License. See `LICENSE`.
