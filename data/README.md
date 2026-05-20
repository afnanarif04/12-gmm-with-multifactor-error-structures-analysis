# Data

This folder contains the raw input series used in the two empirical
applications. Every file is downloaded directly from the World Bank World
Development Indicators database and is redistributed here without
modification, under the World Bank's Creative Commons Attribution 4.0
license.

Each file is a standard World Bank CSV export: four header rows, then one row
per country with one column per year from 1990 to 2019.

## Files and exact sources

| File in this folder | Indicator | World Bank code | Source page |
|---------------------|-----------|-----------------|-------------|
| `wdi_gdp_per_capita.csv` | GDP per capita, constant 2015 US dollars | `NY.GDP.PCAP.KD` | <https://data.worldbank.org/indicator/NY.GDP.PCAP.KD> |
| `wdi_energy_per_capita.csv` | Energy use per capita, kg of oil equivalent | `EG.USE.PCAP.KG.OE` | <https://data.worldbank.org/indicator/EG.USE.PCAP.KG.OE> |
| `wdi_co2_per_capita.csv` | Carbon dioxide emissions per capita (AR5, EDGAR v8.0) | `EN.GHG.CO2.PC.CE.AR5` | <https://data.worldbank.org/indicator/EN.GHG.CO2.PC.CE.AR5> |
| `wdi_investment.csv` | Gross capital formation, percent of GDP | `NE.GDI.TOTL.ZS` | <https://data.worldbank.org/indicator/NE.GDI.TOTL.ZS> |
| `wdi_trade_openness.csv` | Trade openness, percent of GDP | `NE.TRD.GNFS.ZS` | <https://data.worldbank.org/indicator/NE.TRD.GNFS.ZS> |
| `wdi_labour_participation.csv` | Labour force participation rate, percent | `SL.TLF.CACT.ZS` | <https://data.worldbank.org/indicator/SL.TLF.CACT.ZS> |
| `wdi_renewable_share.csv` | Renewable energy share of final consumption, percent | `EG.FEC.RNEW.ZS` | <https://data.worldbank.org/indicator/EG.FEC.RNEW.ZS> |

## Notes on the carbon dioxide series

The World Bank retired the previous carbon dioxide per capita series
(`EN.ATM.CO2E.PC`) at the end of 2024 and replaced it with the AR5 series
(`EN.GHG.CO2.PC.CE.AR5`), which is drawn from the EDGAR v8.0 database. This
repository uses the current AR5 series. The download script in
`code/01_download_data.R` tries the new code first and falls back to the old
code only if the new one is unavailable.

## Reproducing the download

The CSV files are already provided, so no download is needed to run the
analysis. To reproduce the raw download from scratch, run:

```sh
Rscript ../code/01_download_data.R
```

This requires the `WDI`, `dplyr`, and `tidyr` packages. Because the World Bank
updates its database periodically, a fresh download may differ slightly from
the files stored here if any series has been revised since these files were
retrieved.
