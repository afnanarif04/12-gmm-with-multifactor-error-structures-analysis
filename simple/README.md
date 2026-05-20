# Simple scripts for quick replication

This folder contains three small, self-contained R scripts. They are written
for anyone who wants to understand the core idea, reproduce the central
simulation result, or apply the estimator to their own panel data, without
reading the full estimator library in the `code/` folder.

Each script runs on its own with base R. The Excel example also needs the
`readxl` package.

## The three scripts

### 1. `simple_01_monte_carlo.R`

A teaching Monte Carlo. It generates a dynamic panel that contains one
unobserved common factor, then compares two estimators:

- a benchmark first-difference instrumental variable estimator that ignores
  the factor, and
- a factor-aware estimator that adds cross-sectional averages as controls.

The script prints the median bias and median absolute error of each
estimator and saves a small results file. The central finding is that the
benchmark estimator overstates the persistence parameter, while the
factor-aware estimator removes most of that bias.

```bash
Rscript simple_01_monte_carlo.R
```

### 2. `simple_02_estimate_from_excel.R`

Shows how to apply both estimators to a panel stored in an Excel file. It
reads the included `example_panel_data.xlsx`, which holds a balanced
country-by-period panel with a dependent variable and two regressors.

To use your own data, edit the three settings lines near the top of the
script (the file name, the dependent variable, and the regressors), and make
sure your Excel file has the same column layout: one row per country-period,
with columns `iso3c`, `period`, and one column per variable.

```bash
Rscript simple_02_estimate_from_excel.R
```

### 3. `simple_03_figure.R`

Reads the results file written by the first script and draws one clean
black-and-white figure comparing the absolute bias of the two estimators for
both parameters. Run the Monte Carlo script first.

```bash
Rscript simple_01_monte_carlo.R
Rscript simple_03_figure.R
```

## Example data

`example_panel_data.xlsx` is a balanced panel of countries observed over six
five-year periods, with log GDP per capita as the dependent variable and log
energy use and log investment as regressors. It is provided only as a
worked example of the expected input format; replace it with your own file
to estimate your own model.

## A note on the simplified estimators

These scripts use a compact version of the estimator that keeps the code
short and readable. The full library in the `code/` folder uses the complete
instrument set, two-step efficient weighting, and proper group
classification. For any serious application, use the full library. The simple
scripts are for learning and quick checks only.
