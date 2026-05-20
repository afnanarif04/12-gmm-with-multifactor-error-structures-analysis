###########################################################################
# Simple figure (single self-contained file)
#
# Produces one clean black-and-white figure summarising the Monte Carlo
# comparison between the benchmark estimator and the factor-aware estimator.
# It reads the small results file written by simple_01_monte_carlo.R.
#
# Run simple_01_monte_carlo.R first, then run this file with:
#   Rscript simple_03_figure.R
# Dependencies: base R graphics only.
###########################################################################

# =================================================================
# 1. READ THE MONTE CARLO RESULTS
# =================================================================
results_file <- "simple_monte_carlo_results.csv"
if (!file.exists(results_file)) {
  stop("Run simple_01_monte_carlo.R first to create ", results_file)
}
res <- read.csv(results_file, stringsAsFactors = FALSE)

# Absolute bias for plotting (smaller is better)
gamma_abs_bias <- abs(res$gamma_bias)
beta_abs_bias  <- abs(res$beta_bias)
estimators     <- res$estimator

# =================================================================
# 2. DRAW THE FIGURE
# =================================================================
# A grouped bar-free comparison: we use a dot-and-line (Cleveland) style,
# which is cleaner than bars and matches the manuscript figure aesthetic.

png("simple_figure_bias_comparison.png", width = 1600, height = 950, res = 200)
par(mfrow = c(1, 2), mar = c(4, 6, 4, 1), oma = c(0, 0, 3, 0), family = "serif")

draw_panel <- function(values, labels, title_text) {
  n <- length(values)
  plot(values, seq_len(n), type = "n",
       xlim = c(0, max(values) * 1.25), ylim = c(0.5, n + 0.5),
       yaxt = "n", xlab = "Absolute bias", ylab = "", main = title_text)
  axis(2, at = seq_len(n), labels = labels, las = 1)
  for (i in seq_len(n)) {
    segments(0, i, values[i], i, lwd = 1.5, col = "grey40")
    points(values[i], i, pch = 19, cex = 1.6, col = "black")
    text(values[i], i, sprintf("%.3f", values[i]), pos = 4, cex = 0.9)
  }
  abline(v = 0, col = "grey70", lty = 3)
}

draw_panel(gamma_abs_bias, estimators, "Persistence parameter")
draw_panel(beta_abs_bias,  estimators, "Slope parameter")

mtext("Absolute bias of the benchmark and factor-aware estimators",
      outer = TRUE, line = 0.5, cex = 1.1, font = 2, family = "serif")

dev.off()
cat("Figure saved to simple_figure_bias_comparison.png\n")
cat("Done.\n")
