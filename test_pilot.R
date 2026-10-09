# =====================================================================
# Sanity tests for the pilot. They check the data and the building blocks,
# not the results.
# Run from the repo folder (the one with A.rds and rolle_hourly.csv):
#   source("tests/test_pilot.R")
# Prints "All 5 tests passed" or stops at the first failure.
# =====================================================================

cat("Running pilot sanity tests\n")
A <- readRDS("A.rds")
Y <- read.csv("rolle_hourly.csv")
M <- as.matrix(Y[, -1])
upper_names <- rownames(A); bottom_names <- colnames(A)

# 1. The hierarchy matrix has the right shape:
#    7 upper series x 24 meters; the total holds 24 meters, each substation 12,
#    each feeder 6; and every meter appears exactly 3 times (feeder, substation, total).
stopifnot(identical(dim(A), c(7L, 24L)))
stopifnot(all(rowSums(A) == c(24, 12, 12, 6, 6, 6, 6)))
stopifnot(all(colSums(A) == 3))
cat("  1. A has the right shape and structure\n")

# 2. The hourly data are complete: 8,928 hours, time + 31 series, no missing values.
stopifnot(nrow(Y) == 8928, ncol(Y) == 32, !anyNA(M))
stopifnot(all(c(upper_names, bottom_names) %in% colnames(M)))
cat("  2. rolle_hourly.csv is complete (8,928 hours, 31 series, no gaps)\n")

# 3. The data are coherent: every upper series equals A x meters, every hour.
gap <- max(abs(M[, upper_names] - M[, bottom_names] %*% t(A)))
stopifnot(gap < 1e-6)
cat("  3. Upper series = A x meters (largest gap", signif(gap, 2), "kW)\n")

# 4. Injected growth keeps the data coherent.
#    Same apply_shift as in the pilot scripts: +20% in meters m01-m05.
apply_shift <- function(W, s, from_row, shift_meters = bottom_names[1:5]) {
  if (s == 0) return(W)
  rows <- from_row:nrow(W)
  W[rows, shift_meters] <- W[rows, shift_meters] * (1 + s)
  W[rows, upper_names]  <- W[rows, bottom_names, drop = FALSE] %*% t(A)
  W
}
W  <- apply_shift(M[1:100, c(upper_names, bottom_names)], 0.20, from_row = 89)
stopifnot(max(abs(W[, upper_names] - W[, bottom_names] %*% t(A))) < 1e-6)
stopifnot(isTRUE(all.equal(W[100, "m01"], 1.2 * M[100, "m01"])))   # grown
stopifnot(isTRUE(all.equal(W[88, "m01"], M[88, "m01"])))            # before the shift: unchanged
cat("  4. Injected growth stays coherent and starts at the right hour\n")

# 5. The interval score gives the textbook values on toy cases (95% interval, alpha = 0.05):
#    inside -> width only; outside -> width + (2 / 0.05) x distance = width + 40 x distance.
interval_score <- function(lo, hi, y, alpha) {
  (hi - lo) + (2 / alpha) * pmax(lo - y, 0) + (2 / alpha) * pmax(y - hi, 0)
}
stopifnot(isTRUE(all.equal(interval_score(0, 10, 5, 0.05), 10)))    # inside
stopifnot(isTRUE(all.equal(interval_score(0, 10, 12, 0.05), 90)))   # 2 above: 10 + 80
stopifnot(isTRUE(all.equal(interval_score(0, 10, -1, 0.05), 50)))   # 1 below: 10 + 40
cat("  5. Interval score matches hand-computed values\n")

cat("All 5 tests passed\n")
