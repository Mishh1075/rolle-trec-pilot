# Build the aggregation matrix A for the Rolle hierarchy and check it on the real data.
# In RStudio: File > Open File > build_A.R, then click "Source" (top right of the script).
# First set the working directory to your rolle folder:
#   Session > Set Working Directory > Choose Directory... > rolle

Y <- read.csv("rolle_hourly.csv")
Y <- as.matrix(Y[, -1])            # drop the time column; 31 numeric columns remain
upper_names  <- c("all", "S1", "S2", "S11", "S12", "S21", "S22")
bottom_names <- sprintf("m%02d", 1:24)

# A: 7 rows (upper series) x 24 columns (meters); 1 = meter belongs to that series
A <- matrix(0, nrow = 7, ncol = 24, dimnames = list(upper_names, bottom_names))
A["S11", 1:6]   <- 1
A["S12", 7:12]  <- 1
A["S21", 13:18] <- 1
A["S22", 19:24] <- 1
A["S1", ]  <- A["S11", ] + A["S12", ]
A["S2", ]  <- A["S21", ] + A["S22", ]
A["all", ] <- A["S1", ] + A["S2", ]

print(A[, c(1, 6, 7, 12, 13, 18, 19, 24)])   # a readable slice of A
cat("Row sums (meters per series):", rowSums(A), "\n")

# Check: upper = A x bottom, for every hour
bottom <- Y[, bottom_names]
upper  <- Y[, upper_names]
upper_from_A <- bottom %*% t(A)
gap <- max(abs(upper - upper_from_A))
cat("Largest gap between real upper series and A x bottom:", signif(gap, 3), "kW\n")
if (gap < 0.01) {
  cat("A is correct: the hierarchy is coherent.\n")
} else {
  cat("Gap is larger than expected: send this output to Claude.\n")
}

saveRDS(A, "A.rds")   # saved for the pilot
