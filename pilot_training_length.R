# =====================================================================
# Follow-up to pilot_core.R: does t-Rec help when the training window is SHORT?
# The t-Rec paper's gains are largest with short training (n = 25) and shrink
# as training grows (n = 110). Run 1 used 672 hours, far into the "long" regime.
# Same data, same origins, same shifts, same seasonal-AR base model; only L changes.
# Needs rolle_hourly.csv and A.rds. Run in RStudio with Source (about 2 minutes).
# =====================================================================

library(bayesRecon)
A <- readRDS("A.rds")
upper_names <- rownames(A); bottom_names <- colnames(A)
Y <- as.matrix(read.csv("rolle_hourly.csv")[, -1])[, c(upper_names, bottom_names)]
n_u <- nrow(A); n <- n_u + ncol(A)

L_values     <- c(48, 96, 168, 672)          # training window in hours: 2, 4, 7, 28 days
season       <- 24
shift_sizes  <- c(0, 0.05, 0.10, 0.20)
shift_meters <- bottom_names[1:5]
shift_lead   <- 12
affected     <- c("all", "S1", "S11")
n_origins    <- 100; step <- 73; level <- 0.95
z <- qnorm(1 - (1 - level) / 2)

# Same origins as pilot_core.R, so every training length is scored on the same hours
origins <- (max(L_values) + season + 1) + (0:(n_origins - 1)) * step

apply_shift <- function(W, s, from_row) {
  if (s == 0) return(W)
  rows <- from_row:nrow(W)
  W[rows, shift_meters] <- W[rows, shift_meters] * (1 + s)
  W[rows, upper_names]  <- W[rows, bottom_names, drop = FALSE] %*% t(A)
  W
}

base_sar <- function(train, L) {
  T <- nrow(train); fc <- numeric(n); res <- matrix(NA, L, n)
  idx <- (season + 1):T
  for (j in 1:n) {
    X <- cbind(1, train[idx - 1, j], train[idx - season, j])
    fit <- .lm.fit(X, train[idx, j])
    res[, j] <- fit$residuals
    fc[j] <- sum(fit$coefficients * c(1, train[T, j], train[T + 1 - season, j]))
  }
  names(fc) <- colnames(train)
  list(fc = fc, res = res)
}

rows <- list(); fails <- 0
for (L in L_values) {
  for (t0 in origins) {
    W0 <- Y[(t0 - L - season):t0, ]
    for (s in shift_sizes) {
      W <- apply_shift(W0, s, nrow(W0) - shift_lead)
      train <- W[-nrow(W), ]; truth <- W[nrow(W), upper_names]
      y_tr <- train[(season + 1):nrow(train), ]
      b <- base_sar(train, L)
      out <- tryCatch(list(
        g  = reconc_gaussian(A, b$fc, residuals = b$res, return_upper = TRUE),
        tr = reconc_t(A, b$fc, y_train = y_tr, residuals = b$res,
                      return_upper = TRUE, freq = season)),
        error = function(e) NULL)
      if (is.null(out)) { fails <- fails + 1; next }
      hw_g <- z * sqrt(diag(out$g$upper_rec_cov))
      hw_t <- qt(1 - (1 - level) / 2, out$tr$upper_rec_df) *
              sqrt(diag(out$tr$upper_rec_scale_matrix))
      rows[[length(rows) + 1]] <- data.frame(
        L = L, shift = s, series = upper_names,
        MinT_hit = abs(truth - out$g$upper_rec_mean) <= hw_g,
        tRec_hit = abs(truth - out$tr$upper_rec_mean) <= hw_t,
        MinT_width = 2 * hw_g, tRec_width = 2 * hw_t,
        tRec_df = out$tr$upper_rec_df)
    }
  }
  cat("Training length", L, "h done\n")
}
d <- do.call(rbind, rows)
d <- d[d$series %in% affected, ]

tab <- aggregate(cbind(MinT_hit, tRec_hit, MinT_width, tRec_width, tRec_df) ~ L + shift,
                 data = d, FUN = mean)
tab <- tab[order(tab$L, tab$shift), ]
tab$shift <- paste0("+", 100 * tab$shift, "%")
names(tab) <- c("train_hours", "shift", "MinT_coverage", "tRec_coverage",
                "MinT_width_kW", "tRec_width_kW", "tRec_df")
tab[, 3:4] <- round(tab[, 3:4], 3); tab[, 5:7] <- round(tab[, 5:7], 1)
cat("\n== Seasonal-AR base | affected series (all, S1, S11) | 95% coverage by training length ==\n")
print(tab, row.names = FALSE)
cat("\nFailed reconciliations:", fails, "\n")
write.csv(tab, "results_training_length.csv", row.names = FALSE)
cat("Saved results_training_length.csv\n")
