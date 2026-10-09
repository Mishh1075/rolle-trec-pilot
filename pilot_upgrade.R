# =====================================================================
# Pilot upgrade: fuller scores and the slide 8 plot.
# Same data, origins, growth and seasonal-AR base model as pilot_core.R,
# for two training windows (48 h and 672 h). Adds:
#   - coverage at 80, 90 and 95%
#   - interval score (MIS) at 95% (lower is better)
#   - results per hierarchy level (total, substation, feeder) and for the
#     affected series (all, S1, S11)
#   - coverage_vs_growth.png for slide 8
# Needs rolle_hourly.csv and A.rds. Run in RStudio with Source (about 2 minutes).
# =====================================================================

library(bayesRecon)
A <- readRDS("A.rds")
upper_names <- rownames(A); bottom_names <- colnames(A)
Y <- as.matrix(read.csv("rolle_hourly.csv")[, -1])[, c(upper_names, bottom_names)]
n_u <- nrow(A); n <- n_u + ncol(A)

L_values     <- c(48, 672)
season       <- 24
shift_sizes  <- c(0, 0.05, 0.10, 0.20)
shift_meters <- bottom_names[1:5]
shift_lead   <- 12
affected     <- c("all", "S1", "S11")
level_of     <- c(all = "total", S1 = "substation", S2 = "substation",
                  S11 = "feeder", S12 = "feeder", S21 = "feeder", S22 = "feeder")
levels_cov   <- c(0.80, 0.90, 0.95)
n_origins    <- 100; step <- 73
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

# Interval score for a central (1 - alpha) interval [lo, hi] and truth y
interval_score <- function(lo, hi, y, alpha) {
  (hi - lo) + (2 / alpha) * pmax(lo - y, 0) + (2 / alpha) * pmax(y - hi, 0)
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
      mu <- list(Base = b$fc[upper_names],
                 MinT = out$g$upper_rec_mean,
                 tRec = out$tr$upper_rec_mean)
      sd <- list(Base = apply(b$res[, 1:n_u], 2, sd),
                 MinT = sqrt(diag(out$g$upper_rec_cov)),
                 tRec = sqrt(diag(out$tr$upper_rec_scale_matrix)))
      for (m in names(mu)) {
        for (lv in levels_cov) {
          a <- 1 - lv
          q <- if (m == "tRec") qt(1 - a / 2, out$tr$upper_rec_df) else qnorm(1 - a / 2)
          lo <- mu[[m]] - q * sd[[m]]; hi <- mu[[m]] + q * sd[[m]]
          rows[[length(rows) + 1]] <- data.frame(
            L = L, shift = s, origin = t0, series = upper_names,
            level = unname(level_of[upper_names]), method = m, nominal = lv,
            hit = (truth >= lo) & (truth <= hi), width = hi - lo,
            mis = interval_score(lo, hi, truth, a))
        }
      }
    }
  }
  cat("Training length", L, "h done\n")
}
d <- do.call(rbind, rows)
write.csv(d, "results_upgrade_detail.csv", row.names = FALSE)

# ---------- Tables ----------
summarise <- function(x) {
  agg <- aggregate(cbind(hit, width, mis) ~ L + shift + method + nominal, data = x, FUN = mean)
  agg[order(agg$L, agg$nominal, agg$method, agg$shift), ]
}
aff <- summarise(d[d$series %in% affected, ])
names(aff)[5:7] <- c("coverage", "width_kW", "interval_score")
aff[, 5:7] <- round(aff[, 5:7], 3)
cat("\n== Affected series (all, S1, S11): coverage, width and interval score ==\n")
print(aff[aff$nominal == 0.95, ], row.names = FALSE)
write.csv(aff, "results_upgrade_affected.csv", row.names = FALSE)

by_level <- aggregate(cbind(hit, mis) ~ L + shift + method + level,
                      data = d[d$nominal == 0.95, ], FUN = mean)
names(by_level)[5:6] <- c("coverage95", "interval_score95")
by_level[, 5:6] <- round(by_level[, 5:6], 3)
by_level <- by_level[order(by_level$L, by_level$level, by_level$method, by_level$shift), ]
cat("\n== 95% coverage and interval score by hierarchy level (all 7 series) ==\n")
print(by_level, row.names = FALSE)
write.csv(by_level, "results_upgrade_by_level.csv", row.names = FALSE)

# ---------- Slide 8 plot ----------
cols <- c(Base = "#8a8f98", MinT = "#d9822b", tRec = "#2b6cb0")
png("coverage_vs_growth.png", width = 1600, height = 800, res = 160)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
for (L in rev(L_values)) {
  p <- aff[aff$L == L & aff$nominal == 0.95, ]
  plot(NA, xlim = c(0, 20), ylim = c(0.75, 1), xaxt = "n",
       xlab = "Load growth in 5 meters (%)", ylab = "95% coverage",
       main = paste0("Training window: ", L, " h"))
  axis(1, at = c(0, 5, 10, 20))
  abline(h = 0.95, lty = 2, col = "#555555")
  for (m in names(cols)) {
    pm <- p[p$method == m, ]
    lines(100 * pm$shift, pm$coverage, col = cols[m], lwd = 3, type = "b", pch = 19)
  }
  legend("bottomleft", legend = c("No reconciliation", "MinT", "t-Rec", "Nominal 0.95"),
         col = c(cols, "#555555"), lwd = c(3, 3, 3, 1), lty = c(1, 1, 1, 2), bty = "n")
}
dev.off()

cat("\nFailed reconciliations:", fails, "\n")
cat("Saved results_upgrade_affected.csv, results_upgrade_by_level.csv,",
    "results_upgrade_detail.csv and coverage_vs_growth.png\n")
