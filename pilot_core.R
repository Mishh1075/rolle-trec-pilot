# =====================================================================
# Pilot: do t-Rec intervals keep their coverage better than MinT
#        when load grows in part of the grid?
# Data: Rolle distribution grid (Nespoli, Medici et al.), hourly.
# Needs: rolle_hourly.csv and A.rds in the working directory
#        (made by make_hourly.py and build_A.R).
# Run in RStudio: Session > Set Working Directory > rolle, then Source.
# Takes a few minutes. Prints the results tables and saves two CSVs.
# =====================================================================

library(bayesRecon)
start_time <- Sys.time()

# ---------- 1. Data and hierarchy ----------
A <- readRDS("A.rds")                       # 7 x 24 aggregation matrix
upper_names  <- rownames(A)                 # all, S1, S2, S11, S12, S21, S22
bottom_names <- colnames(A)                 # m01 ... m24
Y <- as.matrix(read.csv("rolle_hourly.csv")[, -1])
Y <- Y[, c(upper_names, bottom_names)]      # upper first, then bottom (bayesRecon order)
n_u <- nrow(A); n_b <- ncol(A); n <- n_u + n_b

# ---------- 2. Settings (all fixed before looking at any result) ----------
L            <- 24 * 28                     # training window: 28 days of hourly residuals
season       <- 24                          # daily seasonality in hourly data
shift_sizes  <- c(0, 0.05, 0.10, 0.20)      # load growth: none, +5%, +10%, +20%
shift_meters <- bottom_names[1:5]           # 5 meters that grow (all in feeder S11)
shift_lead   <- 12                          # growth started 12 hours before the forecast hour
affected     <- c("all", "S1", "S11")       # upper series that contain the growing meters
n_origins    <- 100                         # number of rolling forecast origins
step         <- 73                          # hours between origins (73, so the hour of day rotates)
level        <- 0.95
z            <- qnorm(1 - (1 - level) / 2)

first_origin <- L + season + 1              # first row we can forecast
origins <- first_origin + (0:(n_origins - 1)) * step
stopifnot(max(origins) <= nrow(Y))

# ---------- 3. Helper: inject load growth into a data window ----------
# Multiplies the 5 meters by (1 + s) from row `from_row` onward,
# then recomputes the upper series as A x bottom so the data stay coherent.
apply_shift <- function(W, s, from_row) {
  if (s == 0) return(W)
  rows <- from_row:nrow(W)
  W[rows, shift_meters] <- W[rows, shift_meters] * (1 + s)
  W[rows, upper_names]  <- W[rows, bottom_names, drop = FALSE] %*% t(A)
  W
}

# ---------- 4. Base forecasts: two simple models ----------
# Each returns the one-step-ahead forecast for all 31 series (fc)
# and L in-sample residuals per series (res).
# `train` has L + season rows; the target hour is the next row.

# (a) Seasonal naive: forecast = same hour yesterday.
base_snaive <- function(train) {
  T <- nrow(train)
  fc  <- train[T + 1 - season, ]
  res <- train[(season + 1):T, ] - train[1:(T - season), ]
  list(fc = fc, res = res)
}

# (b) Seasonal AR, fitted separately per series:
#     y(t) = a + b * y(t-1) + c * y(t-24), by least squares.
#     Each series gets its own coefficients, so the 31 forecasts no longer add up.
base_sar <- function(train) {
  T <- nrow(train)
  fc  <- numeric(n)
  res <- matrix(NA, nrow = L, ncol = n)
  idx <- (season + 1):T                      # rows with both lags available
  for (j in 1:n) {
    y <- train[idx, j]
    X <- cbind(1, train[idx - 1, j], train[idx - season, j])
    fit <- .lm.fit(X, y)
    res[, j] <- fit$residuals
    fc[j] <- sum(fit$coefficients * c(1, train[T, j], train[T + 1 - season, j]))
  }
  names(fc) <- colnames(train)
  list(fc = fc, res = res)
}

base_models <- list(snaive = base_snaive, sar = base_sar)

# ---------- 5. Main loop ----------
rows <- list()
fails <- 0
for (k in seq_along(origins)) {
  t0 <- origins[k]
  W0 <- Y[(t0 - L - season):t0, ]           # training rows + target row (last)
  for (s in shift_sizes) {
    W <- apply_shift(W0, s, from_row = nrow(W0) - shift_lead)
    train <- W[-nrow(W), ]
    truth <- W[nrow(W), upper_names]
    y_tr  <- train[(season + 1):nrow(train), ]   # last L rows, same as residuals
    for (bm in names(base_models)) {
      b <- base_models[[bm]](train)
      out <- tryCatch({
        g  <- reconc_gaussian(A, b$fc, residuals = b$res, return_upper = TRUE)
        tr <- reconc_t(A, b$fc, y_train = y_tr, residuals = b$res,
                       return_upper = TRUE, freq = season)
        list(g = g, tr = tr)
      }, error = function(e) NULL)
      if (is.null(out)) { fails <- fails + 1; next }

      # Base (no reconciliation): Gaussian interval from each series' own residuals
      mu_base <- b$fc[upper_names]
      hw_base <- z * apply(b$res[, 1:n_u], 2, sd)
      # MinT-style Gaussian reconciliation
      mu_g <- out$g$upper_rec_mean
      hw_g <- z * sqrt(diag(out$g$upper_rec_cov))
      # t-Rec: marginal of a multivariate t is a univariate t
      mu_t <- out$tr$upper_rec_mean
      hw_t <- qt(1 - (1 - level) / 2, df = out$tr$upper_rec_df) *
              sqrt(diag(out$tr$upper_rec_scale_matrix))

      for (m in c("Base", "MinT", "tRec")) {
        mu <- switch(m, Base = mu_base, MinT = mu_g, tRec = mu_t)
        hw <- switch(m, Base = hw_base, MinT = hw_g, tRec = hw_t)
        rows[[length(rows) + 1]] <- data.frame(
          base_model = bm, shift = s, origin = t0, series = upper_names,
          method = m, hit = abs(truth - mu) <= hw, width = 2 * hw)
      }
    }
  }
  if (k %% 10 == 0) cat("Origins done:", k, "of", n_origins, "\n")
}
res_all <- do.call(rbind, rows)
write.csv(res_all, "results_detail.csv", row.names = FALSE)

# ---------- 6. Results tables ----------
make_table <- function(d, label) {
  cov <- aggregate(hit ~ shift + method, data = d, FUN = mean)
  wid <- aggregate(width ~ shift + method, data = d, FUN = mean)
  tab <- reshape(merge(cov, wid), idvar = "shift", timevar = "method", direction = "wide")
  tab <- tab[order(tab$shift), ]
  out <- data.frame(shift = paste0("+", 100 * tab$shift, "%"))
  for (m in c("Base", "MinT", "tRec")) {
    out[[paste0(m, "_coverage")]] <- round(tab[[paste0("hit.", m)]], 3)
  }
  for (m in c("Base", "MinT", "tRec")) {
    out[[paste0(m, "_width_kW")]] <- round(tab[[paste0("width.", m)]], 1)
  }
  cat("\n", label, "\n", sep = "")
  print(out, row.names = FALSE)
  out$table <- label
  out
}

tables <- list()
for (bm in names(base_models)) {
  d <- res_all[res_all$base_model == bm, ]
  tables[[length(tables) + 1]] <- make_table(
    d, paste0("== Base model: ", bm, " | all 7 upper series | 95% coverage and mean width =="))
  tables[[length(tables) + 1]] <- make_table(
    d[d$series %in% affected, ],
    paste0("== Base model: ", bm, " | affected series (all, S1, S11) only =="))
}
write.csv(do.call(rbind, tables), "results_core.csv", row.names = FALSE)

cat("\nSetup: Rolle hourly data, ", n_origins, " rolling origins (every ", step,
    " h), one-step-ahead, training window ", L, " h, growth in ",
    length(shift_meters), " meters starting ", shift_lead, " h before the forecast hour.\n", sep = "")
cat("Failed reconciliations:", fails, "\n")
cat("Run time:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), "minutes\n")
cat("Saved results_core.csv and results_detail.csv\n")
