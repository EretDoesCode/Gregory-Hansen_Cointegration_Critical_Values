# =====================================================================
# Gregory-Hansen (1996) ADF* cointegration test with simulated critical
# values for ANY number of regressors m (built for m = 5).
#
# Models (Gregory & Hansen 1996):
#   "C"  : level shift
#   "CT" : level shift with trend
#   "CS" : regime shift (intercept and all slopes shift)
#
# Only base R is needed.
# 
# This code was written and tested only to get critical values for the
# Gregory-Hansen cointegration test with more than 4 regressors. The codes
# outputs were tested against the tables given in the original article.
# In the future I could write additional code to also conduct the test on
# real data, but right now I am very low on time and I need to get this done ASAP.
# If anyone would like to use the functions written here in their r scripts, feel
# free to do so. 
# =====================================================================

set.seed(12345)

# ---------------- SETTINGS ----------------
REPS         <- 10000     # use 5000-10000 for final critical values
N_OBS        <- 1000      # sample size
M_REG        <- 2         # number of regressors (excluding constant/trend/dummies)
TRIM         <- 0.15      # standard trimming for candidate break dates
PMAX         <- 4         # max lags in the ADF regression (AIC-selected)
MODELS       <- c("C", "CT", "CS")
STATS        <- c("ADF", "Za", "Zt")

# ---------------- ADF t-statistic on residuals ----------------
adf_stat <- function(e, pmax = PMAX) {
  n  <- length(e)
  de <- diff(e)
  y  <- de[(pmax + 1):(n - 1)]
  elag <- e[(pmax + 1):(n - 1)]
  N  <- length(y)
  best_aic <- Inf; best_t <- NA_real_
  for (p in 0:pmax) {
    X <- matrix(elag, ncol = 1)
    if (p > 0) for (j in 1:p) X <- cbind(X, de[(pmax + 1 - j):(n - 1 - j)])
    XtXi <- tryCatch(solve(crossprod(X)), error = function(err) NULL)
    if (is.null(XtXi)) next
    b   <- XtXi %*% crossprod(X, y)
    rss <- sum((y - X %*% b)^2)
    aic <- N * log(rss / N) + 2 * (p + 1)
    if (aic < best_aic) {
      best_aic <- aic
      best_t <- b[1] / sqrt(rss / (N - ncol(X)) * XtXi[1, 1])
    }
  }
  best_t
}

# ---------------- Phillips Za and Zt on residuals ----------------
# e_t = rho * e_{t-1} + u_t (no constant); Bartlett long-run variance.
pp_stats <- function(e) {
  n <- length(e); N <- n - 1
  y <- e[2:n]; x <- e[1:(n - 1)]
  sxx <- sum(x * x)
  rho <- sum(x * y) / sxx
  u   <- y - rho * x
  s2  <- sum(u^2) / (N - 1)
  se  <- sqrt(s2 / sxx)
  q   <- floor(4 * (N / 100)^(2 / 9))
  gam <- function(j) sum(u[(j + 1):N] * u[1:(N - j)]) / N
  g0  <- gam(0)
  lam2 <- g0
  if (q > 0) for (j in 1:q) lam2 <- lam2 + 2 * (1 - j / (q + 1)) * gam(j)
  tstat <- (rho - 1) / se
  Za <- N * (rho - 1) - 0.5 * (N^2 * se^2 / g0) * (lam2 - g0)
  Zt <- sqrt(g0 / lam2) * tstat - 0.5 * (lam2 - g0) / sqrt(lam2) * (N * se / sqrt(g0))
  c(Za = Za, Zt = Zt)
}

# ---------------- Gregory-Hansen: all three statistics ----------------
# Returns, for each statistic, the minimum over candidate break dates
# and the break observation where that minimum occurs.
gh_all <- function(y, X, model = c("C", "CT", "CS"), trim = TRIM, pmax = PMAX) {
  model <- match.arg(model)
  X <- as.matrix(X)
  n <- length(y)
  trend <- seq_len(n)
  cand  <- ceiling(trim * n):floor((1 - trim) * n)
  mat <- matrix(NA_real_, length(cand), 3, dimnames = list(NULL, STATS))
  for (k in seq_along(cand)) {
    phi <- as.numeric(trend > cand[k])
    Z <- switch(model,
                C  = cbind(1, phi, X),
                CT = cbind(1, phi, trend, X),
                CS = cbind(1, phi, X, X * phi))
    e <- lm.fit(Z, y)$residuals
    mat[k, "ADF"] <- adf_stat(e, pmax)
    pp <- pp_stats(e)
    mat[k, "Za"] <- pp["Za"]
    mat[k, "Zt"] <- pp["Zt"]
  }
  list(stat = apply(mat, 2, min, na.rm = TRUE),
       break_index = cand[apply(mat, 2, which.min)])
}

# ---------------- Monte Carlo critical values ----------------
# Returns a 3 x 3 matrix: rows = ADF/Za/Zt, columns = 1%, 5%, 10%.
sim_cv <- function(n, m, model, reps) {
  out <- matrix(NA_real_, reps, 3, dimnames = list(NULL, STATS))
  for (r in seq_len(reps)) {
    X <- apply(matrix(rnorm(n * m), n, m), 2, cumsum)
    y <- cumsum(rnorm(n))
    out[r, ] <- gh_all(y, X, model)$stat
  }
  t(apply(out, 2, quantile, probs = c(0.01, 0.05, 0.10), names = FALSE))
}

# ---------------- RUN: critical values ----------------
m_use <- M_REG
cvs <- list()
for (mod in MODELS) {
  cat("Simulating model", mod, "with m =", m_use, "...\n")
  cv <- sim_cv(N_OBS, m_use, mod, REPS)
  dimnames(cv) <- list(STATS, c("1%", "5%", "10%"))
  cvs[[mod]] <- cv
  cat("\nModel", mod, "| n =", N_OBS, "| m =", m_use, "| reps =", REPS, "\n")
  print(round(cv, 3))
  cat("\n")
}
