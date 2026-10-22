## Calibration of many-pair tests (Section 6.6 and 7 of the paper).
##
## With many pairs, nuisance-estimation errors are shared by all pairs that involve a feature, which
## fattens the tail of the null distribution of the z-statistics while leaving its centre close to N(0, 1)
## (study S8b). Two calibrations are provided:
##   empnull()   : Efron's empirical null, N(delta, sigma^2) fitted to the centre of the z-statistics
##                 (median and interquartile range); it corrects a shifted or rescaled centre, not the tails.
##   perm_null() : a residual-permutation null. Each feature's conditional mean and scale are fitted on all
##                 rows, the rows of the matrix of standardised residuals are permuted, and data are rebuilt
##                 with the fitted means and scales. Means and scales keep their dependence on the
##                 covariates and the residual correlation structure is kept marginally, but every pair's
##                 correlation is constant in the covariates (true slope 0). Re-running the analysis on such
##                 data re-estimates the nuisances, so the null z-statistics carry the same shared nuisance
##                 errors as the real ones.
##   perm_pvalues(): two-sided p-values of observed z-statistics against the pooled permutation null.

#' Efron's empirical null by central matching (robust version)
#' @return list with delta, sigma and two-sided p-values under N(delta, sigma^2).
empnull <- function(z) {
  delta <- stats::median(z); sigma <- stats::IQR(z) / (2 * stats::qnorm(0.75))
  list(delta = delta, sigma = sigma, p = 2 * stats::pnorm(-abs(z - delta) / sigma))
}

#' Residual-permutation null z-statistics for pvcor_many.
#' @param Y N x p feature matrix; data, pairs, nuis_formula, w_formula, K, strategy as in pvcor_many.
#' @param B number of permuted data sets.
#' @param coef name of the coefficient of W whose z-statistics are returned.
#' @param comparators also return null z-statistics of pooled LIONESS regression.
#' @param seed optional seed (set before the permutations).
#' @return list with `orth` (pairs x B), `pooled` (pairs x B or NULL) and the fitted means and scales.
perm_null <- function(Y, data, pairs, nuis_formula, w_formula = ~ x, B = 5, K = 5, coef = "x",
                      comparators = FALSE, strategy = "auto", seed = NULL) {
  Y <- as.matrix(Y); N <- nrow(Y); rhs <- .rhs(nuis_formula); all_rows <- seq_len(N)
  Mu <- sapply(seq_len(ncol(Y)), function(j) .fit_mean(Y[, j], data, all_rows, "gam", rhs)$pred)
  S <- sapply(seq_len(ncol(Y)), function(j) .fit_var(Y[, j] - Mu[, j], data, all_rows, "gam", rhs))
  Z <- (Y - Mu) / S
  if (!is.null(seed)) set.seed(seed)
  runs <- lapply(seq_len(B), function(b) {
    Ys <- Mu + S * Z[sample(N), , drop = FALSE]
    tab <- pvcor_many(Ys, data, pairs = pairs, w_formula = w_formula, nuis_formula = nuis_formula, K = K,
                      comparators = comparators, strategy = strategy)
    list(orth = tab[[paste0("b_", coef)]] / tab[[paste0("se_", coef)]],
         pooled = if (comparators) tab[[paste0("pooled_b_", coef)]] / tab[[paste0("pooled_se_", coef)]] else NULL)
  })
  list(orth = sapply(runs, `[[`, "orth"),
       pooled = if (comparators) sapply(runs, `[[`, "pooled") else NULL)
}

#' Two-sided permutation p-values: (1 + #{|z_null| >= |z|}) / (1 + M), with the null pooled over pairs
#' and permutations.
perm_pvalues <- function(z, znull) {
  a <- sort(abs(as.vector(znull))); M <- length(a)
  (1 + M - findInterval(abs(z), a, left.open = TRUE)) / (1 + M)
}
