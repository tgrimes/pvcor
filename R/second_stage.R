## Second-stage regressions of (pseudo-)outcomes on a working design W, with sandwich inference.

.hc_vcov <- function(W, e, h = NULL, type = c("HC3", "HC0")) {
  type <- match.arg(type)
  B <- solve(crossprod(W))
  if (type == "HC3") {
    if (is.null(h)) h <- rowSums((W %*% B) * W)
    e <- e / (1 - h)
  }
  B %*% crossprod(W * e) %*% B
}

#' Best-linear-projection (OLS) fit with HC0/HC3 standard errors.
#' @param phi outcome vector; @param W design matrix (include the intercept column).
blp_fit <- function(phi, W, se = c("HC3", "HC0")) {
  se <- match.arg(se)
  W <- as.matrix(W)
  B <- solve(crossprod(W)); coef <- drop(B %*% crossprod(W, phi))
  e <- drop(phi - W %*% coef)
  V <- .hc_vcov(W, e, type = se)
  names(coef) <- colnames(W); dimnames(V) <- list(colnames(W), colnames(W))
  list(coef = coef, vcov = V, se = sqrt(diag(V)), resid = e, type = "blp")
}

#' tanh-link estimating equation sum_i W_i {phi_i - tanh(W_i' beta)} = 0, solved by damped
#' Newton steps on the convex criterion sum_i {log cosh(W_i' beta) - phi_i W_i' beta}.
#' Sandwich variance A^{-1} B A^{-1}; "HC3" rescales residuals by the leverages of the
#' linearised fit.
tanh_fit <- function(phi, W, se = c("HC3", "HC0"), maxit = 100, tol = 1e-10) {
  se <- match.arg(se)
  W <- as.matrix(W)
  crit <- function(b) { u <- drop(W %*% b); sum(log(cosh(u)) - phi * u) }
  b <- rep(0, ncol(W))
  for (it in seq_len(maxit)) {
    u <- drop(W %*% b); g <- crossprod(W, phi - tanh(u)); Hm <- crossprod(W * (1 / cosh(u))^2, W)
    step <- solve(Hm, g); f0 <- crit(b); a <- 1
    while (crit(b + a * step) > f0 + 1e-12 && a > 1e-8) a <- a / 2
    b <- drop(b + a * step)
    if (max(abs(a * step)) < tol) break
  }
  u <- drop(W %*% b); d <- (1 / cosh(u))^2; e <- phi - tanh(u)
  A <- crossprod(W * d, W); Ai <- solve(A)
  if (se == "HC3") { h <- d * rowSums((W %*% Ai) * W); e2 <- e / (1 - h) } else e2 <- e
  V <- Ai %*% crossprod(W * e2) %*% Ai
  names(b) <- colnames(W); dimnames(V) <- list(colnames(W), colnames(W))
  list(coef = b, vcov = V, se = sqrt(diag(V)), resid = e, type = "tanh", iterations = it)
}

#' Wild-bootstrap test of H0: beta_j = b0 in an OLS fit (restricted residuals, Rademacher
#' weights, HC3-studentised statistic). Returns the observed t and the bootstrap p-value.
wild_boot_test <- function(phi, W, j = 2, b0 = 0, B = 999) {
  W <- as.matrix(W)
  fit <- blp_fit(phi, W, "HC3"); t_obs <- (fit$coef[j] - b0) / fit$se[j]
  y0 <- phi - b0 * W[, j]; W0 <- W[, -j, drop = FALSE]
  b_r <- solve(crossprod(W0), crossprod(W0, y0)); fit_r <- drop(W0 %*% b_r); e_r <- y0 - fit_r
  BWt <- solve(crossprod(W)); h <- rowSums((W %*% BWt) * W)
  tb <- replicate(B, {
    ys <- fit_r + b0 * W[, j] + e_r * sample(c(-1, 1), length(e_r), replace = TRUE)
    cb <- drop(BWt %*% crossprod(W, ys)); eb <- drop(ys - W %*% cb)
    V <- BWt %*% crossprod(W * (eb / (1 - h))) %*% BWt
    (cb[j] - b0) / sqrt(V[j, j])
  })
  list(t = unname(t_obs), p = mean(abs(tb) >= abs(t_obs)))
}
