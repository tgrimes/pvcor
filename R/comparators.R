## Comparators: jackknife (LIONESS-type) Pearson pseudo-values, liquid association, and a
## Gaussian covariance-regression model fitted by maximum likelihood.

#' Exact jackknife pseudo-values of the Pearson correlation for one pair (pooled; LIONESS edges).
#' Uses the closed form r_(-s) = (r - z1 z2 / m) / sqrt((1 - z1^2 / m)(1 - z2^2 / m)), m = N - 1,
#' with pooled z-scores computed with divisor N.
pv_pearson <- function(y1, y2) {
  N <- length(y1); m <- N - 1
  z1 <- (y1 - mean(y1)) / sqrt(mean((y1 - mean(y1))^2)); z2 <- (y2 - mean(y2)) / sqrt(mean((y2 - mean(y2))^2))
  r <- mean(z1 * z2)
  N * r - m * (r - z1 * z2 / m) / sqrt((1 - z1^2 / m) * (1 - z2^2 / m))
}

#' The same for all pairs of columns of Z (z-scored internally); returns an N x P matrix whose
#' columns follow `pairs` (a two-column index matrix).
pv_pearson_many <- function(Y, pairs) {
  N <- nrow(Y); m <- N - 1
  Z <- scale(Y) * sqrt(N / (N - 1))                      # divisor-N z-scores
  R <- crossprod(Z) / N
  a <- Z[, pairs[, 1], drop = FALSE]; b <- Z[, pairs[, 2], drop = FALSE]
  r <- R[pairs]
  rr <- matrix(r, N, length(r), byrow = TRUE)
  N * rr - m * (rr - a * b / m) / sqrt((1 - a^2 / m) * (1 - b^2 / m))
}

#' Liquid association of Li (2002): mean of the product of normal scores, with a z-test.
liquid_association <- function(y1, y2, x, normal_scores = TRUE) {
  ns <- function(v) if (normal_scores) qnorm((rank(v) - 0.5) / length(v)) else (v - mean(v)) / sd(v)
  p <- ns(y1) * ns(y2) * ns(x)
  c(LA = mean(p), se = sd(p) / sqrt(length(p)))
}

#' Gaussian covariance regression with natural-spline means and log-standard deviations and a
#' tanh-linear correlation, rho(x) = tanh(W beta). Fitted by maximum likelihood; standard errors
#' are model based (inverse observed information), as a likelihood user would report them.
#' @param B basis matrix for the means and log-sds (include an intercept column).
#' @param W design for atanh(rho).
covreg_gauss <- function(y1, y2, B, W) {
  B <- as.matrix(B); W <- as.matrix(W); pb <- ncol(B); pw <- ncol(W)
  idx <- list(a1 = 1:pb, a2 = pb + 1:pb, b1 = 2 * pb + 1:pb, b2 = 3 * pb + 1:pb, c = 4 * pb + 1:pw)
  nll <- function(th) {
    m1 <- B %*% th[idx$a1]; m2 <- B %*% th[idx$a2]; ls1 <- B %*% th[idx$b1]; ls2 <- B %*% th[idx$b2]
    rho <- tanh(W %*% th[idx$c]); z1 <- (y1 - m1) / exp(ls1); z2 <- (y2 - m2) / exp(ls2); q <- 1 - rho^2
    sum(ls1 + ls2 + 0.5 * log(q) + (z1^2 - 2 * rho * z1 * z2 + z2^2) / (2 * q))
  }
  gr <- function(th) {
    m1 <- drop(B %*% th[idx$a1]); m2 <- drop(B %*% th[idx$a2]); s1 <- exp(drop(B %*% th[idx$b1])); s2 <- exp(drop(B %*% th[idx$b2]))
    eta <- drop(W %*% th[idx$c]); rho <- tanh(eta); q <- 1 - rho^2
    z1 <- (y1 - m1) / s1; z2 <- (y2 - m2) / s2
    dz1 <- (z1 - rho * z2) / q; dz2 <- (z2 - rho * z1) / q            # d(quad/2)/dz
    g_a1 <- -crossprod(B, dz1 / s1); g_a2 <- -crossprod(B, dz2 / s2)
    g_b1 <- crossprod(B, 1 - dz1 * z1); g_b2 <- crossprod(B, 1 - dz2 * z2)
    ## d/drho of 0.5 log q + (z1^2 - 2 rho z1 z2 + z2^2) / (2q)
    Qd <- z1^2 - 2 * rho * z1 * z2 + z2^2
    drho <- -rho / q + (-z1 * z2) / q + Qd * rho / q^2
    g_c <- crossprod(W, drho * q)                                     # d rho / d eta = 1 - rho^2
    c(g_a1, g_a2, g_b1, g_b2, g_c)
  }
  st <- c(qr.solve(B, y1), qr.solve(B, y2))
  r1 <- y1 - B %*% st[1:pb]; r2 <- y2 - B %*% st[pb + 1:pb]
  st <- c(st, qr.solve(B, log(abs(r1) + 1e-8) + 0.635), qr.solve(B, log(abs(r2) + 1e-8) + 0.635),
          atanh(max(min(cor(r1, r2), 0.95), -0.95)), rep(0, pw - 1))
  op <- optim(st, nll, gr, method = "BFGS", control = list(maxit = 1000, reltol = 1e-12))
  H <- optimHess(op$par, nll, gr)
  V <- tryCatch(solve(H), error = function(e) matrix(NA, length(op$par), length(op$par)))
  cc <- op$par[idx$c]; Vc <- V[idx$c, idx$c, drop = FALSE]
  names(cc) <- colnames(W)
  list(coef = cc, se = sqrt(pmax(diag(Vc), 0)), vcov = Vc, convergence = op$convergence, par = op$par)
}
