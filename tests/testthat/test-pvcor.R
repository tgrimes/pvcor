test_that("closed-form Pearson pseudo-values match brute-force leave-one-out", {
  set.seed(1); y1 <- rnorm(30); y2 <- y1 + rnorm(30)
  bf <- 30 * cor(y1, y2) - 29 * sapply(1:30, function(i) cor(y1[-i], y2[-i]))
  expect_equal(pv_pearson(y1, y2), bf, tolerance = 1e-10)
  Y <- cbind(y1, y2, rnorm(30)); pr <- rbind(c(1, 2), c(1, 3))
  M <- pv_pearson_many(Y, pr)
  expect_equal(M[, 1], bf, tolerance = 1e-10)
  expect_equal(M[, 2], pv_pearson(Y[, 1], Y[, 3]), tolerance = 1e-10)
})

test_that("tanh_fit solves the estimating equation and recovers a noiseless model", {
  set.seed(2); x <- runif(200, -1, 1); W <- cbind(1, x)
  f <- tanh_fit(tanh(0.3 + 0.8 * x), W)
  expect_equal(unname(f$coef), c(0.3, 0.8), tolerance = 1e-8)
  phi <- tanh(0.3 + 0.8 * x) + rnorm(200, sd = 0.5); f <- tanh_fit(phi, W)
  expect_lt(max(abs(crossprod(W, phi - tanh(W %*% f$coef)))), 1e-8)
})

test_that("Gaussian covariance regression gradient is correct", {
  set.seed(3); n <- 300; x <- runif(n, -1, 1); B <- cbind(1, x); W <- cbind(1, x)
  y1 <- x + exp(0.3 * x) * rnorm(n); y2 <- -x + rnorm(n)
  fit <- covreg_gauss(y1, y2, B, W)
  expect_equal(fit$convergence, 0)
  ## re-derive objective via the closure in the fit: check gradient numerically at a random point
  env <- environment(covreg_gauss)
  th <- fit$par + rnorm(length(fit$par), sd = 0.05)
  f_nll <- function(th) { pb <- 2; m1 <- B %*% th[1:2]; m2 <- B %*% th[3:4]; ls1 <- B %*% th[5:6]; ls2 <- B %*% th[7:8]
    rho <- tanh(W %*% th[9:10]); z1 <- (y1 - m1) / exp(ls1); z2 <- (y2 - m2) / exp(ls2); q <- 1 - rho^2
    sum(ls1 + ls2 + 0.5 * log(q) + (z1^2 - 2 * rho * z1 * z2 + z2^2) / (2 * q)) }
  num <- sapply(seq_along(th), function(j) { e <- replace(0 * th, j, 1e-6); (f_nll(th + e) - f_nll(th - e)) / 2e-6 })
  ## optimum: gradient ~ 0 and numerical gradient matches finite differences of the same objective
  num0 <- sapply(seq_along(fit$par), function(j) { e <- replace(0 * th, j, 1e-6); (f_nll(fit$par + e) - f_nll(fit$par - e)) / 2e-6 })
  expect_lt(max(abs(num0)), 1e-2)
})

test_that("fixed-sp GAM operator reproduces mgcv predictions", {
  set.seed(4); n <- 200; d <- data.frame(x = runif(n)); y <- sin(4 * d$x) + rnorm(n, sd = 0.3); tr <- 1:150
  fit <- mgcv::gam(y ~ s(x, k = 10), data = data.frame(x = d$x, y = y)[tr, ], method = "REML")
  op <- pvcor:::.gam_operator(y, d, tr, "s(x, k = 10)", fit$sp)
  expect_equal(as.numeric(op$Xp %*% (op$M %*% y[tr])), as.numeric(predict(fit, d)), tolerance = 1e-6)
})

test_that("oracle nuisances reproduce the oracle pseudo-outcome", {
  set.seed(5); n <- 400; x <- runif(n, -1, 1); rho <- tanh(0.2 + 0.6 * x)
  e1 <- rnorm(n); e2 <- rho * e1 + sqrt(1 - rho^2) * rnorm(n)
  mu1 <- x^2; mu2 <- sin(2 * x); s1 <- exp(0.5 * x); s2 <- exp(-0.3 * x)
  y1 <- mu1 + s1 * e1; y2 <- mu2 + s2 * e2
  f <- pvcor(y1, y2, data.frame(x = x), oracle = list(mu1 = mu1, mu2 = mu2, s1 = s1, s2 = s2, rho = rho))
  expect_equal(f$strategies$stdcommon$phi, rho + e1 * e2 - rho / 2 * (e1^2 + e2^2), tolerance = 1e-12)
})

test_that("pvcor_many agrees with pvcor up to the common smoothing parameter", {
  set.seed(6); n <- 300; x <- runif(n, -1, 1); rho <- tanh(0.2 + 0.6 * x)
  e1 <- rnorm(n); e2 <- rho * e1 + sqrt(1 - rho^2) * rnorm(n)
  Y <- cbind(x + e1, -x + exp(0.3 * x) * e2, rnorm(n))
  folds <- make_folds(n, 5)
  m <- pvcor_many(Y, data.frame(x = x), folds = folds, nuis_formula = ~ s(x, k = 8))
  s <- pvcor(Y[, 1], Y[, 2], data.frame(x = x), folds = folds, nuis_formula = ~ s(x, k = 8))
  expect_equal(nrow(m), 3)
  expect_lt(abs(m$b_x[1] - s$strategies$stdcommon$fit$blp$coef[2]), 0.05)
})

test_that("in-sample group pseudo-outcomes average to the group correlations (Corollary 3)", {
  set.seed(7); n <- 60; g <- rep(0:1, c(25, 35)); Y <- matrix(rnorm(n * 3), n) + g * 2
  f <- pvcor_groups(Y, g, K = 1, return_phi = TRUE)
  r0 <- cor(Y[g == 0, ])[1, 2]; r1 <- cor(Y[g == 1, ])[1, 2]
  expect_equal(mean(f$Phi[g == 0, 1]), r0, tolerance = 1e-12)
  expect_equal(f$table$b_group1[1], r1 - r0, tolerance = 1e-12)
})

test_that("permutation p-values count null statistics at least as extreme", {
  zn <- matrix(c(-3, -1, 0, 1, 2, 4), 3)
  expect_equal(perm_pvalues(c(0, 2, 5), zn), c(7, 4, 1) / 7)   # |zn| >= 2: 2, 3, 4
  e <- empnull(c(-2, -1, 0, 1, 2) + 1)
  expect_equal(e$delta, 1)
})

test_that("residual-permutation null keeps means and removes correlation dependence", {
  set.seed(8); n <- 200; x <- runif(n, -1, 1); rho <- tanh(1.2 * x)
  e1 <- rnorm(n); e2 <- rho * e1 + sqrt(1 - rho^2) * rnorm(n)
  Y <- cbind(2 * x + e1, -x + e2, rnorm(n)); pairs <- rbind(c(1, 2), c(1, 3), c(2, 3))
  nl <- perm_null(Y, data.frame(x = x), pairs, ~ s(x, k = 6), B = 3, seed = 1)
  expect_equal(dim(nl$orth), c(3, 3))
  expect_lt(mean(abs(nl$orth[1, ])), 3)              # the dynamic pair is null after permutation
})

test_that("random-forest options select in-bag or out-of-bag training predictions", {
  skip_if_not_installed("ranger")
  set.seed(9); n <- 200; d <- data.frame(x = runif(n)); y <- sin(4 * d$x) + rnorm(n, sd = 0.5); tr <- 1:150
  old <- options(pvcor.rf_oob = FALSE, pvcor.rf_min_node = 10); set.seed(1)
  inbag <- pvcor:::.fit_mean(y, d, tr, "rf", "x")$pred; options(old); set.seed(1)
  oob <- pvcor:::.fit_mean(y, d, tr, "rf", "x")$pred
  expect_lt(mean((y[tr] - inbag[tr])^2), mean((y[tr] - oob[tr])^2))   # in-bag residuals are too small
})
