## Orthogonal pseudo-outcomes for the conditional PARTIAL correlation (extension of Algorithm 1 from a
## pair to a covariance functional g(Sigma(x)) with p fixed and small).
##
## Per fold: feature-level means and scales (as for a pair), standardised residuals z, and an initial
## correlation MATRIX R(x) obtained by smoothing the products z_j z_k pair by pair. For each held-out row
## the matrix is clipped to a valid correlation matrix (minimum eigenvalue floor), Omega = R^{-1},
## T = Omega z / sqrt(diag Omega), and for every requested pair (j, k) the pseudo-outcome is
##   phi = rho + T_j T_k + (rho / 2) (T_j^2 + T_k^2),   rho = -Omega_jk / sqrt(Omega_jj Omega_kk),
## which is g(Sigma_hat) + D g_{Sigma_hat}[(y - mu_hat)(y - mu_hat)' - Sigma_hat] for g = partial
## correlation (scale invariance lets us work with the standardised residuals).
.make_pd_corr <- function(R, eig_floor = 0.05) {
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) >= eig_floor) return(R)
  lam <- (eig_floor - min(ev)) / (1 - min(ev))   # shrink toward the identity
  (1 - lam) * R + lam * diag(nrow(R))
}

#' @param Y N x p matrix. @param data data.frame of covariates. @param pairs two-column index matrix.
#' @param return_phi return the N x P matrix of pseudo-outcomes as well.
pvcor_pcor <- function(Y, data, w_formula = ~ x, nuis_formula = ~ s(x, k = 10), K = 5, folds = NULL,
                       se = "HC3", rho_clip = 0.95, var_floor = 1e-3, eig_floor = 0.05, strata = NULL,
                       pairs = NULL, learner = "gam", return_phi = FALSE) {
  Y <- as.matrix(Y); N <- nrow(Y); p <- ncol(Y)
  if (p < 3) stop("partial correlation needs at least three features")
  if (is.null(pairs)) pairs <- which(upper.tri(diag(p)), arr.ind = TRUE)[, c("row", "col"), drop = FALSE]
  pairs <- as.matrix(pairs); P <- nrow(pairs)
  if (is.null(folds)) folds <- make_folds(N, K, strata)
  ks <- sort(unique(folds)); rhs <- .rhs(nuis_formula)
  allp <- which(upper.tri(diag(p)), arr.ind = TRUE)[, c("row", "col"), drop = FALSE]
  Phi <- matrix(NA_real_, N, P)
  for (a in seq_along(ks)) {
    te <- which(folds == ks[a]); tr <- if (length(ks) == 1) te else which(folds != ks[a])
    Mu <- sapply(seq_len(p), function(j) .fit_mean(Y[, j], data, tr, learner, rhs)$pred)
    S <- sapply(seq_len(p), function(j) .fit_var(Y[, j] - Mu[, j], data, tr, learner, rhs, var_floor))
    Z <- (Y - Mu) / S
    ## initial correlations: smooth the products pair by pair, predict at the held-out rows
    Rh <- array(0, c(length(te), p, p))
    for (m in seq_len(nrow(allp))) {
      j <- allp[m, 1]; k <- allp[m, 2]
      r <- .fit_mean(Z[, j] * Z[, k], data, tr, learner, rhs)$pred[te]
      r <- pmax(pmin(r, rho_clip), -rho_clip)
      Rh[, j, k] <- r; Rh[, k, j] <- r
    }
    for (u in seq_along(te)) {
      R <- Rh[u, , ]; diag(R) <- 1
      R <- .make_pd_corr(R, eig_floor)
      Om <- solve(R); Tt <- drop(Om %*% Z[te[u], ]) / sqrt(diag(Om))
      rho <- -Om[pairs[, 1] + (pairs[, 2] - 1) * p] / sqrt(diag(Om)[pairs[, 1]] * diag(Om)[pairs[, 2]])
      Phi[te[u], ] <- rho + Tt[pairs[, 1]] * Tt[pairs[, 2]] + (rho / 2) * (Tt[pairs[, 1]]^2 + Tt[pairs[, 2]]^2)
    }
  }
  W <- model.matrix(w_formula, data)
  fits <- lapply(seq_len(P), function(m) blp_fit(Phi[, m], W, se))
  out <- data.frame(i = pairs[, 1], j = pairs[, 2])
  for (nm in colnames(W)) {
    out[[paste0("b_", nm)]] <- sapply(fits, function(f) f$coef[[nm]])
    out[[paste0("se_", nm)]] <- sapply(fits, function(f) f$se[[nm]])
  }
  if (return_phi) attr(out, "phi") <- Phi
  out
}

#' Exact jackknife (LIONESS) pseudo-values of the partial correlation for given pairs, by direct
#' leave-one-out recomputation (p small).
pv_pcor <- function(Y, pairs) {
  Y <- as.matrix(Y); N <- nrow(Y); p <- ncol(Y); pairs <- as.matrix(pairs)
  pc <- function(Yx) { O <- solve(cov(Yx)); -O[pairs[, 1] + (pairs[, 2] - 1) * p] / sqrt(diag(O)[pairs[, 1]] * diag(O)[pairs[, 2]]) }
  th <- pc(Y)
  out <- matrix(NA_real_, N, nrow(pairs))
  for (i in seq_len(N)) out[i, ] <- N * th - (N - 1) * pc(Y[-i, , drop = FALSE])
  out
}
