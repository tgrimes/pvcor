## Orthogonal pseudo-outcome regression for all pairs of many features (Algorithm 1 at scale).
##
## Feature-level nuisances are fitted once per fold. With the "stdcommon" strategy every feature
## uses the same smoothing parameter on the standardised scale, so the mean smoother is one linear
## operator shared by all features; the initial correlations are likewise obtained by applying one
## fixed linear smoother to the products of standardised residuals. The per-pair work is then
## matrix arithmetic: O(p) model fits plus O(p^2 N) operations.

## Linear operator of a gaussian GAM with fixed smoothing parameters: fitted coefficients are
## M %*% y_train, predictions Xp %*% M %*% y_train.
.gam_operator <- function(template_y, data, train, rhs, sp) {
  d <- data; d$.y <- template_y
  fit <- mgcv::gam(as.formula(paste(".y ~", rhs)), data = d[train, , drop = FALSE], sp = sp)
  X <- predict(fit, d[train, , drop = FALSE], type = "lpmatrix")
  M <- (fit$Vp / fit$sig2) %*% t(X)
  list(M = M, Xp = predict(fit, d, type = "lpmatrix"))
}

.geo_sp <- function(Ymat, data, train, rhs, max_fit = 40) {
  cols <- if (ncol(Ymat) > max_fit) sample(ncol(Ymat), max_fit) else seq_len(ncol(Ymat))
  sps <- lapply(cols, function(j) .fit_mean(Ymat[, j], data, train, "gam", rhs)$sp)
  if (length(sps[[1]]) == 0) return(numeric(0))
  sps <- matrix(unlist(sps), ncol = length(cols))
  exp(rowMeans(log(pmax(sps, 1e-10))))
}

#' @param Y N x p matrix of features.
#' @param data data.frame of covariates (nuisance formula and W).
#' @param pairs two-column matrix of column indices (default: all pairs).
#' @param return_phi attach the N x P matrix of pseudo-outcomes as attribute `phi` (use for modest P).
#' @param comparators logical: also compute pooled LIONESS (jackknife) and liquid association.
#' @param la_var name of the covariate for liquid association.
#' @return data.frame with one row per pair: orthogonal slope/se for each W coefficient (columns
#'   `b_<name>`, `se_<name>`), and the comparators.
pvcor_many <- function(Y, data, pairs = NULL, w_formula = ~ x, nuis_formula = ~ s(x, k = 10), K = 5,
                       folds = NULL, se = "HC3", rho_clip = 0.99, var_floor = 1e-3, chunk = 20000,
                       comparators = TRUE, la_var = "x", strata = NULL, max_sp_fit = 40, verbose = FALSE,
                       strategy = c("auto", "stdcommon", "separate"), return_phi = FALSE) {
  strategy <- match.arg(strategy); if (strategy == "auto") strategy <- resolve_strategy(nuis_formula)
  Y <- as.matrix(Y); N <- nrow(Y); p <- ncol(Y)
  if (is.null(pairs)) pairs <- which(upper.tri(diag(p)), arr.ind = TRUE)[, c("row", "col"), drop = FALSE]
  pairs <- as.matrix(pairs); P <- nrow(pairs)
  if (is.null(folds)) folds <- make_folds(N, K, strata)
  ks <- sort(unique(folds)); rhs <- .rhs(nuis_formula)
  W <- model.matrix(w_formula, data); Bw <- solve(crossprod(W), t(W)); h <- rowSums(W * t(Bw))
  ## ---- feature-level nuisances, per fold --------------------------------------------------------
  Zf <- vector("list", length(ks)); Rop <- vector("list", length(ks))
  for (a in seq_along(ks)) {
    te <- which(folds == ks[a]); tr <- if (length(ks) == 1) te else which(folds != ks[a])
    if (verbose) message("fold ", a, ": feature-level fits")
    Mu <- sapply(seq_len(p), function(j) .fit_mean(Y[, j], data, tr, "gam", rhs)$pred)
    S <- sapply(seq_len(p), function(j) .fit_var(Y[, j] - Mu[, j], data, tr, "gam", rhs, var_floor))
    if (strategy == "stdcommon") {
      Tm <- Y / S
      sp <- .geo_sp(Tm, data, tr, rhs, max_sp_fit)
      op <- .gam_operator(Tm[, 1], data, tr, rhs, if (length(sp)) sp else NULL)
      Mu <- S * (op$Xp %*% (op$M %*% Tm[tr, , drop = FALSE]))
      S <- sapply(seq_len(p), function(j) .fit_var(Y[, j] - Mu[, j], data, tr, "gam", rhs, var_floor))
    }
    Z <- (Y - Mu) / S
    ## common smoother for the products
    smp <- pairs[sample(P, min(P, max_sp_fit)), , drop = FALSE]
    sp_r <- .geo_sp(Z[, smp[, 1], drop = FALSE] * Z[, smp[, 2], drop = FALSE], data, tr, rhs, max_sp_fit)
    Rop[[a]] <- .gam_operator(Z[, 1] * Z[, 2], data, tr, rhs, if (length(sp_r)) sp_r else NULL)
    Rop[[a]]$tr <- tr; Rop[[a]]$te <- te
    Zf[[a]] <- Z
  }
  ## ---- pairs, in chunks ---------------------------------------------------------------------------
  wn <- colnames(W); out <- vector("list", ceiling(P / chunk)); phis <- vector("list", ceiling(P / chunk))
  if (comparators) {
    nsc <- function(v) qnorm((rank(v) - 0.5) / length(v))
    NS <- apply(Y, 2, nsc); nsx <- nsc(data[[la_var]])
  }
  for (cidx in seq_along(out)) {
    rows <- ((cidx - 1) * chunk + 1):min(P, cidx * chunk); pr <- pairs[rows, , drop = FALSE]
    Phi <- matrix(0, N, length(rows))
    for (a in seq_along(ks)) {
      Z <- Zf[[a]]; tr <- Rop[[a]]$tr; te <- Rop[[a]]$te
      Ptr <- Z[tr, pr[, 1], drop = FALSE] * Z[tr, pr[, 2], drop = FALSE]
      rho <- Rop[[a]]$Xp[te, , drop = FALSE] %*% (Rop[[a]]$M %*% Ptr)
      rho <- pmin(pmax(rho, -rho_clip), rho_clip)
      z1 <- Z[te, pr[, 1], drop = FALSE]; z2 <- Z[te, pr[, 2], drop = FALSE]
      Phi[te, ] <- rho + z1 * z2 - rho / 2 * (z1^2 + z2^2)
    }
    if (return_phi) phis[[cidx]] <- Phi
    res <- .ols_many(Phi, W, Bw, h, se)
    df <- data.frame(i = pr[, 1], j = pr[, 2])
    for (l in seq_along(wn)) { df[[paste0("b_", wn[l])]] <- res$coef[l, ]; df[[paste0("se_", wn[l])]] <- res$se[l, ] }
    if (comparators) {
      PV <- pv_pearson_many(Y, pr); rp <- .ols_many(PV, W, Bw, h, se)
      for (l in seq_along(wn)) { df[[paste0("pooled_b_", wn[l])]] <- rp$coef[l, ]; df[[paste0("pooled_se_", wn[l])]] <- rp$se[l, ] }
      LP <- NS[, pr[, 1], drop = FALSE] * NS[, pr[, 2], drop = FALSE] * nsx
      df$la <- colMeans(LP); df$la_se <- sqrt(colSums(sweep(LP, 2, df$la)^2) / (N - 1) / N)
    }
    out[[cidx]] <- df
    if (verbose) message("pairs ", max(rows), " / ", P)
  }
  res <- do.call(rbind, out); attr(res, "strategy") <- strategy
  if (return_phi) attr(res, "phi") <- do.call(cbind, phis)
  res
}

.ols_many <- function(Phi, W, Bw, h, se = "HC3") {
  coef <- Bw %*% Phi; E <- Phi - W %*% coef
  if (se == "const") {                                   # classical OLS standard errors
    s2 <- colSums(E^2) / (nrow(W) - ncol(W))
    return(list(coef = coef, se = sqrt(outer(diag(solve(crossprod(W))), s2))))
  }
  if (se == "HC3") E <- E / (1 - h)
  V <- (Bw^2) %*% (E^2)
  list(coef = coef, se = sqrt(V))
}
