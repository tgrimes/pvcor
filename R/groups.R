#' Orthogonal pseudo-outcomes when the covariate is a group label (all pairs of features).
#'
#' Nuisances are the group means, standard deviations (divisor n) and correlations, estimated on
#' the training folds (stratified by group). With K = 1 the nuisances are in sample, and the group
#' averages of the pseudo-outcomes are then exactly the group sample correlations (Corollary 3 of
#' the paper): this is the group-standardised LIONESS arm to first order.
#' @param Y N x p matrix. @param group factor (or vector) of group labels.
#' @param se "const" (classical OLS), "HC0" or "HC3".
#' @return list(table = per-pair coefficients and SEs, Phi = pseudo-outcomes if requested, W).
pvcor_groups <- function(Y, group, pairs = NULL, K = 5, folds = NULL, se = "HC3", rho_clip = 0.99,
                         return_phi = FALSE) {
  Y <- as.matrix(Y); N <- nrow(Y); p <- ncol(Y); group <- factor(group)
  if (is.null(pairs)) pairs <- which(upper.tri(diag(p)), arr.ind = TRUE)[, c("row", "col"), drop = FALSE]
  pairs <- as.matrix(pairs)
  if (is.null(folds)) folds <- make_folds(N, K, strata = group)
  ks <- sort(unique(folds)); Phi <- matrix(0, N, nrow(pairs))
  for (k in ks) {
    te <- which(folds == k); tr <- if (length(ks) == 1) te else which(folds != k)
    Z <- matrix(0, N, p); R <- list()
    for (g in levels(group)) {
      trg <- tr[group[tr] == g]; all_g <- which(group == g)
      m <- colMeans(Y[trg, , drop = FALSE]); s <- sqrt(colMeans(sweep(Y[trg, , drop = FALSE], 2, m)^2))
      Z[all_g, ] <- sweep(sweep(Y[all_g, , drop = FALSE], 2, m), 2, s, "/")
      Rg <- crossprod(Z[trg, , drop = FALSE]) / length(trg)
      R[[g]] <- pmin(pmax(Rg[pairs], -rho_clip), rho_clip)
    }
    for (g in levels(group)) {
      teg <- te[group[te] == g]; if (!length(teg)) next
      rho <- matrix(R[[g]], length(teg), nrow(pairs), byrow = TRUE)
      z1 <- Z[teg, pairs[, 1], drop = FALSE]; z2 <- Z[teg, pairs[, 2], drop = FALSE]
      Phi[teg, ] <- rho + z1 * z2 - rho / 2 * (z1^2 + z2^2)
    }
  }
  W <- model.matrix(~ group); Bw <- solve(crossprod(W), t(W)); h <- rowSums(W * t(Bw))
  res <- .ols_many(Phi, W, Bw, h, se)
  tab <- data.frame(i = pairs[, 1], j = pairs[, 2])
  for (l in seq_len(ncol(W))) { tab[[paste0("b_", colnames(W)[l])]] <- res$coef[l, ]; tab[[paste0("se_", colnames(W)[l])]] <- res$se[l, ] }
  list(table = tab, Phi = if (return_phi) Phi else NULL, W = W, folds = folds)
}
