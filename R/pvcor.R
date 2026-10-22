#' Default nuisance strategy: shared standardised means ("stdcommon") when the nuisance formula
#' involves a single covariate, separate smoothers otherwise.
resolve_strategy <- function(nuis_formula) if (length(all.vars(nuis_formula)) == 1) "stdcommon" else "separate"

#' Cross-fitted orthogonal pseudo-outcome regression for one feature pair (Algorithm 1).
#'
#' @param y1,y2 numeric vectors (the two features).
#' @param data data.frame with the covariates used in `nuis_formula` and `w_formula`.
#' @param w_formula one-sided formula for the working design W of the second stage.
#' @param nuis_formula one-sided formula for the nuisance regressions.
#' @param K number of cross-fitting folds; `K = 1` fits the nuisances in sample (no cross-fitting).
#' @param folds optional fold labels (overrides K).
#' @param strategies nuisance strategies to compute ("stdcommon", "separate"); the first is the
#'   primary one.
#' @param learner "gam", "lm" or "rf".
#' @param link "identity" (best linear projection) and/or "tanh".
#' @param se "HC3" or "HC0".
#' @param oracle optional list of true nuisances (see fit_nuisances).
#' @return list with, for each strategy, the pseudo-outcomes `phi`, the standardised residuals,
#'   the initial correlation, and second-stage fits.
pvcor <- function(y1, y2, data, w_formula = ~ x, nuis_formula = ~ s(x, k = 10), K = 5, folds = NULL,
                  strategies = "stdcommon", learner = "gam", link = "identity", se = "HC3",
                  oracle = NULL, rho_clip = 0.99, strata = NULL, rho_formula = NULL) {
  n <- length(y1)
  if (is.null(folds)) folds <- make_folds(n, K, strata)
  Kf <- length(unique(folds))
  res <- setNames(lapply(strategies, function(s) list(phi = numeric(n), z1 = numeric(n), z2 = numeric(n),
                                                      rho = numeric(n), mu1 = numeric(n), mu2 = numeric(n),
                                                      s1 = numeric(n), s2 = numeric(n))), strategies)
  for (k in sort(unique(folds))) {
    te <- which(folds == k); tr <- if (Kf == 1) te else which(folds != k)
    or_k <- oracle
    nu <- fit_nuisances(y1, y2, data, tr, nuis_formula, learner, strategies, oracle, rho_clip, rho_formula = rho_formula)
    for (s in strategies) for (nm in c("z1", "z2", "rho", "mu1", "mu2", "s1", "s2")) res[[s]][[nm]][te] <- nu[[s]][[nm]][te]
  }
  W <- model.matrix(w_formula, data)
  for (s in strategies) {
    r <- res[[s]]
    r$phi <- phi_from_nuisances(r$z1, r$z2, r$rho)
    r$fit <- list()
    if ("identity" %in% link) r$fit$blp <- blp_fit(r$phi, W, se)
    if ("tanh" %in% link) r$fit$tanh <- tanh_fit(r$phi, W, se)
    res[[s]] <- r
  }
  structure(list(strategies = res, W = W, folds = folds, call = match.call()), class = "pvcor")
}

#' One-step "iterated" update: replace the initial correlation by the second-stage fit
#' (tanh link, fitted on all pseudo-outcomes) and recompute the pseudo-outcomes from the same
#' standardised residuals. Cheap; no nuisance refits.
pvcor_iterate <- function(obj, strategy = names(obj$strategies)[1], se = "HC3") {
  r <- obj$strategies[[strategy]]
  tf <- if (!is.null(r$fit$tanh)) r$fit$tanh else tanh_fit(r$phi, obj$W, se)
  rho1 <- tanh(drop(obj$W %*% tf$coef))
  phi1 <- phi_from_nuisances(r$z1, r$z2, rho1)
  list(phi = phi1, blp = blp_fit(phi1, obj$W, se), tanh = tanh_fit(phi1, obj$W, se))
}
