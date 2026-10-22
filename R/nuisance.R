## Nuisance estimation for the correlation pseudo-outcome (Algorithm 1, step 1).
##
## A nuisance fit for one feature pair returns, for the requested rows, the conditional means
## mu1, mu2, standard deviations s1, s2 and the initial correlation rho, together with the
## standardised residuals z1, z2.
##
## Learners
##   "gam" : mgcv GAMs. Means: gaussian; variances: Gamma(log) GAM of squared residuals;
##           correlation: gaussian GAM of products of standardised residuals.
##   "lm"  : parametric versions of the same three regressions (glm Gamma(log) for variances).
##   "rf"  : random forests (ranger, one thread) for means, squared residuals and products; training
##           rows get out-of-bag predictions, so that the variance and correlation fits see honest
##           residuals.
## Strategies
##   "separate"  : each feature's mean fitted with its own smoothing parameter.
##   "stdcommon" : after a first pass, both means are refitted on the standardised scale
##                 y_j / sigma_j(x) with ONE common smoothing parameter (the geometric mean of
##                 the two REML choices), back-transformed, and the variances refitted.
##                 This makes the noise part of the mean-estimation remainder T3 vanish
##                 (Proposition 6 of the paper). Only meaningful for learner = "gam".

#' K-fold assignment (balanced within strata when `strata` is given)
make_folds <- function(n, K = 5, strata = NULL) {
  if (K <= 1) return(rep(1L, n))
  if (is.null(strata)) return(sample(rep_len(seq_len(K), n)))
  f <- integer(n)
  for (s in unique(strata)) {
    idx <- which(strata == s)
    f[idx] <- sample(rep_len(sample(seq_len(K)), length(idx)))
  }
  f
}

## ---- learner primitives ---------------------------------------------------------------------
.rhs <- function(f) paste(deparse(f[[length(f)]], width.cutoff = 500L), collapse = " ")

.fit_mean <- function(y, data, train, learner, rhs, sp = NULL) {
  d <- data; d$.y <- y
  switch(learner,
    gam = {
      fit <- mgcv::gam(as.formula(paste(".y ~", rhs)), data = d[train, , drop = FALSE], sp = sp, method = "REML")
      list(pred = as.numeric(predict(fit, d)), sp = fit$sp)
    },
    lm = {
      fit <- lm(as.formula(paste(".y ~", rhs)), data = d[train, , drop = FALSE])
      list(pred = as.numeric(predict(fit, d)), sp = NULL)
    },
    rf = {
      ## options (for study S4c): pvcor.rf_oob (default TRUE) and pvcor.rf_min_node ("sqrt", the default,
      ## meaning max(10, sqrt(n_train)), or a number). FALSE and 10 reproduce the forest of the first draft.
      vars <- all.vars(as.formula(paste("~", rhs)))
      mns <- getOption("pvcor.rf_min_node", "sqrt")
      mns <- if (identical(mns, "sqrt")) max(10, round(sqrt(length(train)))) else as.numeric(mns)
      fit <- ranger::ranger(x = d[train, vars, drop = FALSE], y = y[train], num.trees = 300,
                            min.node.size = mns, num.threads = 1,
                            seed = sample.int(.Machine$integer.max, 1))
      pred <- as.numeric(predict(fit, d[, vars, drop = FALSE], num.threads = 1)$predictions)
      if (isTRUE(getOption("pvcor.rf_oob", TRUE)))
        pred[train] <- fit$predictions        # out-of-bag predictions on the training rows: in-bag
      list(pred = pred, sp = NULL)            # residuals of an interpolating learner are far too small
    },
    stop("unknown learner"))
}

.fit_var <- function(r, data, train, learner, rhs, floor_rel = 1e-3) {
  r2 <- r^2
  fl <- floor_rel * stats::median(r2[train])
  d <- data; d$.r2 <- pmax(r2, 1e-12 * mean(r2[train]))
  pred <- switch(learner,
    gam = {
      fit <- mgcv::gam(as.formula(paste(".r2 ~", rhs)), family = Gamma(link = "log"),
                       data = d[train, , drop = FALSE], method = "REML")
      as.numeric(predict(fit, d, type = "response"))
    },
    lm = {
      fit <- glm(as.formula(paste(".r2 ~", rhs)), family = Gamma(link = "log"), data = d[train, , drop = FALSE],
                 control = list(maxit = 100))
      as.numeric(predict(fit, d, type = "response"))
    },
    rf = .fit_mean(d$.r2, data, train, "rf", rhs)$pred,
    stop("unknown learner"))
  sqrt(pmax(pred, fl))
}

## ---- one pair ---------------------------------------------------------------------------------
#' Fit the nuisances for one feature pair on rows `train` and return them at all rows.
#'
#' @param y1,y2 numeric responses.
#' @param data data.frame of covariates used by `nuis_formula`.
#' @param train logical or integer index of training rows.
#' @param nuis_formula one-sided formula for all nuisance regressions, e.g. `~ s(x, k = 10)`.
#' @param learner "gam", "lm" or "rf".
#' @param strategies character vector, any of "separate", "stdcommon".
#' @param oracle optional list with any of mu1, mu2, s1, s2, rho (vectors over all rows) that
#'   replace the corresponding estimates.
#' @param rho_clip clip the initial correlation to [-rho_clip, rho_clip].
#' @param rho_formula optional formula for the initial correlation (default: nuis_formula).
#' Note: "stdcommon" refits each mean on the standardised scale y / sigma_hat with the same model
#'   class. That is appropriate when the class is rich enough for mu / sigma (e.g. one smooth
#'   covariate). With further covariates entering means and log-scales additively, mu / sigma is not
#'   additive and the refit is misspecified; use "separate" (see `resolve_strategy`).
#' @return named list (one element per strategy) of lists with mu1, mu2, s1, s2, rho, z1, z2.
fit_nuisances <- function(y1, y2, data, train, nuis_formula = ~ s(x, k = 10), learner = "gam",
                          strategies = "stdcommon", oracle = NULL, rho_clip = 0.99, var_floor = 1e-3,
                          rho_formula = NULL) {
  if (is.logical(train)) train <- which(train)
  rhs <- .rhs(nuis_formula); rhs_rho <- if (is.null(rho_formula)) rhs else .rhs(rho_formula)
  ## first pass: separate means, then variances
  m1 <- .fit_mean(y1, data, train, learner, rhs); m2 <- .fit_mean(y2, data, train, learner, rhs)
  mu1 <- m1$pred; mu2 <- m2$pred
  s1 <- .fit_var(y1 - mu1, data, train, learner, rhs, var_floor)
  s2 <- .fit_var(y2 - mu2, data, train, learner, rhs, var_floor)
  base <- list(mu1 = mu1, mu2 = mu2, s1 = s1, s2 = s2)
  out <- list()
  for (st in strategies) {
    nu <- base
    if (st == "stdcommon" && learner == "gam" && length(m1$sp) > 0) {
      t1 <- y1 / s1; t2 <- y2 / s2
      sp1 <- .fit_mean(t1, data, train, "gam", rhs)$sp; sp2 <- .fit_mean(t2, data, train, "gam", rhs)$sp
      sp <- exp((log(pmax(sp1, 1e-10)) + log(pmax(sp2, 1e-10))) / 2)
      nu$mu1 <- s1 * .fit_mean(t1, data, train, "gam", rhs, sp = sp)$pred
      nu$mu2 <- s2 * .fit_mean(t2, data, train, "gam", rhs, sp = sp)$pred
      nu$s1 <- .fit_var(y1 - nu$mu1, data, train, learner, rhs, var_floor)
      nu$s2 <- .fit_var(y2 - nu$mu2, data, train, learner, rhs, var_floor)
    }
    for (nm in intersect(names(oracle), c("mu1", "mu2", "s1", "s2"))) nu[[nm]] <- oracle[[nm]]
    nu$z1 <- (y1 - nu$mu1) / nu$s1; nu$z2 <- (y2 - nu$mu2) / nu$s2
    if (!is.null(oracle$rho)) {
      nu$rho <- oracle$rho
    } else {
      nu$rho <- .fit_mean(nu$z1 * nu$z2, data, train, learner, rhs_rho)$pred
    }
    nu$rho <- pmin(pmax(nu$rho, -rho_clip), rho_clip)
    out[[st]] <- nu
  }
  out
}

#' The correlation pseudo-outcome (equation (phihat) of the paper)
phi_from_nuisances <- function(z1, z2, rho) rho + z1 * z2 - rho / 2 * (z1^2 + z2^2)
