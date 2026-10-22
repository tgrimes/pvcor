# pvcor 0.1.0

First public version. It collects the code written for the accompanying manuscript on pseudo-value
regression for sample-specific co-expression networks (in preparation).

## Orthogonal pseudo-outcomes

* `pvcor()` fits the cross-fitted, Neyman-orthogonal pseudo-outcome regression of the conditional
  correlation of one feature pair on covariates. Second stages are the best linear projection and a
  tanh link, with HC0/HC3 standard errors. Nuisance strategies are `"stdcommon"` (shared smoother on
  the standardized scale) and `"separate"`; learners are `"gam"`, `"lm"` and `"rf"`.
* `pvcor_many()` does the same for all pairs of many features at O(p) model fits plus O(p^2 N)
  arithmetic, optionally with the comparators below.
* `pvcor_groups()` handles a group label; with `K = 1` it reproduces the group-standardized LIONESS
  arm.
* `pvcor_pcor()` extends the estimator to the conditional partial correlation, for a small number of
  features (it needs more samples than features).
* `pvcor_iterate()` gives a one-step update of the initial correlation from the tanh-link fit.

## Inference and calibration

* `blp_fit()`, `tanh_fit()` and `wild_boot_test()` are the second-stage fits and tests.
* `empnull()`, `perm_null()` and `perm_pvalues()` calibrate tests over many pairs: Efron's empirical
  null and a residual-permutation null.

## Comparators

* `pv_pearson()`, `pv_pearson_many()` and `pv_pcor()` give the exact jackknife (LIONESS-type)
  pseudo-values of the Pearson and partial correlation.
* `liquid_association()` and `covreg_gauss()` give liquid association and Gaussian covariance
  regression.

## Housekeeping

* Added `vignettes/pvcor.Rmd`, `README.md`, `NEWS.md` and `.gitignore`.
* Added `tests/testthat.R`, so `R CMD check` runs the existing tests.
