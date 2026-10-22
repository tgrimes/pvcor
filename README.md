# pvcor

**What this package is for.** Suppose you want one co-expression network per sample, and you want to ask
whether the connection between two genes changes with a phenotype, a group, or a confounder such as cell
composition. The usual route computes a "pseudo-value" for each sample (the leave-one-out edge of
LIONESS) and regresses it on the covariate. That regression does not estimate the conditional
correlation you care about: it estimates a straight-line approximation taken around the pooled data, and
anything that moves the genes' means or variances moves the answer. `pvcor` fixes this by centering each
sample's score at its own conditional law, using cross-fitted nuisance estimates, so that what you
regress on the covariate is a pseudo-outcome whose conditional mean is the conditional correlation, up
to a second-order error.

In technical terms, `pvcor` implements cross-fitted, Neyman-orthogonal pseudo-outcome regression for the
**conditional correlation** (and partial correlation) of features given covariates, with
best-linear-projection and tanh-link second stages and heteroscedasticity-consistent inference. It also
contains the jackknife (LIONESS-type) pseudo-values and the comparators used in the accompanying
manuscript (in preparation).

## Installation

```r
# install.packages("remotes")
remotes::install_github("tgrimes/pvcor", build_vignettes = TRUE)
```

## A first example

```r
library(pvcor)
set.seed(1)

## two genes whose means rise with a skewed covariate x; their conditional correlation is 0.3 for everyone
N  <- 500
x  <- as.numeric(scale(rexp(N)))
e1 <- rnorm(N); e2 <- 0.3 * e1 + sqrt(1 - 0.3^2) * rnorm(N)
y1 <- 0.8 * x + e1; y2 <- 0.8 * x + e2

fit <- pvcor(y1, y2, data = data.frame(x = x),
             w_formula = ~ x,                 # second-stage working design
             nuis_formula = ~ s(x, k = 8),    # nuisance regressions (mgcv GAMs by default)
             K = 5)                           # cross-fitting folds
fit$strategies$stdcommon$fit$blp$coef         # intercept and slope on x; the true slope is 0
```

Regressing the ordinary jackknife pseudo-values (`pv_pearson(y1, y2)`) on `x` instead gives a clearly
nonzero slope here, because the two genes share a mean trend and `x` is skewed. See the vignette
(`vignette("pvcor")`) for the full walk-through.

## Functions

| Function | Purpose |
|---|---|
| `pvcor()` | one feature pair; nuisance strategies `stdcommon` / `separate`; learners `gam`, `lm`, `rf`; oracle nuisances for simulation |
| `pvcor_many()` | all pairs of p features: feature-level fits, one common smoother, O(p) fits + O(p^2 N) arithmetic; optional comparators |
| `pvcor_groups()` | group covariate; `K = 1` reproduces the group-standardized LIONESS arm |
| `pvcor_pcor()` | conditional partial correlation for a small set of features (needs N > p) |
| `pvcor_iterate()` | one-step update of the initial correlation from the tanh-link fit |
| `blp_fit()`, `tanh_fit()`, `wild_boot_test()` | second stages and inference (HC0/HC3, wild bootstrap) |
| `empnull()`, `perm_null()`, `perm_pvalues()` | calibration over many pairs: empirical null and residual-permutation null |
| `pv_pearson()`, `pv_pearson_many()`, `pv_pcor()` | exact jackknife pseudo-values of Pearson and partial correlation (LIONESS edges) |
| `liquid_association()`, `covreg_gauss()` | comparators: liquid association; Gaussian covariance regression (ML) |

## Tests

```r
devtools::test()   # or testthat::test_local()
```

## Status

Version 0.1.0 collects the code behind the accompanying manuscript. Planned work is listed in
[TODO.md](TODO.md); changes are in [NEWS.md](NEWS.md).
