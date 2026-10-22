# To do

Not part of the installed package (listed in `.Rbuildignore`).

## Gene-level summaries

Return, for each gene, a score built from its edges, with a test of how the score changes with the
covariate. This is the "which genes are differentially connected" question; the manuscript treats it
as a remark, so the package is where it should be usable.

For gene $k$ and a smooth transformation $h$, let $S_k(x)=\sum_{j\neq k}h\{\rho_{kj}(x)\}$. By the
chain rule the one-step pseudo-outcome of $h(\rho_{kj})$ is

$$\hat\phi^{h}_{kj}=h(\hat\rho_{kj})+h'(\hat\rho_{kj})\big(\hat\phi_{kj}-\hat\rho_{kj}\big),$$

and the gene's score is $\hat\phi^{h}_k=\sum_{j\neq k}\hat\phi^{h}_{kj}$.

* **Inputs:** the matrix of pseudo-outcomes (`pvcor_many(..., return_phi = TRUE)` already attaches it
  as attribute `"phi"`) and the initial correlations $\hat\rho_{kj}$. The initial correlations are not
  returned yet, so `pvcor_many()` needs to return them as well.
* **Transformation:** `h = function(r) r^2` (the per-edge score is then $2\hat\rho\hat\phi-\hat\rho^2$)
  and the Gaussian mutual information $h(\rho)=-\tfrac12\log(1-\rho^2)$ as options.
* **Test:** regress each gene's score on the working design with a robust variance. The edges of one
  gene share data, so their errors are dependent and the edge standard errors must not be added.
* **Calibration:** the permutation null (`perm_null()`) should be applied to the gene scores. Each
  edge's leftover is second order, but a gene's own mean and scale errors enter all of its $p-1$ edges
  together, so the summed leftover may grow like $p-1$ times one edge's. Check size at $p=20$ and at
  the stress case $p=100$, $N\approx100$.

## Package housekeeping and known issues

* Help pages: the exported functions have comments but no `man/` pages, and `NAMESPACE` is written by
  hand. Convert the comments to roxygen2 tags, generate `NAMESPACE` and `man/`, and add `@importFrom`
  entries for the base functions used (`setNames`, `cov`, `IQR`).
* `tanh_fit()` is fragile when the covariate has an extreme tail: with a standardized right-skewed
  covariate (exponential, maximum about 5.7 standard deviations) the Newton step hit a singular
  Hessian in 3 of 60 simulated data sets (`pvcor(..., link = "tanh")` stops with a `solve` error), and
  the reported standard error can be inflated in others. With a symmetric covariate it was stable
  (0 of 60). Add step control or a ridge, and a clear error. Similarly, `pvcor_pcor()` gave a null pair
  a slope with z of about -4.6 under the same covariate, but behaved well with a milder skew (gamma
  with shape 4: null z-statistics had mean near 0 and standard deviation near 1 over 20 data sets). A
  likely cause, not verified, is that the initial correlation is clipped (`rho_clip = 0.95`) and the
  correlation matrix is shrunk toward the identity (`eig_floor`) where the pair's correlation is near 1.
  Check both before recommending either for covariates with heavy tails.
* Decide the copyright holder in `LICENSE` (currently "pvcor authors"). The author is at a university,
  so check whether the institution has a claim before choosing.
