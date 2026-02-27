# CRAN submission comments

## R CMD check results

0 errors | 0 warnings | 0 notes

Tested on:
- macOS aarch64 (local, R 4.5.1)

## Notes for reviewers

- The `\donttest{}` examples include calls to `kpop()` and `kpop_summary()`
  which run entropy balancing and cross-validated ridge regression and take
  approximately 60 seconds on the test machine. They are wrapped in
  `\donttest{}` for this reason.
- `glmnet` is a hard dependency (in Imports), required for the linearized
  variance estimator in `calc_linearized_var()`.
