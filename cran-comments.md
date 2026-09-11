## Submission

This is a new submission to CRAN.

## R CMD check results

Checked locally with `R CMD check --as-cran` (R 4.5.3, x86_64-conda-linux-gnu,
Red Hat Enterprise Linux 8.10) and on win-builder (R-devel and R-release).

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

This NOTE is expected for a first submission and does not indicate a problem.

## Downstream dependency notes

* `Imports` is limited to base packages (`stats`, `utils`); there are no
  CRAN package dependencies beyond base R.
* `Suggests` lists `hierfstat`, `inbreedR`, `knitr`, `pegas`, `rmarkdown`,
  and `testthat`. Every use of `hierfstat`, `inbreedR`, and `pegas` in the
  package and its tests is guarded with `requireNamespace()` (in package
  code) or `testthat::skip_if_not_installed()` (in tests), so the package
  and its test suite run correctly when any of these optional packages is
  unavailable.
* There is no compiled code and no use of external system dependencies.

## Other notes

* `inst/NOTICE` documents that one internal (non-exported) helper,
  `.jost_hsht()` in `R/differentiation_stats.R`, adapts the expression of a
  published formula from the MIT-licensed `mmod` package, with full
  attribution; `mmod` is not a runtime dependency.
