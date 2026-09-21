## Submission

This is a new submission to CRAN.

## R CMD check results

Checked locally with `R CMD check --as-cran --no-manual` (R 4.5.3,
x86_64-conda-linux-gnu, Red Hat Enterprise Linux 8.10; there is no LaTeX on
that machine, so the PDF manual was not built).

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Nathan Whelan <nathan.whelan@gmail.com>'
  New submission

This NOTE is expected for a first submission and does not indicate a problem.

An earlier win-builder run caught a test failure: the Wilcoxon p-values
(`p_wilcox`/`p_wilcox_BH`) in `het_between_pops()` differed from the golden
values by ~0.004-0.006 on R 4.6.1, while every other value in the same rows
matched. R 4.6.0 changed `wilcox.test()`'s default for tied values from the
normal approximation to exact conditional inference. The package now passes
`exact` explicitly (the rule R used before 4.6.0), so its p-values no longer
depend on the R version, and the regression tests compare them at the same
1e-6 tolerance as every other column.

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
  attribution and mmod's MIT copyright and permission notice. DESCRIPTION's
  `Copyright` field points to that file. `mmod` is not a runtime dependency.
* Functions that use random numbers do not set a seed unless the user passes
  `seed`; when they do, the user's random-number state is restored on exit.
  Progress messages use `message()` and are silenced by `verbose = FALSE`.
