## Submission

This is a new submission to CRAN.

## R CMD check results

Checked locally with `R CMD check --as-cran` (R 4.5.3, x86_64-conda-linux-gnu,
Red Hat Enterprise Linux 8.10) and on win-builder (R-devel and R-release,
R 4.6.1, Windows Server 2022).

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission

  Possibly misspelled words in DESCRIPTION (Chesser, Cockerham, FIS, FST,
  FSTAT, Genepop, Jost's, Nei, PLINK, RADpainter, VCF, heterozygosity):
  these are all correct -- author surnames, standard population-genetics
  statistic abbreviations, and software/format names.

This NOTE is expected for a first submission and does not indicate a problem.

An earlier win-builder run caught a real, now-fixed test failure: 4 of the
package's golden-value regression tests (`p_wilcox`/`p_wilcox_BH` in
`het_between_pops()`'s pairwise test table) differed from the recorded
golden values by ~0.004-0.006 on R 4.6.1 vs. the R 4.5.3 values the golden
files were checked against, while every other value in the same rows
(means, difference, Welch's p, Hedges' g) matched exactly. That isolates
the cause to `stats::wilcox.test()`'s own tie/exact-p handling differing
slightly between R versions, not a computation change in this package, so
the regression test's tolerance for those two columns was loosened
(1e-6 -> 0.02) rather than pinning golden files to one R version's
`wilcox.test()` output.

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
