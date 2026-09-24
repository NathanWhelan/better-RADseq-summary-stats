# Contributing to RADdiversity

Thank you for helping. Reports of anything that looks wrong are welcome,
even if you are not sure it is a bug.

## Reporting a problem or asking a question

Open an issue at
<https://github.com/NathanWhelan/better-RADseq-summary-stats/issues>.
Please include:

- what you ran (the R code) and what you expected;
- the message or output you got;
- `sessionInfo()` and `packageVersion("RADdiversity")`;
- if you can, a small VCF and popmap that show the problem. The toy files in
  `system.file("extdata", package = "RADdiversity")` are a good starting
  point. Please do not post unpublished data you are not free to share.

## Suggesting a change to a method

Every statistic in the package cites its source, and every number in the
documentation can be re-run (`vignette("rationale")`, Appendix C). A
suggested change is easiest to judge with the same: the publication or
derivation behind it, and, where it matters, a simulation showing what
changes.

## Changing the code

1. Fork the repository and make a branch.
2. Run the tests: `devtools::test()`. Tests that compare with hierfstat,
   pegas, inbreedR, vegan or adegenet run only when those packages are
   installed.
3. Run `R CMD check --as-cran` on the built package.
4. If a change moves a reported number, the regression tests in
   `tests/testthat/test-regression.R` fail. Say in the pull request which
   numbers changed, and why.
5. The simulations behind the documentation are in `inst/sims/`; each
   script's header says how to run it and what it found.

Code comments are written for a reader who knows population genetics but
not necessarily R. Please keep new code the same way.
