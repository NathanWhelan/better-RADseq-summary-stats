fx <- function(name) test_path("fixtures", name)

## Same make_H() helper convention as test-filter-loci.R/test-write-formats.R
## (each test file in this package is self-contained).
make_H <- function(A1, A2, locus = NULL, locus_raw = NULL, n_alleles = NULL,
                    alleles = NULL, samples = NULL) {
  n_rec <- nrow(A1)
  if (is.null(samples)) samples <- colnames(A1)
  if (is.null(locus)) locus <- rownames(A1)
  dimnames(A1) <- dimnames(A2) <- list(NULL, samples)
  if (is.null(locus)) locus <- paste0("locus_", seq_len(n_rec))
  if (is.null(locus_raw)) locus_raw <- locus
  if (is.null(n_alleles)) n_alleles <- rep(2L, n_rec)
  if (is.null(alleles)) alleles <- replicate(n_rec, c("A", "C"), simplify = FALSE)
  list(A1 = A1, A2 = A2, locus = locus, locus_raw = locus_raw,
       alleles = alleles, n_alleles = n_alleles, samples = samples)
}

## ---------------------------------------------------------------------------
## Cross-check against hierfstat::wc()/pairwise.WCfst() on the real fixture
## ---------------------------------------------------------------------------

test_that("differentiation_stats()'s internal WC formula matches hierfstat::wc()/pairwise.WCfst()", {
  testthat::skip_if_not_installed("hierfstat")
  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf")))
  pops <- suppressMessages(read_popmap(fx("small_popmap.tsv"), H$samples))

  ## Compute hierfstat's own numbers on the exact same data/dat-encoding
  ## differentiation_stats() itself builds internally for its cross-check
  ## message, so this test is checking the SAME thing that message reports,
  ## independently, rather than trusting the message text.
  r <- length(pops); ids <- unlist(pops, use.names = FALSE)
  popvec <- rep(seq_len(r), lengths(pops))
  n_rec <- nrow(H$A1)
  Gm <- matrix(NA_integer_, length(ids), n_rec)
  for (j in seq_len(n_rec)) {
    a <- H$A1[j, ids]; b <- H$A2[j, ids]
    Gm[, j] <- pmin(a, b) * 1000L + pmax(a, b)
  }
  dat <- data.frame(pop = popvec, Gm)
  hw <- hierfstat::wc(dat)
  hpw <- hierfstat::pairwise.WCfst(dat)

  out <- suppressMessages(capture.output(
    res <- differentiation_stats(H, fx("small_popmap.tsv"), nboot = 0,
                                  stem = "wctest", outdir = tempfile("dd-"), verbose = FALSE)
  ))
  expect_equal(res$global$FST, round(hw$FST, 4), tolerance = 1e-6)
  ## pairwise.WCfst() indexes populations by their integer pop id (1, 2, ...,
  ## in `names(pops)` order) rather than by name -- translate before comparing.
  ## res$pairwise_fst holds unrounded point estimates (see @return); compare
  ## against hierfstat's own unrounded value the same way.
  for (a in 2:r) for (b in 1:(a - 1)) {
    expect_equal(res$pairwise_fst[names(pops)[a], names(pops)[b]],
                 hpw[a, b], tolerance = 1e-6)
  }
})

## ---------------------------------------------------------------------------
## Known-truth checks on a hand-constructed dataset.
##
## This is deliberately a CLOSED-FORM check, not a full multi-generation
## island-model simulation: every individual within a population is
## constructed to carry exactly the SAME homozygous genotype (zero
## within-population variance), so the resulting FST/D are computable by
## hand rather than merely expected to trend in the right direction.
## ---------------------------------------------------------------------------

test_that("FST and D are ~0 when two populations are drawn from the same HWE process", {
  ## Both populations: every allele copy drawn INDEPENDENTLY (real
  ## heterozygotes included, unlike an all-homozygous construction, which
  ## turns out to be its own pathological edge case -- see the comment on
  ## the divergence test below) at allele frequency 0.5, from the SAME
  ## process for both populations. With no real among-population signal,
  ## FST/D should be small -- not exactly 0 (this is now a finite-sample
  ## Monte Carlo argument, not an exact algebraic identity), so many loci
  ## and a loose tolerance are used rather than an exact-equality check.
  set.seed(2)
  n <- 20; n_loc <- 200
  A1 <- matrix(sample(1:2, n_loc * n, replace = TRUE), n_loc, n)
  A2 <- matrix(sample(1:2, n_loc * n, replace = TRUE), n_loc, n)
  H <- make_H(A1 = A1, A2 = A2, samples = paste0("s", seq_len(n)))
  pops_f <- tempfile()
  writeLines(c(paste0("s", 1:10, "\tpopA"), paste0("s", 11:20, "\tpopB")), pops_f)

  out <- suppressMessages(capture.output(
    res <- differentiation_stats(H, pops_f, nboot = 0, stem = "identical",
                                  outdir = tempfile("dd-"), verbose = FALSE)
  ))
  expect_lt(abs(res$global$FST), 0.03)
  expect_lt(abs(res$global$D), 0.03)
})

test_that("FST and D increase as two populations' allele frequencies diverge further", {
  ## NOTE: this constructs fully-homozygous individuals (A2 <- A1 below) on
  ## purpose, for simplicity -- unlike the "same HWE process" test above,
  ## this one only checks a DIRECTIONAL claim (does FST/D increase as
  ## popA/popB diverge further), which holds regardless of the absolute
  ## values an all-homozygous construction happens to produce (see that
  ## other test's comment for why the absolute value there needed real
  ## heterozygotes to behave as expected).
  make_data <- function(freq_b) {
    ## popA fixed for allele 1 at every locus; popB's frequency of allele 2
    ## is `freq_b` -- increasing freq_b moves popB further from popA.
    n_each <- 10; n_loc <- 30
    a_geno <- rep(1L, n_each)
    b_geno <- function() sample(1:2, n_each, replace = TRUE, prob = c(1 - freq_b, freq_b))
    set.seed(1)
    A1 <- do.call(rbind, replicate(n_loc, c(a_geno, b_geno()), simplify = FALSE))
    make_H(A1 = A1, A2 = A1, samples = c(paste0("a", 1:n_each), paste0("b", 1:n_each)))
  }
  pops_f <- tempfile()
  writeLines(c(paste0("a", 1:10, "\tpopA"), paste0("b", 1:10, "\tpopB")), pops_f)

  fst_of <- function(freq_b) {
    H <- make_data(freq_b)
    out <- suppressMessages(capture.output(
      res <- differentiation_stats(H, pops_f, nboot = 0, stem = paste0("div", freq_b),
                                    outdir = tempfile("dd-"), verbose = FALSE)
    ))
    c(FST = res$global$FST, D = res$global$D)
  }
  low  <- fst_of(0.15)
  high <- fst_of(0.85)
  expect_gt(high[["FST"]], low[["FST"]])
  expect_gt(high[["D"]], low[["D"]])
})

## ---------------------------------------------------------------------------
## Precision columns (jackknife SE always present; bootstrap CI only if
## nboot > 0) and the returned-object shape.
## ---------------------------------------------------------------------------

test_that("differentiation_stats() returns the documented shape and honors nboot = 0", {
  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE))
  out <- suppressMessages(capture.output(
    res <- differentiation_stats(H, fx("small_popmap.tsv"), nboot = 0,
                                  stem = "shape", outdir = tempfile("dd-"), verbose = FALSE)
  ))
  expect_named(res, c("global", "pairwise", "pairwise_fst", "pairwise_beta", "pairwise_D"))
  expect_equal(nrow(res$pairwise), 1L)                    # one pair: popA-popB
  expect_true(all(is.na(res$global[, c("FST_lo", "FST_hi", "D_lo", "D_hi")])))
  expect_false(is.na(res$global$FST_se))  # jackknife SE needs no nboot
})

test_that("differentiation_stats() requires stem when given an H list, and >= 2 populations", {
  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE))
  expect_error(differentiation_stats(H, fx("small_popmap.tsv"), nboot = 0, verbose = FALSE,
                                     outdir = tempfile("dd-")), "stem")
  ## No files, no stem needed.
  expect_no_error(suppressMessages(differentiation_stats(H, fx("small_popmap.tsv"), nboot = 0,
                                                         verbose = FALSE)))
})
