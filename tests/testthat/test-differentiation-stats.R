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
  ## hierfstat is given only the records typed in every population: at a
  ## record typed in just one population, hierfstat keeps the within-population
  ## variance and this package deliberately skips the record (see
  ## ?differentiation_stats). On every other record the two must agree.
  r <- length(pops); ids <- unlist(pops, use.names = FALSE)
  typed_everywhere <- rowSums(RADdiversity:::.typed_by_pop(H, pops) > 0) == r
  H_typed <- RADdiversity:::.subset_H(H, typed_everywhere)
  popvec <- rep(seq_len(r), lengths(pops))
  n_rec <- nrow(H_typed$A1)
  Gm <- matrix(NA_integer_, length(ids), n_rec)
  for (j in seq_len(n_rec)) {
    a <- H_typed$A1[j, ids]; b <- H_typed$A2[j, ids]
    Gm[, j] <- pmin(a, b) * 1000L + pmax(a, b)
  }
  dat <- data.frame(pop = popvec, Gm)
  hw <- hierfstat::wc(dat)
  hpw <- hierfstat::pairwise.WCfst(dat)

  out <- suppressMessages(capture.output(
    res <- differentiation_stats(H, fx("small_popmap.tsv"), nboot = 0,
                                  stem = "wctest", outdir = tempfile("dd-"), verbose = FALSE)
  ))
  expect_equal(res$global$FST, hw$FST, tolerance = 1e-6)
  ## pairwise.WCfst() indexes populations by their integer pop id (1, 2, ...,
  ## in `names(pops)` order) rather than by name -- translate before comparing.
  ## res$pairwise_fst holds unrounded point estimates (see @return); compare
  ## against hierfstat's own unrounded value the same way.
  for (a in 2:r) for (b in 1:(a - 1)) {
    expect_equal(res$pairwise_fst[names(pops)[a], names(pops)[b]],
                 hpw[a, b], tolerance = 1e-6)
  }
})

test_that("FST ignores records typed in only one population (no dilution toward 0)", {
  ## Two populations differentiated at FST ~ 0.1. Adding records at which
  ## population Y has no genotyped individual -- as Stacks' -r leaves them
  ## with 3+ populations and a low -p -- must not change FST or D.
  set.seed(3)
  L <- 1500; n <- 15
  p_anc <- stats::rbeta(L, 2, 2)
  X <- sim_genotypes(sim_diverged_freqs(p_anc, 0.1), n)
  Y <- sim_genotypes(sim_diverged_freqs(p_anc, 0.1), n)
  ids <- c(paste0("x", 1:n), paste0("y", 1:n))
  A1 <- cbind(X$A1, Y$A1); A2 <- cbind(X$A2, Y$A2)
  colnames(A1) <- colnames(A2) <- ids
  pops <- list(X = ids[1:n], Y = ids[n + 1:n])

  extra <- sim_genotypes(stats::rbeta(600, 2, 2), 2 * n)
  B1 <- extra$A1; B2 <- extra$A2
  B1[, n + 1:n] <- NA; B2[, n + 1:n] <- NA
  colnames(B1) <- colnames(B2) <- ids

  both <- differentiation_stats(sim_H(A1, A2), pops, nboot = 0, verbose = FALSE)
  with_extra <- differentiation_stats(sim_H(rbind(A1, B1), rbind(A2, B2)), pops,
                                      nboot = 0, verbose = FALSE)
  expect_equal(with_extra$global$FST, both$global$FST, tolerance = 1e-12)
  expect_equal(with_extra$pairwise$FST, both$pairwise$FST, tolerance = 1e-12)
  expect_equal(with_extra$global$D, both$global$D, tolerance = 1e-12)
  expect_equal(unname(with_extra$settings$fst_records_skipped), c(600, 600))
  expect_equal(unname(both$settings$fst_records_skipped), c(0, 0))

  ## hierfstat::wc() keeps those records, which pulls its FST down.
  skip_if_not_installed("hierfstat")
  H_extra <- sim_H(rbind(A1, B1), rbind(A2, B2))
  hw <- hierfstat::wc(RADdiversity:::.to_hierfstat_df(H_extra, pops))
  expect_lt(hw$FST, 0.9 * with_extra$global$FST)
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
## Which records the global D and FST use
## ---------------------------------------------------------------------------

## Three populations, each fixed for its own allele at every record (true
## D = 1, FST = 1). Population c is untyped at the records in `untyped_c`.
three_fixed_pops <- function(L = 400, n = 10, untyped_c = seq_len(L / 2)) {
  ids <- paste0(rep(c("a", "b", "c"), each = n), seq_len(n))
  A <- matrix(rep(rep(1:3, each = n), L), L, byrow = TRUE, dimnames = list(NULL, ids))
  A[untyped_c, ids[2 * n + seq_len(n)]] <- NA
  H <- structure(list(A1 = A, A2 = A, locus = paste0("r", seq_len(L)),
                      locus_raw = paste0("L", seq_len(L)),
                      alleles = rep(list(c("A", "C", "G")), L), n_alleles = rep(3L, L),
                      samples = ids), class = "raddiv_vcf")
  list(H = H, pops = list(a = ids[1:n], b = ids[n + 1:n], c = ids[2 * n + 1:n]))
}

test_that("global D uses only records typed in every population, and says so", {
  d <- three_fixed_pops()
  res <- differentiation_stats(d$H, d$pops, nboot = 0, beta = FALSE, verbose = FALSE)
  ## Averaging in the records without population c gave 0.875 before this fix.
  expect_equal(res$global$D, 1)
  expect_equal(res$pairwise$D, c(1, 1, 1))
  complete <- differentiation_stats(RADdiversity:::.subset_H(d$H, 201:400), d$pops,
                                    nboot = 0, beta = FALSE, verbose = FALSE)
  expect_equal(res$global$D, complete$global$D)
  expect_equal(res$global$D_records, 200L)
  expect_equal(res$settings$d_records_skipped_global, 200L)
  printed <- capture.output(print(res))
  expect_true(any(grepl("global D uses only the 200 of 400 records typed in all 3 populations",
                        printed, fixed = TRUE)))
  expect_true(any(grepl("global D uses only", capture.output(print(summary(res))))))
  msgs <- capture_messages(differentiation_stats(d$H, d$pops, nboot = 0, beta = FALSE))
  expect_true(any(grepl("global D uses only", msgs)))
})

test_that("global D is NA, with a warning, when no record is typed in every population", {
  d <- three_fixed_pops(L = 20, untyped_c = 1:20)
  d$H$A1[, d$pops$c] <- 1L; d$H$A2[, d$pops$c] <- 1L        # c typed everywhere ...
  d$H$A1[1:20, d$pops$b] <- NA; d$H$A2[1:20, d$pops$b] <- NA # ... but b nowhere
  d$H$A1[1:10, d$pops$b] <- 2L; d$H$A2[1:10, d$pops$b] <- 2L # b back at 1-10
  d$H$A1[1:10, d$pops$a] <- NA; d$H$A2[1:10, d$pops$a] <- NA # a missing at 1-10
  expect_warning(res <- differentiation_stats(d$H, d$pops, nboot = 0, beta = FALSE,
                                              verbose = FALSE), "Global D is NA")
  expect_true(is.na(res$global$D))
  expect_equal(res$global$D_records, 0L)
  ## a and b are never typed together; a-c and b-c each have 10 records.
  expect_equal(is.na(res$pairwise$D), c(TRUE, FALSE, FALSE))
})

test_that("FST and FIS skip records with a single genotyped individual per population", {
  set.seed(8)
  L <- 800; n <- 12
  p_anc <- stats::rbeta(L, 2, 2)
  X <- sim_genotypes(sim_diverged_freqs(p_anc, 0.1), n)
  Y <- sim_genotypes(sim_diverged_freqs(p_anc, 0.1), n)
  ids <- c(paste0("x", 1:n), paste0("y", 1:n))
  A1 <- cbind(X$A1, Y$A1); A2 <- cbind(X$A2, Y$A2)
  colnames(A1) <- colnames(A2) <- ids
  pops <- list(X = ids[1:n], Y = ids[n + 1:n])
  ## 300 extra records: one heterozygous individual typed in each population.
  B1 <- matrix(NA_integer_, 300, 2 * n, dimnames = list(NULL, ids)); B2 <- B1
  B1[, c(1, n + 1)] <- 1L; B2[, c(1, n + 1)] <- 2L
  base <- differentiation_stats(sim_H(A1, A2), pops, nboot = 0, beta = FALSE, verbose = FALSE)
  more <- differentiation_stats(sim_H(rbind(A1, B1), rbind(A2, B2)), pops, nboot = 0,
                                beta = FALSE, verbose = FALSE)
  expect_equal(more$global$FST, base$global$FST, tolerance = 1e-12)
  expect_equal(more$global$FIS, base$global$FIS, tolerance = 1e-12)
  expect_equal(unname(more$settings$fst_records_skipped), c(300, 300))
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
  expect_named(res, c("global", "pairwise", "pairwise_fst", "pairwise_beta", "pairwise_D",
                      "settings"))
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

## hierfstat guesses how many digits each allele has from the genotype numbers.
## With 3 digits per allele and every number below 10000 (alleles 10 and above
## present, but never two of them in one individual), it read 1012 as alleles
## 10 and 12, which gave a wrong beta. .to_hierfstat_df() now writes 2 digits
## per allele over renumbered alleles, which hierfstat always reads correctly.
test_that("genotypes reach hierfstat as numbers it decodes correctly when alleles reach 10+", {
  skip_if_not_installed("hierfstat")
  set.seed(8)
  n_rec <- 150
  draw <- function(n, probs) matrix(sample.int(length(probs), n_rec * n, replace = TRUE, prob = probs),
                                    n_rec, n)
  pa <- c(.3, .2, .15, .1, .05, .05, .05, .03, .03, .02, .01, .01)
  A1 <- cbind(draw(12, pa), draw(12, rev(pa)))
  A2 <- cbind(draw(12, pa), draw(12, rev(pa)))
  both_high <- A1 >= 10 & A2 >= 10
  A1[both_high] <- 1L                        # every old 3-digit number stays below 10000
  samples <- c(paste0("a", 1:12), paste0("b", 1:12))
  H <- make_H(A1, A2, n_alleles = rep(12L, n_rec),
              alleles = replicate(n_rec, paste0("AC", LETTERS[1:12]), simplify = FALSE),
              samples = samples)
  pops <- list(A = samples[1:12], B = samples[13:24])

  ## hierfstat recovers every locus's allele counts (labels aside).
  dat <- RADdiversity:::.to_hierfstat_df(H, pops)
  decoded <- hierfstat::getal(dat)
  for (j in c(1, 50, 150)) {
    original <- sort(as.vector(table(c(A1[j, ], A2[j, ]))))
    expect_equal(sort(as.vector(table(decoded[[paste0("L", j)]]))), original)
  }

  ## beta equals hierfstat's own on unambiguous 2-digit numbers of the original
  ## alleles, and the wc() cross-check no longer disagrees.
  reference <- data.frame(pop = rep(1:2, each = 12), t(pmin(A1, A2) * 100L + pmax(A1, A2)))
  expected_beta <- hierfstat::pairwise.betas(reference)[2, 1]
  expect_no_warning(res <- differentiation_stats(H, pops, nboot = 0, hierfstat_check = TRUE,
                                                 verbose = FALSE))
  expect_equal(res$pairwise$beta, expected_beta, tolerance = 1e-12)
  expect_equal(res$global$FST, hierfstat::wc(reference)$FST, tolerance = 1e-8)
})

test_that(".to_hierfstat_df() writes 2 digits per allele unless a record has more than 99 alleles", {
  A1 <- matrix(c(1L, 3L, 3L, 1L), 1, dimnames = list(NULL, c("a", "b", "c", "d")))
  A2 <- matrix(c(3L, 3L, NA, 1L), 1, dimnames = list(NULL, c("a", "b", "c", "d")))
  H <- make_H(A1, A2, n_alleles = 3L, alleles = list(c("A", "C", "G")))
  dat <- RADdiversity:::.to_hierfstat_df(H, list(p = c("a", "b"), q = c("c", "d")))
  ## alleles 1 and 3 are carried, so they become 1 and 2
  expect_equal(unname(unlist(dat[, 2])), c(102L, 202L, NA, 101L))
  expect_true(RADdiversity:::.hierfstat_reads_3_digits(c(1001L, 12100L)))
  expect_false(RADdiversity:::.hierfstat_reads_3_digits(c(1001L, 1123L)))
})

test_that("pairs whose population names would give the same label keep their own values", {
  ## "a" + "b__c" and "a__b" + "c" would both be labelled "a__b__c".
  set.seed(21)
  ns <- c(a = 6, a__b = 6, b__c = 6, c = 6)
  freqs <- lapply(seq_along(ns), function(i) stats::runif(120, 0.05 * i, 0.2 * i))
  gs <- Map(sim_genotypes, freqs, ns)
  A1 <- do.call(cbind, lapply(gs, `[[`, "A1")); A2 <- do.call(cbind, lapply(gs, `[[`, "A2"))
  ids <- unlist(Map(function(p, n) paste0(p, "_", seq_len(n)), names(ns), ns))
  colnames(A1) <- colnames(A2) <- ids
  pops <- split(ids, factor(rep(names(ns), ns), levels = names(ns)))
  H <- sim_H(A1, A2)
  res <- differentiation_stats(H, pops, nboot = 0, beta = FALSE, verbose = FALSE)
  one_pair <- function(p1, p2)
    differentiation_stats(H, pops[c(p1, p2)], nboot = 0, beta = FALSE, verbose = FALSE)$global
  expect_equal(res$pairwise_fst["a", "b__c"], one_pair("a", "b__c")$FST)
  expect_equal(res$pairwise_fst["a__b", "c"], one_pair("a__b", "c")$FST)
  expect_false(isTRUE(all.equal(res$pairwise_fst["a", "b__c"], res$pairwise_fst["a__b", "c"])))
  expect_equal(anyDuplicated(names(res$settings$fst_records_skipped)), 0L)

  vcf <- tempfile(fileext = ".vcf")
  H$fields <- cbind(CHROM = "1", POS = seq_len(nrow(A1)), ID = ".", REF = "A", ALT = "C",
                    QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT")
  write_vcf(H, vcf, verbose = FALSE)
  pi_res <- pi_allsites(vcf, pops, locus_from = "window", window_bp = 0, nboot = 0, verbose = FALSE)
  one_dxy <- function(p1, p2)
    pi_allsites(vcf, pops[c(p1, p2)], locus_from = "window", window_bp = 0, nboot = 0,
                verbose = FALSE)$dxy$dxy
  expect_equal(pi_res$dxy$dxy[pi_res$dxy$pop1 == "a" & pi_res$dxy$pop2 == "b__c"], one_dxy("a", "b__c"))
  expect_equal(pi_res$dxy$dxy[pi_res$dxy$pop1 == "a__b" & pi_res$dxy$pop2 == "c"], one_dxy("a__b", "c"))
})
