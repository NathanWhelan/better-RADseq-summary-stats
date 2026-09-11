hdr4 <- c("##fileformat=VCFv4.2",
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta1\ta2\tb1\tb2")
site <- function(chrom, pos, gts)
  paste(c(chrom, pos, ".", "A", "C", ".", "PASS", ".", "GT", gts), collapse = "\t")
write_tmp <- function(lines, ext = ".vcf") { f <- tempfile(fileext = ext); writeLines(lines, f); f }
pm4 <- write_tmp(c("a1\tpopA", "a2\tpopA", "b1\tpopB", "b2\tpopB"), ".tsv")

test_that("pi, dxy and per-individual heterozygosity match a hand calculation", {
  ## site 1: popA 0/0 0/1 -> copies 3 A, 1 C: comparisons 6, differences 3
  ##         popB 1/1 1/1 -> 0 of 6;  between: 16 comparisons, 16 - (3*0 + 1*4) = 12
  ## site 2: invariant everywhere -> popA 0/6, popB 0/6, between 0/16
  ## site 3: popA 0/1 ./. -> 1 of 1;   popB 0/0 0/1 -> 3 of 6;   between 8 - (1*3 + 1*1) = 4
  vcf <- write_tmp(c(hdr4, site("1", 1, c("0/0", "0/1", "1/1", "1/1")),
                           site("1", 2, c("0/0", "0/0", "0/0", "0/0")),
                           site("2", 1, c("0/1", "./.", "0/0", "0/1"))))
  res <- suppressMessages(pi_allsites(vcf, pm4, nboot = 0))
  expect_equal(res$pi$pi, signif(c(4 / 13, 3 / 18), 4))
  expect_equal(res$dxy$dxy, signif(16 / 40, 4))
  expect_equal(res$dxy$da, signif(16 / 40 - (4 / 13 + 3 / 18) / 2, 4))
  expect_equal(res$individual$het_sites, c(1, 1, 0, 1))
  expect_equal(res$individual$called_sites, c(3, 2, 3, 3))
})

## Simulated all-sites data: L sites, a fraction of them variable.
sim_allsites <- function(L, n_each, frac_var = 0.3, miss = 0, seed = 1,
                         p_range = c(0.05, 0.95)) {
  set.seed(seed)
  p <- ifelse(stats::runif(L) < frac_var, stats::runif(L, p_range[1], p_range[2]), 0)
  n <- 2 * n_each
  a1 <- matrix(stats::runif(L * n) < p, L) * 1L
  a2 <- matrix(stats::runif(L * n) < p, L) * 1L
  gt <- matrix(paste0(a1, "/", a2), L)
  gt[matrix(stats::runif(L * n) < miss, L)] <- "./."
  locus <- rep(seq_len(ceiling(L / 50)), each = 50)[seq_len(L)]
  lines <- vapply(seq_len(L), function(i)
    paste(c(locus[i], i, ".", "A", "C", ".", "PASS", ".", "GT", gt[i, ]), collapse = "\t"), "")
  samp <- c(paste0("a", seq_len(n_each)), paste0("b", seq_len(n_each)))
  list(lines = c("##fileformat=VCFv4.2",
                 paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                         "FORMAT", samp), collapse = "\t"), lines),
       popmap = write_tmp(paste0(samp, "\t", rep(c("popA", "popB"), each = n_each)), ".tsv"),
       p = p, a1 = a1, a2 = a2)
}

test_that("with no missing data, pi is the per-site 2n-corrected gene diversity averaged over every site", {
  d <- sim_allsites(400, 6)
  res <- suppressMessages(pi_allsites(write_tmp(d$lines), d$popmap, nboot = 0))
  per_site <- vapply(seq_len(400), function(i)
    gene_div_2n_counts(tabulate(c(d$a1[i, 1:6], d$a2[i, 1:6]) + 1L, nbins = 2)), 0)
  expect_equal(res$pi$pi[1], signif(mean(per_site), 4))
})

test_that("with missing data pixy's estimator stays unbiased; missing-as-reference does not", {
  ## ALT alleles mostly rare (frequency 0.02-0.3), as in real site-frequency
  ## spectra; 30% of genotypes missing at random.
  d <- sim_allsites(6000, 10, frac_var = 0.3, miss = 0.3, seed = 2, p_range = c(0.02, 0.3))
  truth <- mean(2 * d$p * (1 - d$p))
  res <- suppressMessages(pi_allsites(write_tmp(d$lines), d$popmap, nboot = 0))
  expect_lt(abs(mean(res$pi$pi) / truth - 1), 0.05)
  ## The common shortcut -- call every missing genotype homozygous reference --
  ## pulls each ALT frequency p down to about 0.7p, so with mostly-rare ALT
  ## alleles it underestimates pi (with common ALT alleles it would
  ## overestimate: the bias follows the allele frequencies, not a fixed sign).
  naive <- write_tmp(gsub("\\./\\.", "0/0", d$lines))
  res_naive <- suppressMessages(pi_allsites(naive, d$popmap, nboot = 0))
  expect_lt(mean(res_naive$pi$pi) / truth - 1, -0.15)
})

test_that("results do not depend on chunk size, and .gz input works", {
  d <- sim_allsites(300, 5, miss = 0.1, seed = 3)
  f <- write_tmp(d$lines)
  big   <- suppressMessages(pi_allsites(f, d$popmap, nboot = 50, chunk_lines = 1e5))
  small <- suppressMessages(pi_allsites(f, d$popmap, nboot = 50, chunk_lines = 7))
  expect_equal(small$pi, big$pi)
  expect_equal(small$dxy, big$dxy)
  gz <- tempfile(fileext = ".vcf.gz"); con <- gzfile(gz, "w"); writeLines(d$lines, con); close(con)
  expect_equal(suppressMessages(pi_allsites(gz, d$popmap, nboot = 50))$pi, big$pi)
  expect_equal(big$pi$sites[1], 300L)
})

test_that("complete_sites uses only sites typed in every individual of the population", {
  vcf <- write_tmp(c(hdr4, site("1", 1, c("0/1", "./.", "0/0", "0/1")),
                           site("1", 2, c("0/1", "0/0", "0/0", "0/0"))))
  res <- suppressMessages(pi_allsites(vcf, pm4, nboot = 0, complete_sites = TRUE))
  expect_equal(res$pi$sites, c(1L, 2L))           # popA loses site 1
  expect_equal(res$pi$pi[1], signif(3 / 6, 4))    # site 2 only: copies 3 A, 1 C
})
