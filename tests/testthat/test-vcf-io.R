fx <- function(name) test_path("fixtures", name)

test_that("read_stacks_vcf() parses the small fixture", {
  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf")))
  expect_setequal(H$samples,
                   c(paste0("popA_", 1:4), paste0("popB_", 1:3)))
  expect_equal(nrow(H$A1), 80L)
  expect_equal(ncol(H$A1), 7L)
  expect_true(all(H$n_alleles %in% c(2L, 3L)))
})

test_that("reading in chunks gives exactly the same object as reading in one pass", {
  for (file in c("small.haps.vcf", "small_ad.haps.vcf", "small.snps.vcf")) {
    whole <- read_stacks_vcf(fx(file), verbose = FALSE)
    chunked <- read_stacks_vcf(fx(file), chunk_lines = 3, verbose = FALSE)
    expect_identical(chunked, whole, info = file)
  }
})

test_that("read_stacks_vcf() stores depths and allele depths instead of genotype text", {
  f <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC",
    "un\t1\tl1\tA\tC,T\t.\tPASS\t.\tGT:DP:AD\t1/2:9:2,3,4\t1/1:.:1,5,0\t0/0:4",
    "un\t2\tl2\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t./."), f)
  H <- read_stacks_vcf(f, chunk_lines = 1, verbose = FALSE)
  expect_equal(colnames(H$fields), c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                                     "INFO", "FORMAT"))
  ## A: DP 9; B: DP "." so sum of AD = 6; C: DP 4, AD dropped.
  expect_equal(unname(H$depth), rbind(c(9L, 6L, 4L), c(NA, NA, NA)))
  expect_equal(unname(H$ad_ref), rbind(c(2L, 1L, NA), c(NA, NA, NA)))
  ## A is 1/2: reads for ALT 1 and ALT 2 (3 + 4); B is 1/1: ALT 1 once (5).
  expect_equal(unname(H$ad_alt), rbind(c(7L, 5L, NA), c(NA, NA, NA)))
  ## Subsetting records keeps the depth matrices lined up.
  expect_equal(nrow(filter_call_rate(H, 0.5, verbose = FALSE)$depth), 2L)
})

test_that("read_popmap() splits samples by population", {
  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE))
  pops <- suppressMessages(read_popmap(fx("small_popmap.tsv"), H$samples))
  expect_equal(names(pops), c("popA", "popB"))
  expect_equal(lengths(pops), c(popA = 4L, popB = 3L))
})

test_that("read_popmap() follows Stacks' rules and explains malformed lines", {
  write_pm <- function(lines) { f <- tempfile(fileext = ".tsv"); writeLines(lines, f); f }
  ## A third (group) column, a comment, a blank line, CRLF endings and IDs
  ## that look numeric are all fine.
  ok <- write_pm(c("# my popmap", "001\tnorth\tgroup1\r", "", "002\tnorth\tgroup1",
                   "003\tsouth\tgroup2"))
  suppressMessages(expect_message(pops <- read_popmap(ok), "third \\(group\\) column"))
  expect_equal(pops, list(north = c("001", "002"), south = "003"))
  ## A line with no population column names the line.
  expect_error(read_popmap(write_pm(c("s1\tpopA", "s2")), verbose = FALSE),
               "line 2: \"s2\" \\(no population column\\)")
  ## Space-separated columns get a hint.
  expect_error(read_popmap(write_pm(c("s1 popA", "s2 popB")), verbose = FALSE),
               "separated by spaces")
  ## Four columns.
  expect_error(read_popmap(write_pm("s1\tpopA\tg\textra"), verbose = FALSE), "4 columns")
  ## An empty population.
  expect_error(read_popmap(write_pm(c("s1\t\tg")), verbose = FALSE), "empty sample or population")
})

test_that("printing the data object gives a short summary, not the matrices", {
  H <- read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE)
  expect_s3_class(H, "raddiv_vcf")
  out <- capture.output(print(H))
  expect_lt(length(out), 20)
  expect_true(any(grepl("80 records on 80 RAD loci; 7 samples", out, fixed = TRUE)))
})

test_that("read_stacks_vcf() explains a VCF with no sample columns", {
  f <- tempfile(fileext = ".vcf")
  writeLines(c("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
               "1\t1\tx\tA\tC\t.\tPASS\t.\tGT"), f)
  expect_error(read_stacks_vcf(f, verbose = FALSE), "no sample columns")
})

test_that("read_popmap() keeps every sample when `samples` is not given", {
  pops <- read_popmap(fx("small_popmap.tsv"), verbose = FALSE)
  expect_equal(lengths(pops), c(popA = 4L, popB = 3L))
})

test_that("a popmap can be given as a path or as the list from read_popmap(), with the same result", {
  H <- read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE)
  pops <- read_popmap(fx("small_popmap.tsv"), verbose = FALSE)
  by_path <- diversity_stats(H, fx("small_popmap.tsv"), g = 4, nboot = 0, verbose = FALSE)
  by_list <- diversity_stats(H, pops, g = 4, nboot = 0, verbose = FALSE)
  expect_identical(by_path$per_population, by_list$per_population)
  expect_error(diversity_stats(H, list(c("popA_1", "popA_2")), g = 4, verbose = FALSE),
               "named list")
})

test_that("an H list needs a stem only when files are written", {
  H <- read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE)
  expect_no_error(suppressMessages(diversity_stats(H, fx("small_popmap.tsv"), g = 4, nboot = 0)))
})

test_that("read_popmap() keeps IDs as text (leading zeros and T/F survive)", {
  ## Regression: read.delim()'s type guessing turned "001" into 1 (no longer
  ## matching the VCF) and a population named "T" into TRUE.
  pm <- tempfile(fileext = ".tsv")
  on.exit(unlink(pm), add = TRUE)
  writeLines(c("001\tT", "002\tT", "010\tF"), pm)
  pops <- read_popmap(pm, c("001", "002", "010"), verbose = FALSE)
  expect_identical(pops, list(T = c("001", "002"), F = "010"))
})

## ---------------------------------------------------------------------------
## Which RAD locus each record belongs to (locus_from). A wrong answer is
## silent -- the block bootstrap would treat linked SNPs as independent -- so
## every format this has to handle gets a test.
## ---------------------------------------------------------------------------
vcf_from <- function(recs) {
  f <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
               recs), f)
  f
}
rec <- function(chrom, pos, id)
  paste(chrom, pos, id, "A", "C", ".", "PASS", ".", "GT", "0/1", "0/0", sep = "\t")
loci_of <- function(recs, ...)
  read_stacks_vcf(vcf_from(recs), verbose = FALSE, ...)$locus_raw

test_that("locus_from = 'auto' groups SNPs into RAD loci for each common VCF layout", {
  ## Stacks with locus:column:strand IDs (CHROM = "un" in Stacks 1, the locus
  ## number in Stacks 2)
  expect_equal(loci_of(c(rec("un", 1014, "1:14:+"), rec("un", 1060, "1:60:+"),
                         rec("un", 2022, "2:22:+"))), c("1", "1", "2"))
  ## Stacks 2 builds that leave ID empty: the locus is in CHROM, POS is the column
  expect_equal(length(unique(loci_of(c(rec("1", 14, "."), rec("1", 60, "."),
                                       rec("2", 22, "."))))), 2L)
  ## ipyrad: CHROM = "RAD_<n>", ID empty
  expect_equal(loci_of(c(rec("RAD_0", 5, "."), rec("RAD_0", 90, "."),
                         rec("RAD_1", 12, "."))), c("RAD_0:5", "RAD_0:5", "RAD_1:12"))
  ## Reference-aligned: two clusters of SNPs 50 kb apart on one chromosome
  expect_equal(loci_of(c(rec("chr1", 100, "."), rec("chr1", 150, "."),
                         rec("chr1", 50000, "."), rec("chr1", 50100, "."))),
               c("chr1:100", "chr1:100", "chr1:50000", "chr1:50000"))
})

test_that("locus_from can be set explicitly, and says what is wrong when it cannot apply", {
  ref <- c(rec("chr1", 100, "."), rec("chr1", 50000, "."), rec("chr2", 100, "."))
  expect_equal(loci_of(ref, locus_from = "CHROM"), c("chr1", "chr1", "chr2"))
  expect_equal(length(unique(loci_of(ref, locus_from = "window", window_bp = 1e5))), 2L)
  expect_error(loci_of(ref, locus_from = "ID"), "some records have no ID")
  expect_error(loci_of(ref, locus_from = "tag"), "must be one of")
})

test_that(".resolve_H() passes a path through to read_stacks_vcf() and an H list straight through", {
  H_from_path <- suppressMessages(RADdiversity:::.resolve_H(fx("small.haps.vcf"), verbose = FALSE))
  expect_equal(H_from_path$samples, suppressMessages(read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE))$samples)

  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf"), verbose = FALSE))
  expect_identical(RADdiversity:::.resolve_H(H, verbose = FALSE), H)
})

test_that(".resolve_H() rejects something that isn't a path or a proper H list", {
  expect_error(RADdiversity:::.resolve_H(42), "path to a VCF file")
  expect_error(RADdiversity:::.resolve_H(list(A1 = matrix(1))), "missing element")
  bad <- list(A1 = matrix(1), A2 = matrix(1:2), locus = "x", locus_raw = "x",
              alleles = list("A"), n_alleles = 1L, samples = "s1")
  expect_error(RADdiversity:::.resolve_H(bad), "same size")
})

test_that("diversity_stats()/het_between_pops() accept a pre-parsed H list in place of a path", {
  outdir <- tempfile("raddiversity-test-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE), add = TRUE)
  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf")))

  expect_error(
    suppressMessages(diversity_stats(H, fx("small_popmap.tsv"), g = 4, nboot = 10, outdir = outdir)),
    "stem"
  )
  res <- suppressMessages(capture.output(
    result <- diversity_stats(H, fx("small_popmap.tsv"), g = 4, nboot = 10, outdir = outdir, stem = "test")
  ))
  expect_s3_class(result, "raddiv_diversity")
  expect_true(file.exists(file.path(outdir, "diversity_per_population.test.tsv")))

  expect_error(
    suppressMessages(het_between_pops(H, fx("small_popmap.tsv"), min_call = 0.5, outdir = outdir)),
    "stem"
  )
  het_res <- suppressMessages(capture.output(
    het_result <- het_between_pops(H, fx("small_popmap.tsv"), min_call = 0.5,
                                   outdir = outdir, stem = "test")
  ))
  expect_s3_class(het_result, "raddiv_het")
  expect_true(file.exists(file.path(outdir, "individual_heterozygosity.test.tsv")))
})

test_that("a Stacks 2.68 de novo haplotype VCF (ID '.', POS 0) is read as one locus per record", {
  ## The layout VcfHapsExport::write_batch() writes: CHROM = locus ID, POS = 0,
  ## ID = ".", haplotype alleles joined from the SNP columns, INFO snp_columns.
  ## Locus 12 has one SNP (single-base alleles), the others several.
  f <- tempfile(fileext = ".vcf")
  set.seed(4)
  samp <- c(paste0("n", 1:5), paste0("s", 1:5))
  rec <- function(chrom, ref, alt, cols) {
    k <- length(strsplit(alt, ",")[[1]])
    gts <- vapply(1:10, function(i) paste(sort(sample(0:k, 2, TRUE)), collapse = "/"), "")
    gts[3] <- "./."
    paste(c(chrom, "0", ".", ref, alt, ".", "PASS", paste0("snp_columns=", cols), "GT", gts),
          collapse = "\t")
  }
  writeLines(c("##fileformat=VCFv4.2", "##source=\"Stacks v2.68\"",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                       samp), collapse = "\t"),
               rec("3", "ACG", "TCG,ATG", "12,40,77"),
               rec("7", "GT", "AT,GC,AC", "5,61"),
               rec("12", "C", "T", "33"),
               rec("15", "TTA", "CTA", "8,9,90")), f)
  msgs <- capture_messages(H <- read_stacks_vcf(f))
  expect_true(any(grepl("RAD loci from CHROM \\+ POS", msgs)))
  expect_equal(length(unique(H$locus_raw)), 4L)
  expect_true(RADdiversity:::.is_haplotype_H(H))
  expect_true(all(is.na(H$A1[, "n3"])))
  pm <- list(north = samp[1:5], south = samp[6:10])
  res <- suppressWarnings(diversity_stats(H, pm, g = 4, nboot = 0, verbose = FALSE))
  expect_true(res$settings$is_haplotype)
  expect_equal(res$settings$n_loci, 4L)
})

test_that("symbolic ALT alleles do not make a SNP VCF look like haplotype data, and are reported", {
  lines <- c("##fileformat=VCFv4.2",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti1\ti2",
             "chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0",
             "chr1\t20\t.\tG\tT,*\t.\tPASS\t.\tGT\t0/1\t1/2",
             "chr1\t30\t.\tC\tA,<NON_REF>\t.\tPASS\t.\tGT\t0/1\t0/0")
  vcf <- tempfile(fileext = ".vcf")
  writeLines(lines, vcf)
  msgs <- capture_messages(H <- read_stacks_vcf(vcf))
  expect_false(RADdiversity:::.is_haplotype_H(H))
  expect_true(any(grepl("2 records have a symbolic ALT allele", msgs)))
  ## A real multi-base allele still counts.
  H$alleles[[1]] <- c("AG", "CT")
  expect_true(RADdiversity:::.is_haplotype_H(H))
})

test_that("popmap samples missing from the VCF are named; no match at all shows both name sets", {
  H <- read_stacks_vcf(fx("small.snps.vcf"), verbose = FALSE)
  pops <- list(popA = c("popA_1", "popA_2", "popA_3", "popA_4", "popA_5_typo"),
               popB = c("popB_1", "popB_2", "popB_3"))
  msgs <- capture_messages(resolved <- RADdiversity:::.resolve_pops(pops, H$samples))
  expect_true(any(grepl("1 popmap sample(s) not in the VCF, so left out: popA_5_typo", msgs,
                        fixed = TRUE)))
  expect_equal(lengths(resolved), c(popA = 4L, popB = 3L))
  expect_no_message(RADdiversity:::.resolve_pops(pops, H$samples, verbose = FALSE))

  ## The same from a popmap file.
  pm <- tempfile(fileext = ".tsv")
  writeLines(c(readLines(fx("small_popmap.tsv")), "popB_9\tpopB"), pm)
  msgs <- capture_messages(read_popmap(pm, samples = H$samples))
  expect_true(any(grepl("left out: popB_9", msgs, fixed = TRUE)))

  ## Nothing matches: both sets of names are shown.
  expect_error(RADdiversity:::.resolve_pops(list(a = c("x1", "x2"), b = "y1"), H$samples,
                                            verbose = FALSE),
               "VCF samples: +popA_1.*popmap samples: x1, x2, y1")
})
