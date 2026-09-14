fx <- function(name) test_path("fixtures", name)

## Same helper as test-filter-loci.R (each test file is self-contained, per
## the existing convention of one fx() per file rather than a shared setup
## file).
make_H <- function(A1, A2, locus = NULL, locus_raw = NULL, n_alleles = NULL,
                    alleles = NULL, samples = NULL, fields = NULL) {
  n_rec <- nrow(A1)
  if (is.null(samples)) samples <- colnames(A1)
  if (is.null(locus)) locus <- rownames(A1)
  dimnames(A1) <- dimnames(A2) <- list(NULL, samples)
  if (is.null(locus)) locus <- paste0("locus_", seq_len(n_rec))
  if (is.null(locus_raw)) locus_raw <- locus
  if (is.null(n_alleles)) n_alleles <- rep(2L, n_rec)
  if (is.null(alleles)) alleles <- replicate(n_rec, c("A", "C"), simplify = FALSE)
  list(A1 = A1, A2 = A2, locus = locus, locus_raw = locus_raw,
       alleles = alleles, n_alleles = n_alleles, samples = samples, fields = fields)
}

## ---------------------------------------------------------------------------
## write_vcf()
## ---------------------------------------------------------------------------

test_that("write_vcf() round-trips A1/A2, including masked cells, and FORMAT is GT-only", {
  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf")))
  ## Mask one real cell first, so the round-trip check also covers a
  ## genuinely-changed genotype, not just an unmodified one.
  H$A1[1, "popA_1"] <- NA; H$A2[1, "popA_1"] <- NA

  out <- tempfile(fileext = ".vcf")
  write_vcf(H, out, verbose = FALSE)
  back <- suppressMessages(read_stacks_vcf(out))

  expect_equal(back$A1, H$A1)
  expect_equal(back$A2, H$A2)
  expect_true(all(grepl("^GT$", back$fields[, "FORMAT"])))
})

test_that("write_vcf() supports gzip output", {
  H <- suppressMessages(read_stacks_vcf(fx("small.haps.vcf")))
  out <- tempfile(fileext = ".vcf.gz")
  write_vcf(H, out, verbose = FALSE)
  back <- suppressMessages(read_stacks_vcf(out))
  expect_equal(back$A1, H$A1)
})

## ---------------------------------------------------------------------------
## write_plink()
## ---------------------------------------------------------------------------

test_that("write_plink() writes one .map row per locus and one .ped row per individual", {
  H <- make_H(
    A1 = rbind(c(1L, 1L, 2L), c(1L, 2L, 1L)), A2 = rbind(c(1L, 2L, 2L), c(1L, 2L, 1L)),
    alleles = list(c("A", "C"), c("A", "C")), samples = c("s1", "s2", "s3")
  )
  prefix <- tempfile()
  write_plink(H, prefix, verbose = FALSE)
  map_lines <- readLines(paste0(prefix, ".map"))
  ped_lines <- readLines(paste0(prefix, ".ped"))
  expect_length(map_lines, 2L)
  expect_length(ped_lines, 3L)
  expect_length(strsplit(ped_lines[1], " ")[[1]], 6L + 2L * 2L)  # 6 pedigree cols + 2 loci x 2 alleles
})

test_that("write_plink() errors on a multiallelic locus unless drop_multiallelic = TRUE", {
  H <- make_H(
    A1 = rbind(c(1L, 2L, 3L)), A2 = rbind(c(1L, 2L, 3L)),
    alleles = list(c("A", "C", "T")), n_alleles = 3L, samples = c("s1", "s2", "s3")
  )
  prefix <- tempfile()
  expect_error(write_plink(H, prefix, verbose = FALSE), "biallelic")
  suppressMessages(expect_message(write_plink(H, prefix, drop_multiallelic = TRUE, verbose = TRUE),
                                  "excluding"))
  expect_length(readLines(paste0(prefix, ".map")), 0L)
})

test_that("write_plink() writes PLINK's missing-allele code for a missing genotype", {
  H <- make_H(A1 = rbind(c(NA_integer_, 1L)), A2 = rbind(c(NA_integer_, 1L)),
              samples = c("s1", "s2"))
  prefix <- tempfile()
  write_plink(H, prefix, verbose = FALSE)
  ped_lines <- readLines(paste0(prefix, ".ped"))
  fields <- strsplit(ped_lines[1], " ")[[1]]  # s1's row: missing genotype
  expect_equal(tail(fields, 2), c("0", "0"))
})

## ---------------------------------------------------------------------------
## write_structure()
## ---------------------------------------------------------------------------

test_that("write_structure() writes 2 rows per individual and STRUCTURE's -9 missing code", {
  H <- make_H(A1 = rbind(c(1L, NA_integer_)), A2 = rbind(c(2L, NA_integer_)), samples = c("s1", "s2"))
  out <- tempfile()
  write_structure(H, out, verbose = FALSE)
  lines <- readLines(out)
  expect_length(lines, 1L + 2L * 2L)  # 1 header + 2 rows per individual
  s2_rows <- lines[grepl("^s2", lines)]
  expect_length(s2_rows, 2L)
  expect_true(all(grepl("-9$", s2_rows)))
})

test_that("write_structure() emits a placeholder population column when pops is NULL", {
  H <- make_H(A1 = rbind(c(1L, 1L)), A2 = rbind(c(1L, 1L)), samples = c("s1", "s2"))
  out <- tempfile()
  suppressMessages(expect_message(write_structure(H, out, verbose = TRUE), "placeholder"))
  lines <- readLines(out)
  expect_true(all(grepl("^s\\d\\t1\\t", lines[-1])))  # every data row has population code "1"
})

## ---------------------------------------------------------------------------
## write_genepop()
## ---------------------------------------------------------------------------

test_that("write_genepop() requires popmap and lays out POP blocks correctly", {
  H <- make_H(A1 = rbind(c(1L, 2L)), A2 = rbind(c(1L, 2L)), samples = c("a1", "b1"))
  expect_error(write_genepop(H, tempfile()), "needs `popmap`")          # not given at all
  expect_error(write_genepop(H, tempfile(), popmap = NULL), "needs `popmap`")
  expect_error(write_genepop(H, tempfile(), popmap = list(popA = "a1")), "not in `popmap`")

  pops <- list(popA = "a1", popB = "b1")
  out <- tempfile()
  write_genepop(H, out, popmap = pops, verbose = FALSE)
  lines <- readLines(out)
  expect_equal(sum(lines == "POP"), 2L)          # one per population
  expect_true(any(grepl("^a1 ,", lines)))
  expect_true(any(grepl("^b1 ,", lines)))
})

test_that("write_genepop() genotype codes decode back to the same (sorted) allele pair", {
  H <- make_H(A1 = rbind(c(1L, 3L)), A2 = rbind(c(3L, 3L)), n_alleles = 3L, samples = c("a1", "b1"))
  pops <- list(popA = "a1", popB = "b1")
  out <- tempfile()
  write_genepop(H, out, popmap = pops, verbose = FALSE)
  lines <- readLines(out)
  a1_line <- lines[grepl("^a1 ,", lines)]
  geno <- trimws(strsplit(a1_line, ",")[[1]][2])
  expect_equal(geno, "0103")  # width 2 (n_alleles <= 99): sorted alleles 1 and 3
})

test_that("write_genepop() writes the all-zero code for a missing genotype", {
  H <- make_H(A1 = rbind(c(NA_integer_, 1L)), A2 = rbind(c(NA_integer_, 1L)), samples = c("a1", "b1"))
  pops <- list(popA = "a1", popB = "b1")
  out <- tempfile()
  write_genepop(H, out, popmap = pops, verbose = FALSE)
  a1_line <- readLines(out)[grepl("^a1 ,", readLines(out))]
  expect_equal(trimws(strsplit(a1_line, ",")[[1]][2]), "0000")
})

## ---------------------------------------------------------------------------
## write_fstat()
## ---------------------------------------------------------------------------

test_that("write_fstat() requires popmap and its header line matches the data", {
  H <- make_H(A1 = rbind(c(1L, 2L), c(1L, 1L)), A2 = rbind(c(1L, 2L), c(1L, 1L)), samples = c("a1", "b1"))
  expect_error(write_fstat(H, tempfile()), "needs `popmap`")

  pops <- list(popA = "a1", popB = "b1")
  out <- tempfile()
  write_fstat(H, out, popmap = pops, verbose = FALSE)
  header <- as.integer(strsplit(readLines(out)[1], " ")[[1]])
  expect_equal(header[1:2], c(2L, 2L))  # 2 populations, 2 loci
})

test_that("write_fstat() writes the same file whether or not hierfstat is installed", {
  H <- make_H(A1 = rbind(c(1L, 2L), c(1L, 1L)), A2 = rbind(c(1L, 2L), c(1L, 1L)), samples = c("a1", "b1"))
  pops <- list(popA = "a1", popB = "b1")
  out <- tempfile()
  write_fstat(H, out, popmap = pops, verbose = FALSE)
  lines <- readLines(out)
  header <- as.integer(strsplit(lines[1], " ")[[1]])
  expect_equal(header, c(2L, 2L, 2L, 2L))  # 2 pops, 2 loci, max 2 alleles seen, 2-digit codes
  expect_equal(lines[2:3], c("locus_1", "locus_2"))
  expect_equal(lines[4], "1 0101 0101")  # a1 is (1,1) at both loci -> "0101" each
})

test_that("write_fstat()'s header gives the highest allele number used, not the count observed", {
  ## A 3-allele record where only alleles 2 and 3 are carried: two alleles are
  ## observed, but the file contains allele code 03, so the header must say 3.
  H <- make_H(A1 = rbind(c(2L, 3L)), A2 = rbind(c(2L, 3L)), samples = c("a1", "b1"),
              n_alleles = 3L, alleles = list(c("A", "C", "G")))
  out <- tempfile()
  write_fstat(H, out, popmap = list(popA = "a1", popB = "b1"), verbose = FALSE)
  lines <- readLines(out)
  header <- as.integer(strsplit(lines[1], " ")[[1]])
  expect_equal(header[3], 3L)
  expect_equal(lines[3:4], c("1 0202", "2 0303"))
})

## ---------------------------------------------------------------------------
## write_radpainter()
## ---------------------------------------------------------------------------

test_that("write_radpainter() writes tab-separated haplotype pairs with no leading blank column", {
  H <- make_H(A1 = rbind(c(1L, 2L)), A2 = rbind(c(1L, 2L)),
              alleles = list(c("AACGT", "AACGG")), samples = c("s1", "s2"))
  out <- tempfile()
  write_radpainter(H, out, verbose = FALSE)
  lines <- readLines(out)
  expect_equal(lines[1], "s1\ts2")             # header: no leading blank column
  expect_equal(lines[2], "AACGT/AACGT\tAACGG/AACGG")
})

test_that("write_radpainter() writes a completely empty field for a missing genotype", {
  H <- make_H(A1 = rbind(c(NA_integer_, 1L)), A2 = rbind(c(NA_integer_, 1L)),
              alleles = list(c("A", "C")), samples = c("s1", "s2"))
  out <- tempfile()
  write_radpainter(H, out, verbose = FALSE)
  cells <- strsplit(readLines(out)[2], "\t")[[1]]
  expect_equal(cells, c("", "A/A"))
})

test_that("write_radpainter() warns when H looks like per-site SNP data", {
  H <- make_H(A1 = rbind(c(1L, 2L)), A2 = rbind(c(1L, 2L)), samples = c("s1", "s2"))  # single-letter alleles
  suppressMessages(expect_message(write_radpainter(H, tempfile(), verbose = TRUE),
                                  "looks like SNP data"))
})

test_that("writers stop on sample names their format cannot hold", {
  A <- matrix(c(1L, 2L, 1L, 2L), 2, 2)
  H <- make_H(A, A, locus = c("l1", "l2"), samples = c("ind 1", "ind,2"))
  pops <- list(p1 = c("ind 1", "ind,2"))
  ## STRUCTURE: white space separates columns.
  expect_error(write_structure(H, tempfile(), verbose = FALSE), "ind 1")
  ## Genepop: a comma ends the name, but a space is fine.
  expect_error(write_genepop(H, tempfile(), popmap = pops, verbose = FALSE), "ind,2")
  H_space <- make_H(A, A, locus = c("l1", "l2"), samples = c("ind 1", "ind 2"))
  expect_no_error(write_genepop(H_space, tempfile(), popmap = list(p1 = c("ind 1", "ind 2")),
                                verbose = FALSE))
  ## PLINK: white space, in sample or population names.
  expect_error(write_plink(H_space, tempfile(), verbose = FALSE), "white space")
})
