## as_genind() and as_genlight(): the data object as adegenet classes.

ex <- function(f) system.file("extdata", f, package = "RADdiversity")

test_that("as_genind() holds the same genotypes as H, for multi-allelic haplotype loci", {
  skip_if_not_installed("adegenet")
  H <- read_stacks_vcf(ex("small.haps.vcf"), verbose = FALSE)
  gi <- as_genind(H, ex("small_popmap.tsv"))
  expect_s4_class(gi, "genind")
  expect_equal(adegenet::nInd(gi), length(H$samples))
  expect_equal(adegenet::nLoc(gi), nrow(H$A1))
  expect_true(any(adegenet::nAll(gi) > 2L))                   # multi-allelic loci kept
  ## genotypes, as allele sequences, match H (order within a genotype aside)
  df <- adegenet::genind2df(gi, sep = "/", usepop = FALSE)
  for (j in seq_len(nrow(H$A1))) {
    a <- H$alleles[[j]]
    want <- ifelse(is.na(H$A1[j, ]), NA,
                   paste(a[pmin(H$A1[j, ], H$A2[j, ])], a[pmax(H$A1[j, ], H$A2[j, ])], sep = "/"))
    got <- vapply(strsplit(df[[j]], "/"), function(x)
      if (anyNA(x)) NA_character_ else paste(sort(x), collapse = "/"), character(1))
    want <- vapply(strsplit(want, "/"), function(x)
      if (anyNA(x)) NA_character_ else paste(sort(x), collapse = "/"), character(1))
    expect_identical(unname(got), unname(want), info = j)
  }
  ## populations in popmap order
  expect_identical(levels(adegenet::pop(gi)), c("popA", "popB"))
  expect_identical(as.character(adegenet::pop(gi)), rep(c("popA", "popB"), c(4, 3)))
})

test_that("as_genind() replaces dots in locus names and keeps the originals", {
  skip_if_not_installed("adegenet")
  H <- read_stacks_vcf(ex("example.snps.vcf.gz"), verbose = FALSE)
  expect_true(any(grepl(".", H$locus, fixed = TRUE)))        # e.g. "6.1": several SNPs on a tag
  gi <- as_genind(H)
  expect_false(any(grepl(".", adegenet::locNames(gi), fixed = TRUE)))
  expect_identical(unname(adegenet::other(gi)$locus), H$locus)
  expect_identical(adegenet::other(gi)$position, H$fields[, "POS"])
  expect_null(adegenet::pop(gi))
})

test_that("as_genlight() stores the ALT dosage, positions and SNP alleles", {
  skip_if_not_installed("adegenet")
  H <- read_stacks_vcf(ex("example.snps.vcf.gz"), verbose = FALSE)
  gl <- as_genlight(H, ex("example_popmap.tsv"))
  expect_s4_class(gl, "genlight")
  dosage <- t((H$A1 == 2L) + (H$A2 == 2L))
  expect_equal(unname(as.matrix(gl)), unname(dosage))        # NA where missing
  expect_identical(adegenet::locNames(gl), H$locus)
  expect_identical(adegenet::position(gl), as.integer(H$fields[, "POS"]))
  expect_identical(adegenet::alleles(gl)[1], paste(H$alleles[[1]], collapse = "/"))
  expect_identical(levels(adegenet::pop(gl)), c("north", "south"))
})

test_that("as_genlight() stops on multi-allelic records unless told to drop them", {
  skip_if_not_installed("adegenet")
  H <- read_stacks_vcf(ex("small.haps.vcf"), verbose = FALSE)
  expect_error(as_genlight(H), "more than 2 observed alleles")
  gl <- as_genlight(H, drop_multiallelic = TRUE, verbose = FALSE)
  n_biallelic <- sum(locus_allele_stats(H)$n_observed_alleles <= 2L)
  expect_equal(adegenet::nLoc(gl), n_biallelic)
  ## haplotype alleles are not single bases, so they go in other() instead
  expect_true(all(grepl("/", adegenet::other(gl)$alleles)))
})

test_that("the adegenet converters follow the writers' popmap rule and say what is missing", {
  skip_if_not_installed("adegenet")
  H <- read_stacks_vcf(ex("small.snps.vcf"), verbose = FALSE)
  part <- list(popA = c("popA_1", "popA_2", "popA_3", "popA_4"), popB = c("popB_1", "popB_2"))
  expect_error(as_genind(H, part), "filter_samples")
  expect_error(as_genlight(H, part), "filter_samples")
  gi <- as_genind(filter_samples(H, part, verbose = FALSE), part)
  expect_equal(adegenet::nInd(gi), 6L)
  ## a record nobody is genotyped at is left out, with a message
  H$A1[1, ] <- NA; H$A2[1, ] <- NA
  expect_message(gi <- as_genind(H), "no genotyped individual")
  expect_equal(adegenet::nLoc(gi), nrow(H$A1) - 1L)
  ## without adegenet, a clear message
  local_mocked_bindings(.adegenet_available = function() FALSE)
  expect_error(as_genind(H), "needs the adegenet package")
  expect_error(as_genlight(H), "needs the adegenet package")
})

test_that("as_genlight() keeps REF/ALT for a SNP fixed for ALT, so subsets agree", {
  skip_if_not_installed("adegenet")
  H <- read_stacks_vcf(ex("example.snps.vcf.gz"), verbose = FALSE)
  pm <- read_popmap(ex("example_popmap.tsv"), verbose = FALSE)
  south <- filter_samples(H, pm["south"], verbose = FALSE)
  typed <- rowSums(!is.na(south$A1)) > 0L
  south <- RADdiversity:::.subset_H(south, typed)       # as as_genlight() does
  ## records where every typed south individual is ALT/ALT
  fixed_alt <- which(rowSums(south$A1 == 2L & south$A2 == 2L, na.rm = TRUE) ==
                       rowSums(!is.na(south$A1)))
  expect_gt(length(fixed_alt), 0L)
  gl <- as_genlight(south, verbose = FALSE)
  j <- fixed_alt[1]
  expect_identical(adegenet::alleles(gl)[j], paste(south$alleles[[j]], collapse = "/"))
  expect_true(all(as.matrix(gl)[!is.na(south$A1[j, ]), j] == 2))
  ## every record is coded as copies of its ALT allele
  expect_equal(unname(as.matrix(gl)), unname(t((south$A1 == 2L) + (south$A2 == 2L))))
  ## and the same SNP has the same alleles in another subset
  gl_north <- as_genlight(filter_samples(H, pm["north"], verbose = FALSE), verbose = FALSE)
  common <- intersect(adegenet::locNames(gl), adegenet::locNames(gl_north))
  expect_identical(adegenet::alleles(gl)[match(common, adegenet::locNames(gl))],
                   adegenet::alleles(gl_north)[match(common, adegenet::locNames(gl_north))])
})
