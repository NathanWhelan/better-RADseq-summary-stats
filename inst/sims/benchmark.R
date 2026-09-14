###############################################################################
#
#  inst/sims/benchmark.R -- time and memory on realistically sized data.
#
#  Writes a synthetic Stacks-like SNP VCF (FORMAT GT:DP:AD:GQ:GL, 1-4 SNPs per
#  RAD tag, 10% missing genotypes, three populations) and times one package
#  function on it. Peak memory is best measured from outside R, one function
#  per R process:
#
#    Rscript inst/sims/benchmark.R make 100000 200 /tmp/bench
#    for fn in read filter_depth diversity het differentiation kinship hwe; do
#      /usr/bin/time -f "$fn %e s, peak %M KB" \
#        Rscript inst/sims/benchmark.R run $fn /tmp/bench
#    done
#
#  `make <records> <samples> <stem>` writes <stem>.snps.vcf.gz and
#  <stem>.popmap.tsv; `run <function> <stem>` reads them and prints the
#  elapsed time of the call itself (reading included, since every function
#  starts by reading). Needs the package installed, or run under
#  devtools::load_all().
#
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: see the header of inst/sims/benchmark.R")

make_vcf <- function(n_rec, n_samp, stem, chunk = 20000L) {
  set.seed(1)
  pops <- rep(c("popA", "popB", "popC"), length.out = n_samp)
  samples <- sprintf("ind%04d", seq_len(n_samp))
  writeLines(paste(samples, pops, sep = "\t"), paste0(stem, ".popmap.tsv"))
  con <- gzfile(paste0(stem, ".snps.vcf.gz"), "w")
  on.exit(close(con))
  writeLines(c("##fileformat=VCFv4.2", "##source=\"RADdiversity benchmark\"",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                       "FORMAT", samples), collapse = "\t")), con)
  written <- 0L
  tag <- 0L
  while (written < n_rec) {
    m <- min(chunk, n_rec - written)
    ## RAD tags of 1-4 SNPs.
    snps_per_tag <- sample(1:4, m, replace = TRUE, prob = c(0.45, 0.3, 0.17, 0.08))
    tag_of <- rep(tag + seq_along(snps_per_tag), snps_per_tag)[seq_len(m)]
    tag <- max(tag_of)
    col <- stats::ave(tag_of, tag_of, FUN = seq_along) * 7L
    p <- pmin(pmax(stats::rbeta(m, 0.4, 0.4), 0.02), 0.98)
    a1 <- matrix(stats::runif(m * n_samp) < p, m)
    a2 <- matrix(stats::runif(m * n_samp) < p, m)
    dosage <- a1 + a2
    depth <- matrix(stats::rpois(m * n_samp, 20) + 1L, m)
    alt_reads <- ifelse(dosage == 0, 0L, ifelse(dosage == 2, depth,
                        stats::rbinom(m * n_samp, depth, 0.5)))
    gt <- c("0/0", "0/1", "1/1")[dosage + 1L]
    cell <- paste0(gt, ":", depth, ":", depth - alt_reads, ",", alt_reads, ":",
                   sample(10:40, m * n_samp, TRUE), ":",
                   sprintf("%.2f,%.2f,%.2f", -stats::runif(m * n_samp, 0, 20),
                           -stats::runif(m * n_samp, 0, 2), -stats::runif(m * n_samp, 0, 20)))
    cell[stats::runif(m * n_samp) < 0.10] <- "./."
    cell <- matrix(cell, m)
    fixed <- paste(tag_of, col, paste0(tag_of, ":", col, ":+"), "A", "G", ".", "PASS",
                   ".", "GT:DP:AD:GQ:GL", sep = "\t")
    writeLines(paste(fixed, do.call(paste, c(asplit(cell, 2), sep = "\t")), sep = "\t"), con)
    written <- written + m
  }
  invisible(NULL)
}

run_one <- function(fn, stem) {
  suppressMessages(library(RADdiversity))
  vcf <- paste0(stem, ".snps.vcf.gz")
  popmap <- paste0(stem, ".popmap.tsv")
  q <- function(expr) suppressMessages(suppressWarnings(expr))
  elapsed <- system.time(switch(fn,
    read            = q(read_stacks_vcf(vcf, verbose = FALSE)),
    filter_depth    = q(filter_genotype_depth(read_stacks_vcf(vcf, verbose = FALSE), min_dp = 6,
                                              verbose = FALSE)),
    diversity       = q(diversity_stats(vcf, popmap, g = 20, nboot = 1000, seed = 1,
                                        verbose = FALSE)),
    het             = q(het_between_pops(vcf, popmap, seed = 1, verbose = FALSE)),
    differentiation = q(differentiation_stats(vcf, popmap, nboot = 1000, seed = 1,
                                              verbose = FALSE)),
    differentiation_nobeta = q(differentiation_stats(vcf, popmap, nboot = 1000, beta = FALSE,
                                                     seed = 1, verbose = FALSE)),
    kinship         = q(kinship_check(vcf, popmap, verbose = FALSE)),
    hwe             = q(hwe_test(vcf, popmap, seed = 1, verbose = FALSE)),
    stop("Unknown function: ", fn)))[["elapsed"]]
  cat(sprintf("%s: %.1f s in R\n", fn, elapsed))
}

if (args[1] == "make") {
  make_vcf(as.integer(args[2]), as.integer(args[3]), args[4])
} else if (args[1] == "run") {
  run_one(args[2], args[3])
} else {
  stop("First argument must be `make` or `run`.")
}
