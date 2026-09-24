###############################################################################
#
#  inst/sims/make_ibd_example_data.R -- builds the isolation-by-distance
#  example in inst/extdata.
#
#  SIMULATED data, shaped like a Stacks 2 `populations` SNP VCF, for the
#  isolation_by_distance() examples, README and workflow vignette. The main
#  example data (make_example_data.R) have only two populations, and a Mantel
#  test needs several.
#
#    ibd_example.snps.vcf.gz      one SNP per RAD locus (FORMAT GT)
#    ibd_example_popmap.tsv       sample <TAB> population
#    ibd_example_distances.csv    river km between sites: a square matrix with
#                                 only the upper half filled in (the lower half
#                                 is left blank, which read_distances() allows)
#
#  What is built in:
#    * six sites along one river, at km 0, 15, 27, 48, 60 and 85, with 8
#      individuals each;
#    * isolation by distance along the river (a 1-D habitat): moving
#      downstream, each site's allele frequencies are drawn from the previous
#      site's (Balding & Nichols 1995), with a drift F that grows with the
#      river distance between the two sites (0.0012 per km). So sites farther
#      apart differ more;
#    * random mating within sites, and 5% of genotypes missing at random.
#
#  Run from the package source directory:  Rscript inst/sims/make_ibd_example_data.R
#
###############################################################################

set.seed(20260923)
out_dir <- file.path("inst", "extdata")
sites <- paste0("site", 1:6)
river_km <- c(0, 15, 27, 48, 60, 85)
n_per_site <- 8L
n_loci <- 600L                      # before dropping records with no variation
f_per_km <- 0.0012
missing_rate <- 0.05
bases <- c("A", "C", "G", "T")

samples <- unlist(lapply(sites, function(s) sprintf("%s_%02d", s, seq_len(n_per_site))))
pop <- rep(sites, each = n_per_site)
n <- length(samples)

## Allele frequencies at each site: a starting frequency at the top of the
## river, then one Balding-Nichols draw per step downstream.
step_down <- function(p, F) stats::rbeta(length(p), p * (1 - F) / F, (1 - p) * (1 - F) / F)
freq <- matrix(NA_real_, n_loci, length(sites))
freq[, 1] <- stats::runif(n_loci, 0.05, 0.95)
for (k in 2:length(sites))
  freq[, k] <- step_down(freq[, k - 1], f_per_km * (river_km[k] - river_km[k - 1]))

## Genotypes: two independent draws from the site's frequency (random mating).
site_of <- match(pop, sites)
alt_copies <- matrix(0L, n_loci, n)
for (i in seq_len(n))
  alt_copies[, i] <- stats::rbinom(n_loci, 2L, freq[, site_of[i]])
gt <- ifelse(alt_copies == 0L, "0/0", ifelse(alt_copies == 1L, "0/1", "1/1"))
gt[matrix(stats::runif(n_loci * n) < missing_rate, n_loci)] <- "./."

## As Stacks writes: only records variable among the typed genotypes; CHROM is
## the RAD locus, ID "locus:position:+".
typed <- gt != "./."
alt_typed <- ifelse(typed, alt_copies, 0L)
variable <- rowSums(alt_typed) > 0L & rowSums(alt_typed) < 2L * rowSums(typed)
rows <- character(0)
for (t in which(variable)) {
  pos <- sample(10:75, 1L)
  ref <- sample(bases, 1L)
  alt <- sample(setdiff(bases, ref), 1L)
  locus <- length(rows) + 1L
  rows[locus] <- paste(locus, pos, paste0(locus, ":", pos, ":+"), ref, alt, ".", "PASS", ".", "GT",
                       paste(gt[t, ], collapse = "\t"), sep = "\t")
}

header <- c(
  "##fileformat=VCFv4.2",
  "##source=\"RADdiversity isolation-by-distance example (simulated; inst/sims/make_ibd_example_data.R)\"",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples),
        collapse = "\t"))
con <- gzfile(file.path(out_dir, "ibd_example.snps.vcf.gz"), "w")
writeLines(c(header, rows), con)
close(con)

writeLines(paste(samples, pop, sep = "\t"), file.path(out_dir, "ibd_example_popmap.tsv"))

## River km between sites: upper half only, 0 on the diagonal.
km <- abs(outer(river_km, river_km, "-"))
matrix_rows <- vapply(seq_along(sites), function(i) {
  cells <- ifelse(seq_along(sites) < i, "", format(km[i, ], trim = TRUE))
  paste(c(sites[i], cells), collapse = ",")
}, character(1))
writeLines(c(paste(c("", sites), collapse = ","), matrix_rows),
           file.path(out_dir, "ibd_example_distances.csv"))

cat(sprintf("Wrote %d SNP records for %d individuals in %d sites to %s\n",
            length(rows), n, length(sites), out_dir))
