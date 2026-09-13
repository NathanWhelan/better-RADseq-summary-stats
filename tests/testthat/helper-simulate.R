## Helpers for tests that need genotypes with a known structure, built in
## memory rather than read from a fixture file.

## A raddiv_vcf object from two allele matrices (records x samples, 1 = REF,
## 2 = ALT, NA = missing), with `records_per_locus` biallelic SNP records on
## each RAD locus.
sim_H <- function(A1, A2, records_per_locus = 1L) {
  n_rec <- nrow(A1)
  structure(list(A1 = A1, A2 = A2, locus = paste0("r", seq_len(n_rec)),
                 locus_raw = paste0("L", ceiling(seq_len(n_rec) / records_per_locus)),
                 alleles = rep(list(c("A", "C")), n_rec), n_alleles = rep(2L, n_rec),
                 samples = colnames(A1)),
            class = "raddiv_vcf")
}

## Genotypes of `n` individuals drawn under Hardy-Weinberg from the per-record
## ALT allele frequencies `p`: list(A1, A2), records x individuals.
sim_genotypes <- function(p, n) {
  L <- length(p)
  list(A1 = 1L + matrix(stats::rbinom(L * n, 1, p), L, n),
       A2 = 1L + matrix(stats::rbinom(L * n, 1, p), L, n))
}

## Population allele frequencies around an ancestral `p_anc`, differentiated
## by `fst` (Balding-Nichols).
sim_diverged_freqs <- function(p_anc, fst) {
  stats::rbeta(length(p_anc), p_anc * (1 - fst) / fst, (1 - p_anc) * (1 - fst) / fst)
}
