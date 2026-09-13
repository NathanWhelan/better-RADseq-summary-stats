#' RADdiversity: Population Genetic Statistics from RAD-seq Data
#'
#' Summary statistics for RAD-seq data, designed around Stacks 2
#' `populations` output:
#'
#' * [diversity_stats()]: Ho, He (Nei & Chesser 1983), FIS as a ratio of sums,
#'   rarefied allelic and private allelic richness, per-sequenced-site Ho/He.
#' * [het_between_pops()]: do populations differ in heterozygosity, with the
#'   individual (not the locus) as the unit of replication.
#' * [differentiation_stats()]: Weir & Cockerham FST, Jost's D and Weir &
#'   Goudet's beta.
#' * [pi_allsites()]: nucleotide diversity and divergence from an all-sites VCF.
#' * [individual_inbreeding()] and [identity_disequilibrium()]: inbreeding of
#'   each individual, and whether individuals differ in inbreeding.
#' * [kinship_check()] and [hwe_test()]: a relatedness screen and an exact
#'   Hardy-Weinberg test (a report, never a filter).
#' * [read_stacks_vcf()], [read_popmap()], the `filter_*()` functions and the
#'   `write_*()` functions: read, filter and export.
#'
#' Every analysis function takes `vcf` (a file path or the object from
#' [read_stacks_vcf()]) and, where populations matter, `popmap` (a file path
#' or the list from [read_popmap()]). Printing a result shows its main
#' tables; `summary()` shows the full report. Standard errors and confidence
#' intervals resample whole RAD loci. The command-line scripts are in
#' `system.file("scripts", package = "RADdiversity")`; the statistical
#' rationale is in `vignette("rationale")`.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats cor cor.test median p.adjust quantile rnorm runif setNames sd t.test var wilcox.test
#' @importFrom utils head read.delim write.table
## usethis namespace: end
NULL
