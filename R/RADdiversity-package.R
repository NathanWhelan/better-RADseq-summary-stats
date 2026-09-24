#' RADdiversity: Population Genetic Statistics from RAD-seq Data
#'
#' Summary statistics for RAD-seq data, designed around Stacks 2
#' `populations` output:
#'
#' * [diversity_stats()]: Ho, He (Nei & Chesser 1983), FIS as a ratio of sums,
#'   rarefied allelic richness and rarefied private alleles, per-sequenced-site
#'   Ho/He.
#' * [het_between_pops()]: do populations differ in heterozygosity? One value
#'   per individual, and a test that counts both individuals and RAD loci.
#' * [diversity_table()]: the SNP and haplotype results as one table for a
#'   manuscript, with a caption.
#' * [differentiation_stats()]: Weir & Cockerham FST, Jost's D and Weir &
#'   Goudet's beta.
#' * [isolation_by_distance()]: do populations farther apart differ more? A
#'   Mantel test of pairwise FST or D against distance, with distances read by
#'   [read_distances()].
#' * [pi_allsites()]: nucleotide diversity and divergence from an all-sites VCF.
#' * [individual_inbreeding()] and [identity_disequilibrium()]: inbreeding of
#'   each individual, and whether individuals differ in inbreeding.
#' * [kinship_check()] and [hwe_test()]: a relatedness screen and an exact
#'   Hardy-Weinberg test (a report, never a filter).
#' * [read_stacks_vcf()], [read_popmap()], the `filter_*()` functions (with
#'   [filter_samples()] to keep only the popmap's individuals), the
#'   `write_*()` functions and [as_genind()] / [as_genlight()]: read, filter
#'   and export.
#'
#' Every analysis function takes `vcf` (a file path or the object from
#' [read_stacks_vcf()]) and, where populations matter, `popmap` (a file path
#' or the list from [read_popmap()]). Printing a result shows its main
#' tables; `summary()` is a short summary and `summary(details = TRUE)` the full
#' report. Standard errors and confidence
#' intervals resample whole RAD loci. The command-line scripts are in
#' `system.file("scripts", package = "RADdiversity")`; the statistical
#' rationale is in `vignette("rationale")`.
#'
#' Genotypes are assumed to be diploid. Each function's help page lists the
#' methods it implements; please cite those alongside the package
#' (`citation("RADdiversity")`).
#'
#' @references
#' Nei, M. & Chesser, R.K. (1983) Estimation of fixation indices and gene
#' diversities. *Annals of Human Genetics* 47:253-259.
#'
#' Weir, B.S. & Cockerham, C.C. (1984) Estimating F-statistics for the
#' analysis of population structure. *Evolution* 38:1358-1370.
#'
#' Kalinowski, S.T. (2004) Counting alleles with rarefaction: private alleles
#' and hierarchical sampling designs. *Conservation Genetics* 5:539-543.
#'
#' Jost, L. (2008) GST and its relatives do not measure differentiation.
#' *Molecular Ecology* 17:4015-4026.
#'
#' Mantel, N. (1967) The detection of disease clustering and a generalized
#' regression approach. *Cancer Research* 27:209-220.
#'
#' Rousset, F. (1997) Genetic differentiation and estimation of gene flow from
#' F-statistics under isolation by distance. *Genetics* 145:1219-1228.
#'
#' Schmidt, T.L., Jasper, M.-E., Weeks, A.R. & Hoffmann, A.A. (2021) Unbiased
#' population heterozygosity estimates from genome-wide sequence data.
#' *Methods in Ecology and Evolution* 12:1888-1898.
#'
#' Van Dongen, S. (1995) How should we bootstrap allozyme data? *Heredity*
#' 74:445-447.
#'
#' Rochette, N.C., Rivera-Colon, A.G. & Catchen, J.M. (2019) Stacks 2:
#' analytical methods for paired-end sequencing improve RADseq-based
#' population genomics. *Molecular Ecology* 28:4737-4754.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats cor cor.test median p.adjust quantile rnorm runif setNames sd t.test var wilcox.test
#' @importFrom utils head write.table
## usethis namespace: end
NULL
