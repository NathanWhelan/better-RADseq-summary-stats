# RADdiversity 0.1.0 (development version)

First release as an R package. Notes for anyone who used the earlier
stand-alone scripts:

## Changes you will notice

* `diversity_stats()`, `het_between_pops()` and `differentiation_stats()`
  return result objects; printing one shows the familiar report. Files are
  written only when `outdir` is given. The command-line scripts still print the
  report and write to the working directory.
* The command-line scripts are installed with the package, in
  `system.file("scripts", package = "RADdiversity")`, instead of living at the
  top of the repository. `het_between_pops.R` also accepts `--min-call=` and
  `--outdir=`.
* `het_between_pops()` output files now include the input's stem
  (`individual_heterozygosity.snps.tsv`, `het_between_pops_tests.snps.tsv`), so
  the SNP-VCF and haplotype-VCF runs no longer overwrite each other.
* `read_haps_vcf()` is now `read_stacks_vcf()`; the old name remains as an
  alias for this release.
* hierfstat is optional. Every statistic is computed by the package itself;
  `hierfstat_check = TRUE` compares the results against hierfstat.
* `differentiation_stats()` also returns the per-pair table (`$pairwise`).

None of these change any reported number: on the reference datasets the
results agree with the previous version to within 5e-13.

## New

* `individual_inbreeding()`, and a test on individual F in
  `het_between_pops()` (`pairwise_F_tests`): whether populations differ in
  inbreeding, next to whether they differ in diversity.
* `identity_disequilibrium()`: g2 per population, also printed by
  `het_between_pops()`. `diversity_stats(se_individuals = TRUE)` adds
  delete-one-individual jackknife SEs, for when individuals differ in
  inbreeding (`inst/sims/coverage_se.R` shows why locus-based intervals are
  then too narrow).
* `pi_allsites()`: pi, dxy, Da and per-individual heterozygosity per
  sequenced site from an all-sites VCF (Stacks `populations --vcf-all`),
  handling missing data site by site.
* `read_stacks_vcf(locus_from = ...)`: assigns records to RAD loci for VCFs
  from other pipelines and for Stacks builds that leave the ID column empty.
* Three vignettes: `workflow` (one analysis start to finish), `rationale`
  (most of the former README) and `reviewer-faq`.

## Bug fixes

* A popmap with all-numeric sample IDs such as `001` was read as numbers, so
  the IDs no longer matched the VCF.
* `--sites=N` on the command line failed with "file not found".
* `hwe_test()`: a Monte Carlo p-value could be exactly 0 (now
  `(ge + 1) / (n + 1)`), and draws that tied the observed table were not
  counted as "at least as extreme".
* A VCF with a header but no records, and a VCF in which every genotype is
  missing, stopped with misleading error messages.
* `kinship_check()` returned `NaN` instead of `NA` for a pair with no
  heterozygous locus.

## Speed

* The jackknife and the bootstrap over RAD loci are computed from per-locus
  sums, and hierfstat is no longer called by default. On 5,000 RAD loci and 60
  individuals (`nboot = 1000`): `diversity_stats()` 43.5 s -> 1.8 s,
  `differentiation_stats()` 40.0 s -> 13.5 s.
* `hwe_test()` draws its Monte Carlo samples in vectorised batches and stops
  early once a p-value is clearly not small (Besag & Clifford 1991).
