# RADdiversity 0.1.0

First release as an R package. Notes for anyone who used the earlier
stand-alone scripts:

## Changes you will notice

* Every analysis function takes the same two first arguments: `vcf`, a file
  path or the object returned by `read_stacks_vcf()` (optionally filtered),
  and `popmap`, a file path or the list returned by `read_popmap()`.

  | function | old argument | new argument |
  |---|---|---|
  | `diversity_stats()`, `het_between_pops()`, `differentiation_stats()`, `pi_allsites()` | `vcf_file`, `popmap_f` | `vcf`, `popmap` |
  | `individual_inbreeding()`, `identity_disequilibrium()`, `hwe_test()` | `H`, `pops` | `vcf`, `popmap` |
  | `kinship_check()` | `vcf_file` | `vcf` |
  | `filter_call_rate()`, `write_plink()`, `write_structure()`, `write_genepop()`, `write_fstat()` | `pops` | `popmap` |
  | `filter_maf()`, `filter_mac()` | `stats` | `allele_stats` |
  | `read_stacks_vcf()` | `read_haps_vcf()` (old function name) | removed; use `read_stacks_vcf()` |

* Printing a result shows its main tables; `summary(result)` shows the full
  report that printing used to show. The command-line scripts still print the
  full report and write to the working directory.
* Result tables keep full precision; printing and the TSV files round them
  as before. Tables that used to be hidden in the `"report"` attribute are now
  list elements: `het_between_pops()` returns `population_summary`,
  `overdispersion`, `g2` and `missingness_confound`, and `diversity_stats()`
  returns `estimator_comparison` and `he_difference`. Run settings are in
  `$settings`.
* `seed` defaults to `NULL`: results use R's random-number stream, so
  `set.seed()` before a call makes it reproducible, as in base R. Passing a
  number still gives a reproducible result without changing your session's
  random numbers. The command-line scripts use `--seed=2024` by default.
* `verbose = FALSE` now silences every progress message.
* `filter_low_conf_alt()` returns the filtered data, so it chains with the
  other filters. The tables it used to return come from the new
  `low_conf_alt_calls()`.
* `read_stacks_vcf()` returns an object of class `raddiv_vcf`, which prints
  as a short summary instead of every matrix.
* `kinship_check()` returns its result visibly, and `threshold = NULL` lists
  every pair with a kinship value.
* Former fixed limits are now arguments: `het_between_pops(min_loci = 50,
  nboot_g2 = 1000)` and `kinship_check(min_shared_loci = 30,
  low_confidence_loci = 200)`. `nboot_g2` defaults to 1000, the same as
  `identity_disequilibrium(nboot = )`; 200 replicates gave a noisy 95%
  interval.
* `diversity_stats(sites = NULL)` is the new default meaning "no per-site
  conversion" (`0` still works).
* `write_fstat()` always writes the file itself, so the output no longer
  depends on whether hierfstat is installed.
* `diversity_core_selftest()` and `het_between_pops_selftest()` gain
  `verbose` and return their results as a data frame.
* The command-line scripts are installed with the package, in
  `system.file("scripts", package = "RADdiversity")`, instead of living at the
  top of the repository. `het_between_pops.R` also accepts `--min-call=` and
  `--outdir=`, and both analysis scripts accept `--seed=`.
* `het_between_pops()` output files include the input's stem
  (`individual_heterozygosity.snps.tsv`, `het_between_pops_tests.snps.tsv`), so
  the SNP-VCF and haplotype-VCF runs no longer overwrite each other.
* hierfstat is optional. Every statistic is computed by the package itself;
  `hierfstat_check = TRUE` compares the results against hierfstat.
* `differentiation_stats()` also returns the per-pair table (`$pairwise`).

None of these change any reported number: on the reference datasets the
results agree with the previous version (golden-value tests in
`tests/testthat/test-regression.R`).

## Scientific corrections (these change reported numbers)

Found in a review of every estimator against the literature and the Stacks
2.68 source, each confirmed by simulation before being changed.

* `diversity_stats()`: `Ho_autosomal`/`He_autosomal` scaled every population
  by the total number of VCF records. Stacks' `Sites` counts only the sites
  where that population has a genotyped individual, so a population absent
  from some records (Stacks' `-r` with 3 or more populations and a low `-p`)
  got an inflated value: 1.43× when absent from 30% of records. Each
  population is now scaled by its own count of variant records (the file's
  `Variant_Sites` when `sites` is `populations.sumstats_summary.tsv`), shown
  in the new `variant_records` column. A warning is given when records were
  removed in R before the per-site conversion. Datasets where every
  population is typed at every record are unchanged.
* `differentiation_stats()`: F<sub>ST</sub> used records typed in only one
  of the compared populations, adding their within-population variance to
  the denominator and pulling F<sub>ST</sub> toward 0 (0.102 → 0.076 with one
  population absent from 30% of records). Such records are now skipped, as
  in VCFtools and Stacks' own pairwise F<sub>ST</sub>; `hierfstat::wc()` keeps
  them, so the two now differ on such data. The count skipped is in
  `settings$fst_records_skipped`. Jost's D was already computed this way.
* `kinship_check(method = "beta")` computed beta against the whole pooled
  sample, so with differentiated populations most unrelated same-population
  pairs were flagged (277 of 378 in a simulation at F<sub>ST</sub> = 0.1). The
  new `popmap` argument computes beta within each population; without it a
  message explains the problem. KING is unaffected.
* Recommended Stacks flags: `--min-mac 3` is removed from the diversity run.
  Stacks turns a SNP failing it into a fixed site, removing about 4/(N − 1)
  of π for N gene copies (21% at 10 diploids, 10% at 20).

## New

* `het_between_pops()` returns `omnibus` for 3 or more populations: one
  overall test of whether any population differs in heterozygosity or F
  (Welch's one-way ANOVA and Kruskal-Wallis), to report before the pairwise
  tests. Written to `het_between_pops_omnibus.<stem>.tsv` with `outdir`.
* `diversity_stats()` returns `fis_by_call_rate`: FIS and He within groups of
  records by each population's call rate. FIS rising as call rate falls is
  the signature of null alleles or allele dropout (Gautier et al. 2013)
  rather than inbreeding.
* `diversity_stats()`'s `summary()` describes what the data look filtered
  at (the rarest allele count left, the lowest record call rates, the
  highest record Ho), for data filtered before they reached the package.
* `differentiation_stats(beta = )`: `beta = FALSE` skips
  `hierfstat::pairwise.betas()`, which takes about 95% of the run time on
  large datasets.
* `read_stacks_vcf()` reads the VCF in chunks (`chunk_lines`) and keeps
  depths and allele depths as integer matrices (`$depth`, `$ad_ref`,
  `$ad_alt`) instead of the genotype text, so `$fields` now holds only the
  nine fixed columns. Memory use on large VCFs drops several-fold (see
  `inst/sims/benchmark.R`); every result is unchanged.
* `filter_genotype_depth()`: masks genotype calls below a minimum (and above
  a maximum) read depth, whatever the genotype. `filter_low_conf_alt()` masks
  only calls carrying an ALT allele, which removes heterozygotes
  preferentially; its help page and message now say so.
* `pi_allsites()` adds `pi_nc`: pi among gene copies of different individuals
  only, equal at each site to Nei & Chesser's He and unbiased when individuals
  are inbred (the pixy estimator, `pi`, runs low by FIS/(2n − 1)).
* `hwe_test()` adds a per-locus `F` column (sign of any departure: deficit or
  excess) and a message when all samples are tested as one group.
* `read_stacks_vcf()` records `n_records_read`, which filters leave unchanged.
* `individual_inbreeding()`, and a test on individual F in
  `het_between_pops()` (`pairwise_F_tests`): whether populations differ in
  inbreeding, next to whether they differ in diversity.
* `identity_disequilibrium()`: g2 per population, also reported by
  `het_between_pops()`. `diversity_stats(se_individuals = TRUE)` adds
  delete-one-individual jackknife SEs, for when individuals differ in
  inbreeding (`inst/sims/coverage_se.R` shows why locus-based intervals are
  then too narrow).
* `pi_allsites()`: pi, dxy, Da and per-individual heterozygosity per
  sequenced site from an all-sites VCF (Stacks `populations --vcf-all`),
  handling missing data site by site.
* `read_stacks_vcf(locus_from = ...)`: assigns records to RAD loci for VCFs
  from other pipelines and for Stacks builds that leave the ID column empty.
* `low_conf_alt_calls()`: the genotype calls `filter_low_conf_alt()` would
  flag.
* Three vignettes: `workflow` (one analysis start to finish), `rationale`
  (most of the former README) and `reviewer-faq`.

## Bug fixes

* `kinship_check()`: at a biallelic record whose two observed alleles are not
  numbered 1 and 2 (a haplotype record declaring three alleles of which only
  the 2nd and 3rd are carried), KING missed opposite homozygotes and beta
  miscounted allele dosage, so unrelated individuals read as related (mean
  KING kinship 0.25 instead of 0.03 in a test). Observed alleles are now
  renumbered 1 and 2 first.
* `write_fstat()`: the header's third number is now the highest allele number
  used, as FSTAT requires, not the number of distinct alleles observed.
* `het_between_pops()`: the Wilcoxon p-value now uses the same rule on every
  R version (exact without ties and with fewer than 50 per group, otherwise
  the normal approximation). R 4.6.0 changed `wilcox.test()`'s default for
  tied values, which made `p_wilcox` depend on the R version.
* `write_plink()` stops with a clear message when a sample or population
  name contains a space, which would shift the `.ped` columns.
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
* Arguments of the wrong length or type (e.g. `g = c(10, 20)`, `g = 4.9`)
  stopped with "the condition has length > 1" or were silently truncated;
  they now stop with a message naming the argument.
* `write_genepop()` and `write_fstat()` without `popmap` stopped with
  R's generic "argument is missing" error instead of their own explanation.
* A VCF with no sample columns, a `populations.sumstats_summary.tsv` block
  with no population rows, and `filter_low_conf_alt()` on data without the raw
  VCF columns all stopped with unhelpful errors.
* The `diversity_stats()` hierfstat cross-check skipped a failed
  `basic.stats()` call silently; it now warns, as it does for
  `allelic.richness()`.

## Speed

* The jackknife and the bootstrap over RAD loci are computed from per-locus
  sums, allele counts are tabulated for all records at once, and hierfstat is
  no longer called by default. On 5,000 RAD loci and 60 individuals in 3
  populations (`nboot = 1000`), `diversity_stats()` runs in about 0.3 s.
* `hwe_test()` draws its Monte Carlo samples in vectorised batches and stops
  early once a p-value is clearly not small (Besag & Clifford 1991).
