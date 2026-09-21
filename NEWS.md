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

* Printing a result shows its main tables. `summary(result)` is now a short
  view, and `summary(result, details = TRUE)` is the full report that
  `summary()` used to print. This holds for `diversity_stats()`,
  `differentiation_stats()`, `het_between_pops()` and `pi_allsites()`:
  * The short view leads with the results as `estimate (SE)`, then what to
    take from the run, then a list of checks marked `ok`, `info` or `look`.
    Every column of `estimate (SE)` cells has "(SE)" in its header, and a line
    under the table says which standard error it is and what each column means.
  * Only the standard error the package recommends is shown, one per
    statistic: `_se_combined` for Ho, He, FIS and the per-site Ho and He; `_se`
    for Ar, privAr, FST, D, pi, `dxy` and `da_nc`; the combined test for
    `het_between_pops()`. The standard errors of FST, D, pi and dxy count RAD
    loci only and have not been checked by simulation; the summaries say so. The
    other uncertainty columns are in the tables and in `details = TRUE`.
  * The full report gives each statistic its own small table, so no row of
    numbers is left without its population name (the wide tables used to wrap
    into three blocks at the default console width). Its explanations are
    rewritten in shorter sentences, and "TAKE FROM THIS RUN" comes first.
  * The command-line scripts print the short summary; `--details` prints the
    full report. They still write to the working directory.
* `diversity_stats(se_individuals = )` now defaults to `TRUE`, because
  `_se_combined` is the recommended standard error for Ho, He and FIS. It adds
  about 30% to the run time (300,000 records and 60 individuals: 3.9 seconds
  became 5.2). `se_individuals = FALSE` skips it. On populations with very few
  individuals it can warn that the individual standard errors are probably too
  large; the help-page examples that use the tiny toy data set it to `FALSE`. In
  the shell script `diversity_stats.R`, `--se-individuals` (never in a release)
  is replaced by `--no-se-individuals`, which turns them off.
* The per-site values `Ho_autosomal` and `He_autosomal` have the same
  uncertainty columns as Ho and He (`_se`, `_se_ind`, `_se_combined`, `_lo`,
  `_hi`). A per-site value is Ho or He times `variant_records / sites_used`, so
  each is the matching Ho or He column times that multiplier. This treats the
  multiplier as fixed: it leaves out uncertainty in how many SNPs there are per
  sequenced site.
* `summary()` now returns the tables its short view prints, in `summary(x)$tables`
  (`results` and `per_site` for `diversity_stats()`; `results` for
  `differentiation_stats()`; `results`, `between` and `individual_het` for
  `pi_allsites()`; `heterozygosity`, `inbreeding` and `populations` for
  `het_between_pops()`). They hold every population or pair, where the printed
  view cuts long lists of pairs. A table for a paper is exactly these tables:
  the README ("Make a table for a manuscript"), the workflow vignette ("A table
  for your manuscript") and the help pages ("Tables for a paper") show how to
  build it, what stays out of it and why, and when to add more.
* `print()` of a diversity or pi result repeats the population name on every
  block of a wide table, as the summaries do. The per-site table has 15 columns
  now, and R's own wrapping left its later blocks without a population name.
* The short `diversity_stats()` summary says why a population has no per-site
  value (a `sites` value smaller than its variant records is ignored for that
  population), instead of showing `NA (NA)` without a reason.
* The heterozygosity missingness check no longer says "every individual was
  called at every locus" when an individual with no genotypes is present. It
  says call rate does not vary among the individuals checked.
* Small fixes to the wording of the reports: "1 copies" is now "1 copy",
  "absent from all 1 others" is now "absent from the other population", and the
  pointer to "the per-population lines printed while it ran" says the lines are
  shown unless `verbose = FALSE`.
* Result tables keep full precision; printing and the TSV files round them
  as before. Tables that used to be hidden in the `"report"` attribute are now
  list elements: `het_between_pops()` returns `population_summary`, `g2` and
  `missingness_confound`, and `diversity_stats()`
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
  population is now scaled by its own count of variant records in the data,
  shown in the new `variant_records` column. A warning is given when records were
  removed in R before the per-site conversion. Datasets where every
  population is typed at every record are unchanged.
* `differentiation_stats()`: F<sub>ST</sub> used records typed in only one
  of the compared populations, adding their within-population variance to
  the denominator and pulling F<sub>ST</sub> toward 0 (0.102 → 0.076 with one
  population absent from 30% of records). Such records are now skipped, as
  in VCFtools and Stacks' own pairwise F<sub>ST</sub>; `hierfstat::wc()` keeps
  them, so the two now differ on such data. The count skipped is in
  `settings$fst_records_skipped`. Pairwise D was already computed this way;
  for the global D see the next section.
* `kinship_check(method = "beta")` computed beta against the whole pooled
  sample, so with differentiated populations most unrelated same-population
  pairs were flagged (277 of 378 in a simulation at F<sub>ST</sub> = 0.1). The
  new `popmap` argument computes beta within each population; without it a
  message explains the problem. KING is unaffected.
* Recommended Stacks flags: `--min-mac 3` is removed from the diversity run.
  Stacks turns a SNP failing it into a fixed site, removing about 4/(N − 1)
  of π for N gene copies (21% at 10 diploids, 10% at 20).

## Standard errors, per-site values and other changes from the last review

Found by checking earlier advice against simulations with a known truth. On
the shipped example data every previously reported number is unchanged; the
changes are new columns, removed columns, and what the documentation tells
you to report.

* **Which standard error to report.** A population's value is uncertain
  because loci are a sample of the genome and because individuals are a
  sample of the population. The package used to advise the jackknife over
  individuals *alone* when g2 > 0. The simulation behind that advice reused
  the same loci in every repeat, so it never included the uncertainty from
  which loci were typed. With new loci and individuals in every repeat
  (`inst/sims/uncertainty_sources.R`), 95% intervals contained the true value:
  * from the locus SE alone: He 93–97%, but Ho and FIS 30–65% in most
    settings once individuals differed in inbreeding;
  * from the individual SE alone: Ho and FIS 75–94%, but He 47–86%;
  * from the combined SE, `sqrt(_se^2 + _se_ind^2)`: Ho, He and FIS 90–99%.
  `diversity_stats()` now adds `Ho_se_combined`,
  `He_se_combined` and `Fis_se_combined`; report those. A bootstrap over
  individuals was also checked and is not a substitute: its He intervals
  contained the truth 6–8% of the time. How each SE is calculated is
  explained in `?diversity_stats` ("How standard errors are calculated") and
  `vignette("rationale")`, section 4.
* **Ar and privAr no longer have individual SEs** (`Ar_se_ind` and
  `privAr_se_ind` are removed, with their warnings). Their locus-based SE and
  bootstrap interval covered the truth 88–96% of the time; their individual
  SEs were too wide for Ar (1.4–1.5 times the true SE) and far too small for
  privAr.
* **The individual jackknife uses the same records in every replicate.** A
  record with exactly `min_n` genotyped individuals used to drop out whenever
  one of them was left out, which inflated the SEs: in a one-off simulation with 10 individuals per
  population, 30% missing genotypes and `min_n = 5`, `He_se_ind` was 1.64
  times its true value (1.35 after the change) and `Ho_se_ind` 1.23 (1.01).
  Results at the default `min_n = 2` are unchanged. A new warning flags
  populations where more than 5% of records have only 2 genotyped
  individuals, the one case this cannot fix (He `_se_ind` 1.32 times its true
  value at 6% such records, 1.10 at none).
* **Per-site values (`Ho_autosomal`, `He_autosomal`) after filtering.** When
  `sites` was `populations.sumstats_summary.tsv`, He was multiplied by the
  file's `Variant_Sites`, which also counts SNPs removed after Stacks. That
  treated each removed SNP as having the average He of those kept; after
  `filter_mac(min_mac = 3)` it inflated per-site He by 41% in a simulation
  with 20 diploids. He is now always multiplied by the SNP records in the
  data, so removed SNPs add nothing (the same loss Stacks' own `--min-mac`
  causes). Unfiltered data give the same result as before. New warnings say
  when records or whole RAD loci were removed after reading (Stacks' `Sites`
  still counts the sites of removed loci) and when the VCF's SNP count does not
  match the Stacks file (for example a VCF thinned to one SNP per locus).
  `?diversity_stats` now explains what `sites` must count.
* **`het_between_pops()` no longer reports the overdispersion factor.** It
  equals 1 + g2 × (mean heterozygosity)² ÷ (binomial variance), so it only
  rescaled g2 by the number of loci, and missing genotypes inflated it (1.15
  at 10% missing when individuals did not differ). g2, which handles missing
  data, remains.
* Haplotype allelic richness and private allelic richness depend on the
  number of SNPs per RAD locus. As SNPs per locus rose from 1.3 to 3.3 in a
  simulation (`inst/sims/vignette_sims.R`, part `snp_density`), haplotype Ar
  rose by 51–70% and the ratio between two populations grew from 1.12 to
  1.26. The printed result and the documentation now say to compare them only
  among populations analyzed together (same run, filters and g), and that
  which population is richer is reliable but the size of the difference is
  not.
* `write_structure()` stops on sample names containing white space, and
  `write_genepop()` on sample names containing a comma, which would corrupt
  those files.
* Symbolic ALT alleles (`*`, `<*>`, `<NON_REF>`) no longer make a SNP VCF
  count as haplotype data; `read_stacks_vcf()` counts records that have them.
* `inst/sims/uncertainty_sources.R` replaces `inst/sims/coverage_se.R`,
  `inst/sims/jackknife_richness.R` and the `boot_modes` part of
  `inst/sims/vignette_sims.R`.

## Production-readiness review

None of these change a reported number on the shipped example data.

* **R 4.0.0 or later is required** (was R 3.5). Before R 4.0, `data.frame()`
  turned text into factors, and `het_between_pops()` would then have given
  individuals the wrong `F` and `kinship_check()` the wrong population labels,
  without any error.
* **The data object records its filters.** Every `filter_*()` function adds
  a row to `H$filter_log` (filter, settings, records and RAD loci removed,
  calls masked, samples removed). `print(H)` lists the steps, and
  `summary()` of a `diversity_stats()` result shows them, ready for a methods
  section.
* **No false per-site warning after a MAC or MAF filter.** After
  `filter_mac(H, 3)` on the example data, `diversity_stats(sites = ...)` warned
  that 48 "whole RAD loci were removed" and that per-site values were too low.
  Those loci only lost their SNPs to the filter; their sites were still
  sequenced, as with Stacks' `--min-mac`, so Stacks' `Sites` is right.
  Following the old advice would have inflated `He_autosomal`. The loci
  warning now fires only for loci removed by other filters (it names them),
  and gives both possibilities for loci removed outside a `filter_*()`
  function. The records-removed warning names the filters too, and after
  `filter_thin_one_snp()` says to use the unthinned data.
* **Populations that cannot be analyzed stop with a message** instead of
  giving a row of `NA`: in `diversity_stats()`, a population with no record at
  `min_n` (the message gives how far to lower it); in `individual_inbreeding()`
  and `identity_disequilibrium()`, a population with fewer than 2 individuals
  or no locus at `min_call`, as in `het_between_pops()`. These two functions
  and `hwe_test()` also check that the popmap file exists before reading the
  VCF.
* **Popmap samples missing from the VCF are named** in a message (they used
  to be dropped silently, so a typo made a population smaller). When no name
  matches, the error shows a few names from both files.
* `kinship_check(outdir = )` writes `kinship_pairwise.king.tsv` or
  `kinship_pairwise.beta.tsv`, so runs with the two methods no longer
  overwrite each other.
* `print()` of a `differentiation_stats()` result says why beta is missing
  (`beta = FALSE`, hierfstat not installed, more than 99 alleles, or
  `pairwise.betas()` failed); it used to blame a missing hierfstat every time.
  The reason is in `settings$beta_status`.
* A haplotype VCF given `sites` is reported as needing the SNP VCF, whatever
  the value of `sites`.
* `?individual_inbreeding`: the note on failed libraries had ended up inside
  the description of `min_call`.

## Corrections from the final pre-release review

Each was confirmed by simulation or a constructed example before the change.
On the shipped example data, the only reported number that changes is
`fis_by_call_rate`'s `Fis_se`; Ho, He, FIS, Ar, privAr, F<sub>ST</sub>, D and
π are identical.

* `differentiation_stats()`: **the global Jost's D now uses only records
  genotyped in every population.** A record missing some populations measures
  differentiation among a different set of populations, and averaging it in
  biased the global D low: 0.875 instead of 1 for three completely
  differentiated populations with one untyped at half the records. The
  number of records used is in the new `D_records` column of `$global`, and
  `print()`, `summary()` and the progress messages state the rule. Pairwise D
  (which already required both populations) and datasets where every
  population is typed at every record are unchanged. If no record is typed
  in every population, the global D is `NA` with a warning.
* `differentiation_stats()`: F<sub>ST</sub> and the global F<sub>IS</sub> also
  skip records where every typed population has a single genotyped
  individual (no variance among individuals can be estimated there, but its
  heterozygosity was added to the denominator). Counted in
  `settings$fst_records_skipped`.
* `diversity_stats()`: `fis_by_call_rate`'s `Fis_se` counted loci with no
  record in the call-rate group as jackknife blocks, inflating the SE
  (by 12% for a group with 5 loci). Only the group's own loci are used now.
* Haplotype detection: one multi-allelic SNP (e.g. `A` -> `C,T`, as GATK or
  freebayes emit) made a whole SNP VCF count as haplotype data, silently
  dropping `Ho_autosomal`/`He_autosomal`. A file is now haplotype data only
  when some allele is longer than one base. The toy `small.haps.vcf` now has
  multi-base haplotype alleles like real Stacks output; its genotypes, and
  every number computed from it, are unchanged.
* `diversity_stats()`'s summary of what the data look filtered at is now
  computed over the popmap's individuals only, like every statistic.
* `read_popmap()` follows Stacks' popmap rules: an optional third (group)
  column is accepted, and a line with one column, more than three, or spaces
  instead of a TAB stops with the line number and text.
* `hwe_test(stop_after = )` must be a whole number or `Inf`.
* `differentiation_stats()`: **Weir & Goudet's beta could be wrong on
  haplotype data with 10 or more alleles at a locus.** Genotypes were passed
  to hierfstat with 3 digits per allele (alleles 1 and 12 as `1012`), and
  hierfstat, which guesses the number of digits from the data, read them with
  2 digits (`1012` as alleles 10 and 12) whenever every number was below
  10000. That happens when no individual carries two alleles numbered 10 or
  more. In a constructed example beta was 0.049 instead of 0.065, and
  `hierfstat_check = TRUE` falsely reported that this package's F<sub>ST</sub>
  (which was correct) disagreed with `hierfstat::wc()`. Alleles are now
  renumbered within each record and written with 2 digits, which hierfstat
  always reads correctly. The `hierfstat_check` comparisons of
  `diversity_stats()` had the same problem. Data with fewer than 10 alleles
  per locus, including the shipped examples, are unaffected.
* `het_between_pops()` **no longer removes individuals.** It removed any
  individual genotyped at fewer than `min_loci` of the loci that pass a
  POOLED call-rate filter. When populations are absent from many records
  (Stacks' `-r` with several populations), few loci pass that filter, and in
  one example every individual was removed although each population had 180
  loci of its own. A failed library belongs to assembly and filtering, and
  one still in the popmap was probably kept on purpose. So each individual is
  now kept, and one genotyped at fewer than `min_loci` records (when the rest
  of its population has that many), or at under half of the records where
  the rest of its population is genotyped, is named
  in a warning that says what keeping it does (also in
  `settings$flagged_individuals`, `print()` and `summary()`).
  `individual_inbreeding()` and `identity_disequilibrium()` give the same
  warning. An individual not genotyped at the loci a table or test uses is
  left out of that table or test only. Its row stays in
  `individual_heterozygosity` with `NA`.
* `het_between_pops()` stopped with "missing value where TRUE/FALSE needed"
  for a population of one individual or a population with no locus at
  `min_call`, and gave `NaN` in `population_summary` for an individual with no
  call at its population's loci. Now it stops with a message naming the
  population and its highest call rate, and leaves such individuals out of
  the summary (`n` counts individuals with a value).
* `differentiation_stats()` and `pi_allsites()`: two pairs whose population
  names join to the same label (e.g. `a` + `b__c` and `a__b` + `c`) were
  both given one pair's values. Statistics are now looked up by pair number.
* `diversity_stats(sites = )`: a population given two values used the first
  silently; it now stops. The note about names that match no population
  respects `verbose = FALSE`.
* A popmap given as a data frame, or as a list of factors, stops with a
  message saying how to convert it. `read_sumstats_summary()` errors no
  longer print an internal function call.

## Scientific review before journal submission

These change reported numbers.

* **`het_between_pops()` has a new primary test.** Welch's *t* on one value
  per individual counts how much individuals vary but treats the loci as
  fixed. When populations are differentiated, each one's heterozygosity at
  the typed loci differs by chance from its genome-wide value, and Welch's *t*
  treats that as real. In simulations with identical true heterozygosity and
  each population drawing its own allele frequencies
  (`inst/sims/het_test_null.R`), Welch's *t* rejected up to 14% of the time
  at FST 0.05 and up to 30% at FST 0.2, at a nominal 5%, when individuals
  were alike in inbreeding. The new combined test adds a jackknife over RAD
  loci to Welch's variance, without counting genotype noise twice, and stayed
  between 3% and 7.5% in every setting. `pairwise_tests` and
  `pairwise_F_tests` gain `se_combined`, `df`, `p_combined` and
  `p_combined_BH`; `ci_lo` and `ci_hi` are now the combined interval (on the
  package's golden test data they widened slightly, e.g. -0.0190 to 0.0153
  became -0.0201 to 0.0164). `p_welch` and `p_wilcox` are unchanged and kept
  for comparison. For F, the locus part is usually close to 0, so the combined
  test is close to Welch's. The overall test for 3 or more populations is
  still Welch's ANOVA and Kruskal-Wallis, with a note that it uses
  individuals only.
* `het_between_pops_selftest()` now gives each population its own allele
  frequencies and reports the combined test, Welch's *t*, Wilcoxon and a locus
  bootstrap in three settings; its result gains an `fst` column.
* **Jost's D uses Nei & Chesser's Hs and Ht**, which subtract observed
  heterozygosity, as `hierfstat::basic.stats()` does for `Dest`. The
  2N/(2N - 1) correction used before (from mmod) assumes random mating, so D
  came out too high when individuals were inbred: 0.017, 29 standard errors
  above the true 0, for two samples of 8 from one population with F = 0.4. On
  the example data D is about 7% lower (haplotype VCF 0.0554 to 0.0516, SNP
  VCF 0.0407 to 0.0379). FST is unchanged.
* `pi_allsites()` adds `da_nc`, net divergence built on `pi_nc`, which
  inbreeding does not push up (`da` is kept).
* `diversity_stats()` notes when one population's Ar rests on under 90% of
  the records of another (populations are then compared over partly
  different loci).
* Documentation: `het_between_pops()` tests **observed** heterozygosity,
  He × (1 − F), which mixes diversity and inbreeding (it had been described
  as a test of diversity); KING-robust reads low when individuals are inbred
  (about −F/(1 − F) for unrelated pairs); sex-linked loci and linkage within
  RAD loci (for g2) are discussed; a sentence that described Petit & Pons
  (1998) wrongly was removed; and the rationale's claim that loci with fewer typed
  individuals "contribute less" to FIS was corrected (every record counts
  equally).
* `fis_by_call_rate` is now checked by simulation (`inst/sims/vignette_sims.R`,
  part `dropout`): false alarms 2–4%, strong dropout always found, weak
  dropout (null allele at 10% of loci, frequency 0.1) found 23% of the time.

## New

* `filter_samples(H, popmap)` removes the samples that are not in the popmap.
  The `filter_*()` functions and `write_*()` functions use every sample in
  `H`, including ones left out of the popmap (an outgroup, say); running
  `filter_samples()` first bases filtering, and exported files, on the
  analyzed individuals only.
* `?read_stacks_vcf` gains "What Stacks writes": the CHROM, POS and ID
  columns of Stacks 2.68's SNP and haplotype VCFs, de novo and
  reference-aligned, and how each is grouped into RAD loci.

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
* `read_stacks_vcf()` records `n_records_read` and `n_loci_read`, which
  filters leave unchanged.
* `individual_inbreeding()`, and a test on individual F in
  `het_between_pops()` (`pairwise_F_tests`): whether populations differ in
  inbreeding, next to whether they differ in diversity.
* `identity_disequilibrium()`: g2 per population, also reported by
  `het_between_pops()`.
* `diversity_stats()` adds, for Ho, He and FIS, a
  jackknife SE over individuals (`_se_ind`) and `_se_combined` (on by default;
  `se_individuals = FALSE` skips them), which puts
  the uncertainty from which loci and which individuals were sampled
  together. See "Standard errors" below.
* `pi_allsites()`: pi, dxy, Da and per-individual heterozygosity per
  sequenced site from an all-sites VCF (Stacks `populations --vcf-all`),
  handling missing data site by site.
* `read_stacks_vcf(locus_from = ...)`: assigns records to RAD loci for VCFs
  from other pipelines and for Stacks builds that leave the ID column empty.
* `low_conf_alt_calls()`: the genotype calls `filter_low_conf_alt()` would
  flag.
* Two vignettes: `workflow` (one analysis start to finish) and `rationale`
  (most of the former README).

## Bug fixes

* `summary()` of `diversity_stats()` no longer says that Ar is "capped at 2
  because these are biallelic SNPs" (and to run the haplotype VCF) for a
  haplotype VCF. The note was triggered by a mean Ar of 2 or less, which is
  common on haplotype data with few SNPs per RAD locus. It now appears only
  for SNP VCFs.
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
