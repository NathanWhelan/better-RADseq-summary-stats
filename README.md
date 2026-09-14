# RADdiversity

This package provides several methods for calculating population-genetic
statistics from RAD-seq data. It was built with output from Stacks 2 in mind,
but VCFs from other assembly pipelines work too. Genotypes are assumed to be
diploid.

The motivation for this package was a desire to calculate statistics more
robustly than Stacks does. For instance, Stacks calculates F<sub>IS</sub> as
a mean of per-site ratios rather than a ratio of sums (1 −
ΣH<sub>o</sub>/ΣH<sub>e</sub>), so every site counts equally however little
information it carries, including sites typed in only a few individuals.

## Why use it

Three things routine RAD-seq summaries tend to get wrong:

1. **F<sub>IS</sub>.** Stacks averages per-site F<sub>IS</sub> ratios and
   divides by an expected heterozygosity that assumes F<sub>IS</sub> = 0.
   Here, H<sub>e</sub> is Nei & Chesser's (1983) estimator, unbiased at any
   F<sub>IS</sub>, and F<sub>IS</sub> = 1 − ΣH<sub>o</sub>/ΣH<sub>e</sub>
   over loci.
2. **Comparing populations.** Bootstrapping or testing over loci treats loci
   as independent replicates of a population mean. They are repeated
   measurements on the same animals: in simulations of two populations with
   identical diversity, a locus bootstrap rejected the null about half the
   time when individuals differed moderately in inbreeding, and about 80% of
   the time when they differed more, at a nominal 5%
   (`het_between_pops_selftest()`). Here, populations are compared on one
   number per individual -- heterozygosity, or individual F.
3. **Uncertainty.** SNPs on one RAD tag are linked, so standard errors and
   intervals resample whole RAD loci, never SNP rows. And when individuals
   differ in inbreeding, even that is too optimistic for population-level
   inference: the package measures this (g2) and also gives standard errors
   over individuals.

## Install

```r
install.packages("RADdiversity")        # once the package is on CRAN
# the development version:
# install.packages("remotes")
remotes::install_github("NathanWhelan/better-RADseq-summary-stats")
```

R ≥ 4.0; base R is enough. Optional: `hierfstat` (Weir & Goudet's beta in
`differentiation_stats()`, and `hierfstat_check = TRUE` cross-checks).

## Quick start

On the simulated example dataset that ships with the package (two
populations of 15 and 10 individuals, 1,000 RAD loci):

```r
library(RADdiversity)
haps     <- system.file("extdata", "example.haps.vcf.gz", package = "RADdiversity")
snps     <- system.file("extdata", "example.snps.vcf.gz", package = "RADdiversity")
popmap   <- system.file("extdata", "example_popmap.tsv", package = "RADdiversity")
sumstats <- system.file("extdata", "example.sumstats_summary.tsv", package = "RADdiversity")

set.seed(2024)                              # reproducible bootstrap intervals
div_haps <- diversity_stats(haps, popmap, g = 16)                    # FIS, Ar, privAr
div_snps <- diversity_stats(snps, popmap, g = 16, sites = sumstats)  # Ho, He, per-site values
div_haps                                    # the main tables
summary(div_haps)                           # the full report, with notes on interpretation
het <- het_between_pops(snps, popmap)       # do populations differ?
het$pairwise_tests
```

Every result is a list of tables (`div_haps$per_population`,
`div_haps$richness`, ...) kept at full precision; printing rounds them. Add
`outdir = "results"` to also write the tables as TSV files, and
`verbose = FALSE` to silence progress messages. Every analysis function
accepts the data as a file path or as the object from `read_stacks_vcf()`,
and the popmap as a file path or as the list from `read_popmap()`. On your own
data, `g` is the rarefaction size in **gene copies** (10 diploids = 20), at
most twice your smallest population.

From the shell, the same analyses run with the scripts installed alongside
the package:

```bash
RD=$(Rscript -e 'cat(system.file("scripts", package = "RADdiversity"))')
Rscript $RD/diversity_stats.R out/populations.haps.vcf popmap.tsv --g=20
Rscript $RD/diversity_stats.R out/populations.snps.vcf popmap.tsv --g=20 \
        --sites=out/populations.sumstats_summary.tsv
Rscript $RD/het_between_pops.R out/populations.snps.vcf popmap.tsv --min-call=0.9 --outdir=het_out
```

Run any script with no arguments for its options. `vignette("workflow",
package = "RADdiversity")` walks through a whole analysis.

## Before you start: STACKS `populations` flags to consider

```bash
populations --in-path ./stacks_out --popmap popmap.tsv -O ./out \
            -p 2 -r 0.8 --min-gt-depth 6 --max-obs-het 0.70 \
            --vcf --vcf-all -t 8
```

* `--min-gt-depth 6` (Stacks ≥ 2.67): a heterozygote seen in only a few
  reads is sometimes called a confident homozygote. This does not mean the
  allele call is wrong, but it could be. `--min-gt-depth` makes such
  genotypes missing data instead. Values of 6 or 10 could be justified. In the
  author's tests, not filtering by genotype depth made little difference to
  the final statistics, so this is a reasonable precaution rather than a
  requirement.
* Use both `-r` and `-p`, not a global `-R`. A global `-R` can leave a locus
  mostly missing in one population while it passes on the others' coverage.
  A missing-data filter applied within each population (`-r`, with `-p` set
  to your number of populations) ensures each locus that passes is well
  genotyped in every population. This matters most when sampling is uneven:
  `-R` lets the larger population's coverage carry the smaller one, which
  then carries more missing data -- an artificial difference between exactly
  the two samples you mean to compare.
* Minor allele count or frequency filters (`--min-mac`, `--min-maf`) are
  common, and fine for many analyses, but they remove rare variants, which
  are real diversity: under a neutral site frequency spectrum, sites with a
  minor allele count of 2 or less hold about 21% of π at 10 diploids and 10%
  at 20. If you use one, report the threshold and compare diversity values
  only with data filtered the same way; if you can, a separate run without it
  for diversity statistics avoids the issue.
* `--vcf-all` (Stacks ≥ 2.62) writes an all-sites VCF, which `pi_allsites()`
  uses to compute nucleotide diversity directly.
* Per-site values read `Sites` from `populations.sumstats_summary.tsv`. That
  count fits only data not filtered after `populations`; if you remove whole
  RAD loci later, give `sites` for the loci you kept (`?diversity_stats`).

Data filtered elsewhere, or with other settings, work too:
`vignette("workflow")` has a table of what common upstream filters change and
what to report.

## What each function answers

| question | function |
|---|---|
| Ho, He, F<sub>IS</sub>, % polymorphic, allelic and private allelic richness | `diversity_stats()` |
| nucleotide diversity π, d<sub>xy</sub> and per-individual heterozygosity per sequenced site | `pi_allsites()` (all-sites VCF), or `diversity_stats(sites = ...)` (approximation) |
| do populations differ in diversity? in inbreeding? | `het_between_pops()` |
| each individual's inbreeding coefficient | `individual_inbreeding()` |
| do individuals differ in inbreeding (g2)? | `identity_disequilibrium()` |
| F<sub>ST</sub>, Jost's D, Weir & Goudet's beta (global D uses only records typed in every population) | `differentiation_stats()` |
| are there close relatives in the sample? | `kinship_check()` |
| Hardy–Weinberg departures (a report, never a filter) | `hwe_test()` |
| read and filter | `read_stacks_vcf()`, `read_popmap()`, `filter_samples()`, `filter_call_rate()`, `filter_genotype_depth()`, `filter_maf()`, `filter_mac()`, `filter_max_het()`, `filter_thin_one_snp()`, `filter_low_conf_alt()`, `locus_allele_stats()`, `low_conf_alt_calls()`, `low_conf_alt_sensitivity()` |
| export | `write_vcf()`, `write_plink()`, `write_structure()`, `write_genepop()`, `write_fstat()`, `write_radpainter()` |
| read Stacks' `sumstats_summary.tsv` | `read_sumstats_summary()` |
| the individual estimators | `?estimators` |
| check the installation | `diversity_core_selftest()`, `het_between_pops_selftest()` |

## Which file for which statistic

Stacks writes one record per SNP (`populations.snps.vcf`) and one per RAD
locus (`populations.haps.vcf`). They answer different questions:

| statistic | file | why |
|---|---|---|
| H<sub>o</sub>, H<sub>e</sub>, per-site values (π) | `.snps` | per-site values are on a scale other studies can compare |
| F<sub>IS</sub> | `.haps` | a ratio, so the scale cancels; multi-allelic loci are more precise |
| allelic richness, private allelic richness | `.haps` | on biallelic SNPs richness can only be 1 or 2; depends on SNPs per locus, so compare only populations analyzed together (same run, filters and *g*): which is richer is reliable, how much richer is not |
| between-population tests | either -- say which | |

H<sub>o</sub> and H<sub>e</sub> are **not** comparable between the two files
(haplotype H<sub>o</sub> asks "is this tag heterozygous?", per-site
H<sub>o</sub> "is this site?"); F<sub>IS</sub> is. Each printed result says
which of its numbers to take.

## Which standard error to report

A population's value is uncertain for two reasons: you typed some of the
genome's RAD loci, and you caught some of the population's individuals.

* `_se`, `_lo`, `_hi` resample **RAD loci** and hold the individuals fixed.
* `_se_ind` (`diversity_stats(se_individuals = TRUE)`) is a jackknife over
  **individuals** that holds the loci fixed.
* `_se_combined` puts both together: `sqrt(_se^2 + _se_ind^2)`.

**Report `_se_combined` for H<sub>o</sub>, H<sub>e</sub> and F<sub>IS</sub>,
and `_se` (or `_lo`, `_hi`) for A<sub>r</sub> and privA<sub>r</sub>.** In
simulations with a known truth (`inst/sims/uncertainty_sources.R`; 2 × 15
individuals, 1,000 loci, new loci and individuals in every repeat), 95%
intervals contained the true value this often:

| | locus SE | individual SE | combined SE |
|---|---|---|---|
| H<sub>o</sub> | 41–96% | 88–91% | 93–98% |
| H<sub>e</sub> | 95–96% | 72–77% | 98–99% |
| F<sub>IS</sub> | 31–94% | 92–93% | 93–98% |

The low locus-SE values are for populations whose individuals differ in
inbreeding. `vignette("rationale")`, section 4, explains how each SE is
calculated and why A<sub>r</sub> and privA<sub>r</sub> use the locus SE.
Never compare populations by overlapping intervals -- use
`het_between_pops()`.

## Glossary

* **H<sub>o</sub>** -- observed heterozygosity: the fraction of typed
  individuals that are heterozygous, averaged over loci.
* **H<sub>e</sub> (H<sub>s</sub>, gene diversity)** -- the chance that two
  gene copies drawn from the population differ (Nei & Chesser 1983).
* **π (nucleotide diversity)** -- the same quantity averaged over every
  sequenced site, invariant ones included.
* **F<sub>IS</sub>** -- 1 − H<sub>o</sub>/H<sub>e</sub> for a population,
  summed over loci; **F** -- the same for one individual.
* **gene copies** -- a diploid carries two: 10 individuals are 20 gene copies.
* **RAD locus vs record** -- a locus is one RAD tag; a SNP VCF has one record
  per SNP, so several records can share a locus.
* **A<sub>r</sub>, privA<sub>r</sub>** -- allelic and private allelic
  richness rarefied to the same number of gene copies (Kalinowski 2004).
* **g2** -- identity disequilibrium: how much more often an individual is
  heterozygous at two loci at once than two different individuals are; 0 when
  individuals do not differ in inbreeding (David et al. 2007).
* **jackknife / bootstrap** -- recompute a statistic with one unit left out /
  with units resampled; here the unit is a RAD locus (or, for `_se_ind`, an
  individual).
* **standard error (SE)** -- how much a number would change if the study were
  repeated; estimate ± 1.96 SE is a 95% confidence interval.

## What this package does not do

| left out | why |
|---|---|
| HWE filtering | a heterozygote deficit is the signal. `hwe_test()` reports departures; nothing removes loci on them. |
| null-allele correction | restriction-site null alleles cannot be removed by depth filtering; `diversity_stats()` reports `fis_by_call_rate` to detect them, and `vignette("reviewer-faq")` has text for the methods. |
| paralog detection beyond an excess-heterozygosity screen | `filter_max_het()` is that screen; dedicated tools (e.g. HDplot; McKinney et al. 2017) go further. |
| N<sub>e</sub>, AMOVA, neutrality tests | different questions. |

## How it relates to other tools

* **Stacks** reports H<sub>o</sub>, π and F<sub>IS</sub> per site and per
  population. This package recomputes F<sub>IS</sub> (ratio of sums,
  Nei–Chesser H<sub>e</sub>), adds rarefied and private allelic richness,
  standard errors over RAD loci and individuals, and individual-level tests.
  `vignette("rationale")` maps every Stacks column onto its counterpart here.
* **hierfstat** (Goudet 2005) uses the same H<sub>s</sub> and Weir &
  Cockerham estimators; the package's own implementations are checked
  against it in the tests. Its bootstrap resamples loci as independent rows,
  ignoring linkage within RAD tags. One deliberate difference: at a record
  where only one of the compared populations is genotyped, `hierfstat::wc()`
  keeps the within-population variance, which pulls F<sub>ST</sub> toward 0;
  `differentiation_stats()` skips such records, as VCFtools and Stacks' own
  pairwise F<sub>ST</sub> do.
* **dartR** can correct heterozygosity for invariant sites (Schmidt et al.
  2021); **snpR** (Hemstrom & Jones 2023) computes many SNP statistics across
  categorical metadata.
* **pixy** (Korunes & Samuk 2021) computes π and d<sub>xy</sub> from
  all-sites VCFs; `pi_allsites()` uses the same estimator, with standard
  errors over RAD loci, and adds `pi_nc`, which leaves out the comparison
  between the two gene copies inside each individual and so stays unbiased
  when individuals are inbred (Nei & Chesser's H<sub>e</sub> per site).
* **inbreedR** (Stoffel et al. 2016) computes g2; `identity_disequilibrium()`
  estimates it from its definition, using exactly the locus pairs typed when
  data are missing.

## Citation

`citation("RADdiversity")`. Please also cite the methods you report -- the
references are in each function's help page and in `vignette("rationale")`.

## Further reading

* `vignette("workflow", package = "RADdiversity")` -- one analysis from start
  to finish.
* `vignette("rationale", package = "RADdiversity")` -- why these flags,
  estimators and tests; formulas; the Stacks column mapping.
* `vignette("reviewer-faq", package = "RADdiversity")` -- methods text to
  adapt, and answers to the usual reviewer objections.

## Key references

David, P., Pujol, B., Viard, F., Castella, V. & Goudet, J. (2007) Reliable
selfing rate estimates from imperfect population genetic data. *Molecular
Ecology* 16:2474–2487.

Goudet, J. (2005) HIERFSTAT, a package for R to compute and test hierarchical
F-statistics. *Molecular Ecology Notes* 5:184–186.

Hemstrom, W. & Jones, M. (2023) snpR: user friendly population genomics for
SNP data sets with categorical metadata. *Molecular Ecology Resources*
23:962–973.

Kalinowski, S.T. (2004) Counting alleles with rarefaction: private alleles
and hierarchical sampling designs. *Conservation Genetics* 5:539–543.

McKinney, G.J., Waples, R.K., Seeb, L.W. & Seeb, J.E. (2017) Paralogs are
revealed by proportion of heterozygotes and deviations in read ratios in
genotyping-by-sequencing data from natural populations. *Molecular Ecology
Resources* 17:656–669.

Korunes, K.L. & Samuk, K. (2021) pixy: Unbiased estimation of nucleotide
diversity and divergence in the presence of missing data. *Molecular Ecology
Resources* 21:1359–1368.

Nei, M. & Chesser, R.K. (1983) Estimation of fixation indices and gene
diversities. *Annals of Human Genetics* 47:253–259.

Rochette, N.C., Rivera-Colón, A.G. & Catchen, J.M. (2019) Stacks 2.
*Molecular Ecology* 28:4737–4754.

Schmidt, T.L., Jasper, M.-E., Weeks, A.R. & Hoffmann, A.A. (2021) Unbiased
population heterozygosity estimates from genome-wide sequence data. *Methods
in Ecology and Evolution* 12:1888–1898.

Stoffel, M.A., Esser, M., Kardos, M., Humble, E., Nichols, H., David, P. &
Hoffman, J.I. (2016) inbreedR: an R package for the analysis of inbreeding
based on genetic markers. *Methods in Ecology and Evolution* 7:1331–1339.

Van Dongen, S. (1995) How should we bootstrap allozyme data? *Heredity*
74:445–447.

Weir, B.S. & Cockerham, C.C. (1984) Estimating F-statistics for the analysis
of population structure. *Evolution* 38:1358–1370.
