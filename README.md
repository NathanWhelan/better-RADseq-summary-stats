# RADdiversity

## PLEASE NOTE: This was mostly written by Claude AI, but it has been checked and edited. If something seems off, please post an issue. If you notice something that's horribly wrong, I would love to hear feedback. --N. Whelan (Sept 8, 2026).

This package provides several methods for calculating population-genetic
statistics from RAD-seq data. It was built with output from Stacks 2 in mind,
but VCFs from other assembly pipelines work too. Genotypes are assumed to be
diploid.

The motivation for this package was a desire to calculate statistics more
robustly than Stacks does. For instance, Stacks calculates F<sub>IS</sub> as
a mean of per-site ratios rather than a ratio of sums (1 −
ΣH<sub>o</sub>/ΣH<sub>e</sub>), so every site counts equally regardless of how little
information it carries, including sites typed in only a few individuals.

## Why use this package

Three things routine RAD-seq summaries tend to get wrong:

1. **F<sub>IS</sub>.** Stacks averages per-site F<sub>IS</sub> ratios and
   divides by an expected heterozygosity that assumes F<sub>IS</sub> = 0.
   Here, H<sub>e</sub> is Nei & Chesser's (1983) estimator, unbiased at any
   F<sub>IS</sub>, and F<sub>IS</sub> = 1 − ΣH<sub>o</sub>/ΣH<sub>e</sub>
   over loci.
2. **Comparing populations.** Bootstrapping or testing over loci treats loci
   as independent replicates of a population mean. However, they are repeated
   measurements on the same animals. In simulations of two populations with
   identical diversity, a locus bootstrap rejected the null about half the
   time when individuals differed moderately in inbreeding, and about 80% of
   the time when they differed considerably in inbreeding
   (`het_between_pops_selftest()`). Here, populations are compared on one
   number per individual -- heterozygosity, or individual F.
3. **Uncertainty.** SNPs on one RAD tag are linked, so standard errors and
   intervals resample whole RAD loci, never SNP rows. And when individuals
   differ in inbreeding, even that is too optimistic for population-level
   inference. The package measures this with the g2 statistic and also gives
   standard errors over individuals.

## Install

```r
install.packages("remotes")
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
summary(div_haps)                           # short: each estimate (SE), what to take, checks
het <- het_between_pops(snps, popmap)       # do populations differ in heterozygosity?
summary(het)
dif <- differentiation_stats(haps, popmap)  # how different are they? FST and Jost's D
summary(dif)                                # on haplotype data, report D (or beta) with FST

## A table for your paper: the columns summary() shows, and no more
snp    <- summary(div_snps)$tables          # Ho, He, pct_poly, and per-site Ho, He
hap    <- summary(div_haps)$tables          # Fis, Ar, privAr
parts  <- list(snp$results, snp$per_site, hap$results[names(hap$results) != "n"])
table1 <- Reduce(function(a, b) merge(a, b, by = "population", sort = FALSE),
                 Filter(Negate(is.null), parts))
table1 <- table1[match(snp$results$population, table1$population), ]   # popmap order
table1
table2 <- summary(dif)$tables$results       # FST (SE), D (SE), beta
```

The standard errors to report are computed by default: `_se_combined` for Ho,
He and FIS (see "Which standard error to report" below). Every result is a
list of tables (`div_haps$per_population`, `div_haps$richness`, ...) kept at
full precision; printing rounds them. `verbose = FALSE` silences progress
messages. The table for your paper is exactly what the short `summary()`
shows, and no more; "Make a table for a manuscript" below explains what it
holds and what stays out. To read a `summary()` and save the tables, see
"Printing and saving results" below.

On a full dataset the slow parts are the bootstrap (`nboot`, default 10,000)
and the `beta` step of `differentiation_stats()` (it calls hierfstat). For a
fast first look, use `nboot = 1000` and `beta = FALSE`, and run again with the
defaults for the results you report. The standard errors over individuals
(`se_individuals`, on by default) add only about 30% to the run time, about a
second per 100,000 records with 60 individuals, so leave them on.

### Using your own files

The four `system.file()` lines are the only place the example files appear.
To analyze your own data, replace them with paths to your files.

| variable | your file | what it is |
|---|---|---|
| `snps` | `populations.snps.vcf` | Stacks output: one record per SNP. `.vcf.gz` works too. |
| `haps` | `populations.haps.vcf` | Stacks output: one record per RAD locus (haplotypes). `.vcf.gz` works too. |
| `popmap` | `popmap.tsv` | The popmap you gave Stacks. No header. Each line is the sample name, a TAB, then the population. Sample names must match the VCF's. |
| `sumstats` | `populations.sumstats_summary.tsv` | Stacks output. Optional: only needed for per-site values (`sites =`). |

```r
snps     <- "out/populations.snps.vcf"
haps     <- "out/populations.haps.vcf"
popmap   <- "popmap.tsv"
sumstats <- "out/populations.sumstats_summary.tsv"   # optional

g <- 20    # gene copies: at most 2 x the individuals in your smallest population
div_haps <- diversity_stats(haps, popmap, g = g)
div_snps <- diversity_stats(snps, popmap, g = g, sites = sumstats)
dif      <- differentiation_stats(haps, popmap)
```

`g` is the rarefaction size in **gene copies** (10 diploids = 20). Every
analysis function accepts file paths as they are. `pi_allsites()` takes only a
path (to an all-sites VCF). The others also accept the R objects from the
readers below, which lets you look at a file, filter it first, or set
`locus_from`:

```r
H_snps <- read_stacks_vcf(snps)                # a VCF -> R object
H_haps <- read_stacks_vcf(haps)                # same, for the haplotype VCF
pops   <- read_popmap(popmap, H_snps$samples)  # popmap -> sample names per population, checked against the VCF
ss     <- read_sumstats_summary(sumstats)      # Stacks' summary table -> R list
ss$all_positions[, c("population", "sites")]   # sequenced sites per population

div_haps <- diversity_stats(H_haps, pops, g = g)
div_snps <- diversity_stats(H_snps, pops, g = g,
                            sites = sumstats)   # `sites` still takes the path
het      <- het_between_pops(H_snps, pops)
dif      <- differentiation_stats(H_haps, pops)
```

VCF from another pipeline (ipyrad, dDocent, ...)? Use it as `snps` and skip
`haps`. If its RAD loci are not recognized, read it with
`read_stacks_vcf(snps, locus_from = "CHROM")` (see `?read_stacks_vcf`).

### Printing and saving results

**Look at the results.** Printing a result shows its main tables, every
column, rounded. `summary()` is short. It shows the main results as
`estimate (SE)`, what to take from this run, and a list of checks marked `ok`,
`info` or `look`. `summary(x, details = TRUE)` is the full report, with every
uncertainty column and a note on how to read each number. Both work for
`diversity_stats()`, `differentiation_stats()`, `het_between_pops()` and
`pi_allsites()`. Each table can also be pulled out by name:

```r
div_haps                          # main tables, every column
summary(div_haps)                 # short: each estimate (SE), what to take, checks
summary(div_haps, details = TRUE) # the full report
div_haps$per_population           # one table, at full precision
dif$global                        # FST, FIS and D over all populations
dif$pairwise                      # one row per pair of populations
```

In a short `summary()`, every column of `estimate (SE)` cells has "(SE)" in its
header, and the lines under the table say which standard error it is and what
each column means. Only the recommended standard error is shown for each
statistic; the other columns are in `details = TRUE` and in the tables. The
names that are not obvious:

* `Ar`, allelic richness: the number of different alleles a population shows
  at a locus when you look at `g` gene copies, averaged over loci. Every
  population is cut down to the same `g`, so a bigger sample does not look
  richer.
* `privAr`, private allelic richness: the same, counting only alleles found in
  this population and in no other.
* `Ar_n`, `privAr_n` (in the tables and in `details`): how many loci went into
  each. A locus counts for a population's `Ar` if that population has at least
  `g` gene copies typed there, and for `privAr` only if every population has,
  so `privAr_n` is never larger than `Ar_n`.

**Save the working tables.** Add `outdir = "results"` to a call and its tables
are written there as TSV files, rounded as printed. The folder is created if
needed.

| call | files written |
|---|---|
| `diversity_stats()` | `diversity_per_population.<stem>.tsv`, `diversity_richness.<stem>.tsv`, and `diversity_autosomal.<stem>.tsv` (only when `sites` gave per-site values) |
| `differentiation_stats()` | `differentiation_global.<stem>.tsv`, `differentiation_pairwise.<stem>.tsv` |
| `het_between_pops()` | `individual_heterozygosity.<stem>.tsv`, `het_between_pops_tests.<stem>.tsv`, `het_between_pops_F_tests.<stem>.tsv` |

`<stem>` is `haps` or `snps`, taken from the VCF's file name, so the two runs
can share one folder. If you pass an object from `read_stacks_vcf()`, give
`stem = "haps"` yourself. The other tables of a diversity result
(`estimator_comparison`, `he_difference`, `fis_by_call_rate`) are not written;
`summary(details = TRUE)` prints them.

**Make a table for a manuscript.** Use what the short `summary()` reports, and
nothing more. For each statistic it shows the one standard error the package
found most accurate, so the table already holds what you should report.
`summary(x)$tables` gives those tables as data frames, with the same
`estimate (SE)` text the summary prints and one row for every population or
pair (the printed summary cuts long lists of pairs):

```r
snp <- summary(div_snps)$tables   # results: Ho, He, pct_poly.  per_site: Ho, He per site (needs `sites`)
hap <- summary(div_haps)$tables   # results: Fis, Ar, privAr

## Ho and He come from the SNP VCF, Fis, Ar and privAr from the haplotype VCF.
## `n` is in both tables, so keep one copy. Join by population name, never by row order.
parts  <- list(snp$results, snp$per_site, hap$results[names(hap$results) != "n"])
table1 <- Reduce(function(a, b) merge(a, b, by = "population", sort = FALSE),
                 Filter(Negate(is.null), parts))          # per_site is NULL without `sites`
table1 <- table1[match(snp$results$population, table1$population), ]   # popmap order
table1

table2 <- summary(dif)$tables$results   # FST (SE), D (SE) and beta for each pair

dir.create("results", showWarnings = FALSE)
write.csv(table1, "results/table1_diversity.csv", row.names = FALSE)
write.csv(table2, "results/table2_differentiation.csv", row.names = FALSE)
```

CSV files open in Excel, and you can paste from there into Word. In R
Markdown, `knitr::kable(table1)` prints the table. Rename columns to suit the
journal (for example, `Ho_autosomal (SE)` to `Ho per site (SE)`).

The other functions work the same way. `summary(het)$tables` holds
`heterozygosity` and `inbreeding` (the difference with its SE, 95% interval, p
and Hedges' g for each pair) and `populations`. `summary(pi_res)$tables` holds
`results` (`pi_nc` and `pi`), `between` (`dxy` and `da_nc`) and
`individual_het`.

**What the table holds.** Population, `n`, `Ho (SE)`, `He (SE)`, `pct_poly (%)`,
`Ho_autosomal (SE)` and `He_autosomal (SE)` (if you gave `sites`), `Fis (SE)`,
`Ar (SE)` and `privAr (SE)`. Every column that says "(SE)" holds the estimate
and its standard error. The lines under the table in `summary()` say which
standard error each one is: `_se_combined` (RAD loci and individuals) for Ho,
He, Fis and the per-site values, and `_se` (RAD loci only) for Ar and
privAr. `pct_poly` and `n` have no SE.

**What stays out, and where it goes.** Each column the summary leaves out is
left out for a reason.

| left out | why | where it goes |
|---|---|---|
| `_se`, `_se_ind`, `_lo`, `_hi` | The SE in the table is the most accurate of them. In the simulations, the locus-only versions held the true value as rarely as 41% (Ho) and 31% (FIS) of the time when individuals differed in inbreeding. Two uncertainty measures side by side invite readers to pick the narrower one. | the supplement (`outdir` files, or `summary(x, details = TRUE)`) |
| `Ar_n`, `privAr_n`, `priv_total`, `sites_used`, `variant_records` | These are counts for you to check, not results. | one sentence in the methods ("Ar was averaged over 739 loci") |
| Ho and He from the haplotype VCF, Fis from the SNP VCF | Haplotype Ho and He cannot be compared with other studies, and Fis is more precise from the haplotype VCF. | nowhere |
| The Welch and Wilcoxon p-values, the `se` of a population's mean heterozygosity | They are for comparison. The mean's `se` invites comparing populations by overlapping intervals, which the combined test replaces. | the supplement |

**What goes in the caption and methods** (text, not columns): the rarefaction
size `g` (Ar and privAr change with it), which VCF each column came from, which
standard error each column is, the filters and thresholds you used (call rate,
minor allele count), and the sample sizes. A caption that fits:
*"Ho, He, the percent of polymorphic SNPs and the per-site values are from the
SNP VCF; Fis, Ar and privAr are from the haplotype VCF, and Ar and privAr are
rarefied to g = 16 gene copies. Values are estimate (SE). SE for Ho, He, Fis and
the per-site values combine variation among RAD loci and among individuals; SE
for Ar and privAr are over RAD loci."*

**When to report more.** Usually never. Add something only when a check said
`look` and readers need the number to judge the table. For example, if `Ar` is
averaged over different loci in different populations, give the counts
(`Ar_n`) in a table footnote. If a journal asks for intervals for Ar and
privAr, replace `(SE)` by the `_lo` to `_hi` interval; do not show both. If the
data look filtered by minor allele count, give the threshold in the caption.

The per-site standard errors are the Ho and He ones times the same records /
sites multiplier, so they leave out uncertainty in how many SNPs there are per
sequenced site. The FST and D standard errors count RAD loci only, and no
simulation has checked them. If hierfstat is installed, `beta` is also in
`summary(dif)$tables$results`, as a point estimate with no SE.

### From the shell

The same analyses run from the shell with the scripts installed alongside
the package:

```bash
RD=$(Rscript -e 'cat(system.file("scripts", package = "RADdiversity"))')
Rscript $RD/diversity_stats.R out/populations.haps.vcf popmap.tsv --g=20 \
        --outdir=results
Rscript $RD/diversity_stats.R out/populations.snps.vcf popmap.tsv --g=20 \
        --outdir=results --sites=out/populations.sumstats_summary.tsv
Rscript $RD/het_between_pops.R out/populations.snps.vcf popmap.tsv --min-call=0.9 --outdir=het_out
```

The scripts print the short summary. Add `--details` for the full report. The
standard errors over individuals (`_se_combined`) are computed by default;
`--no-se-individuals` skips them. The scripts write their tables to the
working directory unless you give `--outdir`. There is no script for
`differentiation_stats()`; use it from R. Run any script with no arguments for
its options. `vignette("workflow", package = "RADdiversity")`
walks through a whole analysis.

## Before you start: Stacks `populations` flags to consider

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
* Use both `-r` and `-p`, not a global `-R`.
  A missing-data filter applied within each population (`-r`, with `-p` set
  to your number of populations) ensures each locus that passes is well
  genotyped in every population. This matters most when sampling is uneven:
  the global `-R` filter lets the larger population's coverage carry the smaller one, which
  then carries more missing data -- an artificial difference between exactly
  the two samples you mean to compare.
* Minor allele count or frequency filters (`--min-mac`, `--min-maf`) are
  common, and fine for many analyses, but they remove rare variants, which
  are real diversity: under a neutral site frequency spectrum, sites with a
  minor allele count of 2 or less hold about 21% of π at 10 diploids and 10%
  at 20. If you use one, report the threshold and compare diversity values
  only with data filtered the same way.
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
| Ho, He, F<sub>IS</sub>, % polymorphic, allelic richness, and private allelic richness | `diversity_stats()` |
| nucleotide diversity (π), d<sub>xy</sub>, and per-individual heterozygosity per sequenced site | `pi_allsites()` (all-sites VCF), or `diversity_stats(sites = ...)` (approximation) |
| do populations differ in observed heterozygosity? in inbreeding? | `het_between_pops()` |
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

When allowing multiple SNPs per locus, <i>populations</i> writes two VCF files: one with each SNP as a record (`populations.snps.vcf`) and one with each locus as a record (`populations.haps.vcf`). They should be used to answer different questions:

| statistic | file | why |
|---|---|---|
| H<sub>o</sub>, H<sub>e</sub>, per-site values (π) | `.snps` | per-site values are on a scale other studies can compare |
| F<sub>IS</sub> | `.haps` | a ratio, so the scale cancels; multi-allelic loci are more precise |
| allelic richness, private allelic richness | `.haps` | on biallelic SNPs richness can only be 1 or 2. These calculations depend on the number of SNPs per locus, so comparisons are most valid on populations analyzed together (same run, filters and *g*): which population has a greater number of private alleles and allelic richness is reasonable to compare, but saying how much richer a population is (e.g., 50% more) is not robust. |
| between-population tests | either -- say which | |

H<sub>o</sub> and H<sub>e</sub> are **not** comparable between the two files
(haplotype H<sub>o</sub> asks "are there any SNPs on this locus that are heterozygous?", per-site
H<sub>o</sub> "is this site heterozygous?"); F<sub>IS</sub> is. Each printed result says
which of its numbers to take.

## Which standard error to report

A population's value is uncertain for two reasons: you typed some of the
genome's RAD loci, and you caught some of the population's individuals.

* `_se`, `_lo`, `_hi` resample **RAD loci** and hold the individuals fixed.
* `_se_ind` is a jackknife over **individuals** that holds the loci fixed.
  `diversity_stats()` computes it by default (`se_individuals = FALSE` skips
  it).
* `_se_combined` puts both together: `sqrt(_se^2 + _se_ind^2)`.

**Report `_se_combined` for H<sub>o</sub>, H<sub>e</sub>, F<sub>IS</sub> and
the per-site H<sub>o</sub> and H<sub>e</sub>, and `_se` for A<sub>r</sub> and
privA<sub>r</sub>.** These are the standard errors a short `summary()` shows.
In simulations with a known truth (`inst/sims/uncertainty_sources.R`; 2 × 15
individuals, 1,000 loci, new loci and individuals in every repeat), 95%
intervals contained the true value this often (95% is perfect):

| | locus SE | individual SE | combined SE |
|---|---|---|---|
| H<sub>o</sub> | 41–96% | 88–91% | **93–98%** |
| H<sub>e</sub> | 95–96% | 72–77% | **98–99%** |
| F<sub>IS</sub> | 31–94% | 92–93% | **93–98%** |
| A<sub>r</sub> | **93–94%** | 98–100% | 100% (too wide) |
| privA<sub>r</sub> | **92–93%** | 31–36% | 92–94% |

The low locus-SE values are for populations whose individuals differ in
inbreeding. Over more settings (many rare variants, 8 or 30 individuals, 20%
missing genotypes, haplotype loci) the ranges widen to 90–99% for the combined
SE of H<sub>o</sub>, H<sub>e</sub> and F<sub>IS</sub>, and 88–96% for the
locus SE of A<sub>r</sub> and privA<sub>r</sub>. The bootstrap intervals over
loci (`_lo`, `_hi`) held the true value about as often as the locus SE, or
slightly less often, in every row (for A<sub>r</sub>, 92–93% against
93–94%), so the summaries show standard errors. For H<sub>e</sub> alone the locus SE was also close;
the combined SE is used for all three so that no choice depends on how much
individuals differ in inbreeding.

No simulation has checked the standard errors of F<sub>ST</sub>, D, beta, π,
d<sub>xy</sub> and d<sub>a</sub>; they count RAD loci only, and the summaries
say so. The combined heterozygosity test has its own simulation: with no true
difference it wrongly rejected 3–7.5% of the time, where Welch's t alone
rejected up to 30%. `vignette("rationale")`, section 4, explains how each SE
is calculated and why A<sub>r</sub> and privA<sub>r</sub> use the locus SE.
Never compare populations by overlapping intervals -- use
`het_between_pops()`.

## Glossary

* **H<sub>o</sub>** -- observed heterozygosity: the fraction of typed
  individuals that are heterozygous, averaged over loci.
* **H<sub>e</sub> (H<sub>s</sub>, gene diversity)** -- the chance that two
  gene copies drawn from the population differ (Nei & Chesser 1983).
* **π (nucleotide diversity)** -- H<sub>e</sub> averaged over every
  sequenced site, invariant sites included. Invariant sites contribute 0, so
  π is much smaller than H<sub>e</sub> per variant site. `pi_allsites()`
  reports two versions. `pi_nc` uses Nei & Chesser's H<sub>e</sub> at each
  site, and stays unbiased whether or not individuals are inbred. `pi` also
  compares the two gene copies inside each individual, which is what Stacks'
  `Pi`, pixy and VCFtools do. Inbreeding makes those two gene copies alike,
  so `pi` runs low in an inbred population. Report `pi_nc` as the estimate,
  and report `pi` when comparing with those programs.
* **F<sub>IS</sub>** -- 1 − H<sub>o</sub>/H<sub>e</sub> for a population,
  summed over loci; **F** -- the same for one individual.
* **gene copies** -- a diploid carries two: 10 individuals are 20 gene copies.
* **RAD locus vs record** -- a locus is one RAD tag (similar to an assembled contig). A SNP VCF has one record
  per SNP, so several records can be on the same locus unless filtered at some point to only allow for one SNP per locus.
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
| HWE filtering | a heterozygote deficit is the signal. `hwe_test()` reports departures. This package could be used to create a blacklist of SNPs or loci outside HWE that could then be used for filtering with a different program (e.g., <i>populations</i>). |
| null-allele correction | restriction-site null alleles cannot be removed by depth filtering; `diversity_stats()` reports `fis_by_call_rate` to detect them. |
| paralog detection beyond an excess-heterozygosity screen | `filter_max_het()` is that screen. Dedicated tools (e.g. HDplot; McKinney et al. 2017) go further and users are encouraged to use them if they think this could be an issue in their data. |
| N<sub>e</sub>, AMOVA, neutrality tests | different questions. |

## How it relates to other tools

* **Stacks** reports H<sub>o</sub>, π, and F<sub>IS</sub> per site and per
  population. This package recomputes F<sub>IS</sub> (ratio of sums,
  Nei–Chesser H<sub>e</sub>), adds rarefied allelic richness and rarefied private alleles,
  standard errors over RAD loci and individuals, and individual-level tests.
  Stacks measures π at a site as the chance that two of the gene copies
  typed there differ: `Pi` = 1 − Σ C(*n<sub>i</sub>*, 2)/C(*n*, 2), for *n*
  gene copies carrying allele counts *n<sub>i</sub>* (Hohenlohe et al. 2010).
  That formula is the same as (*n*/(*n* − 1))(1 − Σ*p<sub>i</sub>*<sup>2</sup>),
  which is H<sub>e</sub> with the 2*N* correction rather than Nei & Chesser's
  H<sub>e</sub>. The "All positions (variant and fixed)" block of
  `populations.sumstats_summary.tsv` averages `Pi` over every sequenced site,
  so that block alone is π per sequenced site; the "Variant positions" block
  is per variant site.
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

## Key references

David, P., Pujol, B., Viard, F., Castella, V. & Goudet, J. (2007) Reliable
selfing rate estimates from imperfect population genetic data. *Molecular
Ecology* 16:2474–2487.

Goudet, J. (2005) HIERFSTAT, a package for R to compute and test hierarchical
F-statistics. *Molecular Ecology Notes* 5:184–186.

Hemstrom, W. & Jones, M. (2023) snpR: user friendly population genomics for
SNP data sets with categorical metadata. *Molecular Ecology Resources*
23:962–973.

Hohenlohe, P.A., Bassham, S., Etter, P.D., Stiffler, N., Johnson, E.A. &
Cresko, W.A. (2010) Population genomics of parallel adaptation in threespine
stickleback using sequenced RAD tags. *PLoS Genetics* 6:e1000862.

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
