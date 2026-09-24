# RADdiversity 0.1.0

First release as an R package. The full development history, including every
review and the numbers each change moved, is in the git log of the GitHub
repository.

## What the package does

* `diversity_stats()`: Ho, He (Nei & Chesser 1983), FIS as a ratio of sums,
  rarefied allelic richness and rarefied private alleles, and Ho and He per
  sequenced site. Standard errors over RAD loci and, by default, over individuals
  (`_se_combined`, the one to report for Ho, He and FIS).
* `diversity_table()`: the SNP and haplotype results as one table for a
  manuscript, with a caption that says where each column came from and which
  standard error it holds.
* `het_between_pops()`: do populations differ in heterozygosity or in
  inbreeding? One value per individual, and a test that counts both
  individuals and RAD loci (3–7.5% false positives in simulation, where
  Welch's t alone reached 30%).
* `differentiation_stats()`: Weir & Cockerham FST, Jost's D (with Nei &
  Chesser's Hs and Ht) and, with hierfstat, Weir & Goudet's beta.
* `pi_allsites()`: nucleotide diversity (`pi_nc` and pixy's `pi`), dxy and
  net divergence from an all-sites VCF.
* `isolation_by_distance()`: do populations farther apart differ more? A
  one-sided Mantel test (exact with 7 or fewer populations) of pairwise
  FST/(1 - FST) (Rousset 1997) or Jost's D against distance (`habitat = "1D"`)
  or log distance (`"2D"`), with the slope, `summary()` and `plot()`.
  `read_distances()` reads the distances as a square matrix (one half may be
  blank) or a list of pairs, and matches populations by name. A simulated
  six-site river dataset (`ibd_example.*` in `inst/extdata`) shows it.
* `individual_inbreeding()`, `identity_disequilibrium()` (g2, leaving out
  pairs of SNPs on the same RAD locus), `kinship_check()` (KING-robust or
  Goudet's beta) and `hwe_test()` (an exact test, a report and never a
  filter, with Benjamini-Hochberg p-values).
* Reading, filtering and export: `read_stacks_vcf()`, `read_popmap()`, the
  `filter_*()` functions (each recorded in `H$filter_log`), the `write_*()`
  functions (VCF, PLINK, STRUCTURE, Genepop, FSTAT, RADpainter), and
  `as_genind()` / `as_genlight()` for adegenet and dartR. Genepop and
  STRUCTURE files write a `.` in a locus name as `_` (`6.1` becomes `6_1`),
  so `read.genepop()` and `read.structure()` in adegenet can read them; the
  tests read every exported format back with adegenet or hierfstat.
* Every result prints its main tables; `summary()` is a short view (each
  estimate with its recommended SE, what to take from the run, and checks
  marked `ok`, `info` or `look`); `summary(details = TRUE)` is the full report.

## For users of the earlier stand-alone scripts

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

* Numbers that changed from the scripts, and why:
  * Jost's D uses Nei & Chesser's Hs and Ht, which stay unbiased when
    individuals are inbred (about 7% lower on the example data).
  * The global D uses only records typed in every population.
  * FST skips records typed in only one of the compared populations, as
    VCFtools and Stacks do.
  * Per-site values multiply by the data's own SNP count, not Stacks'
    `Variant_Sites`.
  * g2 leaves out pairs of SNPs on the same RAD locus. Counting them made the
    test find g2 above 0 in 54-56% of simulated data sets with no inbreeding;
    now 2-4% (`inst/sims/g2_linked_snps.R`).
  * The short summaries and tables report rarefied private alleles as a
    count (`priv_total`, with SE = `privAr_se` × `privAr_n`) rather than per
    locus (`privAr`, still in the results), and leave `pct_poly` to the full
    tables.
* `seed` defaults to `NULL` everywhere: call `set.seed()` first for
  reproducible intervals, or pass `seed`, which leaves your session's random
  numbers as they were. The command-line scripts use `--seed=2024` by default.
* The package needs R 4.0 or later, and only base R.
