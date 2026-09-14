###############################################################################
#
#  R/report_text.R -- the fixed paragraphs of the full reports printed by
#  summary() on a result (and by the command-line scripts).
#
#  The report functions (.diversity_report() in R/diversity_stats.R,
#  .het_report() in R/het_between_pops.R, and so on) print the numbers and
#  decide which paragraphs apply; the wording lives here, so it can be read
#  and edited in one place, separately from the code. Each element is a
#  character vector with one element per printed line (the report indents
#  each line by two spaces).
#
###############################################################################

.report_text <- list(

  ## ---- diversity_stats() ----------------------------------------------------

  diversity_complete_case_dropped = c(
    "^ most records were dropped for missing data. complete_case keeps a",
    "record only if genotyped in EVERY individual (Schmidt et al. 2021), so a",
    "few poor libraries can cost most of the dataset. Check the per-individual",
    "table from het_between_pops() (loci_called column) before accepting",
    "this -- or drop complete_case."),

  diversity_available_data = c(
    "Each population's Ho/He/Fis is a ratio of sums over whichever loci it",
    "cleared min_n at -- see the per-population lines printed while it ran."),

  diversity_se_columns = c(
    "_se columns: delete-one-block jackknife over RAD loci (Weir 1996) --",
    "always available (no nboot needed), on the same loci axis as _lo/_hi."),

  diversity_se_ind_columns = c(
    "_se_ind columns: delete-one-INDIVIDUAL jackknife -- the uncertainty from",
    "which individuals were sampled, which the loci-based _se/_lo/_hi hold fixed.",
    "Report _se_ind when individuals differ in inbreeding (g2 > 0; see",
    "identity_disequilibrium() or the het_between_pops() report)."),

  diversity_ar_n = c(
    "(Ar needs only THAT population at >= g gene copies; privAr needs EVERY",
    "population at >= g simultaneously, so privAr_n <= Ar_n, often far less",
    "under available data -- see vignette(\"rationale\") if privAr_n is small)."),

  diversity_ar_capped = c(
    "NOTE: Ar is capped at 2 because these are biallelic SNPs, so rarefaction",
    "has almost nothing to correct. Run this on the HAPLOTYPE VCF for a",
    "richness worth reporting; private allelic richness is still informative."),

  diversity_two_estimators = c(
    "He_NeiChesser  (n/(n-1))(1 - sum p^2 - Ho/2n)   unbiased at any FIS",
    "He_2n_corr     (2n/(2n-1))(1 - sum p^2)         unbiased only if FIS = 0",
    "               -- this is the estimator behind Stacks' `Pi` column."),

  diversity_fis_sensitivity = c(
    "He barely moves; FIS moves ~(1-F)/F times as much, because FIS divides by",
    "the small quantity the two estimators disagree about. So Stacks' Pi is safe",
    "to quote and its Fis is not (and its Fis carries a second, larger error:",
    "it is a mean of per-locus ratios)."),

  diversity_looks_mac_filtered = c(
    "No variable record has an allele seen only once or twice, so these data look",
    "filtered by minor allele count (a common RAD-seq step). State the threshold",
    "in your methods: removing rare alleles lowers He per sequenced site,",
    "pct_poly, Ar and privAr, by more in smaller samples, so values are comparable",
    "only between datasets filtered the same way."),

  diversity_prior_filters = c(
    "These describe the data as given (an HWE filter cannot be detected). See",
    "vignette(\"workflow\"), \"If your data were already filtered\"."),

  diversity_fis_by_call_rate = c(
    "Each population's records grouped by ITS OWN call rate there. Inbreeding",
    "raises FIS at every record alike; null alleles (restriction-site mutations)",
    "and allele dropout make heterozygotes look homozygous or go missing, so they",
    "show up as FIS RISING as call rate falls. If so, report FIS from the",
    "well-typed records too (complete_case = TRUE, or a higher min_n) and treat",
    "differences in FIS between populations with care: dropout is more common",
    "in more diverse populations (Gautier et al. 2013)."),

  diversity_sites_on_haplotypes = c(
    "`sites` was supplied but this is a HAPLOTYPE VCF, so the autosomal",
    "conversion was SKIPPED: He here is per-tag gene diversity, not per-site,",
    "and dividing by sequenced sites would understate nucleotide diversity.",
    "Re-run on populations.snps.vcf with the same `sites` for He_autosomal."),

  diversity_autosomal_meaning = c(
    "He_autosomal is the SAME per-site quantity as He above, averaged over every",
    "sequenced site instead of over variant records. That denominator is what",
    "makes it nucleotide diversity (pi). There is no column named `pi`.",
    "Compare against Stacks' `Pi` from the All-positions block -- not the",
    "variant-positions block, and not `Exp_Het` (Schmidt et al. 2021)."),

  diversity_complete_case_representative = c(
    "them being representative. If missingness is concentrated in a few poor",
    "libraries it is not -- drop complete_case (the default already uses",
    "available data per population) or drop those individuals and re-run."),

  diversity_sites_implausible = c(
    "this run called, which is not a plausible sequenced-site count -- it was",
    "IGNORED rather than used. Pass the `Sites` column of the All-positions",
    "block of populations.sumstats_summary.tsv, which should be orders of",
    "magnitude larger than the number of variant records."),

  diversity_per_ascertained_record = c(
    "Per ASCERTAINED record: the denominator is however many markers this run",
    "called. Valid for comparing these populations to each other and nothing",
    "else. Pass `sites` (the `Sites` column of the All-positions block of",
    "populations.sumstats_summary.tsv) to get the autosomal value."),

  diversity_take_haplotype = c(
    "TAKE FROM THIS RUN:  Fis, Ar, privAr.",
    "NOT Ho / He: haplotype gene diversity has no fixed ceiling, so it moves",
    "with read length, enzyme and filters and is not comparable across",
    "studies. Run populations.snps.vcf for those. Fis is a ratio, so the",
    "scale cancels -- which is why it belongs here."),

  diversity_take_snp = c(
    "TAKE FROM THIS RUN:  Ho, He, pct_poly, Ho_autosomal, He_autosomal.",
    "Run populations.haps.vcf for Fis, Ar and privAr: on biallelic SNPs",
    "rarefied richness is bounded at 2, and Fis is more precise from",
    "multi-allelic loci."),

  diversity_take_common = c(
    "Fis is comparable between the two files; Ho and He are NOT (haplotype Ho",
    "asks 'is this tag heterozygous', per-site Ho asks 'is this site'). Label",
    "every table row with the file it came from.",
    "",
    "He and Ho above are averaged over ALL retained records including those",
    "monomorphic within a population -- the POOLED denominator, matching",
    "Stacks' Variant-positions block. Report pct_poly beside them:",
    "    He(pooled) = He(polymorphic within pop) x pct_poly/100   (exactly)",
    "Do not switch to the within-population denominator; it conditions on the",
    "outcome and can reverse a real difference."),

  diversity_he_difference_loci = c(
    "CAUTION: this interval treats LOCI as the replicate (boot=\"loci\"). For a",
    "comparison of population MEAN heterozygosity that is pseudoreplication --",
    "the uncertainty is dominated by which INDIVIDUALS you sampled. Use",
    "het_between_pops() for the p-value you report."),

  diversity_he_difference_individuals = c(
    "locus bootstrap -- but it is still not a paired individual-level test. Use",
    "het_between_pops() for the p-value you report."),

  ## ---- het_between_pops() ---------------------------------------------------

  het_confound_positive = c(
    "POSITIVE and significant: individuals with more missing data look LESS",
    "heterozygous, which is the allele-dropout signature. Re-run on a",
    "complete-data locus set (min_call = 1.0) before quoting the test",
    "below; if the difference survives that, it is not a coverage artefact."),

  het_call_rate_gap = c(
    "WARNING: mean call rate differs between populations by more than 2",
    "percentage points. That asymmetry is itself a candidate explanation."),

  het_sd_among_individuals = c(
    "'sd' is the spread AMONG INDIVIDUALS. That spread, divided by sqrt(n), is",
    "the real uncertainty in a population mean -- and it does not shrink at all",
    "as you add loci, which is why locus-based tests can be invalid here."),

  het_overdispersion_undefined = c(
    "Overdispersion is undefined for every population (no variation in either",
    "heterozygosity or per-locus heterozygote frequency) -- too degenerate a",
    "dataset to say whether a locus bootstrap would fail here."),

  het_overdispersion_small = c(
    "Overdispersion is near 1: individuals differ about as much as coin-flip",
    "noise alone predicts. For this dataset a locus bootstrap would be roughly",
    "correctly calibrated, and the two approaches should broadly agree. The",
    "individual-level test below remains the safer default."),

  het_overdispersion_moderate = c(
    "Overdispersion is moderate. Individuals vary by more than coin-flip noise,",
    "so a locus bootstrap understates the standard error by roughly the factor",
    "in the last column and its p-values are too small. Use the test below."),

  het_overdispersion_large = c(
    "Overdispersion is LARGE. Individuals differ far more than coin-flip noise,",
    "so a locus bootstrap would understate the standard error by the factor in",
    "the last column and reject far too often. Use the test below, and inspect",
    "the per-individual table for outliers (a bad library looks like one",
    "individual with unusually LOW heterozygosity)."),

  het_van_dongen = c(
    "This is the point made by Van Dongen (1995, Heredity 74:445-447): the unit",
    "of resampling changes what the bootstrap means, and loci are usually the",
    "wrong unit because they are all measured on the same individuals."),

  het_g2 = c(
    "95% CI: bootstrap over individuals (see ?identity_disequilibrium). A CI above",
    "0 means individuals differ in inbreeding, so diversity_stats()'s locus-based",
    "intervals are too narrow for population-level inference: report its",
    "individual-jackknife SEs (se_individuals = TRUE) instead."),

  het_omnibus = c(
    "Welch's one-way ANOVA (unequal variances allowed) and Kruskal-Wallis ask ONE",
    "question: do any of the populations differ? Report it first. If it is not",
    "significant, do not interpret individual pairs; if it is, the pairwise tests",
    "below (BH-adjusted) say which populations differ."),

  het_F_tests = c(
    "The same tests on individual inbreeding, F = 1 - observed/expected",
    "heterozygosity (expected from each individual's own population; see",
    "?individual_inbreeding). Heterozygosity asks whether the populations differ",
    "in DIVERSITY; F asks whether they differ in INBREEDING."),

  het_interpretation = c(
    "* Power is limited by the NUMBER OF INDIVIDUALS, not the number of loci.",
    "  At n = 15 vs 10, Welch's test detected a difference in mean F of half",
    "  the among-individual SD 17% of the time, and of 0.8 SD 39% of the time",
    "  (inst/sims/vignette_sims.R). A non-significant result is weak evidence",
    "  of no difference, not evidence of no difference.",
    "* If some individuals are relatives, they are not independent units and",
    "  even this test is anti-conservative. Check relatedness first",
    "  (kinship_check()).",
    "* A single poor-quality library shows up as one individual with unusually",
    "  LOW heterozygosity. Inspect the per-individual table before trusting a",
    "  result, especially at small n where one outlier can drive it."),

  ## ---- differentiation_stats() ----------------------------------------------

  differentiation_se_columns = c(
    "_se columns: delete-one-block jackknife over RAD loci (same convention",
    "as diversity_stats()); _lo/_hi: bootstrap 95% CI over the same RAD loci."),

  differentiation_one_pop_records = c(
    "they say nothing about differences between populations, and keeping them",
    "would pull FST toward 0 (hierfstat::wc() keeps them; VCFtools and Stacks'",
    "own Fst skip them). See ?differentiation_stats."),

  differentiation_haplotype = c(
    "This is a HAPLOTYPE VCF: prefer D (or beta) over FST here -- FST's",
    "ceiling shrinks as marker diversity grows, and haplotype loci are",
    "exactly the highly-variable-marker case that distorts (see",
    "?differentiation_stats and Jost 2008)."),

  ## ---- pi_allsites() ---------------------------------------------------------

  pi_definition = c(
    "pi = sum of pairwise differences / sum of pairwise comparisons, site by",
    "site, among the gene copies typed there (pixy; Korunes & Samuk 2021)"),

  pi_columns = c(
    "_se: delete-one-block jackknife over RAD loci; sites: sites with at least",
    "two typed gene copies in that population.",
    "pi_nc: the same, comparing only gene copies from DIFFERENT individuals --",
    "Nei & Chesser's He per site, unbiased when FIS != 0. `pi` (pixy, Stacks' Pi)",
    "also compares the two copies inside each individual and runs low when",
    "individuals are inbred; quote it only to compare with those tools.")
)
