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

## How often 95% intervals built from the recommended SE held the true value in
## the package's own simulations (inst/sims/uncertainty_sources.R, parts `main`
## and `robust`; the coverage tables in vignette("rationale"), section 4). The
## short summaries quote these ranges. They are the wider ranges over all the
## simulated settings, so nothing is overstated; keep them in step with that
## vignette. The het test's number is its false-positive rate when populations
## do not differ (inst/sims/het_test_null.R; ?het_between_pops).
##   combined  Ho, He, Fis with _se_combined (SE over loci AND individuals)
##   loci      Ar, privAr with _se (SE over RAD loci only)
##   het_test  the combined heterozygosity test, nominal 5%
.coverage <- list(combined = "90-99%", loci = "88-96%", het_test = "3-7.5%")

.report_text <- list(

  ## ---- summary() short views, diversity_stats() -------------------------------

  ## What each column of the short RESULTS table is. {se} is filled in with the
  ## kind of SE the table shows.
  diversity_legend_snp = c(
    "Ho, He    observed and expected heterozygosity, averaged over SNPs.",
    "pct_poly  percent of SNPs that vary within the population. It has no SE."),
  diversity_legend_hap = c(
    "Fis       inbreeding: 0 = random mating; above 0 = fewer heterozygotes",
    "          than expected.",
    "Ar        allelic richness: alleles per locus, counted in g = {g} gene",
    "          copies and averaged over loci.",
    "privAr    private allelic richness: the same, counting only alleles found",
    "          in this population and in no other."),

  diversity_take_snp_short = c(
    "Get Fis, Ar and privAr from the haplotype VCF: on SNPs, Ar cannot pass 2.",
    "Never compare Ho or He between the SNP VCF and the haplotype VCF.",
    "Report pct_poly beside Ho and He."),
  diversity_take_haplotype_short = c(
    "Get Ho and He from the SNP VCF; haplotype Ho and He are not comparable.",
    "Compare Ar and privAr only among populations analyzed together (same",
    "run, filters and g)."),

  ## ---- diversity_stats() ----------------------------------------------------

  diversity_complete_case_dropped = c(
    "^ most records were dropped for missing data. complete_case keeps a",
    "record only if genotyped in EVERY individual (Schmidt et al. 2021), so a",
    "few poor libraries can cost most of the dataset. Check the per-individual",
    "table from het_between_pops() (loci_called column) before accepting",
    "this -- or drop complete_case."),

  diversity_available_data = c(
    "Each population's Ho, He and Fis use every locus where that population has at",
    "least min_n genotyped individuals. The number used for each population is",
    "shown while the analysis runs (unless verbose = FALSE)."),

  diversity_se_columns = c(
    "_se: leave one RAD locus out at a time, recompute, and see how much the",
    "answer moves (a jackknife; Weir 1996). It needs no bootstrap and holds the",
    "individuals fixed. _lo and _hi: the middle 95% of the answers when RAD loci",
    "are drawn again at random, many times (a bootstrap over loci). They also",
    "hold the individuals fixed."),

  diversity_se_ind_columns = c(
    "_se_ind columns: leave one INDIVIDUAL out at a time (a jackknife). It holds",
    "the loci fixed. _se_combined = sqrt(_se^2 + _se_ind^2): it counts both which",
    "loci and which individuals were sampled. Report _se_combined for Ho, He and",
    "Fis. For Ar and privAr, report _se (see the note under those tables). The",
    "percentages under each table are how often a 95% interval held the true",
    "value in the package's simulations (2 populations of 15 individuals, 1,000",
    "loci). More in ?diversity_stats, \"How standard errors are calculated\", and",
    "vignette(\"rationale\"), section 4."),

  diversity_se_loci_only = c(
    "These SEs hold the individuals fixed, because se_individuals = FALSE was set.",
    "For Ho, He and Fis, re-run with se_individuals = TRUE (the default) and",
    "report _se_combined (?diversity_stats)."),

  ## %s is the number of loci in all.
  diversity_ar_n_columns = c(
    "Ar_n: the loci that went into that population's Ar (it has at least g gene",
    "copies typed there). privAr_n: the loci that went into privAr (EVERY",
    "population has at least g typed there), so privAr_n is never larger than",
    "Ar_n, and can be much smaller when there is missing data. There are %s loci",
    "in all."),

  ## Under each "one small table per statistic" block of the full report. What
  ## to report, and how often a 95% interval from each measure held the true
  ## value in the main simulation (2 populations of 15 individuals, 1,000 loci;
  ## vignette("rationale") section 4, inst/sims/uncertainty_sources.R).
  diversity_block_notes = list(
    Ho = c("Report Ho_se_combined (loci and individuals): its 95% intervals held the truth",
           "93-98% of the time. Ho_se (loci only): 41-96%. Ho_lo to Ho_hi (bootstrap over",
           "loci): 39-95%. Ho_se_ind (individuals only): 88-91%."),
    He = c("Report He_se_combined: 98-99%. He_se alone: 95-96% (He depends mostly on which",
           "loci were typed). He_lo to He_hi: 93-95%. He_se_ind: 72-77%."),
    Fis = c("Report Fis_se_combined: 93-98%. Fis_se alone: 31-94%. Fis_lo to Fis_hi:",
            "30-94%. Fis_se_ind: 92-93%."),
    Ar = c("Report Ar_se (loci only): 93-94%. Ar_lo to Ar_hi: 92-93%. No individual SE is",
           "computed for Ar: it was too wide (1.4 to 1.5 times the true SE)."),
    privAr = c("Report privAr_se (loci only): 92-93%. privAr_lo to privAr_hi: 91-93%. No",
               "individual SE is computed: it was far too small (31-36%)."),
    Ho_autosomal = "Report Ho_autosomal_se_combined: the Ho SE times records / sites.",
    He_autosomal = "Report He_autosomal_se_combined: the He SE times records / sites."),

  ## The same, when se_individuals = FALSE (there is no _se_combined to report).
  diversity_block_notes_loci_only = list(
    Ho = "Ho_se counts RAD loci only. The SE to report needs se_individuals = TRUE.",
    He = "He_se counts RAD loci only. The SE to report needs se_individuals = TRUE.",
    Fis = "Fis_se counts RAD loci only. The SE to report needs se_individuals = TRUE.",
    Ar = c("Report Ar_se (loci only): 93-94%. Ar_lo to Ar_hi: 92-93%."),
    privAr = c("Report privAr_se (loci only): 92-93%. privAr_lo to privAr_hi: 91-93%.")),

  diversity_ar_capped = c(
    "NOTE: Ar is capped at 2 because these are biallelic SNPs, so rarefaction",
    "has almost nothing to correct. Run this on the HAPLOTYPE VCF for a",
    "richness worth reporting; private allelic richness is still informative."),

  diversity_two_estimators = c(
    "He_NeiChesser  (n/(n-1))(1 - sum p^2 - Ho/2n): unbiased at any FIS.",
    "He_2n_corr     (2n/(2n-1))(1 - sum p^2): unbiased only if FIS = 0. This is the",
    "               estimator behind Stacks' `Pi` column."),

  diversity_fis_sensitivity = c(
    "He barely changes between the two. FIS changes about (1 - F)/F times as much,",
    "because FIS divides by the small number the two estimators disagree about.",
    "So Stacks' Pi is safe to quote, but its Fis is not. (Its Fis also has a",
    "second, larger error: it is a mean of per-locus ratios.)"),

  diversity_looks_mac_filtered = c(
    "No variable record has an allele seen only once or twice, so these data look",
    "filtered by minor allele count (a common RAD-seq step). State the threshold in",
    "your methods. Removing rare alleles lowers He per sequenced site, pct_poly, Ar",
    "and privAr, and lowers them more in smaller samples. So compare values only",
    "between datasets filtered the same way."),

  diversity_prior_filters = c(
    "These describe the popmap's individuals as given; samples in the VCF but not",
    "in the popmap are left out (an HWE filter cannot be detected). See",
    "vignette(\"workflow\"), \"If your data were already filtered\"."),

  diversity_fis_by_call_rate = c(
    "Each population's records are grouped by that population's own call rate",
    "there. Inbreeding raises FIS at every record alike. Null alleles",
    "(restriction-site mutations) and allele dropout make heterozygotes look",
    "homozygous or go missing, so they show up as FIS RISING as call rate falls.",
    "If you see that, also report FIS from the well-typed records (complete_case",
    "= TRUE, or a higher min_n), and be careful comparing FIS between",
    "populations: dropout is more common in more diverse populations (Gautier et",
    "al. 2013)."),

  diversity_sites_on_haplotypes = c(
    "`sites` was supplied but this is a HAPLOTYPE VCF, so the autosomal",
    "conversion was SKIPPED: He here is per-tag gene diversity, not per-site,",
    "and dividing by sequenced sites would understate nucleotide diversity.",
    "Re-run on populations.snps.vcf with the same `sites` for He_autosomal."),

  diversity_autosomal_meaning = c(
    "He_autosomal is the same per-site quantity as He above, averaged over every",
    "sequenced site instead of over variant records. That denominator is what",
    "makes it nucleotide diversity (pi). There is no column named `pi`. Compare",
    "it with Stacks' `Pi` from the All-positions block, not the variant-positions",
    "block and not `Exp_Het` (Schmidt et al. 2021).",
    "The SE and interval of a per-site value are those of Ho or He times the same",
    "records / sites multiplier. They leave out uncertainty in how many SNPs",
    "there are per sequenced site, which the value itself also takes as fixed.",
    "`sites` must count the sequenced sites of the loci in these data. Stacks'",
    "Sites does so only if nothing was filtered after populations",
    "(?diversity_stats)."),

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
    "NOT Ho or He. A haplotype locus's gene diversity has a ceiling that rises",
    "with the number of SNPs in the tag, so it changes with read length, enzyme",
    "and SNP calling and cannot be compared across studies. Run",
    "populations.snps.vcf for those. Fis is a ratio, so the scale cancels: that",
    "is why it belongs here.",
    "Ar and privAr depend on the number of SNPs per locus too (more SNPs, more",
    "possible haplotypes). Compare them only among populations analyzed together",
    "(same run, filters and g), never with other studies. Which population is",
    "richer is reliable; how much richer also depends on SNPs per locus. The loci",
    "averaged over are the ones variable in this run, so adding or dropping",
    "populations from the Stacks run changes them."),

  diversity_take_snp = c(
    "TAKE FROM THIS RUN:  Ho, He, pct_poly, Ho_autosomal, He_autosomal.",
    "Run populations.haps.vcf for Fis, Ar and privAr: on biallelic SNPs, rarefied",
    "richness cannot pass 2, and Fis is more precise from multi-allelic loci."),

  diversity_take_common = c(
    "Fis is comparable between the two files. Ho and He are NOT: haplotype Ho",
    "asks 'is this tag heterozygous?', per-site Ho asks 'is this site",
    "heterozygous?'. Label every table row with the file it came from.",
    "",
    "Ho and He above are averaged over ALL retained records, including records",
    "where the population does not vary (the 'pooled' denominator, as in Stacks'",
    "Variant-positions block). Report pct_poly beside them, because:",
    "    He(pooled) = He(variable in that population) x pct_poly/100   (exactly)",
    "Do not switch to averaging only over records that vary within a population:",
    "that conditions on the outcome and can reverse a real difference."),

  diversity_he_difference_loci = c(
    "CAUTION: this interval treats LOCI as the only replicate (boot=\"loci\") and",
    "holds the individuals fixed. It is a description, not a test. To test a",
    "difference in heterozygosity with both individuals and loci counted, use",
    "het_between_pops()."),

  diversity_he_difference_individuals = c(
    "locus bootstrap -- but it is still not a paired individual-level test. Use",
    "het_between_pops() for the p-value you report."),

  ## ---- het_between_pops() ---------------------------------------------------

  het_failed_individuals = c(
    "Each rests on few loci, yet counts as a full individual in the means and",
    "tests. Heavy missing data usually comes with allele dropout (heterozygosity",
    "biased low, F high). It also lowers its population's call rate, so fewer",
    "loci clear min_call. Remove them from the popmap (or with filter_samples())",
    "if they were not meant to be analyzed."),

  het_confound_positive = c(
    "POSITIVE and significant: individuals with more missing data look LESS",
    "heterozygous, which is the allele-dropout signature. Re-run on a",
    "complete-data locus set (min_call = 1.0) before quoting the test below; if",
    "the difference survives that, it is not a coverage artifact."),

  het_call_rate_gap = c(
    "WARNING: mean call rate differs between populations by more than 2",
    "percentage points. That asymmetry is itself a candidate explanation."),

  het_sd_among_individuals = c(
    "'sd' is the spread AMONG INDIVIDUALS. That spread divided by sqrt(n) is one",
    "part of the uncertainty in a population mean. It does not shrink as you add",
    "loci, which is why tests over loci alone can be invalid here. The other part",
    "is which loci were typed. The combined test counts both."),

  het_van_dongen = c(
    "This is the point made by Van Dongen (1995, Heredity 74:445-447): the unit",
    "of resampling changes what the bootstrap means, and loci alone are the",
    "wrong unit because they are all measured on the same individuals."),

  het_g2 = c(
    "95% CI: bootstrap over individuals (see ?identity_disequilibrium). A CI above",
    "0 means individuals differ in inbreeding, so standard errors over loci alone",
    "are too small for Ho and Fis: report _se_combined from diversity_stats()",
    "(computed by default)."),

  het_omnibus = c(
    "Welch's one-way ANOVA (unequal variances allowed) and Kruskal-Wallis ask ONE",
    "question: do any of the populations differ? Report it first. If it is not",
    "significant, do not interpret single pairs. If it is, the pairwise tests",
    "below (BH-adjusted) say which populations differ. These overall tests use",
    "individuals only, so with differentiated populations they can reject too",
    "readily: confirm a difference with the pairwise p_combined_BH."),

  het_F_tests = c(
    "The same tests on individual inbreeding, F = 1 - observed/expected",
    "heterozygosity (expected from each individual's own population; see",
    "?individual_inbreeding). Heterozygosity (observed, He x (1 - F)) reflects",
    "DIVERSITY and INBREEDING together; F asks whether they differ in INBREEDING."),

  het_interpretation = c(
    "* Power is limited mainly by the NUMBER OF INDIVIDUALS, not the number of",
    "  loci. At n = 15 vs 10, Welch's test detected a difference in mean F of half",
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
    "_se: leave one RAD locus out at a time, recompute, and see how much the answer",
    "moves (a jackknife, as in diversity_stats()). _lo and _hi: the middle 95% of",
    "the answers when RAD loci are drawn again at random (a bootstrap over loci).",
    "Both hold the individuals fixed."),

  ## Under the FST and D tables. No simulation has checked how often these
  ## intervals hold the true value; the diversity_stats() SEs have been checked.
  differentiation_block_note = c(
    "Report the SE. It counts RAD loci only, and no simulation has checked how",
    "often its 95% interval holds the true value."),
  differentiation_fis_note = c(
    "FIS here is Weir & Cockerham's, pooled over populations. The FIS to report is",
    "the Fis from diversity_stats() (with Fis_se_combined)."),

  differentiation_one_pop_records = c(
    "they say nothing about differences between populations (or, with one",
    "individual per population, about variation among individuals), and keeping",
    "them would pull FST toward 0 (hierfstat::wc() keeps them; VCFtools and",
    "Stacks' own Fst skip records typed in one population). See ?differentiation_stats."),

  differentiation_d_global_records = c(
    "Why: Jost's D for k populations measures differentiation among those k",
    "populations. A record typed in only some of them measures it among a",
    "different set, and averaging such records in biased the global D low",
    "(0.875 instead of 1 for three completely differentiated populations, one",
    "of them untyped at half the records). Pairwise D is not affected: each",
    "pair uses every record typed in both of its populations."),

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
    "_se: leave one RAD locus out at a time and see how much the answer moves (a",
    "jackknife). _lo and _hi: the middle 95% of the answers when RAD loci are drawn",
    "again at random (a bootstrap over loci). sites: sites with at least two typed",
    "gene copies in that population.",
    "pi_nc: the same, but comparing only gene copies from DIFFERENT individuals. At",
    "each site it is Nei & Chesser's He, unbiased when FIS != 0. `pi` (as in pixy and",
    "Stacks' Pi) also compares the two copies inside each individual, so it runs low",
    "when individuals are inbred; quote it only to compare with those tools."),

  ## Under each small table of the full report.
  pi_block_notes = list(
    pi_nc = c("Report pi_nc and pi_nc_se. The SE counts RAD loci only, and no simulation has",
              "checked how often its 95% interval holds the true value."),
    pi = c("pi is for comparing with pixy, VCFtools and Stacks. It runs low when",
           "individuals are inbred."),
    dxy = c("dxy compares only copies from different individuals, so inbreeding does not",
            "affect it. Report dxy_se (loci only; not checked by simulation)."),
    da_nc = "Report da_nc. It is dxy minus the mean pi_nc of the two populations.",
    da = c("da is dxy minus the mean pi, so it comes out a little high when individuals",
           "are inbred."))
)
