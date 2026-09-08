#!/usr/bin/env Rscript
# =============================================================================
#  FIS, population heterozygosity, allelic richness, and private alleles
#  -- ROBUST VERSION of the script you wrote.
# =============================================================================
#
#  This is your original script, with five real problems fixed. They're
#  summarised here, in one place, so you can see WHAT changed and WHY without
#  diffing code line by line. Ranked by how much each one could have changed
#  your actual conclusions:
#
#  1. THE POPULATION-COMPARISON TESTS WERE INVALID (your old Step 9: Welch
#     ANOVA, Friedman test, paired pairwise tests -- all run on Ho ~ Locus).
#     Every locus in your data is measured on the SAME set of individuals. So
#     treating "locus" as the unit that gets compared is like weighing the
#     same 20 people 2,000 times and pretending you have 2,000 independent
#     people -- it isn't independence, it's the same people measured
#     repeatedly. The real unit of replication is the INDIVIDUAL.
#     How much does this matter? This was measured directly, not guessed at:
#     with two populations that have IDENTICAL true heterozygosity, treating
#     locus as the replicate gives a false-positive rate of 5% ONLY if every
#     individual happens to have exactly the same heterozygosity as every
#     other individual. As soon as individuals differ from each other even a
#     little -- which real biological data always does -- the false-positive
#     rate jumps to 55-76%. In plain terms: with your old Step 9, more than
#     half of the "significant" population differences it reports could be
#     nothing more than a couple of unusual individuals, not a real
#     population-level effect.
#     FIXED: Step 9 below computes ONE heterozygosity number per INDIVIDUAL,
#     then compares populations using individuals (not loci) as the unit --
#     the design that was actually measured to hold a 5% false-positive rate.
#
#  2. THE BOOTSTRAP RESAMPLED SNPs, NOT WHOLE RAD TAGS. A "RAD tag" is one
#     piece of DNA that got sequenced; often several SNPs sit on the SAME
#     tag. Those SNPs are NOT independent evidence -- they were read off the
#     same DNA fragment, in the same individuals, at the same time. Your
#     original bootstrap picked random SNP ROWS one at a time
#     (`sample(nloc, replace=TRUE)`), which silently treats every SNP as
#     independent even when several came from one tag. That makes your
#     confidence intervals falsely narrow (too precise), which can make a
#     difference LOOK statistically significant when it's really just noise
#     dressed up by overcounted "independent" evidence.
#     FIXED: the bootstrap below resamples whole RAD TAGS -- every SNP on a
#     chosen tag moves together, every time.
#
#  3. THE PRIVATE-ALLELE "RAREFACTION" WAS BIASED. Your original repeatedly
#     drew random INDIVIDUALS (without replacement) and recounted private
#     alleles each time, then averaged. The problem: the two gene copies
#     inside ONE individual are only an independent draw from the population
#     if FIS is exactly 0. Your data's FIS is NOT 0 (it's often clearly
#     negative, meaning MORE heterozygotes than random mating predicts) --
#     so drawing individuals (2 linked copies at a time) instead of drawing
#     gene copies directly gives a biased answer, in a direction and size
#     that depends on FIS.
#     FIXED: private allele richness below uses the standard ANALYTIC formula
#     (the same one the field-standard tools HP-RARE and ADZE use) that
#     computes the exact expected count for a random draw of GENE COPIES
#     directly. No simulation, no bias, and it gives the exact same answer
#     every time you run it (nothing random about it at all).
#
#  4. THE ALLELIC RICHNESS TARGET WAS FRAGILE. hierfstat's automatic choice of
#     how many gene copies to rarefy down to is set by the single WORST
#     locus/population combination anywhere in your ENTIRE dataset -- your
#     own comment on the original noted it came out to just 4 gene copies (2
#     individuals). One unlucky cell can crater the whole analysis.
#     FIXED: below, we choose the rarefaction target ourselves (twice your
#     smallest population's individual count -- a standard, sensible choice)
#     and report exactly how many loci actually reach it per population, so
#     nothing gets silently thrown away without you knowing.
#
#  5. POPULATION LABELS WERE ASSIGNED BY POSITION, WITH NOTHING CHECKING IT.
#     Your original pasted in eleven population names in the order you
#     BELIEVED the samples were stored, with a comment telling yourself to
#     double-check it -- but the code never actually checks anything. Get
#     that order wrong just once and EVERY number this script produces is
#     silently mislabeled: no error, no warning, nothing.
#     FIXED: population labels are now read from the same popmap file (a
#     sample-ID-to-population lookup table) used elsewhere in this project,
#     and matched to your genepop file BY SAMPLE NAME, not position. A
#     mismatch throws a loud error instead of silently corrupting every
#     number downstream.
#
#  Two smaller fixes, worth knowing about even though they're less dramatic:
#    - `colSums(Ho, na.rm=TRUE) / colSums(Hs, na.rm=TRUE)` can, in principle,
#      sum a locus into the numerator (Ho can be defined with just 1 typed
#      individual) while leaving it OUT of the denominator (Hs needs at least
#      2). Fixed by requiring both to be defined before a locus counts at
#      all -- "usable" below.
#    - A "total" computed as (mean across loci) x (total locus count) is only
#      correct if EVERY locus contributed to that mean. Fixed to multiply by
#      each population's own actual contributing-locus count instead.
#
#  WHAT WAS ALREADY RIGHT AND IS UNCHANGED: the core FIS arithmetic itself --
#  summing Ho and Hs across loci FIRST, then dividing ONCE (a "ratio of
#  sums"), rather than averaging a pile of separate per-locus Fis ratios
#  (a "mean of ratios", which is what Stacks itself does and is a well known,
#  separate problem -- not relevant here since you weren't doing that).
# =============================================================================


# =============================================================================
# CONFIG -- this is the ONLY section you should need to edit for a new
# dataset. Everything below reads these variables; nothing else changes.
# =============================================================================
work_dir    <- "/fs/project/whelan.105/virgata/virgata-assembly/post-submission-filtering-tests/virgata_m3M3n3_r80p11_maf025_gt6"
snps_gen    <- "populations.snps.genepop"   # SNP genepop file: Fis, Ho, population comparisons
haps_gen    <- "populations.haps.genepop"   # haplotype genepop file: allelic richness, private alleles
                                             # (Stacks writes ".genepop"; if yours came out ".gen"
                                             #  instead, just change these two lines to match)
popmap_file <- "/fs/project/whelan.105/virgata/virgata-assembly/popmap_virgata.txt"
ncode       <- 2L                       # allele digit-width used when the genepop files were written
nboot       <- 10000L                   # bootstrap replicates (Fis and privAr confidence intervals)
alpha       <- 0.05                     # significance cutoff, AFTER Benjamini-Hochberg correction
seed        <- 1L                       # makes every random draw below reproducible; change or
                                         # remove the set.seed() call for a fresh draw each run

setwd(work_dir)
set.seed(seed)

suppressMessages(library(adegenet))   # reads genepop files
suppressMessages(library(hierfstat))  # Ho/Hs/allelic-richness engine

for (f in c(snps_gen, haps_gen, popmap_file))
  if (!file.exists(f))
    stop("Can't find '", f, "' in ", getwd(), ".\n",
         "  If it's a .gen file, generate it with Stacks' `populations --genepop`\n",
         "  flag from the SAME run as your VCF files. Otherwise fix the path in\n",
         "  the CONFIG block above.")


# =============================================================================
# HELPER 1 (fixes #5): assign population labels BY SAMPLE NAME, not position.
# =============================================================================
# Reads a popmap file (sample_id, population -- no header, tab-separated) and
# returns a factor of population labels, in the SAME ORDER as the samples in
# `gi`, matched by comparing sample ID STRINGS -- not by trusting that the
# genepop file's sample order matches some list you typed in by hand. If any
# sample can't be matched, or if any resulting population ends up with fewer
# than 2 individuals (Fis/He are mathematically undefined at n = 1), this
# stops with a clear error instead of quietly producing wrong numbers.
assign_pops_from_popmap <- function(gi, popmap_file) {
  pm <- read.delim(popmap_file, header = FALSE, sep = "\t",
                    col.names = c("sample", "population"),
                    colClasses = "character", strip.white = TRUE)
  ids <- indNames(gi)
  m <- match(ids, pm$sample)
  missing_ids <- ids[is.na(m)]
  if (length(missing_ids))
    stop("These sample(s) are in the genepop file but have NO matching row in\n",
         "  '", popmap_file, "': ", paste(missing_ids, collapse = ", "), "\n",
         "  Fix the popmap (or check the genepop file's sample names) before\n",
         "  trusting anything this script produces -- a silent mismatch here is\n",
         "  exactly the bug this fix exists to prevent.")
  pop_assigned <- factor(pm$population[m])
  sizes <- table(pop_assigned)
  cat("\nPopulation sizes after matching to the popmap (double-check this looks right):\n")
  print(sizes)
  too_small <- names(sizes)[sizes < 2]
  if (length(too_small))
    stop("Population(s) with fewer than 2 individuals after matching: ",
         paste(too_small, collapse = ", "), "\n",
         "  Fis and He are undefined at n = 1. Fix the popmap or drop that population.")
  pop_assigned
}


# =============================================================================
# HELPER 2 (fixes #2): figure out which RAD TAG each locus/row belongs to.
# =============================================================================
# Several SNPs can sit on the same RAD tag (piece of DNA) and are NOT
# independent measurements of anything -- they were sequenced together, in
# the same individuals, on the same fragment. The bootstraps below need to
# resample whole TAGS, or they will silently treat linked SNPs as if they
# were independent evidence (see fix #2 at the top of this file).
#
# PRIMARY METHOD: Stacks writes every output format (VCF, genepop, ...) from
# ONE pass over the same internal catalog, in the same row order. So the VCF
# sitting right next to your genepop file already tells us, in order, which
# tag each record came from -- its CHROM column literally IS the tag ID. We
# only read that one column (not genotypes), so this is fast even on a big
# VCF, and it was verified against real data this session: on the gt6
# dataset, the VCF's distinct-CHROM count matched the independently-known
# RAD-locus count for that file exactly.
#
# FALLBACK, only used if that VCF isn't found: try to recover the tag ID from
# the genepop LOCUS NAMES themselves, assuming Stacks' usual
# "<tag>_<position>" naming (strip the trailing "_<number>"). This fallback
# was NOT checked against a real Stacks genepop file this session (none
# existed yet to inspect) -- it prints a sanity number (mean SNPs per tag)
# for you to eyeball, and refuses to continue if that number looks
# implausible for a SNP file (should be noticeably more than 1).
get_rad_tags <- function(gen_file, n_records, loc_names) {
  ## Stacks itself writes ".genepop"; some pipelines rename to the shorter
  ## ".gen" adegenet's own docs use -- accept either so this isn't brittle
  ## to that naming choice.
  vcf_guess <- sub("\\.gen(epop)?$", ".vcf", gen_file)
  if (file.exists(vcf_guess)) {
    lines <- readLines(vcf_guess)
    body  <- lines[!startsWith(lines, "#")]
    if (length(body) != n_records)
      stop("Found ", vcf_guess, " but it has ", length(body), " variant record(s); ",
           "'", gen_file, "' has ", n_records, " loci. They must come from the ",
           "SAME `populations` run (same catalog, same filters) for the row-order ",
           "trick above to be valid -- regenerate both files together and re-run.")
    tag <- vapply(strsplit(body, "\t", fixed = TRUE), `[`, character(1), 1)
    message(sprintf(
      "RAD-tag grouping for %s: read from %s (%s distinct tags on %s records, %.2f per tag).",
      gen_file, vcf_guess, format(length(unique(tag)), big.mark = ","),
      format(n_records, big.mark = ","), n_records / length(unique(tag))))
    return(tag)
  }
  message("No sibling VCF ('", vcf_guess, "') found next to '", gen_file, "'.\n",
          "  Falling back to guessing RAD-tag IDs from the genepop locus names\n",
          "  (assuming Stacks' usual '<tag>_<position>' naming). This fallback was\n",
          "  NOT verified against real Stacks output this session -- check the\n",
          "  sanity number below before trusting anything downstream of it.")
  tag <- sub("_[0-9]+$", "", loc_names)
  ratio <- n_records / length(unique(tag))
  message(sprintf("  guessed %s distinct tags on %s records (%.2f per tag)",
                   format(length(unique(tag)), big.mark = ","),
                   format(n_records, big.mark = ","), ratio))
  if (ratio < 1.05)
    stop("That guess puts almost exactly 1 SNP per tag, which is implausible for a\n",
         "  SNP file that should have several SNPs per RAD locus -- the locus-name\n",
         "  guess is almost certainly wrong. Put the matching VCF ('", vcf_guess,
         "')\n  next to '", gen_file, "' so the reliable method above can be used instead.")
  tag
}


# =============================================================================
# HELPER 4: adegenet's read.genepop() insists on a file literally named
# "*.gen" and REFUSES anything else -- confirmed directly this session:
# Stacks itself writes ".genepop", and pointing read.genepop() straight at
# one fails with "File extension .gen expected" before it reads a single
# line. Rather than making you rename or copy your real data files by hand
# every time, this makes a throwaway SYMLINK (same file, just a second name
# ending in ".gen"), reads THROUGH that, and deletes the symlink immediately
# afterward -- your actual file is never touched, copied, or modified. Falls
# back to an actual (temporary) copy only if symlinks aren't supported on
# your filesystem.
# =============================================================================
read_genepop_flexible <- function(path, ncode) {
  if (grepl("\\.gen$", path)) return(read.genepop(path, ncode = ncode))   # already fine as-is
  tmp <- tempfile(fileext = ".gen")
  linked <- tryCatch(file.symlink(normalizePath(path), tmp), error = function(e) FALSE)
  if (!isTRUE(linked)) file.copy(path, tmp)   # fallback: some filesystems don't do symlinks
  on.exit(unlink(tmp), add = TRUE)            # clean up the throwaway link/copy either way
  read.genepop(tmp, ncode = ncode)
}


# =============================================================================
# HELPER 3: make sure `ncode` is actually big enough for this data.
# =============================================================================
# genind2hierfstat() packs a genotype into ONE integer as allele1*10^ncode +
# allele2. If any locus has 100 or more distinct alleles but ncode = 2 (each
# allele given only 2 digits), allele index 100 and allele index 0 become
# indistinguishable in that packed number -- genotypes at that locus would be
# silently misdecoded (wrong answer, no error). Haplotype-based ("haps")
# files are the ones at real risk of this, since one RAD tag can carry many
# distinct haplotypes; plain biallelic SNPs never come close. Check this
# ONCE, right after reading each file, so a bad `ncode` fails loudly here
# instead of quietly corrupting every number downstream.
check_ncode_is_big_enough <- function(gi, ncode) {
  max_alleles <- max(nAll(gi))
  if (max_alleles >= 10 ^ ncode)
    stop("At least one locus has ", max_alleles, " distinct alleles, but ncode = ",
         ncode, " only leaves room for ", 10 ^ ncode - 1, " per allele slot.\n",
         "  Re-export the genepop file with a bigger `ncode` (e.g. 3) and update\n",
         "  the CONFIG block above -- otherwise genotypes at that locus will be\n",
         "  silently misdecoded rather than erroring.")
}


# =============================================================================
# STEP 1-4 (same as your original): read the SNP genepop file, assign
# populations (now via HELPER 1, not by hand), convert to hierfstat's format.
# =============================================================================
gi  <- read_genepop_flexible(snps_gen, ncode)
# If adegenet just printed "Individuals with no scored loci have been
# removed" above, one or more samples had EVERY genotype missing and were
# silently dropped before we ever saw them -- indNames(gi) below will be
# short by that many, and assign_pops_from_popmap() only checks that
# whatever names ARE present have a popmap match, not that nobody's missing.
# Watch for that warning; it will not otherwise raise an error here.
check_ncode_is_big_enough(gi, ncode)
pop(gi) <- assign_pops_from_popmap(gi, popmap_file)
dat <- genind2hierfstat(gi)

# ---- Ho and Hs, with numerator/denominator ALWAYS counting the same loci ----
# hierfstat can leave Ho defined at a locus with just 1 typed individual (you
# can still tell if that ONE person is heterozygous) while leaving Hs
# undefined there (Hs needs at least 2 typed individuals -- its formula
# divides by n - 1). If we summed Ho and Hs separately with na.rm = TRUE, a
# locus like that would silently count in the Ho sum but NOT the Hs sum --
# comparing a total that includes it against a total that doesn't. `usable`
# below fixes this: a locus counts for a population only when BOTH Ho and Hs
# are defined there, matching the same fix built into the RADseq-Claude
# diversity_stats.R this session.
st <- basic.stats(dat)
Ho <- as.matrix(st$Ho)
Hs <- as.matrix(st$Hs)
usable <- is.finite(Ho) & is.finite(Hs)
Ho[!usable] <- NA_real_
Hs[!usable] <- NA_real_
pops <- colnames(Ho)
nloc <- nrow(Ho)


# =============================================================================
# STEP 5-6 (fixed #2): bootstrap over RAD TAGS, not SNP rows.
# =============================================================================
tag       <- get_rad_tags(snps_gen, nloc, locNames(gi))
rows_of_tag <- split(seq_len(nloc), tag)   # tag id -> which SNP rows belong to it
utag      <- unique(tag)
n_tags    <- length(utag)

# Fis from a chosen set of rows: pool Ho and Hs across those rows FIRST, then
# divide ONCE -- a "ratio of sums", never a "mean of ratios". This part of
# your original was already correct; unchanged here.
fis_from_rows <- function(rows) {
  1 - colSums(Ho[rows, , drop = FALSE], na.rm = TRUE) /
      colSums(Hs[rows, , drop = FALSE], na.rm = TRUE)
}
point_fis <- fis_from_rows(seq_len(nloc))

boot_fis <- matrix(NA_real_, nboot, length(pops), dimnames = list(NULL, pops))
for (i in seq_len(nboot)) {
  # pick a random set of TAGS (with replacement) -- every SNP belonging to a
  # chosen tag moves together, so linked SNPs are never split apart and
  # resampled as if they were independent of each other.
  sampled_tags <- sample(utag, n_tags, replace = TRUE)
  rows <- unlist(rows_of_tag[sampled_tags], use.names = FALSE)
  boot_fis[i, ] <- fis_from_rows(rows)
}

ci <- apply(boot_fis, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
summary_tab <- data.frame(Pop = pops, Fis = point_fis, lower = ci[1, ], upper = ci[2, ])
cat("\n--- Per-population Fis with 95% confidence interval (RAD-tag bootstrap) ---\n")
print(summary_tab, row.names = FALSE)


# =============================================================================
# STEP 7-8 (same as your original -- this logic was already correct, it just
# now runs on top of the corrected, RAD-tag-based bootstrap above).
# =============================================================================
pairs <- combn(pops, 2)   # every possible pair of population names

pair_res <- apply(pairs, 2, function(pr) {
  # Comparing the bootstrap difference WITHIN each replicate (rather than
  # comparing two separate confidence intervals) correctly accounts for both
  # populations being resampled from the SAME tags each time.
  d <- boot_fis[, pr[1]] - boot_fis[, pr[2]]
  p_low  <- (sum(d <= 0) + 1) / (nboot + 1)
  p_high <- (sum(d >= 0) + 1) / (nboot + 1)
  p_raw  <- min(1, 2 * min(p_low, p_high))
  data.frame(PopA = pr[1], PopB = pr[2],
             diff  = point_fis[pr[1]] - point_fis[pr[2]],
             p_raw = p_raw)
})
pair_res <- do.call(rbind, pair_res)

# BH correction: we ran one test per pair, so some will look "significant"
# just by chance -- this controls the expected proportion of false positives
# among everything we call significant.
pair_res$p_adj <- p.adjust(pair_res$p_raw, method = "BH")
pair_res$significant <- pair_res$p_adj < alpha
pair_res <- pair_res[order(pair_res$p_adj), ]

cat("\n--- Pairwise Fis comparisons, BH-corrected ---\n")
print(pair_res, row.names = FALSE, digits = 4)
cat("\nSignificant pairs after BH correction:", sum(pair_res$significant),
    "out of", nrow(pair_res), "\n")


# =============================================================================
# STEP 9 (REPLACED -- fixes #1): does heterozygosity differ between
# populations?
# =============================================================================
# Your original tested this treating every LOCUS as its own independent data
# point (Welch ANOVA, Friedman test). That's the problem explained at the top
# of this file: the SAME individuals are measured at every locus, so loci
# aren't independent replicates of anything -- INDIVIDUALS are. This section
# computes ONE heterozygosity number per INDIVIDUAL (the fraction of THAT
# PERSON'S OWN genotyped loci where they're heterozygous), then compares
# populations using individuals as the unit -- the design measured (in this
# project's het_between_pops.R self-test) to hold the nominal 5%
# false-positive rate, versus 55-76% for the locus-based version.
#
# Decoding a genotype code back into its two alleles: genind2hierfstat()
# stores a diploid genotype as ONE integer, allele1 * 10^ncode + allele2
# (checked directly this session: with ncode = 2, "01/02" becomes exactly
# 102, and a fully-missing genotype, code "0000", correctly becomes R's NA --
# not a fake "0/0" homozygote -- also checked directly). Integer division and
# remainder recover the two alleles; if either allele was missing, the whole
# code is NA, and both operations correctly stay NA -- nothing needs to be
# handled specially.
geno      <- as.matrix(dat[, -1, drop = FALSE])   # drop the "pop" column
storage.mode(geno) <- "integer"   # %/% and %% need a number, not text; belt-and-
                                   # braces in case genind2hierfstat ever hands
                                   # back character codes instead of integers
divisor   <- 10 ^ ncode
a1        <- geno %/% divisor
a2        <- geno %%  divisor
het_mat   <- (a1 != a2)          # TRUE = heterozygous at that locus, for that individual
                                 # (NA carries through automatically for missing genotypes)

ind_het   <- rowMeans(het_mat, na.rm = TRUE)   # this person's own heterozygosity
ind_ncall <- rowSums(!is.na(het_mat))          # how many loci this person was actually typed at
ind_pop   <- dat$pop                           # already correctly assigned, via the popmap fix

# Someone typed at almost nothing gives an unreliable heterozygosity number,
# and rowMeans() of an all-missing row is NaN -- which would then silently
# poison every test that touches it. Exclude them LOUDLY instead.
too_thin <- ind_ncall < 50
if (any(too_thin)) {
  message(sprintf("Excluding %d individual(s) typed at fewer than 50 loci: %s",
                   sum(too_thin), paste(indNames(gi)[too_thin], collapse = ", ")))
  ind_het   <- ind_het[!too_thin]
  ind_ncall <- ind_ncall[!too_thin]
  ind_pop   <- droplevels(ind_pop[!too_thin])
}

cat("\n--- Per-individual heterozygosity, summarised by population ---\n")
ind_summary <- do.call(rbind, lapply(levels(ind_pop), function(p) {
  v <- ind_het[ind_pop == p]
  data.frame(Pop = p, n = length(v), mean_Ho = round(mean(v), 4),
             sd = round(sd(v), 4), se = round(sd(v) / sqrt(length(v)), 4))
}))
print(ind_summary, row.names = FALSE)

# Missingness-confound check: if individuals with MORE missing data also look
# LESS heterozygous, you might be measuring library quality, not biology
# (classic allele-dropout signature) -- worth ruling out before trusting the
# pairwise tests below.
if (sd(ind_ncall) > 0) {
  cr <- ind_ncall / ncol(het_mat)   # ncol(het_mat), not nloc: same number here, but this
                                    # keeps the ratio's denominator tied to the matrix it
                                    # actually came from rather than a separately-computed value
  ct <- suppressWarnings(cor.test(cr, ind_het))
  cat(sprintf("\ncor(per-individual call rate, heterozygosity) = %+.3f (p = %.3g)\n",
              unname(ct$estimate), ct$p.value))
  if (is.finite(ct$p.value) && ct$p.value < 0.05 && ct$estimate > 0)
    cat("  POSITIVE and significant: less-genotyped individuals look LESS\n",
        "  heterozygous -- the classic allele-dropout signature. Treat the\n",
        "  pairwise results below with extra caution until you've checked\n",
        "  whether this tracks specific populations.\n", sep = "")
  else
    cat("  No worrying positive association -- heterozygosity doesn't appear\n",
        "  to just be tracking who got sequenced better.\n", sep = "")
}

pairs2 <- combn(levels(ind_pop), 2)
ind_pair_res <- apply(pairs2, 2, function(pr) {
  a <- ind_het[ind_pop == pr[1]]; b <- ind_het[ind_pop == pr[2]]
  tt <- t.test(a, b)                        # Welch's t-test: doesn't assume equal variance
  wt <- suppressWarnings(wilcox.test(a, b))  # distribution-free check, no normality assumption
  data.frame(PopA = pr[1], PopB = pr[2], diff = round(mean(a) - mean(b), 4),
             p_welch = tt$p.value, p_wilcox = wt$p.value)
})
ind_pair_res <- do.call(rbind, ind_pair_res)
ind_pair_res$p_welch_BH  <- p.adjust(ind_pair_res$p_welch,  method = "BH")
ind_pair_res$p_wilcox_BH <- p.adjust(ind_pair_res$p_wilcox, method = "BH")
ind_pair_res$significant <- ind_pair_res$p_welch_BH < alpha
ind_pair_res <- ind_pair_res[order(ind_pair_res$p_welch_BH), ]

cat("\n--- Pairwise Ho comparisons, INDIVIDUAL as the replicate, BH-corrected ---\n")
print(ind_pair_res, row.names = FALSE, digits = 4)
cat("\nSignificant pairs after BH correction:", sum(ind_pair_res$significant),
    "out of", nrow(ind_pair_res), "\n")

# A NOTE ON YOUR ORIGINAL STEP 10 (paired Ho-vs-He per population): that
# question -- "is this population's heterozygosity different from what
# random mating predicts?" -- is exactly what Fis asks, and the corrected Fis
# confidence interval up in Step 5-6 already answers it correctly (and
# accounts for linked SNPs, via the RAD-tag bootstrap, which a plain paired
# test on raw loci does not). Re-answering the same question here with a
# less careful test would just add a second, weaker version of it, so it's
# been removed rather than patched.


# =============================================================================
# PART 3-4: allelic richness and private alleles, from the HAPLOTYPE file
# =============================================================================
gi_haps  <- read_genepop_flexible(haps_gen, ncode)
# Haplotype files are where this actually matters (see HELPER 3 above): one
# RAD tag can carry many distinct haplotypes, unlike a plain biallelic SNP.
check_ncode_is_big_enough(gi_haps, ncode)
pop(gi_haps) <- assign_pops_from_popmap(gi_haps, popmap_file)
dat_haps <- genind2hierfstat(gi_haps)
pops_h     <- levels(pop(gi_haps))
pop_sizes  <- table(pop(gi_haps))


# =============================================================================
# PART 3 (fixed #4): rarefied allelic richness, with an EXPLICIT target.
# =============================================================================
# hierfstat's automatic target is set by the single WORST locus/population
# combination anywhere in your ENTIRE dataset -- one unlucky cell can crater
# it for everyone (your own comment on the original: it came out to just 4
# gene copies, i.e. 2 individuals). Instead, we choose the target ourselves:
# twice your smallest population's individual count -- a standard, sensible
# choice -- and report exactly how many loci reach it per population, so
# nothing gets silently thrown away without you seeing it happen.
g <- 2 * min(pop_sizes)
cat(sprintf("\nRarefying allelic richness to g = %d gene copies (2x smallest population, n=%d)\n",
            g, min(pop_sizes)))
ar     <- allelic.richness(dat_haps, min.n = g)
ar_mat <- as.matrix(ar$Ar)
ar_n   <- colSums(!is.na(ar_mat))
ar_summary <- data.frame(Pop = pops_h, Ar = round(colMeans(ar_mat, na.rm = TRUE), 4),
                          Ar_n = ar_n, Ar_total_loci = nrow(ar_mat))
cat("\n--- Rarefied allelic richness (per locus) ---\n")
print(ar_summary, row.names = FALSE)
cat("Ar_n: how many loci actually had >= g gene copies for THAT population and\n",
    "  contributed to its Ar average (out of Ar_total_loci available).\n", sep = "")
if (min(ar_n) < 0.5 * nrow(ar_mat))
  cat("WARNING: at least one population's Ar rests on fewer than half its loci.\n",
      "  Consider a smaller g.\n", sep = "")


# =============================================================================
# PART 4 (fixed #3): rarefied private allele richness, done ANALYTICALLY --
# no simulation, no bias.
# =============================================================================
# "Private" = an allele found in only ONE of your populations. A population
# with MORE individuals sampled gets more chances to turn up a rare allele
# nobody else has, purely by luck -- raw private-allele COUNTS aren't
# comparable across differently-sized populations for exactly that reason.
# Rarefaction fixes it: mathematically ask "if I'd only drawn g GENE COPIES
# from this population, how many private alleles would I expect to see?"
#
# Your original did this by repeatedly drawing random INDIVIDUALS (without
# replacement) and recounting -- but the two gene copies inside ONE
# individual are only an independent draw from the population if Fis happens
# to be exactly 0. This dataset's Fis is NOT 0 (often clearly negative, i.e.
# MORE heterozygotes than random mating predicts), so drawing individuals
# (2 linked copies at a time) instead of drawing gene copies directly gives a
# biased answer.
#
# The fix: the standard ANALYTIC formula (the same one HP-RARE and ADZE use)
# that computes the EXACT expected count for a random draw of g gene copies
# directly. No simulation, no bias, and it's deterministic -- run it twice,
# get the same answer both times. The formulas (Hurlbert 1971; Kalinowski
# 2004; Szpiech et al. 2008) are short enough to just write out here.

# p_sampled: for each allele (given as a vector of how many gene copies of
# each allele exist in the whole population, e.g. c(30, 12, 2) for 3
# alleles), what's the chance THAT allele shows up at least once if you draw
# g gene copies at random, WITHOUT replacement, from the whole pool? Uses
# lchoose ("log of n-choose-k") rather than choose() directly, because the
# raw counts can get big enough to overflow an ordinary number -- logs keep
# the arithmetic stable no matter how big the population is.
p_sampled <- function(counts, g) {
  N <- sum(counts)
  if (!is.finite(N) || g > N || g < 1) return(rep(NA_real_, length(counts)))
  1 - exp(lchoose(N - counts, g) - lchoose(N, g))
}

# rare_private: for population j at ONE locus, add up, across every allele,
# the chance that allele (a) shows up in population j's random draw of g
# copies, AND (b) does NOT show up in ANY other population's random draw of g
# copies. That's the textbook definition of "private", applied at a fixed,
# fair sample size for everyone being compared.
# count_mat = populations (rows) x alleles (columns) of raw gene-copy counts
# at this one locus.
rare_private <- function(count_mat, j, g) {
  if (!is.matrix(count_mat)) count_mat <- rbind(count_mat)
  ps <- matrix(NA_real_, nrow(count_mat), ncol(count_mat))
  for (rr in seq_len(nrow(count_mat))) ps[rr, ] <- p_sampled(count_mat[rr, ], g)
  # If ANY population lacks g gene copies at this locus, "private at sample
  # size g" isn't well defined for ANYONE here (not just the short
  # population) -- this locus contributes NA for everyone, not a silent
  # zero, so it's correctly excluded from every population's average below.
  if (anyNA(ps)) return(NA_real_)
  term <- ps[j, ]
  for (k in setdiff(seq_len(nrow(count_mat)), j)) term <- term * (1 - ps[k, ])
  sum(term)
}
rare_private_all <- function(count_mat, g)
  vapply(seq_len(nrow(count_mat)), function(j) rare_private(count_mat, j, g), numeric(1))

# ---- build, for every locus, a (population x allele) table of how many gene
# copies of each allele each population actually has. Reuses the same
# genotype-code decoding trick as Step 9 above. ----
geno_h    <- as.matrix(dat_haps[, -1, drop = FALSE])
storage.mode(geno_h) <- "integer"   # same belt-and-braces as the SNP decode above
divisor_h <- 10 ^ ncode
n_loc_h   <- ncol(geno_h)
pop_h     <- dat_haps$pop

cmats <- vector("list", n_loc_h)     # one (population x allele) count matrix per locus
for (j in seq_len(n_loc_h)) {
  A1 <- geno_h[, j] %/% divisor_h
  A2 <- geno_h[, j] %%  divisor_h
  alleles_here <- sort(unique(c(A1, A2)))
  alleles_here <- alleles_here[!is.na(alleles_here)]
  m <- matrix(0L, length(pops_h), length(alleles_here), dimnames = list(pops_h, NULL))
  for (i in seq_along(pops_h)) {
    v <- c(A1[pop_h == pops_h[i]], A2[pop_h == pops_h[i]])
    v <- v[!is.na(v)]
    m[i, ] <- tabulate(match(v, alleles_here), nbins = length(alleles_here))
  }
  cmats[[j]] <- m
}

pr_point <- matrix(NA_real_, n_loc_h, length(pops_h), dimnames = list(NULL, pops_h))
for (j in seq_len(n_loc_h)) pr_point[j, ] <- rare_private_all(cmats[[j]], g)

# RAD-tag block bootstrap for privAr's confidence interval, reusing HELPER 2
# from up top -- in the haplotype file every row is already its own RAD
# locus (one haplotype per tag), so this mostly double-checks that this file
# and its sibling VCF are actually in sync with each other.
tag_h        <- get_rad_tags(haps_gen, n_loc_h, locNames(gi_haps))
rows_of_tag_h <- split(seq_len(n_loc_h), tag_h)
utag_h       <- unique(tag_h)

boot_pr <- matrix(NA_real_, nboot, length(pops_h), dimnames = list(NULL, pops_h))
for (i in seq_len(nboot)) {
  sampled <- sample(utag_h, length(utag_h), replace = TRUE)
  rows    <- unlist(rows_of_tag_h[sampled], use.names = FALSE)
  boot_pr[i, ] <- colMeans(pr_point[rows, , drop = FALSE], na.rm = TRUE)
}
pr_ci <- t(apply(boot_pr, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
pr_n  <- colSums(!is.na(pr_point))

# raw (NOT rarefied) private allele count -- kept from your original, which
# was already fine: a plain count has no "mean of ratios" or subsampling
# issue to fix.
count_private <- function(freqs, pops) {
  pc <- setNames(rep(0L, length(pops)), pops)
  for (locus_mat in freqs) {
    present <- locus_mat > 0
    npres   <- rowSums(present, na.rm = TRUE)
    priv    <- which(npres == 1)
    for (a in priv) {
      pop_with_it <- colnames(present)[which(present[a, ])]
      pc[pop_with_it] <- pc[pop_with_it] + 1
    }
  }
  pc
}
raw_private <- count_private(pop.freq(dat_haps), pops_h)

# priv_total: the dataset-wide total, in gene-copy-rarefied units.
# colMeans(..., na.rm=TRUE) above already divides by however many loci were
# NON-NA for THAT population -- so the correct total is that mean times THAT
# population's own non-NA locus count (pr_n), never the dataset's total
# locus count. Multiplying by the wrong (larger) count silently inflates
# whichever population happens to have more missing loci -- exactly
# backwards from what you'd want.
privAr     <- colMeans(pr_point, na.rm = TRUE)
priv_total <- round(privAr * pr_n, 1)

private_summary <- data.frame(Pop = pops_h, Raw = raw_private,
                                privAr = round(privAr, 4),
                                privAr_lo = round(pr_ci[, 1], 4),
                                privAr_hi = round(pr_ci[, 2], 4),
                                privAr_n = pr_n, priv_total = priv_total)
cat("\n--- Private allele richness: raw count vs analytically rarefied ---\n")
print(private_summary, row.names = FALSE)
cat("privAr_n: how many of this file's ", n_loc_h, " loci had EVERY population at\n",
    "  >= g gene copies simultaneously (private richness needs everyone comparable\n",
    "  at the same sample size at once, so this is often less than Ar_n above --\n",
    "  though the two aren't quite the same statistic: Ar_n comes from hierfstat's\n",
    "  own allelic.richness(), privAr_n from the rare_private() code just above,\n",
    "  so treat this as 'in the same ballpark', not 'exactly comparable').\n", sep = "")
if (min(pr_n) < 0.5 * n_loc_h)
  cat(sprintf(paste0("WARNING: privAr_n is below 50%% of loci for at least one population --\n",
                      "  its rarefied private richness rests on a minority of loci. Consider\n",
                      "  a smaller g (currently %d gene copies).\n"), g))

cat("\nDone.\n")
