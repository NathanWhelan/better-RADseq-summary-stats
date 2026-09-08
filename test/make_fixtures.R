#!/usr/bin/env Rscript
## Build the synthetic VCFs and popmaps the test suite runs against.
## Two populations with a known FIS, tags carrying 1-4 SNPs, plus a set of
## pathological files that every guard in the scripts should reject.
set.seed(4242)
n1 <- 15; n2 <- 10; ntag <- 400; Fis <- 0.10
samples <- c(sprintf("A%02d", 1:n1), sprintf("B%02d", 1:n2))
pm <- data.frame(s = samples, p = rep(c("popA", "popB"), c(n1, n2)))
write.table(pm, "popmap.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
nuc <- c("A", "C", "G", "T")
snp <- hap <- list()

for (t in seq_len(ntag)) {
  k  <- sample(1:4, 1, prob = c(.45, .30, .17, .08))
  hm <- unique(matrix(rbinom(max(2, sample(2:min(2^k, 6), 1)) * k, 1, .5),
                      ncol = k))
  if (nrow(hm) < 2) next
  nh <- nrow(hm); hf <- rgamma(nh, 1); hf <- hf / sum(hf)
  G <- matrix(0L, n1 + n2, k); H1 <- H2 <- integer(n1 + n2)
  for (i in seq_len(n1 + n2)) {
    a <- sample.int(nh, 1, prob = hf)
    b <- if (runif(1) < Fis) a else sample.int(nh, 1, prob = hf)
    G[i, ] <- hm[a, ] + hm[b, ]; H1[i] <- a; H2[i] <- b
  }
  gs <- matrix(c("0/0", "0/1", "1/1")[G + 1L], nrow(G), ncol(G))
  for (s in seq_len(k))
    snp[[length(snp) + 1]] <- paste(c("un", t * 100 + s, paste0(t, ":", s, ":+"),
      "A", "C", ".", "PASS", "NS=25", "GT", gs[, s]), collapse = "\t")
  alle <- apply(hm, 1, function(r) paste(nuc[r + 1L], collapse = ""))
  alle <- make.unique(alle, sep = "")
  hap[[length(hap) + 1]] <- paste(c("un", t, as.character(t), alle[1],
    if (nh > 1) paste(alle[-1], collapse = ",") else ".", ".", "PASS", "NS=25",
    "GT", paste0(pmin(H1, H2) - 1L, "/", pmax(H1, H2) - 1L)), collapse = "\t")
}

hdr <- c("##fileformat=VCFv4.2", "##source=synthetic",
         paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                 "FORMAT", samples), collapse = "\t"))
writeLines(c(hdr, unlist(snp)), "sim.allsnps.vcf")
writeLines(c(hdr, unlist(hap)), "sim.haps.vcf")
writeLines(c(hdr, unlist(snp)[!duplicated(sub("^un\t[0-9]+\t([0-9]+):.*", "\\1",
            unlist(snp)))]), "sim.onesnp.vcf")

## pathological fixtures
rec <- unlist(snp)
sub_gt <- function(r, f) { x <- strsplit(r, "\t")[[1]]
  x[10:length(x)] <- f(x[10:length(x)]); paste(x, collapse = "\t") }
writeLines(rec[1:50],                                        "bad_nohdr.vcf")
writeLines(hdr,                                              "bad_norec.vcf")
writeLines(c(hdr, sapply(rec[1:50], sub_gt, function(g) rep("./.", length(g)))),
           "bad_allmiss.vcf")
writeLines(c(hdr, sapply(rec[1:50], sub_gt, function(g) rep("0/0", length(g)))),
           "bad_mono.vcf")
writeLines(c(hdr, sapply(rec[1:50], sub_gt,
           function(g) { g[sample.int(length(g), 1)] <- "./."; g })),
           "bad_nocomplete.vcf")
writeLines(c(hdr, sapply(rec, sub_gt,
           function(g) ifelse(runif(length(g)) < .10, "./.", g))), "miss10.vcf")
writeLines(c(hdr, sapply(rec, sub_gt, function(g) { g[4] <- "./."; g })),
           "one_dead_ind.vcf")

## LOCUS-DRIVEN missingness, as opposed to individual-driven: 20% of loci each
## drop a random 40% of popA (indices 1:n1) while popB stays fully genotyped,
## and no single individual is disproportionately missing. This is the
## pattern found in real RAD data (per-individual call rates tight, dropout
## concentrated by LOCUS) that --complete-case guts almost entirely -- every
## degraded locus fails "genotyped in every individual" even though most of
## popA and all of popB are still typed there -- but available-data mode
## should retain heavily (popA keeps well over min_n = 2 at every degraded
## locus; popB is untouched).
set.seed(55)
degraded <- runif(length(rec)) < 0.20
writeLines(c(hdr, mapply(function(r, deg) {
              if (!deg) return(r)
              sub_gt(r, function(g) { g[sample.int(n1, floor(0.4 * n1))] <- "./."; g })
            }, rec, degraded, SIMPLIFY = TRUE)),
           "locus_driven_miss.vcf")

w <- function(f, d) write.table(d, f, sep = "\t", quote = FALSE,
                                row.names = FALSE, col.names = FALSE)
w("pm_dup.tsv",       rbind(pm, pm[1, ]))
w("pm_nomatch.tsv",   data.frame(s = c("ZZ1", "ZZ2"), p = c("a", "b")))
w("pm_onepop.tsv",    data.frame(s = pm$s, p = "only"))
w("pm_n1.tsv",        data.frame(s = pm$s, p = c("solo", rep("rest", nrow(pm) - 1))))
w("pm_subset.tsv",    pm[c(1:8, 17:nrow(pm)), ])
w("pm_numeric.tsv",   data.frame(s = pm$s, p = ifelse(pm$p == "popA", 1, 2)))
w("pm_underscore.tsv",data.frame(s = pm$s, p = ifelse(pm$p == "popA",
                                                      "pop_north", "pop_south")))
w("pm_three.tsv",     data.frame(s = pm$s, p = c("g1", "g2", "g3")[
                                   seq_len(nrow(pm)) %% 3 + 1]))
R.utils <- NULL
con <- gzfile("sim.allsnps.vcf.gz", "w")
writeLines(readLines("sim.allsnps.vcf"), con); close(con)
cat("fixtures written\n")
