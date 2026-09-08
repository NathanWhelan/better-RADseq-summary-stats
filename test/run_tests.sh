#!/bin/bash
# Full test suite. Run from anywhere: bash test/run_tests.sh (or cd test &&
# bash run_tests.sh). Fixtures, golden files, and this script live in
# test/; the scripts under test (diversity_stats.R, diversity_core.R,
# haps_common.R, het_between_pops.R) live one directory up, in the repo root.
# Regenerates fixtures, runs both self-tests, the label guard, a numeric
# regression check against known-good output, and every edge case the
# scripts are expected to reject.
set -u
cd "$(dirname "$0")"
R=".."
pass=0; fail=0
run () { local label="$1" expect="$2"; shift 2
  out=$("$@" 2>&1); rc=$?
  if { [ "$expect" = ok ] && [ $rc -eq 0 ]; } || { [ "$expect" = err ] && [ $rc -ne 0 ]; }
  then pass=$((pass+1))
  else echo "  FAIL  $label (rc=$rc)"; echo "$out" | tail -3 | sed 's/^/          /'
       fail=$((fail+1)); fi
}
## Unlike run(), which only checks exit status, this diffs a script's TSV
## output against a checked-in golden copy (check_golden.R) -- catching a
## script that runs to completion and prints a WRONG number, which an
## exit-code check cannot.
run_golden () { local label="$1" actual="$2" golden="$3"
  out=$(Rscript check_golden.R "$actual" "$golden" 2>&1); rc=$?
  if [ $rc -eq 0 ]; then pass=$((pass+1))
  else echo "  FAIL  $label"; echo "$out" | sed 's/^/          /'
       fail=$((fail+1)); fi
}
D="Rscript $R/diversity_stats.R"; H="Rscript $R/het_between_pops.R"

echo "regenerating fixtures"; Rscript make_fixtures.R >/dev/null || exit 1

echo "self-tests"
run "diversity_core --selftest"     ok  Rscript $R/diversity_core.R --selftest
run "het_between_pops --selftest"   ok  Rscript $R/het_between_pops.R --selftest
run "output labels match columns"   ok  Rscript check_labels.R

echo "golden-value regression (pins numeric output; exit code alone can't catch a wrong-but-successful run)"
## --g=20 reproduces the old default (twice the smallest population, 10, in
## popmap.tsv) so the golden files stay valid without regenerating them.
$D sim.allsnps.vcf popmap.tsv --nboot=0 --sites=0 --g=20 >/dev/null 2>&1
run_golden "allsnps per-population vs golden" diversity_per_population.allsnps.tsv golden/diversity_per_population.allsnps.tsv
run_golden "allsnps richness vs golden"       diversity_richness.allsnps.tsv       golden/diversity_richness.allsnps.tsv
$D sim.haps.vcf popmap.tsv --nboot=0 --sites=0 --g=20 >/dev/null 2>&1
run_golden "haps per-population vs golden"    diversity_per_population.haps.tsv    golden/diversity_per_population.haps.tsv
run_golden "haps richness vs golden"          diversity_richness.haps.tsv          golden/diversity_richness.haps.tsv
$D miss10.vcf popmap.tsv --nboot=0 --sites=0 --g=20 >/dev/null 2>&1
run_golden "miss10 per-population vs golden"  diversity_per_population.miss10.tsv  golden/diversity_per_population.miss10.tsv
run_golden "miss10 richness vs golden"        diversity_richness.miss10.tsv        golden/diversity_richness.miss10.tsv
$H sim.allsnps.vcf popmap.tsv 0.9 /tmp/hbp_golden_allsnps >/dev/null 2>&1
run_golden "hbp allsnps individual het vs golden" /tmp/hbp_golden_allsnps/individual_heterozygosity.tsv golden/hbp_allsnps.individual_heterozygosity.tsv
run_golden "hbp allsnps pairwise tests vs golden" /tmp/hbp_golden_allsnps/het_between_pops_tests.tsv    golden/hbp_allsnps.het_between_pops_tests.tsv
$H one_dead_ind.vcf popmap.tsv 0.9 /tmp/hbp_golden_onedead >/dev/null 2>&1
run_golden "hbp one_dead_ind individual het vs golden" /tmp/hbp_golden_onedead/individual_heterozygosity.tsv golden/hbp_onedead.individual_heterozygosity.tsv
run_golden "hbp one_dead_ind pairwise tests vs golden" /tmp/hbp_golden_onedead/het_between_pops_tests.tsv    golden/hbp_onedead.het_between_pops_tests.tsv

echo "arguments"
run "no args"                       err Rscript $R/diversity_stats.R
run "one arg"                       err $D sim.allsnps.vcf
run "missing --g"                   err $D sim.allsnps.vcf popmap.tsv
run "missing vcf"                   err $D nope.vcf popmap.tsv --g=20
run "missing popmap"                err $D sim.allsnps.vcf nope.tsv --g=20
run "unrecognized flag"             err $D sim.allsnps.vcf popmap.tsv --g=20 --bogus=1
run "bare --g (no =value)"          err $D sim.allsnps.vcf popmap.tsv --g 20
run "leftover positional"           err $D sim.allsnps.vcf popmap.tsv --g=20 extra
run "g > smallest population"       err $D sim.allsnps.vcf popmap.tsv --nboot=0 --sites=0 --g=999
run "g = 1"                         err $D sim.allsnps.vcf popmap.tsv --nboot=0 --sites=0 --g=1
run "nboot = 0"                     ok  $D sim.allsnps.vcf popmap.tsv --nboot=0 --g=20
run "sites < n records"             ok  $D sim.allsnps.vcf popmap.tsv --nboot=0 --sites=10 --g=20
run "gzipped vcf"                   ok  $D sim.allsnps.vcf.gz popmap.tsv --nboot=0 --g=20

echo "popmaps"
run "duplicate sample"              err $D sim.allsnps.vcf pm_dup.tsv --nboot=0 --g=20
run "no sample matches"             err $D sim.allsnps.vcf pm_nomatch.tsv --nboot=0 --g=2
run "single population"             err $D sim.allsnps.vcf pm_onepop.tsv --nboot=0 --g=2
run "population of one"             err $D sim.allsnps.vcf pm_n1.tsv --nboot=0 --g=2
run "popmap subsets vcf"            ok  $D sim.allsnps.vcf pm_subset.tsv --nboot=0 --g=16
run "numeric population names"      ok  $D sim.allsnps.vcf pm_numeric.tsv --nboot=0 --g=20
run "underscore in name"            ok  $D sim.allsnps.vcf pm_underscore.tsv --nboot=0 --g=20
run "three populations"             ok  $D sim.allsnps.vcf pm_three.tsv --nboot=0 --g=16

echo "vcf pathologies"
run "no #CHROM line"                err $D bad_nohdr.vcf popmap.tsv --nboot=0 --g=20
run "no records"                    err $D bad_norec.vcf popmap.tsv --nboot=0 --g=20
run "all genotypes missing"         err $D bad_allmiss.vcf popmap.tsv --nboot=0 --g=20
run "nothing survives --complete-case" err $D bad_nocomplete.vcf popmap.tsv --nboot=0 --g=20 --complete-case
run "available-data survives 1-missing-per-record" ok $D bad_nocomplete.vcf popmap.tsv --nboot=0 --g=20
run "monomorphic everywhere"        err $D bad_mono.vcf popmap.tsv --nboot=0 --g=20
run "10% missing"                   ok  $D miss10.vcf popmap.tsv --nboot=0 --g=20
run "10% missing, --complete-case"  ok  $D miss10.vcf popmap.tsv --nboot=0 --g=20 --complete-case

echo "available-data mode"
run "locus-driven missingness, default"        ok  $D locus_driven_miss.vcf popmap.tsv --nboot=0 --g=20
run "locus-driven missingness, --complete-case" ok  $D locus_driven_miss.vcf popmap.tsv --nboot=0 --g=20 --complete-case
run "min_n = 3"                                 ok  $D locus_driven_miss.vcf popmap.tsv --nboot=0 --sites=0 --g=20 --min-n=3
run "min_n = 1 rejected"                        err $D locus_driven_miss.vcf popmap.tsv --nboot=0 --sites=0 --g=20 --min-n=1

echo "file types and the sites path"
run "haps + sites (guard fires)"    ok  $D sim.haps.vcf popmap.tsv --nboot=0 --sites=1000000 --g=20
run "snps + sites"                  ok  $D sim.allsnps.vcf popmap.tsv --nboot=0 --sites=1000000 --g=20
run "one-snp-per-tag typed as SNP"  ok  $D sim.onesnp.vcf popmap.tsv --nboot=0 --sites=1000000 --g=20

echo "het_between_pops"
run "hbp no args"                   err $H
run "hbp min_call out of range"     err $H sim.allsnps.vcf popmap.tsv 5
run "hbp on snps"                   ok  $H sim.allsnps.vcf popmap.tsv 0.9 /tmp/hbp1
run "hbp on haps"                   ok  $H sim.haps.vcf popmap.tsv 0.9 /tmp/hbp2
run "hbp three populations"         ok  $H sim.allsnps.vcf pm_three.tsv 0.9 /tmp/hbp3
run "hbp excludes a dead library"   ok  $H one_dead_ind.vcf popmap.tsv 0.9 /tmp/hbp4

echo
echo "passed $pass, failed $fail"
[ $fail -eq 0 ]
