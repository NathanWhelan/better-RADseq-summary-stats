# Session log — RADseq diversity pipeline: bootstrap CIs and the `--boot` flag

**Date:** 2026-09-07 to 2026-09-08
**Repo:** `/home/whelan.105/RADseq-Claude/` (not a git repo)
**Files touched:** `diversity_stats.R`, `het_between_pops.R`, `MINIMAL_WORKFLOW.md`

This is a written record of a long working session, reconstructed from the
conversation. Early parts are summarized (the original transcript was
compacted mid-session); everything from "Plan Mode, round 2" onward is
reproduced in full detail since it was directly observed.

---

## 1. How it started

**User:** "Does the pipeline in this folder generate confidence intervals
for all metrics like Ho, He?"

Investigation found: Fis, Ar, and privAr already had bootstrap confidence
intervals (`_lo`/`_hi` columns) from a pre-existing block-bootstrap over RAD
loci. Ho and He did not.

**User:** "Add them and incorporate it into screen output and file output.
Add info to readme and workflow documents as appropriate. But don't go
overboard and make things more complicated than they need to be."

→ Added `Ho_lo`/`Ho_hi`/`He_lo`/`He_hi` columns to `diversity_stats.R`,
matching the existing Fis/Ar/privAr pattern exactly (same bootstrap
machinery, same table format). Updated `MINIMAL_WORKFLOW.md` accordingly.

## 2. Stress-testing and doc review

**User:** "Does Ar and private alleles get confidence intervals? If not,
add them like you did for Ho and He. Test all code for bugs. Stress test.
Deeply study md files to identify issues. Then let me know if md files
need editing to describe everything CORRECTLY; do not make edits on md
files, but tell me what changes you recommend for md files"

Found: Ar/privAr already had CIs (see above); `priv_total` (the
dataset-wide private-allele sum) did not. Stress testing turned up two real
bugs:
- **Bug A**: a crash when Ar is undefined for every population (no locus
  reaches the rarefaction target `g`) — the code assumed at least one
  finite Ar value existed.
- **Bug C**: a cosmetic `NaN%`/`Inf%` in a He-comparison table when a
  population is monomorphic everywhere (He = 0).
- A similar NaN bug in `het_between_pops.R`'s overdispersion table for a
  population with zero variance in heterozygosity.
- Several `MINIMAL_WORKFLOW.md` inaccuracies (wrong columns claimed to
  have CIs, a stale filename reference, wrong output-file count, a
  mis-scoped "Denominators" paragraph, a `pi_2n_corr`/`He_2n_corr` naming
  inconsistency).

**User:** "fix but ask me for clarifying questions if you're uncertain"

→ Fixed Bug A, Bug C, the `het_between_pops.R` NaN, and all the doc issues
(two of them via `AskUserQuestion` since they were genuinely ambiguous).

## 3. The `priv_total` question, and the start of Plan Mode

**User:** "patch it and deal with priv_total's missing CI (e.g., does it
even make sense to put a CI on private alleles? Some populations will have
no private alleles, and that makes sense biologically)"

This opened **Plan Mode** and a long back-and-forth (many rounds, each one
catching a real gap):

1. *"Double check that private allele bootstrapping and CI. Does it make
   sense to bootstrap individuals or loci?"* — Prompted research into
   whether resampling loci or individuals is the right axis for a
   bootstrap CI at all.
2. *"Are there any papers or studies you can find that would suggest what
   is more standard in the field...?"* — Literature search surfaced Van
   Dongen (1995), Petit & Pons (1998), `diveRsity::divBasic` (verified
   against its own source).
3. *"OK. recode so bootstrapping to do bootstrapping by individuals for Ar
   and private alleles. Should we reconsider Ho, He and Fis bootstrapping?"*
4. *"I don't want the time or effort to affect testing whether Ho and He
   CIs should be based on loci or individuals... It sort of seems like
   everything should be consistent."*
5. *"I've never heard of bootstrapping two things and weighting. is that
   statistically valid? ... can you keep the current bootstrapping for Ho
   and He as an option so I can compare it... against the field's
   standard, hierfstat?"* — Led to designing a `--boot=loci|individuals|both`
   flag rather than a hard switch.
6. (`AskUserQuestion`) User chose: Ar/privAr hard-switched to individuals,
   only Ho/He/Fis get the flag — **later reversed**, see below.
7. *"Why does hierfstat sample by locus? ... Explain it to me like I'm an
   idiot ... Check with the advisor. Then show me the full plan!"*
8. *"Wait, should Fis default be by loci still? Or should it be both? You
   can check the math in the md file."*
9. *"Uh, it sounds like Ho, He, and Fis should all be bootstrapping on
   both... are you sure you're doing this right? I heard you weren't able
   to fully read Owen & Eckles. I could upload a copy."*
10. *"It's in your folder read to read"* — the user uploaded the actual
    Owen & Eckles (2012) PDF, read in full this session (34 pages).
11. *"There's something in your plan saying you're keeping loci as the
    default... But I thought the default was going to be both for He, Ho,
    and Fis. Also, you sure it shouldn't also be default for Ar and
    private alleles? Consult the advisor."* — caught a stale, contradictory
    bullet left over from round 6. Consulted the advisor, who pointed out
    the original Ar/privAr hard-switch reasoning didn't actually hold up:
    resampling loci captures marker-choice uncertainty independent of
    within-tag SNP linkage, so there was no principled reason to exempt
    Ar/privAr from the same flag. **Unified everything under one
    `--boot=loci|individuals|both` flag, default `both`.**
12. *"There will always be more loci than individuals. Does that influence
    how bootstrapping and weighting should be done. Reference Owen &
    Eckles to make decision. Check advisor."* — Re-read the Owen & Eckles
    PDF specifically for its imbalance diagnostic (`ε`), worked out from
    first principles that `ε ≈ 1/n_individuals` in this project's setting
    (loci ≫ individuals is exactly the "modest number of levels" case the
    paper names) — didn't change the design, but added a precise,
    citable note about it.
13. *"Check the logic and statistical foundation, possible with more web
    searching, of the 'Both' bootstrapping. I need this to be right and
    I'm skeptical of trusting an AI model with this."* — Found an
    independent field precedent: Lipson, Loh, Levin, Reich, Patterson &
    Berger (2013), the MixMapper paper, uses the *exact* same combined
    scheme (locus-block resample + within-group individual resample) —
    verified directly from the paper's own Methods text.

**Plan approved** at this point (default `both`, one flag, covering all
nine bootstrapped quantities).

## 4. Implementation

Implemented in `diversity_stats.R`:
- `stat_from_resampled(sel)` — the shared engine for `individuals`/`both`
  modes: draws one fresh per-population individual weight vector per
  replicate, builds weighted allele/typed counts per locus, feeds them into
  the *existing, unmodified* `diversity_core.R` functions
  (`hs_from_counts`, `gene_div_2n_counts`, `rare_richness`,
  `rare_private_all`).
- `--boot=loci|individuals|both` flag, with validation and usage text.
- One bootstrap dispatch (`switch` on mode) covering Ho, He, Fis, He2n,
  Fis2n, poly, Ar, privAr, and a **new** `priv_total` (a genuine per-replicate
  sum, correctly `NA` — not the old silent `NaN` — when a population has zero
  defined loci at all, as opposed to zero private alleles at defined loci).
- Filename-suffix logic so an explicit non-default `--boot` mode can't
  silently overwrite the default output.
- Fixed a `het_between_pops.R` bug (Fix 1, unrelated but flagged during
  stress-testing): `t.test()`/`wilcox.test()` return a "successful" object
  with `NaN` p-values (instead of erroring) when both groups are constant
  at exactly zero — R's zero-stderr guard is a *relative* check that never
  fires there. Fixed with an explicit `is.nan()` check.

Verified: golden-fixture regression (55/0), exact bit-for-bit match against
recorded pre-session numbers under `--boot=loci`, `lo <= point <= hi`
invariant checks, degenerate stress fixtures (n=2 population), filename
logic, `--complete-case` interaction, self-tests.

## 5. The coverage-simulation discovery (the pivotal finding)

While checking the `lo <= point <= hi` invariant, a real anomaly turned
up: under `--boot=individuals`, the point estimate for He/Fis sometimes
fell *outside* its own 95% CI — always on the same side (point above the
upper bound). Ruled out an implementation bug directly: forcing the
resampling weights to identity (no duplication) reproduced the point
estimate to 13 decimal places.

**User:** "explain the individual bootstrap issue for He and Fis like I'm
an idiot reviewer." — answered in plain terms (see the actual explanation
given in-session; summarized: resampling individuals with replacement
means some are counted more than once, but the `n/(n-1)` correction factor
in the He formula doesn't know that, so it comes out biased low).

**User:** "Or maybe disable He and Fis bootstrapping with just
individuals? Are you sure that doing both is statistically valid given the
He and Fis sample correction? Check online and consult advisor if
necessary. Does what we learned about boostrapping mean we should modify
the statistical tests in het_between_pops.R?"

This triggered the decisive step: an actual **coverage simulation** — not
more theory, but synthetic data with a *known true* He/Fis/Ar, bootstrapped
thousands of times, checking what fraction of 95% CIs actually contained
the truth.

**Result — individual-resampling and combined resampling both fail badly:**

```
                loci      individuals   both
He  (n=10)      93.7%       0.0%        0.0%
He  (n=15)      92.7%       0.0%        0.0%
Fis (n=10)      94.0%       4.0%       13.7%
Fis (n=15)      94.0%       4.7%       22.7%
Ar  (g=2n)       ~95%       0.0%        0.0%
Ho  (control)    n/a       88–92%  (point-in-CI 100%)
```

Mechanism: Nei-Chesser's `He` and the rarefaction formula behind `Ar` both
assume every individual contributing to a replicate is *distinct*. A
with-replacement resample is only ~63% distinct individuals; the built-in
finite-sample correction has no way to detect the duplication, so it comes
out systematically biased. **Critically, this gets worse with more loci,
not better** — the bias is fixed, but the bootstrap's spread shrinks as
locus count grows, so the bias-to-noise ratio grows too. `both` mode
doesn't rescue this since it inherits the same individual-axis bias. Ho
(no correction factor) was unaffected — the clean control confirming the
mechanism.

This also meant the earlier claim (that Owen & Eckles' delta-method
argument extends their guarantee to Fis/Ar) was an overstatement: their
proof covers smooth functions of a resampled *mean*; He/Ar have an
explicit, separate sample-size dependence the argument doesn't reach.

**Conclusion: `--boot=loci` is the default for every statistic.**
`individuals`/`both` remain as explicit, clearly-labeled comparison modes,
documented with the measured undercoverage so nobody mistakes them for a
valid CI.

**On `het_between_pops.R`:** checked explicitly and confirmed no change is
needed — its actual test runs on one raw heterozygosity value per real
individual, with no resampling-with-replacement and no `n/(n-1)`-style
correction anywhere in that path. The bias mechanism found above requires
all three ingredients; none are present there.

## 6. Final state

- `diversity_stats.R`: `--boot=loci|individuals|both`, default `loci`.
  Default run is fast (~4s at `nboot=10000` on the test fixture);
  `individuals`/`both` are expensive (5m37s / 3m35s on the same fixture)
  and exist for comparison only.
- `priv_total_lo`/`priv_total_hi` columns added; `priv_total` is now `NA`
  (not a silent `NaN`) when a population has no defined locus at all.
- `het_between_pops.R`: pairwise NaN fix applied and verified.
- `MINIMAL_WORKFLOW.md`: new "Bootstrap mode" section with the full
  coverage-simulation numbers and mechanism; new References (McCullagh
  2000, Owen & Eckles 2012, Hall 1992, Mammen 1992, Keenan et al. 2013).
- Plan file (`~/.claude/plans/replicated-snacking-wombat.md`) rewritten to
  match what was actually measured, including the corrected, narrowed
  delta-method claim.
- Full test suite: **55/0**.

### Follow-up Q&A after the work was reported

**User:** "What was the Ar bootstrapping default before the edits in this
session? If it changes, are we sure? Please explain to an idiot reviewer."

**Answer given:** Before this session there was no `--boot` flag at all —
resampling RAD loci was the *only* method that existed for Ar/privAr's CI,
same as Ho/He/Fis. That is still exactly the default today. Confirmed by
(1) a byte-for-bit regression check against numbers recorded earlier this
session, (2) the fact that `stat_from()`'s existing Ar/privAr lines were
never edited — only a new field was appended alongside them for the new
`priv_total` CI, and (3) 55/0 on the full test suite. The plan to change
the default (to individual- or combined-resampling) was abandoned *because*
the coverage simulation showed it performed catastrophically (0% coverage)
for Ar — so the "change" that was seriously considered never shipped, and
the end state is behaviorally identical to session start for the default
path.

**User:** "Wait, wasn't something bootstrapping by individuals?"

**Answer given:** Two different things had that name, worth separating.
(1) Mid-session, the plan briefly had Ar/privAr *hard-switched* to
individual-resampling with no flag and no loci fallback at all (your own
choice at one `AskUserQuestion` step) — that design was abandoned before it
was ever implemented, superseded first by a unified `both`-default flag,
then by the full revert to `loci`-default after the coverage simulation.
It never shipped. (2) `--boot=individuals` genuinely exists in the shipped
code today, fully working and tested — it's just not the default; you have
to type it explicitly, and the docs tell you not to trust its interval
when you do.

### Did this session touch rarefaction? And is the pre-existing loci-bootstrap for Ar actually *ideal*?

**User:** "Did you touch rarefaction? Ask the advisor about the Ar
confidence interval bootstrap method before the edits you made. Are you
and the advisor SURE that this approach is ideal?"

**Rarefaction itself: untouched, confirmed by file evidence, not just
memory.** `diversity_core.R`'s file modification time is 19:43, hours
before any of today's edits (23:09–23:50 across the other three files).
`p_sampled()`, `rare_richness()`, `rare_private()`, `rare_private_all()`
are byte-identical to session start. `stat_from_resampled()` (the new code
added today) *calls* these functions; it never modifies them.

**Is the pre-existing loci-bootstrap for Ar "ideal"? No — and the honest
answer required correcting something said earlier in this same file.**
Re-examining the coverage simulation's own numbers (52%/70.7% for `loci`
mode with inbreeding in the generative model, vs. ~95% without, at the
identical g=2n case) with the advisor turned up a real distinction that
had been glossed as a "simulation artifact" above — it isn't one:

- Kalinowski's rarefaction formula is a hypergeometric expectation
  *conditional on the observed 2n gene copies actually sampled*: "if I drew
  g of these SAME copies without replacement, how many distinct alleles
  would I expect."
- The coverage check's target (`sum_k[1-(1-p_k)^g]`) is a *different*
  quantity: "g independent draws from the infinite source population."
- These two coincide only when the observed 2n copies are themselves 2n
  independent draws — which inbreeding breaks, since identical-by-descent
  copies are correlated, not independent. Under inbreeding the two
  quantities genuinely differ, and a CI built around the first will not
  cover the second — neither is "wrong," they answer different questions.

So what was actually established is narrower than "the loci bootstrap is
correct": it achieves ~95% coverage **of its own, sample-conditional
estimand**, in a simplified single-population test. It was never evaluated
against — and makes no claim about — an infinite-population richness
target under real inbreeding. And it still only captures marker-choice
(which loci got sequenced) uncertainty, not individual-sampling (which
animals got caught) uncertainty — the entire question this session set out
to investigate, and the one the individual-resampling experiment showed
cannot be patched in by simply swapping the resampling axis, for the same
finite-sample-correction reason found for He/Fis. A method that captured
both correctly would need a variance formula derived specifically for the
rarefaction estimator, not a different resampling axis — that is methods
development beyond this session's scope.

**The accurate claim, going in the plan file and this log**: `--boot=loci`
is the best available option for the estimand it targets, with that
limitation named plainly. Not a claim of statistical ideality.

**User (separately):** "I mean, you were awfully confident about the plan
for both bootstrapping method."

**Acknowledged plainly, not deflected.** The Owen & Eckles literature work
earlier in this session was individually accurate — the paper was read in
full and quoted correctly, the math (the Fis Jacobian, the ε imbalance
calculation) was worked out correctly, and the MixMapper precedent was
verified directly from its own text. None of that was the load-bearing
check, and presenting it with that much confidence before running one was
the actual mistake. The paper's guarantee is proved for a linear mean;
Fis/Ar are not means, and the "delta-method extends this to smooth
functions of means" citation was invoked without verifying that Fis/Ar
actually satisfy that extension's precondition — they don't, because both
have an explicit, separate dependence on sample size that a smoothness
argument about means doesn't reach. Careful literature work answers "is
this mechanism real and did I represent it correctly"; it does not answer
"does it work for my specific estimators," and the two got conflated here
for several rounds. What actually caught it was the user declining to
accept "verified against the primary text" as a stopping point and asking,
repeatedly, for it to be checked — which is what eventually produced the
coverage simulation in section 5 above.

### The definitive per-metric table, with hierfstat and diveRsity as reference points

**User:** "Provide standard error calculations for each metric. I'm
worried bootstrapping theory for this is unclear... I need a real answer!"

Then, after the table below was given: "Add to that table what hierfstat
does and what diveRsity does as a point of reference." All three tools'
behavior was re-verified directly from their own source code this session
(not from memory or the earlier plan file's claims) — `hierfstat`'s
functions were printed directly since the package is installed locally;
`diveRsity`'s were fetched from its GitHub source since it is not
installed here.

`boot.ppfis()` is `hierfstat`'s *only* bootstrap function touching any of
these statistics, and it returns an Fis interval only (its source returns
`1 - colSums(Ho[x,])/colSums(Hs[x,])` and nothing else) — an error caught
by the advisor on review of the first version of this table, which had
wrongly implied it also covers Ho.

| Metric | Correction factor baked in? | This project's default | `hierfstat` | `diveRsity::divBasic` |
|---|---|---|---|---|
| Ho | No | `--boot=loci` | No Ho bootstrap exists (`boot.ppfis()` returns Fis only) | Resamples **individuals** (`sample(dim(pasub)[1], replace=TRUE)`) |
| He | Yes (Nei-Chesser `n/(n-1)`) | `--boot=loci` | No He bootstrap exists, same reason | Resamples individuals, but on the **uncorrected** `1-sum(p^2)` estimator (no `n/(n-1)` term at all, confirmed in both its point-estimate and bootstrap code) — not vulnerable to the SAME compositional bias found here, because it never had a correction factor to interact badly with duplication. Its numbers are a different, already-known-biased-low quantity regardless of resampling axis — this project's own docs already say to compare against Stacks' `Pi`, never the uncorrected `Exp_Het`-style estimator diveRsity uses here. |
| Fis | Yes (built from He) | `--boot=loci` | `boot.ppfis()` resamples loci, but **per-SNP row**, not per-RAD-block (source: `x <- sample(nloc, replace=TRUE)`) — this is the mechanism behind this project's own "~1.8x too narrow" finding | Resamples individuals; same uncorrected-He caveat as above |
| Ar / privAr / priv_total | Yes (hypergeometric, keyed on total gene copies) | `--boot=loci` | `allelic.richness()` has **no bootstrap at all** — point estimate only | `rarefactor()` computes rarefied richness but has **no bootstrap and no private-allele richness at all** — point estimate only |

Two things worth being explicit about, since they change how to read the
"comparison mode" framing:
- `hierfstat`'s own loci bootstrap (the field's usual comparison point) is
  *itself* not what this project's `--boot=loci` does — hierfstat resamples
  individual SNP rows, ignoring within-RAD-tag linkage, which is exactly
  the ~1.8x-too-narrow problem this project's own block-bootstrap (over
  whole RAD tags, not rows) was built to fix. So "matches hierfstat" means
  matching hierfstat's resampling *axis* (loci), not its resampling *unit*.
- Neither `hierfstat` nor `diveRsity` offers any bootstrap for Ar/privAr at
  all — this project's Ar/privAr confidence interval has no equivalent to
  benchmark against in either tool.

### Standard errors: the jackknife-over-loci addition

Rather than a closed-form SE formula (none exists that correctly handles
within-RAD-tag SNP linkage without effectively becoming a block-resampling
method anyway), added the field-standard alternative to bootstrapping:
**delete-one-block jackknife over RAD loci** (Weir 1996, *Genetic Data
Analysis II* — the same method behind GENEPOP/FSTAT's own Fst/Fis standard
errors). Implemented as a new `jack_se()` function in `diversity_stats.R`,
using the existing, unmodified `stat_from()` (the `--boot=loci` mechanism)
with one RAD-locus block dropped per replicate — never resamples with
replacement, so it cannot have the compositional-duplication bias found
for individual-resampling. Runs regardless of `--boot` mode or `--nboot`
(it needs no random numbers at all — with `nL` RAD loci there are exactly
`nL` delete-one replicates, deterministic run to run).

```
SE_jackknife = sqrt( (nL-1)/nL × Σ_b (θ̂₍₋b₎ − θ̄)² )
```

New `_se` columns added alongside every `_lo`/`_hi` pair in both output
tables (`Ho_se`, `He_se`, `Fis_se`, `Ar_se`, `privAr_se`, `priv_total_se`).

**Cross-checked against the bootstrap's own spread**, at two very
different locus counts (a jackknife's reliability is known to depend on
block count, so checking only one size wouldn't have been enough): on
`sim.haps.vcf` (304 RAD loci, `--nboot=10000`, `--boot=loci`, the default),
`1.96 × SE_jackknife` matched the bootstrap CI's half-width to within about
1% for every metric checked (He, Fis, Ar, privAr, priv_total); re-checked
on a much smaller 60-RAD-locus, n=2-population fixture, the same
comparison held within about 2–6% — still close agreement at the small end
of what this pipeline supports, not a divergence.

**A real, previously-unverified bug turned up while checking this**, and
was fixed rather than left for later: `Ar`/privAr's point estimate printed
the literal text `NaN` (not `NA`) when a population had it undefined at
every retained locus — `colMeans(..., na.rm=TRUE)` on an all-`NA` column is
`0/0`, which R evaluates to `NaN`. `priv_total` already had this guarded
from earlier in the session; `Ar`/privAr did not, and the new jackknife SE
inherited the same `NaN` from them. Fixed in `stat_from()` with the same
explicit guard already used for `priv_total`, verified via a direct
synthetic-override test (not reachable through any existing test fixture),
confirmed the normal-case numbers are bit-for-bit unchanged, and re-ran the
full suite (55/0). This was pre-existing behavior the session did not
introduce; it surfaced only because validating the new SE's NA-handling
required constructing exactly the input that triggers it.

All 55 tests still pass; the three richness/per-population golden files
were regenerated (new `_se` columns only — no value changes to any
existing column, confirmed by `check_golden.R`'s diff before updating).

### On the earlier overconfidence

**User:** "This needs to be right. Don't mess up like with the both
bootstrapping nonsense."

Addressed by re-verifying hierfstat/diveRsity from their actual source
code this session (see table above) rather than trusting the earlier
session's memory of having done so, and by cross-checking the new
jackknife SE numerically against the bootstrap before reporting it as
correct, rather than asserting it from the formula alone.

### How does Stacks itself compute SE, and a bug found along the way

**User:** "How does STACKS calculate SE. Is the method you implemented
similar?"

Verified directly from Stacks' own source (downloaded the 2.68 source
tarball from catchenlab.life.illinois.edu — 2.68 confirmed to be the
current release as of 2024-08-23 via the official changelog page, not a
stale version), specifically `src/populations.cc`,
`SumStatsSummary::final_calculation()`:

- Stacks computes `StdErr = sqrt(sample variance) / sqrt(n)` (classical
  SD/√n), variance via Welford's online algorithm, accumulated over
  **every variant SNP site** — not RAD-locus blocks.
- This is the identical per-row-not-per-block mechanism already documented
  in this project as the reason hierfstat's `boot.ppfis()` gives intervals
  "~1.8x too narrow" — so Stacks' own `StdErr` should be expected to be
  optimistically narrow for the same reason.
- Stacks' `Fis` mean is a running average of the per-site Fis ratio (the
  mean-of-ratios aggregation this project's own docs argue against).
- **Answer to "is my method similar": no**, on both dimensions that
  matter — per-block vs per-site replication unit, and a delete-one
  recompute of the ratio-of-sums estimator vs a classical SD/√n on raw
  per-site ratios.

**A bug found while reading the source, verified rather than assumed**:
the `Fis` `StdErr` output line in that same file reads
`sqrt(_num_indv_var[j]) / _sq_n[j]` — every other statistic's `StdErr` in
the same output block correctly uses its own matching `_..._var[j]`
variable (checked each one, header-to-data-column, 1:1); only `Fis`'s
uses `_num_indv_var[j]` instead of `_fis_var[j]`, a copy-paste error. The
`Var` column for Fis is computed and printed correctly right next to it —
only the adjacent `StdErr` computation uses the wrong variable. Confirmed
this is distinct from an already-fixed 2015 bug (Stacks 1.30: "Fis values
... were incorrect", about the Fis *mean*, not the StdErr) via the
official changelog, so it isn't a rediscovery of old, resolved history.

**User:** "Double check that the STACKS error is real. Write something
about the STACKS 2.68 error that I can send to the developers."

Re-verified the column alignment (header list vs. data-row list) line by
line to rule out a misread before writing anything up. Wrote
`stacks_bug_report_fis_stderr.md` in the project directory: version
affected, exact file/function, the code excerpt with every StdErr line
shown for comparison (not just the wrong one, so the pattern is visible),
expected vs. actual, a one-line suggested fix, how a user could spot it
independently in their own output (Fis StdErr won't track Fis Var the way
every other stat's does), and the note distinguishing it from the 2015
issue. Also found and reported: no GitHub issue tracker exists for Stacks
itself — the `stacks-users` Google Group (where the actual developers,
Rochette and Catchen, are active) is the real channel, so the report was
written in a form suitable for posting there or emailing directly, not
filed anywhere by this session.

---

*This file was written at the user's request to preserve a record of the
session. It is a reconstruction from the conversation, not a raw log
export — quotes from the user are verbatim where shown; narration of
intermediate reasoning is summarized.*
