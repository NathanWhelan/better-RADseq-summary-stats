# Bug report: `Fis` `StdErr` in `populations.sumstats_summary.tsv` reports `Num_Indv`'s standard error instead of `Fis`'s

**Version affected:** Stacks 2.68 (current release as of 2024-08-23, source
tarball from catchenlab.life.illinois.edu)

**File / function:** `src/populations.cc`, `SumStatsSummary::final_calculation()`,
in the block that writes the "# Variant positions" table of
`populations.sumstats_summary.tsv`.

## Summary

The `StdErr` column reported for `Fis` in `populations.sumstats_summary.tsv`
is not `Fis`'s standard error — it's `Num_Indv`'s standard error, printed a
second time by mistake.

## The code

Every other `StdErr` value in this output block is computed from its own
matching variance variable, printed immediately after that variable:

```cpp
<< _num_indv_mean[j]       << "\t"
<< _num_indv_var[j]        << "\t"
<< sqrt(_num_indv_var[j]) / _sq_n[j] << "\t"     // Num_Indv StdErr -- correct
<< _p_mean[j]              << "\t"
<< _p_var[j]               << "\t"
<< sqrt(_p_var[j])         / _sq_n[j] << "\t"     // P StdErr -- correct
...
<< _pi_mean[j]             << "\t"
<< _pi_var[j]              << "\t"
<< sqrt(_pi_var[j])        / _sq_n[j] << "\t"     // Pi StdErr -- correct
<< _fis_mean[j]            << "\t"
<< _fis_var[j]             << "\t"
<< sqrt(_num_indv_var[j])  / _sq_n[j] << "\n";    // Fis StdErr -- WRONG: should be _fis_var[j]
```

`_fis_var[j]` is computed correctly (via the same Welford's-algorithm
accumulation as every other statistic, over `s->nucs[pos].fis()`) and is
printed correctly in the preceding `Var` column — it's simply not the
value used in the adjacent `StdErr` computation. This looks like a
copy-paste of the first `StdErr` line in the block whose variable name was
never updated to match.

## Expected vs. actual

- **Expected:** `Fis` `StdErr` = `sqrt(_fis_var[j]) / sqrt(_n[j])`
- **Actual:** `Fis` `StdErr` = `sqrt(_num_indv_var[j]) / sqrt(_n[j])`

## Suggested fix

```diff
-           << sqrt(_num_indv_var[j])  / _sq_n[j] << "\n";
+           << sqrt(_fis_var[j])       / _sq_n[j] << "\n";
```

## How to confirm this independently

`Num_Indv` variance (variance in the count of individuals genotyped per
site — an integer-scale quantity, typically small) and `Fis` variance (a
ratio/proportion, typically in the 0–1 range or slightly negative) are on
completely different scales, so this is usually visible by inspection: the
`Fis` `StdErr` column will not track the adjacent `Fis` `Var` column the
way every other stat's `StdErr`/`Var` pair does in the same row, and will
instead track the `Num_Indv` `Var` column near the start of the row.

## Note on a related, but different, historical issue

The Stacks changelog records a 1.30 (2015-05-07) fix: *"Fis values in
batch_X.sumstats_summary.tsv were incorrect (although raw values in
batch_X.sumstats.tsv were correct)."* That was about the `Fis` **mean**
being wrong and was fixed at the time. This report is a distinct,
apparently still-present issue in the `Fis` **StdErr** column specifically,
not a regression of the 2015 bug.

## How this was found

Found while comparing Stacks' `populations.sumstats_summary.tsv` `StdErr`
methodology against a different pipeline's confidence-interval approach,
by reading Stacks' own source rather than inferring behavior from output
files. Not found via a crash or a runtime error — this is a silent,
plausible-looking wrong number, which is why it's worth flagging even
without a specific dataset that "looks wrong" on its face.
