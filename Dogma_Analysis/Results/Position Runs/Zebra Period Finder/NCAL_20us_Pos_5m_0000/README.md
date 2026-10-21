# 5 m position run - Zebra Period Finder

**For the corrected full-population ch0 fold, open the
[Ch0 Burst Folding results](<../../Ch0 Burst Folding/NCAL_20us_Pos_5m_0000/README.md>).**
The Step04 fold below is preserved as the earlier support-selected diagnostic,
not the current full-population alternative. Steps00-03 remain the verified
fresh upstream inputs used by the corrected run.

## Start here

**The five step folders contain final results only.** There are no pilot,
full or attempt choices to make inside them.

| Folder | Contents |
|---|---|
| [step00 - Decode](<step00 - Decode/>) | Final decode validation and raw framing checks |
| [step01 - Cleaning](<step01 - Cleaning/>) | Final cleaning summary, reference quality and QC |
| [step02 - Rates](<step02 - Rates/>) | Final all-channel rate plots and report |
| [step03 - Time and ToT](<step03 - Time and ToT/>) | Final unfolded trigger/ch0 timing and ToT plots |
| [step04 - Ch0 period and folding](<step04 - Ch0 period and folding/>) | Final clock diagnostics, folded plots and analysis controls |
| [Archive](Archive/) | Technical history only: pilots, rejected fit, logs, scripts and provenance |

For the comparison you requested:

- [Confirmed millisecond trigger envelope: all-window absolute timestamps](<step00 - Decode/Trigger timing/README.md>)
- [NCAL ch2 - final 3P versus ToT](<step04 - Ch0 period and folding/Folded_RF_3x/ch02_rf_phase3_vs_tot.png>)
- [All final 3P plots](<step04 - Ch0 period and folding/Folded_RF_3x/>)
- [Ch2 - fixed-period and shuffled-ToT controls](<step04 - Ch0 period and folding/Analysis checks/ch02_period_and_tot_controls.png>)
- [Ch2 - rebinned control](<step04 - Ch0 period and folding/Analysis checks/ch02_rebinned_control.png>)
- [Analysis findings and limitations](<step04 - Ch0 period and folding/Analysis checks/ANALYSIS_ISSUES.md>)

## What was run

Steps00-04 completed on padiwa for
`NCAL_Test_Prescaler1024_TrigWindow_0x83e803e8_20us_Pos_5m_20260425_0000`.
The original position-run outputs and previously copied original ROOT files
are unchanged.

The unmodified single-clock estimator rejected irregular acquisition gaps.
That rejected attempt is archived, not presented as a passing result.
The final alternative explicitly uses gap-aware ch0-only periods, with no
detector-feature scoring or detector-derived phase alignment.

**Subsequent context correction:** an unbiased timestamp investigation
confirms approximately 150 Hz bursts with approximately 20% recorded-trigger
occupancy. Ordinary millisecond gaps are a recurring acquisition feature,
not automatically bad/missing windows. The existing regional fold is still
a support-selected diagnostic, not a matched-population comparison. No old
results were changed. See the
[trigger timing report](<step00 - Decode/Trigger timing/README.md>).

| Setting | Value |
|---|---|
| Hardware / TDC ordinal | TDC_001049 / 1 |
| Native coarse quantum / assumed prescaler | 5 ns / 1024 |
| First passing ch0 rise / ToT | -410 to -395 ns / 16.5 to 19.5 ns |
| Cleaning ToT | 1 to 5000 ns |
| Detailed trigger/ch0 axes | -6000 to 6000 ns, 5 ns bins |
| Displayed ToT | 0 to 128 ns, 0.1 ns bins |
| Clock support | At least 32 references per contiguous region; eight-reference boundary margin |
| Clock averaging / consistency spans | Half-span 4 / half-spans 2, 4 and 8 |
| Conditional consistency budget at 6000 ns | 0.1 ns |

Zebra defaults and pulser settings were not replaced.

## Measured result and accounting

The supported local-period mean used for display is **38.850056901754 ns**.
The supported-region endpoint estimate is **38.850132620741 ns**,
corresponding to **25.739937872596 MHz**. These are different averages,
not an uncertainty interval. The original ch2-derived period was
**38.8495 ns**, with phase origin **37.612242 ns**.

| Full-run category | Pulses |
|---|---:|
| Cleaned input | 186,894,978 |
| Supported clock; folded | 95,843,347 |
| Unsupported clock; retained upstream, not folded | 89,779,409 |
| Missing passing reference | 1,272,207 |
| Nonreference ch0 in supported windows | 15 |

The four output categories sum exactly to the cleaned input. Unsupported
windows were not deleted from cleaning, rates or full ToT distributions.
The all-reference constant-period ch2 control contains 31,159,659 pulses;
the support-selected main ch2 fold contains 16,203,826.

The diagonal ch2 structure remains with both local and fixed ch0-only
periods, including the all-reference control. Shuffling ToT removes the
diagonal structure while preserving both one-dimensional marginals.
This constrains an estimator-only explanation; it does not prove beam physics.

## Validation and storage

- 34 analysis regression cases and seven transfer-safety cases passed.
- Independent QC reconstructed support, periods and every saved histogram
  bin. Maximum independent period difference: 4.98e-13 ns.
- 1,573 source artifacts were copied and SHA-256 verified, including
  682 PNG/PDF pairs. Reorganization preserved every copied source byte.
- Subsequent all-trigger timing diagnostics have their own
  [supplementary copy manifest](Archive/Provenance/TRIGGER_TIMING_COPY_MANIFEST.json);
  they are additional products, not a replacement for that original inventory.
- The [current local manifest](Archive/Provenance/COPY_MANIFEST.json) maps
  each remote source path to its simplified local destination. The
  [current hash list](Archive/Provenance/SOURCE_SHA256SUMS.txt) uses paths
  relative to this run folder.
- The [original detailed report](<Archive/Provenance/Original run report.txt>)
  and [original issue ledger](<Archive/Provenance/Original analysis issues.txt>)
  are preserved byte-for-byte as historical snapshots. Paths inside those
  snapshots describe the original layout, not this local browsing layout.

Only the local presentation changed. The padiwa tree, processing scripts,
bulk data and ROOT files remain in their original remote locations:

```text
/data6/Dogma_analysis_by_Dachi/Results/NCAL_20us_Pos_5m_0000/Zebra_Period_Finder_20261015/
```

Archived scripts are frozen records of that remote workflow; do not execute
them from their relocated local archive. DOGMA processing remains padiwa-only.
