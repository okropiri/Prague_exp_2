# 5 m analysis findings and limitations

**Historical support-selected diagnostic.** The
[corrected full-population run](<../../../../Ch0 Burst Folding/NCAL_20us_Pos_5m_0000/README.md>)
retains every passing-reference window and preserves all original windows
and cleaned NCAL pulses in ROOT TTrees. The findings below describe the
earlier diagnostic and are retained for provenance.

These are final-run findings. Tests and rejected attempts are kept in the
[technical archive](<../../Archive/>) and are not alternate final datasets.

## Subsequent timing-context correction

The [unbiased all-trigger diagnostic](<../../step00 - Decode/Trigger timing/README.md>)
establishes a dominant approximately 150 Hz, approximately 20% recorded-trigger
envelope. Ordinary inter-burst millisecond gaps are therefore a recurring
data feature, not automatically data loss, corruption or phase resets.
The precise RF/discriminator/trigger-gating cause remains unconfirmed.
Sparse inter-burst triggers and one 1.2803 s exceptional gap are also present.

The symmetric clock-support exclusions in the existing fold remain
implementation choices, not detector-quality failures. They must not be used
to silently disregard valid burst-boundary windows in a period-estimator
comparison. The timing diagnostic retains every trigger window and changes
no existing cleaned or folded product.

## What the controls show

The [ch2 four-panel control](ch02_period_and_tot_controls.png) compares:

1. Local ch0 periods: 16,203,826 supported ch2 pulses.
2. One fixed ch0-derived period, using exactly the same pulses.
3. The same phases with ToT permuted within channel/ROOT chunk, preserving
   both displayed one-dimensional marginals exactly.
4. The fixed period with **all 31,159,659 referenced ch2 pulses**, including
   the windows excluded by the regional clock mask.

The diagonal low-ToT structures remain in panels 1, 2 and 4, but disappear
after permutation. They also remain in the
[2x2 count-conserving rebin](ch02_rebinned_control.png).
Ch2 period optimization, adaptive-period fitting, regional selection and
the original display bin size are therefore not individually necessary
for the visible diagonal structure.

This is evidence of a phase-ToT joint association in the decoded/cleaned
data, not a determination of its physical or instrumental cause.

## Issues tracked

| Issue | Evidence and treatment |
|---|---|
| Duplicated original export | Full line-by-line comparison confirms two consecutive copies, each 610,789,983 lines and 568,714,285 edges; only process pointer addresses are ignored. Identical duplication doubles statistics but does not by itself invent a new joint distribution. |
| Timestamp restart stitching | Original cleaning stitched one backward jump forward by about 2050.818326 s. Fresh cleaning has no restart. Stitching is not deduplication. |
| Mixed boards in the original analysis | Old cleaning used all boards; the new chain explicitly selects board 1049. Board 1048 contributes 2,557,691 raw edges. Old/new total-count differences are not purely estimator differences. |
| Integer-tick inference across gaps | The original estimator rejected the pilot, with local estimates about 38.72-39.51 ns and about 199 ns maximum half-interval disagreement at 6000 ns. The subsequent unbiased diagnostic confirms a regular burst envelope. Do not infer integer prescaled ticks through inter-burst intervals without evidence of carrier/divider continuity; the regional fold is a selected diagnostic, not a detector-quality rule. |
| Regional support selection | 5,233,538 of 10,450,559 reference windows lack a supported period. Most are boundary/support exclusions; 4,286 fail phase consistency. All excluded pulses remain counted upstream. Panel 4 checks the all-reference population. |
| Detector-based period/phase selection | Removed in the new clock. Old period 38.8495 ns and phase origin 37.612242 ns came from ch2 scoring. New zero is selected ch0, so old/new phase origins are not aligned. |
| First-valid ch0 reference | Selection is first passing hit in file order, not minimum time. 65,652 windows lack a passing reference. The selected ch0 self-hit always folds at zero and cannot independently validate the period. |
| Pairing and ToT acceptance | Fresh cleaning records 3,352,795 orphan rises, 176,947,155 orphan falls and 6,033,344 below-minimum pairs across all channels. Their contribution to ch2 stripes is not isolated; these are not ch2-only failure counts. |
| Shared rise-time variable | ToT is fall minus rise, while phase uses rise minus ch0. Rise-time error, time walk or calibration can create correlations. Shuffling removes both instrumental and physical correlations and cannot distinguish them. |
| Phase-versus-time plots | Folded phase is a function of unfolded time. Wrapping diagonal stripes in those maps are expected by construction and are not independent physics evidence. |
| Binning and color scales | Phase bins are about 0.249 ns and ToT bins 0.1 ns. The rebinned control preserves all counts; log colors can still exaggerate visually weak structures. |
| True 3P folding | Verified from eligible pulses using modulo 3 local cycles. The output is not three pasted copies of a one-period histogram. |
| Native scale and assumed divider | Native coarse quantum 5 ns and prescaler 1024 remain declared assumptions, not independent hardware certification. No Zebra 150 MHz rescaling or pulser divider 4096 is applied. |
| Truncated raw tail | Final record has 292 bytes missing. All 10,521,678 complete events were decoded; the incomplete record was not recovered. |
| Dead time and acceptance | Rates span 2050.818335 s and are accepted counts per elapsed second, not dead-time/efficiency-corrected flux. |
| Display overflow | Fourteen pulses above the 128 ns Step03 display remain in the full ToT distribution, up to 132.464 ns. None belongs to the supported fold. The display limit is not a cleaning cut. |

## Evidence

- [Independent final Step04 verification](<../verification.json>)
- [Clock assessment](<../assessment.json>)
- [Clock diagnostics](<../Clock diagnostics/>)
- [Exact original-export comparison](<../../Archive/Reports/legacy_export_comparison.txt>)
- [Original clock rejection investigation](<../../Archive/Reports/pilot_clock_guard_investigation.json>)
- [Run overview and current layout](../../README.md)

Original position-run products and the two specially copied original ROOT
files were not replaced. A numerical PASS or surviving zebra pattern does
not certify instantaneous RF phase or establish a beam-physics explanation.
