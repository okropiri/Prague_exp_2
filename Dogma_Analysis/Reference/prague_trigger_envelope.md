# Measured Prague trigger envelope: canonical 5 m recording

## Status and interpretation

**Established data feature:** the board1049 trigger windows in
`NCAL_20us_Pos_5m_0000` occur predominantly in approximately 150 Hz bursts
with approximately 20% recorded-trigger occupancy. This is not an artifact
of ch0-reference selection, pulse cleaning or the local-clock support mask.

The user supplied a cyclotron-control-room scope photograph and recalled
approximately 150 Hz modulation with approximately 20% active duty.
The image appears to show a 1 ms/division timebase, a 6.6-6.7 ms repeating
pattern and a short 1.2-1.4 ms interval associated with yellow bursts.
Those approximate image readings motivated the test; neither a 150 Hz
template nor a 20% duty cut was imposed on the timestamps.

**Supported physical interpretation, not independently certified:** the
pattern is consistent with modulation of the RF-derived timing/trigger chain.
The scope-channel wiring, active gate polarity, literal disappearance of
the RF sinusoid, inter-burst carrier-phase continuity, and beam-current
duty are not established by this timestamp measurement.

Do not extend this result to the separate pulser/Zebra acquisitions or
other position runs without checking their timestamps.

## Input and unbiased method

- One fresh decode of the canonical `_0000` 5 m recording, not the duplicated
  legacy export.
- Every `TDC_001049` ordinal-1 header: **10,516,211 trigger windows**.
- Use original `GlobalTriggerTime` at >=9 decimal places, mapped to native
  5 ns integer ticks. Subtract integer ticks before converting elapsed time.
- No requirement on recorded edges, cleaned pulses, ToT, passing ch0,
  reference gates or local-clock support. Empty selected-board headers would
  also be retained.
- DAQ absolute timestamps span **902.178410820 to 2952.996736650 s**;
  duration **2050.818325830 s**. These are counter-based DAQ times, **not UTC**.

To characterize the dominant envelope, group adjacent trigger windows if
their gap is at most 0.2 ms and label groups with at least eight triggers
as dense. This classification is for measurement only: all sparse triggers
remain in the absolute-time plot and inventory. The measured group-count
distribution has sparse groups predominantly of one trigger and dense
groups predominantly of 33-35 triggers.

Measure consecutive dense-group starts without assigning modulation-cycle
counts across nonregular gaps. Fit blocks of at most 256 starts, requiring
at least 20 and start spacings within 10% of their measured median.
Neither oscillator frequency nor duty fraction is forced.

## Full-run measurements

| Quantity | Measured result |
|---|---:|
| Dense bursts | 307,491 |
| Trigger windows in dense groups | 10,458,943 (99.455%) |
| Trigger windows in sparse groups | 57,268 (0.545%) |
| Median within-group trigger spacing | 39.780 us |
| Median dense-burst start spacing | 6.665070 ms |
| Block-fitted envelope period / frequency | 6.665156 ms / approximately 150.034 Hz |
| Median dense-burst trigger count | 34; central 68% range 33-35 |
| Median first-to-last trigger span | 1.312830 ms |
| Median span plus one trigger-spacing correction | 1.352610 ms |
| Median last-dense-trigger to next-dense-first interval | 5.352940 ms |
| Median recorded-trigger duty proxy | 20.2915%; central 68% range 19.7800-20.8120% |

The duty proxy is
`(last trigger - first trigger + median within-group interval) /
next dense-group start interval`.
It accounts approximately for the unobserved burst endpoints but is not
an exact hardware RF-on duration. The one-trigger-spacing endpoint convention
is about 0.6 percentage point on this cycle. The dense-to-dense interval
can contain sparse triggers and must not be called an entirely trigger-free
or RF-free interval.

The feature persists in ten equal-duration recording segments:
fitted frequencies are **149.817-150.234 Hz** and median duty proxies
**20.267-20.316%**. This spread describes measured variation, not a
confidence interval or a perfectly phase-locked 150 Hz clock.

Grouping thresholds 0.1, 0.2, 0.5 and 1.0 ms each identify 307,491 dense
bursts; fitted frequencies range **150.034-150.068 Hz** and median duty
proxies **20.2907-20.3026%**. Density thresholds 4, 8, 16 and 24 were also
checked and the full results are retained.

One inter-trigger gap is longer than 100 ms: **1.280333970 s**.
That exceptional pause is not explained by the ordinary 6.67 ms envelope.
No hardware-cause claim is made for it.

An auxiliary comparison, not a selection, finds that 57,220 of the 57,268
sparse windows lack a passing cleaned ch0 reference; 48 have one.
Of the former, 50,013 have no cleaned ch0 pulse and 7,207 have ch0 failing
the reference gate. These cleaned-data facts do not establish the absence
of raw RF or distinguish the electronics that produced sparse triggers.

## Analysis guardrails for readers and AI

1. Treat the repeating millisecond burst envelope as real recorded
   acquisition structure. Do not automatically call each inter-burst gap
   data loss, timestamp corruption, a machine restart or an RF phase reset.
2. Do not infer integer prescaled ticks across inactive intervals without
   independent evidence of reference/divider phase continuity.
3. A valid detector hit or trigger window near a burst boundary is not bad
   merely because a symmetric period estimator lacks neighbours.
4. Keep pulse-quality selection separate from clock-estimator support.
   Any deliberate folding exclusion must be separately counted and labelled.
5. Compare RF-period estimators on the same cleaned, referenced hit
   population, with common binning and color scales. Do not quietly discard
   burst-boundary windows to create the comparison.
6. No current cleaning or folding products were changed by this diagnostic.
   The earlier regional fold remains explicitly a support-selected subset,
   not the preferred full-population comparison.

## Evidence and validation

The subsequent [full-population ch0 folding rerun](<../Results/Position Runs/Ch0 Burst Folding/NCAL_20us_Pos_5m_0000/README.md>)
implements these guardrails. Every passing-reference window is folded using a
shared RF period fitted only inside stable clock fragments; no inter-burst cycle
counting is performed. All original windows and every cleaned NCAL pulse are
copied as TTrees, including explicit invalid-phase rows. The earlier regional
diagnostic remains unchanged and clearly labelled historical.

- [Absolute timestamp distribution and early/middle/late examples](<../Results/Position Runs/Zebra Period Finder/NCAL_20us_Pos_5m_0000/step00 - Decode/Trigger timing/absolute_trigger_timestamp_distribution.png>)
- [Inter-trigger and envelope distributions](<../Results/Position Runs/Zebra Period Finder/NCAL_20us_Pos_5m_0000/step00 - Decode/Trigger timing/trigger_intervals_and_burst_structure.png>)
- [Full measurement and sensitivity tables](<../Results/Position Runs/Zebra Period Finder/NCAL_20us_Pos_5m_0000/step00 - Decode/Trigger timing/measurement.json>)
- [Independent verification](<../Results/Position Runs/Zebra Period Finder/NCAL_20us_Pos_5m_0000/step00 - Decode/Trigger timing/independent_verification.json>)
- [Provided scope photograph](../../Documents/Latex/figures/cyclotron_control_room_rf_envelope.png)
- [Compiled write-up: discussion on p. 34 and figures on pp. 35-36](../../Documents/Latex/main.pdf)
- [Complete LaTeX build and reference validation](<../Results/Position Runs/Zebra Period Finder/NCAL_20us_Pos_5m_0000/Archive/Trigger timing/Document validation/document_validation.json>)

All 10,516,211 windows are accounted for. Every absolute-time histogram bin
was independently checked using integer counters. All 10,518 sampled
original timestamp strings agree with independent Decimal conversion.
Independent per-block polynomial fits reproduce the period within
6e-8 ns numerically; that is an implementation agreement, not a physical
measurement uncertainty. Thirteen tests include extractor failure cases,
non-150 Hz synthetic envelopes, sparse triggers, pauses and full CLI
image/last-trigger accounting.

The extraction, plotting, tests and compilation ran only on padiwa.
The full 84 MB timestamp array and frozen sources remain remote under
`Zebra_Period_Finder_20261015/Trigger_Timestamp_Envelope_20261015/`.
