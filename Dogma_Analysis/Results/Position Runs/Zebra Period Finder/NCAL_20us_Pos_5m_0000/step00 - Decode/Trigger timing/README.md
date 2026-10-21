# Absolute trigger timestamps - final measurements

**Confirmed recorded-data feature:** approximately 150 Hz trigger bursts
with approximately 20% trigger occupancy in the canonical 5 m `_0000` run.
This is consistent with the supplied control-room scope photograph.
The exact RF/discriminator/gating mechanism and beam-current duty are not
certified by timestamp data alone.

## Open these first

- [All-window absolute distribution, with early/middle/late zooms](absolute_trigger_timestamp_distribution.png)
- [Inter-trigger intervals and burst structure](trigger_intervals_and_burst_structure.png)
- [Stability across the recording](envelope_stability_across_recording.png)
- [Sparse versus dense group inventory](group_trigger_count_inventory.png)
- [Measurements and threshold/density sensitivity](measurement.json)
- [Independent QC and reference association](independent_verification.json)
- [Full context and reader/AI guardrails](../../../../../../Reference/prague_trigger_envelope.md)

## Result

Every one of **10,516,211** native board1049 trigger-window timestamps
is retained; no edge, cleaned-pulse, ToT, ch0-reference or local-clock cut
selects these distributions. Absolute times are DAQ counter times, not UTC.

| Quantity | Measured |
|---|---:|
| DAQ timestamp range | 902.178410820-2952.996736650 s |
| Duration | 2050.818325830 s |
| Dense bursts | 307,491 |
| Triggers in dense / sparse groups | 10,458,943 / 57,268 |
| Typical triggers per dense burst | 33-35; median 34 |
| Median within-group spacing | 39.780 us |
| Median burst-start spacing | 6.665070 ms |
| Aggregate block-fitted repetition | Approximately 150.034 Hz |
| Median first-to-last trigger span | 1.312830 ms |
| Span plus one trigger-spacing correction | 1.352610 ms |
| Median trigger-duty proxy | 20.29% |
| Last dense trigger to next dense first | 5.352940 ms; can contain sparse triggers |
| Exceptional gap >100 ms | One, lasting 1.280333970 s |

The duty proxy estimates trigger occupancy, not an exact hardware RF-on
gate. Repetition varies slightly during the recording: ten segment fits
range from 149.817 to 150.234 Hz. Threshold and density checks are saved,
not hidden; no 150 Hz template or 20% acceptance cut was imposed.

Dense classification uses at least eight triggers per group, with adjacent
gaps <=0.2 ms. This is a **description of the envelope**, not a detector-hit
selection. All sparse triggers remain in the absolute plot and inventory.

Thirteen targeted tests passed, including complete CLI output/accounting.
Independent QC verifies every 0.1 s histogram bin with integer counters,
all 10,518 Decimal timestamp samples, sparse/dense accounting, and the
burst period using separate per-block polynomial fits.

**Ordinary inter-burst gaps are not automatically data loss or RF phase
resets. Valid detector hits near burst boundaries must not be excluded
merely because a symmetric clock estimator lacks neighbours.**
No existing cleaned or folded product was modified by this diagnostic.

The full timestamp binary remains on padiwa. Local source-copy provenance
is recorded in [the supplementary timing manifest](<../../Archive/Provenance/TRIGGER_TIMING_COPY_MANIFEST.json>).
