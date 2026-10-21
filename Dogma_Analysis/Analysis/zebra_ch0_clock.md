# Parallel Zebra ch0 clock calibration

## Full-population burst-aware mode

For the canonical native 5 m run, `--burst-clock --allow-reference-gaps`
replaces the earlier support-selected regional fold as the full-population
alternative. The default Zebra and multi-cadence pulser paths are unchanged.
The [completed canonical 5 m results](<../Results/Position Runs/Ch0 Burst Folding/NCAL_20us_Pos_5m_0000/README.md>)
include the full local NCAL pulse TTrees and independent bin-by-bin QC.

The [measured millisecond envelope](../Reference/prague_trigger_envelope.md)
requires separating **period-calibration evidence** from **folding eligibility**:

- Fit one shared RF period using ch0 only. Each contiguous unit-cadence
  reference region has its own intercept, so inactive gaps never contribute
  elapsed RF cycles or imply carrier-phase continuity.
- Calibration regions require 32 references. Stable-cadence membership is
  consecutive original IDs and an interval within one declared coarse quantum
  of the dominant interval. This selects clock evidence, **not hits/windows**.
- Validate against every region of at least 16 references, including those
  too short to calibrate. The maximum-delay phase displacement from regional
  slope disagreement plus three conditional slope errors must meet 0.1 ns.
  Errors use the greater of observed regression scatter and the declared
  edge precision. Also vary calibration lengths 16/24/32.
- Isolated non-core intervals are not counted as nominal divider intervals.
  A sustained non-core cadence of 15 intervals rejects the shared-period
  model rather than silently hiding another RF frequency.
- Assign the resulting period to **every passing-reference window**, including
  first/last windows and short or isolated regions. Such regions use the
  explicitly stated shared-period assumption, not a fictitious local fit.
- Anchor each window at its guarded first passing ch0 in file order, never
  at the minimum-time ch0 or a detector-derived peak. All ch0 pulses are
  retained in folding, not just the chosen reference pulse.

The stability budget is conditional evidence for a shared-period fold,
not a certification of instantaneous RF phase during startup or between
unobserved cycles. The nominal DAQ time scale and physical divider remain
hardware assumptions. Failure of the shared model stops this mode; it
does not produce a passing result by excluding inconvenient windows.

Windows without a passing ch0 cannot acquire a justified ch0 phase. They
remain in the all-window TTree, and their NCAL pulses remain in pulse TTrees
with `phase_valid=false` and NaN phase. ToT display limits affect plot bins,
not the retained pulse trees. No cleaning or original ch2 products are changed.

[verify_ch0_burst_fold.py](verify_ch0_burst_fold.py) independently reconstructs
the period with interval-weighted regression, compares every 1P/3P profile
and detailed-map bin against the complete cleaned ROOT input, and writes:

- `Trigger_windows.root`: every original board1049 window, including empty,
  unreferenced and boundary windows, with native timestamp and support flags.
- `NCAL_ch02_pulses_part*.root`: every cleaned NCAL pulse, not a sample;
  explicit phase validity, original timing, ToT and 1P/3P phases.
- `All_channel_plot_replay.root`: weighted-bin **TTrees** reproducing every
  channel profile and all detailed maps exactly. These are binned replay data,
  not substitutes for unbinned event trees or arbitrary finer reanalysis.

Pulse parts contain at most five million rows, with a checked 250 MiB/file
transfer limit. An entire NCAL pulse set is copied only if it fits the declared
total limit; never transfer an unlabeled partial set. Full cleaned inputs and
complete generated products remain on padiwa.

## Native five-metre position run

The separate [5 m comparison](<../Results/Position Runs/Zebra Period Finder/NCAL_20us_Pos_5m_0000/>)
uses the canonical `NCAL_20us_Pos_5m_0000` acquisition, with numbered
decode, cleaning, rates, timing/ToT and ch0-clock directories. It does not
replace the original position-run products or the Zebra/pulser settings.
Its [frozen remote launcher](<../Results/Position Runs/Zebra Period Finder/NCAL_20us_Pos_5m_0000/Archive/Run scripts/run_five_m_chain.py>)
supports explicit start/stop steps and requires verified pilot/prerequisite
steps before full processing.

The local browsing layout contains five final `step - description` folders.
Pilots, rejected attempts, scripts and technical logs are grouped under one
`Archive`; the padiwa processing layout is unchanged. Archived scripts retain
remote paths and are not intended to run from their relocated local folders.

Native settings are hardware `TDC_001049`, ordinal 1, coarse quantum 5 ns,
assumed prescaler 1024, first passing ch0 rise -410 to -395 ns and ToT
16.5 to 19.5 ns. Both timing axes are -6000 to 6000 ns; detailed timing
bins are 5 ns; displayed ToT is 0 to 128 ns in 0.1 ns bins. Cleaning
retains 1 to 5000 ns ToT. No detector-derived period scoring or phase
alignment is performed. All 32 channels retain profiles; detailed maps
cover channels 0 to 8.

The original single-clock method rejects the real pilot: averaging across
irregular acquisition gaps produces 38.72 to 39.51 ns local estimates and
about 199 ns maximum half-interval disagreement at a 6000 ns delay.
That attempt and its source snapshot are preserved. Its rejection is not
turned into a passing reconstruction by weakening the original checks.

The subsequent [all-trigger timestamp investigation](../Reference/prague_trigger_envelope.md)
establishes an approximately 150 Hz, approximately 20% recorded-trigger
envelope in this 5 m run. Ordinary inter-burst gaps are a recurring data
feature, not evidence by themselves of lost/corrupt windows or phase resets.
The regional fold below remains a support-selected diagnostic; its symmetric
boundary exclusions are not detector-quality failures and must not define a
full-population estimator comparison. No folded products were changed by
the timestamp investigation.

The explicitly selected `--split-clock-gaps` variant uses the same
neighbor-interval formula only inside supported contiguous regions:

- Split at a missing original reference ID, or when a consecutive interval
  differs from the median interval by more than one declared coarse quantum
  (5 ns here). Do not infer integer tick counts across such gaps.
- Require at least 32 references per region. For half-span 4, require full
  symmetric support through half-span 8, leaving eight references at each
  region boundary uncalibrated.
- Require half-interval disagreement, half-span 2/4/8 sensitivity, and the
  conditional three-sigma error from assumed 0.1 ns edges each to be at
  most 0.1 ns at the maximum 6000 ns folded delay.
- Record every unsupported reference ID and separately account for every
  pulse excluded from folding. Cleaning and upstream timing products retain
  those pulses.

`clock_regions.tsv` records region IDs/support; `assessment.json` records
the exclusion reasons. The endpoint period combines supported regional
endpoints, not the gap-contaminated full wall-clock span. Tick counters
refer only to supported regions; the missing ticks across gaps are
**unknown**, not measured as zero. The default Zebra single-clock path and
the pulser multi-cadence path remain unchanged.

The [independent verifier](<../Results/Position Runs/Zebra Period Finder/NCAL_20us_Pos_5m_0000/Archive/Run scripts/verify_step04.py>)
reconstructs support, periods and every profile/detailed-map bin from the
cleaned ROOT export. Additional controls compare the same supported
population with one fixed ch0-derived period, a deterministic within-channel
ToT permutation, and 2x2 rebinning. A separately labelled constant-period
panel uses **all** referenced pulses to expose selection effects.
Permutation preserves both displayed one-dimensional marginals; it is
not a proof of a physical origin.

The legacy decoded export was independently compared line by line:
it contains exactly two consecutive copies of the fresh export, ignoring
only process pointer addresses in TDC headers. The old cleaner stitched
the timestamp restart forward, and also processed all boards. The new
chain uses one fresh decode and explicitly selects board 1049. Consequently,
the old/new results are **not** a strict one-variable period-estimator
comparison. The raw recording also ends with a record missing 292 bytes;
all 10,521,678 complete events are decoded. Neither numerical QC nor
persistent zebra lines certify instantaneous RF phase or a beam-physics
origin.

This is a separate alternative to the physics-feature RF scan. It does not edit
or invoke the existing Step04 scanner, and does not require a gamma peak or use
detector pulses to estimate a period.

## Timing reconstruction

The frozen Zebra decoder uses triggered analysis, `SetCustomMhz(150.)`, and
floating hit times relative to its DOGMA trigger counter. The actual counter
conversion in `hadaq::TdcProcessor::DoBufferScan` is:

```text
trigger_time = ((epoch << 12) | (coarse << 1)) * 1000 / clock_mhz
relative_hit_time = calibrated_hit_time - trigger_time
```

Thus the trigger anchor lies on a `2000 / clock_mhz` ns lattice. The original
decoded export rounds that anchor to 1 ns. Recovering its integer counter from
this known lattice, then adding the cleaned ch0 fine offset, undoes the relative
representation without introducing the anchor's formatting error as new jitter.
Integer counters are differenced before conversion to avoid subtracting large
absolute floating-point timestamps.

The six-decimal trigger-seconds column in the cleaned table is **not** used for
precise reconstruction. Original decoded metadata supplies the counter anchor
only; selected hits still come from cleaned data. Every selected rising/falling
edge, event index, TDC identity, and rounded cleaned trigger timestamp is checked
against the original export. Lattice mismatches fail explicitly.

The lattice above is specific to the verified Zebra decoder configuration. It
must not be applied to arbitrary exports just because their timestamps look similar.
The supplied clock frequency is a decoder time-scale assumption, not an
independently calibrated hardware frequency. The separately verified native
Prague configuration is described below; it does not use this Zebra lattice.

## Estimation and acceptance

- AE and AF are fitted independently.
- Select the first valid ch0 in file order, rise -392 to -374 ns and ToT
  16.5 to 19.5 ns. Multiple guarded candidates require assessment and are rejected
  by this runner. The delayed second pulse is not assumed to be an RF cycle.
- Infer missing prescaled ticks from integer interval multiples.
- Measure the full-span mean prescaled period directly from elapsed ch0 time
  and the number of prescaled ticks; divide by the explicitly supplied
  prescaler, 1024 for this test.
- For local folding, average intervals spanning four neighboring events on each
  side (eight prescaled intervals, approximately 0.32 ms). At file boundaries use
  the available one-sided span. There is no merit optimization or detector
  dependence.
- Report sensitivity to half/double the averaging span, and disagreement between
  the first and second halves of each interval. These expose source timing
  variation rather than interpreting it as measurement noise.
- Propagate the explicitly assumed 0.1 ns ch0 edge precision through each
  interval. Require this measurement-error component's three-sigma folded-phase
  displacement over an 11 us delay to be below 0.1 ns. This does **not** certify
  instantaneous RF phase between unobserved RF cycles.
- Check the local-period average against the independent full-span mean.
- The 35-45 ns range is a validity guard, not a search objective or forced fit.
  Out-of-range data fail.

Conditional edge-error estimates exclude clock-scale/calibration systematics and
unobserved source phase noise between prescaled pulses. The data show appreciable
prescaled timing variation; do not claim 100 ps instantaneous RF-phase accuracy
from 100 ps hit resolution alone. No uncertain windows or pulses are silently
discarded. Both boards share the
DAQ timing chain; agreement does not constitute independent absolute-frequency
traceability.

## Folding and outputs

The phase origin is the selected ch0 itself. There is no detector-based phase
fitting or segment realignment. For each event:

```text
cycles = (channel_rise_ns - selected_ch0_rise_ns) / local_period_ns
phase_1x = (cycles modulo 1) * board_mean_period_ns
phase_3x = (cycles modulo 3) * board_mean_period_ns
```

The fixed displayed axis is in local cycles scaled by the board mean period.
The 3x product is a real modulo-3P fold, not three pasted copies.

All 32 channels have folded profiles. Detailed phase-ToT, phase-ch0-time, and
phase-trigger-time maps cover channels 0, 1, 2, 3, and 6. Time axes remain
-11000 to 2000 ns; ToT maps display 0-80 ns. Long-ToT pulses remain in profiles
and are explicitly accounted for outside the displayed ToT maps. Ch0 profiles
contain one selected reference per event, matching existing reference semantics.

The runner verifies input/profile totals, sparse bins, ROOT bins, and PNG/PDF
integrity before writing a passing assessment.

## Execution

Edit on the Mac; compile, test, and process data only on padiwa. Use a new output
directory: the runner refuses to overwrite prior outputs.

```bash
python -m unittest discover -s tests -p 'test_zebra_ch0_clock.py' -v

python run_zebra_ch0_clock.py \
  --input /data6/.../Refined_Data/...AE_rawRefined/...AE_rawRefined_pulses.tsv \
  --decoded /data6/.../Step00_Decode/...0115.dld.dat \
  --output /data6/.../Step04_Ch0Clock/Results/...AE \
  --board AE --tdc-clock-mhz 150 --prescaler 1024 \
  --build-dir /data6/.../Step04_Ch0Clock/build
```

Repeat for AF, sequentially. Set single-thread BLAS environment variables for
bounded remote resource use. Keep large reference/window-period tables, ROOT
files, and sparse data products on padiwa.

Small verified artifacts can be copied to:

```text
Dogma_Analysis/Results/Zebra Runs/run_20261009_0115/step04_ch0/
    AE/
    AF/
    Reports/
```

Keep the earlier `step04/` results unchanged. Exclude every ROOT file and each
file over 10 MiB, and verify copied files with source SHA-256 checksums.

## Native Prague external-trigger pulser

The Prague pulser uses the same ch0-only estimator, with explicit native
overrides. The general runner retains the original Zebra gates, axes, and
detailed-channel defaults. The run-specific
[Prague launcher](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/run_pulser_step04.sh>)
also retains the Zebra settings as comments.

| Setting | Zebra baseline | Prague pulser override |
|---|---|---|
| Hardware / selected ordinal | AE / 1, AF / 2 | TDC_001049 / 1 |
| Trigger-anchor quantum | 2000/150 ns | Native 5 ns; no rescaling of decoded hits |
| Prescaler | 1024 | 4096, explicitly assumed by the user |
| First passing ch0 rise | -392 to -374 ns, width 18 ns | -410 to -395 ns, width 15 ns |
| First passing ch0 ToT | 16.5 to 19.5 ns, width 3 ns | 20.0 to 21.5 ns, width 1.5 ns |
| Trigger-relative display | -11000 to 2000 ns, 2600 bins | -700 to 700 ns, 1400 bins |
| Ch0-relative display | -11000 to 2000 ns, 2600 bins | -300 to 1100 ns, 1400 bins |
| ToT display | 0 to 80 ns, 800 bins | 0 to 128 ns, 1280 bins |
| Detailed channels | 0,1,2,3,6 | 0,1,2,3,4,5,6,7,8 |
| Phase-bin-width target | 0.25 ns | 0.25 ns; actual width follows measured period and integer bin count |
| Neighbor half-span | 4 | 4; eight regular intervals span approximately 1.31 ms |
| Assumed ch0 edge precision | 0.1 ns | 0.1 ns; conditional error checked over a 1100 ns delay |
| Reference-gap handling | Strict by default | Explicit `--allow-reference-gaps`, with original IDs retained |

Supply `--coarse-quantum-ns 5` instead of `--tdc-clock-mhz`; these options are
mutually exclusive. Non-Zebra hardware labels require `--tdc-ordinal`. The
acquisition filename still contains `Prescaler_4`: it is not renamed, and does
not independently confirm the working divider assumption.

This export begins with board-0-only events. Cleaning numbers windows by
**selected-board headers**, so the decoded matcher must do the same. Counting
all board-0 headers would incorrectly associate cleaned pulses with earlier
events. Original high-precision decoded trigger anchors, both selected edges,
hardware identity, and cleaned rounded trigger values are checked before a
reference timestamp is accepted. The delayed second ch0 family remains in
cleaned data; only the selected first passing pulse contributes to ch0 folding.

The full cleaned run has 14,374,543 selected-board windows and 14,374,171
references. Its 372 reference gaps are recorded individually: one empty
cleaned window, 321 nonempty windows without accepted ch0, and 50 nonempty
windows without a passing ch0 candidate. No period or phase is fabricated for
these windows. Their pulses remain in prior cleaning, absolute-rate, and
trigger-relative products and are counted separately in folding accounting.
Nonreference ch0 pulses are a separate exclusion category. Missing references,
events without the selected board, and inferred missing prescaled clock ticks
are distinct quantities.

The [independent Prague verifier](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/verify_pulser_step04.py>)
compares exact gap IDs with Step01, streams every cleaned ROOT pulse, checks
original reference IDs, independently reconstructs endpoint frequency and
every 1x/3x profile and detailed 2D bin, and reconciles eligible counts with
Steps02-03. It also checks ROOT axes/bins and plot completeness/integrity.
Regression coverage includes the original Zebra tests and native Prague
clock, gap, first-valid file-order, and modulo-3P tests.

Compile, test, and run only on padiwa, sequentially with bounded threads and
low priority. Keep the frozen sources, builds, large reference/period tables,
and all ROOT products remote. Local small copies use:

```text
Dogma_Analysis/Results/Pulser Runs/run_<original_acquisition>/
    step04_ch0_clock/
        TDC_001049/
        Reports/
```

The full runner records source/input hashes, refuses existing outputs, and
requires a verified real pilot. Preserve the original 24-byte incomplete
raw-tail warning. Neither detector-feature fitting nor detector-based phase
alignment is used. A passing assessment remains conditional on the decoder
time scale and assumed divider; instantaneous RF phase between unobserved
prescaled edges is not certified.

### Full-run cadence constraint

The [Prague single-clock attempt report](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/attempt_01_single_clock/Reports/full/step04_report.txt>)
records the real pilot passing at 40.000185972235 ns
(24.999883768 MHz, conditional on divider 4096). These numbers describe
the **pilot**, not a certified full-run clock. The original attempt is preserved
under `step04_ch0_clock/attempt_01_single_clock/`.

Full extraction verifies every original/cleaned reference and all 372 gap
IDs/classes, but the estimator explicitly rejects the full clock sequence:

| Measured reference cadence | Intervals | Period if each edge represents divider 4096 |
|---|---:|---:|
| 163840.759033 ns | 14,315,609 | Approximately 40 ns |
| 81920.379883 ns | 57,799 | Approximately 20 ns |
| 40960.192017 ns | 681 | Approximately 10 ns |

The faster populations include sustained blocks and adjacent reference-window
IDs. They are not accounted for by the 372 reference gaps. The observed
cadences alone cannot distinguish changes of oscillator frequency from
changes of divider. No new divider values are assumed, no frequency is
forced, and no faster blocks are silently removed. The original single-clock
acceptance rule is retained by default. The explicitly selected fixed-divider
cadence model below handles this acquisition without changing that default.
Pilot products and full diagnostics are copied separately beneath
`step04_ch0_clock/TDC_001049/{pilot,full}/` and
`step04_ch0_clock/Reports/{pilot,full}/`, with source SHA-256 records.

### Explicit fixed-divider cadence model

The Prague launcher now supplies `--cadence-factors 1,2,4`. The physical
divider remains **4096 in every region**. These factors describe measured
frequency/cadence ratios relative to the dominant interval; they do not assign
new hardware dividers or force exactly 25, 50, or 100 MHz.

Intervals must agree with one of the declared cadence ratios within 1%.
Integer longer intervals can represent missing ticks. A longer interval
across missing reference IDs inside a matching faster region is treated in
that region rather than inventing a one-interval slow-frequency block.
Unclassified intervals fail explicitly.

Each contiguous cadence region is measured independently, using the existing
neighbor-interval estimator and its acceptance checks. No averaging crosses
a cadence boundary. Regions must contain at least **32 references**. A
reference shared by different cadence regions, or belonging only to a
shorter unsupported region, has **no assigned period**. Its original ID is
listed in `Clock_Calibration/unresolved_clock_windows.tsv`.

These are explicit folding exclusions, not cleaning changes. All pulses in
such windows are counted under `without_clock_pulses`; missing-reference
pulses and nonreference ch0 pulses remain separate categories. The verifier
reconstructs the boundary/support mask and local periods independently,
then reconciles:

```text
folded + no_reference + unresolved_clock + nonreference_ch0 = cleaned input
folded eligible + eligible unresolved_clock = Step02/03 relative eligibility
```

`cadence_regions.tsv` records original IDs, interval/tick counts and measured
regional endpoints. `cadence_regimes.tsv` reports the measured periods and
conditional frequencies separately. The overall mean period is only the
fixed **mixed-regime display normalization**, not one physical oscillator
period. All 1x/3x profiles still use each resolved window's local period;
3x is actual modulo 3P. Detailed NCAL maps, axes, first-valid gates, source
hashing, raw-tail warning, and native decoded-hit units are unchanged.

The broadened backend period guard is explicit for this mode (8.75-45 ns
for factors 1,2,4); Zebra's default guard remains 35-45 ns. Changing hardware
dividers remains an alternative physical explanation. These outputs are
conditional diagnostic folds, not a certificate of instantaneous RF phase.

### Numerically verified conditional Prague result

The [full report](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/Reports/full/step04_report.txt>)
and [measured regime table](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/TDC_001049/full/Clock_Calibration/cadence_regimes.tsv>)
record the completed conditional model:

| Cadence factor | Measured period (ns) | Conditional frequency (MHz) | Supported references |
|---|---:|---:|---:|
| 1 | 40.000186025 | 24.999883735 | 14,315,688 |
| 2 | 20.000093003 | 49.999767494 | 57,796 |
| 4 | 10.000046467 | 99.999535332 | 652 |

All 14,374,171 original/cleaned references match. Thirty-five transition or
short-region references remain explicitly unresolved. Their 105 pulses remain
upstream; excluding 70 previously relative-eligible pulses yields 38,871,689
folded pulses. Adding those 70 reproduces the earlier 38,871,759 eligible count.
Every cleaned pulse, both fold factors, detailed bins and ROOT axes reconcile.
The 372 missing-reference IDs/classes are unchanged.

The complete-event count is 14,374,728, but there are only 14,374,712 board-0
headers: 16 complete events lack board 0. Complete events, individual TDC headers,
185 events lacking the selected board, 372 reference gaps, 35 unresolved clock
references, and 84 inferred missing prescaled ticks are separate quantities.

Twenty regression cases, the real pilot, full independent cleaned-ROOT checks,
and source/input provenance checks pass. The local copy contains 555 selected
artifacts (32,771,849 bytes): 182 full analysis PNG/PDF pairs, two report-figure
pairs, and 50 separate pilot pairs, plus small tables/logs/provenance records.
Every selected artifact is SHA-256 verified; all ROOT and bulk timing tables stay
on padiwa. The earlier single-clock attempt remains byte-identical in its archive.

### Pulser views and the limit of numerical verification

The earlier PASS certifies bookkeeping, implementation and output consistency;
it does not show that the oscillator changed between 25, 50 and 100 MHz.
Observed ch0 arrival cadences alone cannot distinguish oscillator changes from
divider/trigger/acquisition changes. The existing folded display also uses
39.918408207808 ns mixed normalization, which is not one physical oscillator period.
Zero delay for the selected ch0 is by construction, not measured reference jitter.

The supplementary [pulser diagnostics](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/TDC_001049/full/Pulser_Diagnostics/>)
show measured delays, separate cadence membership and side-by-side comparisons
with the original conditional 3P fold. The original ROOT bins are cropped,
not refolded; full-range originals remain unchanged. The compact ch2 population
already exists before folding, at median delay 60.582 ns and ToT about 17.96 ns.
No referenced ch2-8 pulses occur in either supported faster cadence.

Use [run_cleaned_pulser_diagnostics.py](run_cleaned_pulser_diagnostics.py) only on
padiwa. It requires `--cleaned-root-dir` (all cleaned ROOT parts),
`--clock-result-dir` (the verified Step04 product directory), `--verification`
(the independent full Step04 JSON), and a new `--output` directory. It refuses
overwrites and reconciles every cleaned pulse against verified Step04 categories.
Detailed channels 0-8 must fit the preserved -300..1100 ns delay and 0..128 ns
ToT axes; an out-of-axis pulse fails explicitly.

The 27 PNG/PDF pairs include both full-population zooms and seven main-pulse
views for ch2-8. Main-pulse views use a declared **40-80 ns measured-delay display
gate, inclusive, with no ToT gate**. Late populations remain in full-population
views and upstream inputs. For ch4 and ch8, their combined central-68% halfwidths
near 37 ns reflect separated populations and must not be called timing jitter.
The gated descriptive halfwidths are 0.140 ns and 0.1375 ns, respectively.
No Gaussian fit or electronics/detector-resolution claim is made.

See the [supplementary report](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/Reports/pulser_diagnostics/pulser_diagnostics_report.txt>)
and [independent diagnostic QC](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/Reports/pulser_diagnostics/independent_verification.json>).
The additional 65 small artifacts (4,952,829 bytes) are copied and source-SHA
verified separately from the original 555-file copy. All ROOT, bulk tables and
frozen source snapshots remain remote. These diagnostic improvements do not
certify physical RF phase or authorize downstream physics analysis.
