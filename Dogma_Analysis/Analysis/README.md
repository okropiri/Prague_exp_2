# DOGMA Analysis Scripts

## Working Model

- Edit code locally in the mounted workspace.
- Run full processing on `padiwa` against the real `/data6` filesystem.
- Keep the active chain in step order so the cleaned-data workflow is readable from filenames alone.
- Keep all small local result copies under [Results](../Results/), grouped into [Position Runs](<../Results/Position Runs/>), [Pulser Runs](<../Results/Pulser Runs/>), [Source Runs](<../Results/Source Runs/>), and [Zebra Runs](<../Results/Zebra Runs/>). Remote `/data6` result paths are unchanged.
- The separate [5 m Zebra Period Finder run](<../Results/Position Runs/Zebra Period Finder/NCAL_20us_Pos_5m_0000/>) uses ch0-only timing, not detector-feature scoring. Its [clock method](zebra_ch0_clock.md#native-five-metre-position-run) records the rejected single-clock attempt and the explicit regional support exclusions. Original position-run outputs are preserved.
- Experiment-specific timing context is summarized in `../Reference/experiment_context.md`.
- The canonical 5 m run has a [measured approximately 150 Hz, approximately 20% trigger envelope](../Reference/prague_trigger_envelope.md), verified from all original board1049 headers before reference cuts. Ordinary inter-burst gaps and burst-boundary windows must not automatically be treated as lost/bad data or RF phase resets. This result is run-specific; other acquisitions require their own checks.
- The [full-population `--burst-clock` mode](zebra_ch0_clock.md#full-population-burst-aware-mode) fits a shared ch0-only period inside bursts, with independent regional intercepts and explicit stability checks. It assigns a period to every passing-reference window, retains all ch0 hits, and never infers cycles through inactive gaps. [Independent QC and ROOT TTree export](verify_ch0_burst_fold.py) preserve all windows and every cleaned NCAL pulse, including explicit invalid-phase rows.
- A human-oriented explanation of the NCAL RF scan is in `../Reference/ncal_rf_phase_scan_walkthrough.md`.
- Activate the existing Python environment before running any helper scripts:

```bash
conda activate conda_env
```

## Active Cleaned-Data Chain

### Controlled 5 m RF-method comparison

[run_matched_five_m.py](run_matched_five_m.py) reproduces shared Steps01-03
from the hash-verified, single-copy native decode, then estimates the period
from either ch0 burst timing or ch2 detector structure. Both estimates feed
the **same** unshifted C++ folding implementation. The selected guarded ch0
is the origin; neither a global nor a regional detector-derived offset is
subtracted. The original detector-aligned Step04 remains available unchanged
for historical reproduction. Its new `--estimate-only` mode preserves the
period search but does not produce aligned folded data.

Execution is restricted to padiwa. The numbered browsing steps have shared
upstream products, two labelled Step04 branches and separately labelled
downstream outputs. Step08's unfolded ToT fit is shared. Step00 stays remote
and is reused only after its original-input/decoded-export identity and board
inventory are verified. The strict Step01 `--reject-trigger-restarts` guard
refuses the historical repeated-recording restart rather than stitching it.
This is an input-integrity guard, not hit-value deduplication.

The comparison uses **0-60 ns ToT displays** while retaining the 0-128 ns
stored map axes, full cleaned pulse widths and unchanged analysis gates
(including background sidebands beyond the display). Step03 provides
`--tot-map-display-max-ns` separately from
`--tot-distrib-display-max-ns`; ROOT contents are not cropped.
The legacy backend accounting name `tot_outside_display` refers to the
configured stored ToT-map range (0-128 ns here), not the new rendering crop;
the equal-bin NCAL comparison separately reports pulses at/above 60 ns.

[analysis_provenance.py](analysis_provenance.py) records resolved arguments
(including defaults) when `DOGMA_STEP_REPORT_DIR` is set by the tracked
workflow. Each step retains settings, exact backend commands, logs, exit
status and hashed outputs. The run manifest identifies the actual dirty
source snapshot, compiler, packages and inputs. Existing failed or modified
outputs are never silently resumed as successful results. No historical
settings are reconstructed by guessing current defaults.

Gamma-relative downstream diagnostics remain explicitly distinguished from
the unshifted ch0-referenced RF products. Accepted-cycle normalization still
uses a declared 10000 ns acquisition-window **assumption**, not the plot
width or a 20% duty correction. A 5 m-only run is not a cross-position
neutron-template/background validation; unavailable controls are warned.

The [completed canonical matched 5 m results](<../Results/Position Runs/NCAL_20us_Pos_5m_0000/README.md>)
include both full Steps01-10 branches, the equal-bin comparison, ROOT replay
instructions and verified copy/archive manifests. Both methods fold the same
31,159,659 NCAL pulses; 58,053 missing-reference pulses remain in the shared
31,217,712-entry pulse set with invalid phases. The measured periods are
38.8500713060152 ns (ch0) and 38.84950000000001 ns (ch2), in nominal DAQ units.
The common-bin comparison uses the same 0.25 ns phase and 0.1 ns ToT edges.

- `step01_dogma_raw_data_refiner.cpp`: raw DOGMA export to cleaned pulse table. This is the gatekeeper step; downstream analyses are expected to consume its cleaned pulse output.
- `step01_run_dogma_raw_data_refiner.py`: step-01 wrapper that compiles and runs the cleaner and can hand the cleaned pulse table to the ROOT writer.
- `step01_write_cleaned_root_tree.py`: converts the cleaned pulse table into the split `CleanedHits` ROOT export plus `RunMetadata` and the parts manifest.
- `step02_dogma_cleaned_all_channel_rates.cpp`: cleaned-data all-channel rate analyzer.
  - Absolute-rate timing starts at the earliest accepted pulse across the table, not the first row; channel-grouped file order must not discard earlier pulses in the first window. This requires a preliminary input pass.
- `step02_write_cleaned_all_channel_rates_outputs.py`: writes the cleaned absolute-rate and ch0-referenced rate plots plus ROOT outputs.
- `step02_run_cleaned_all_channel_rates.py`: step-02 runner and shared helper module used by later cleaned-data steps.
- `step03_dogma_cleaned_all_channel_time_tot.cpp`: cleaned-data all-channel time-vs-ToT analyzer in trigger and guarded ch0 frames.
- `step03_write_cleaned_all_channel_time_tot_outputs.py`: writes step-03 trigger-ref and ch0-ref time-vs-ToT plots plus ROOT outputs.
- `step03_run_cleaned_all_channel_time_tot.py`: step-03 runner.
  - `--tot-min-ns`, `--tot-max-ns` and `--tot-bins` set the 2D ToT axis; out-of-range pulses are counted in the summary. Whole-run ToT tables/ROOT retain the observed upper tail. Use `--tot-distrib-display-max-ns` to bound only their PNG/PDF display, with excluded-bin counts and percentages shown. Prague defaults keep 0-128 ns. For Zebra 0115, use `--tot-min-ns 0 --tot-max-ns 80 --tot-bins 800 --tot-distrib-display-max-ns 80` for useful 0-80 ns plots without changing cleaned data.
- `step04_dogma_cleaned_all_channel_rf_period.cpp`: cleaned-data RF-period deduction and all-channel RF folding analyzer.
- `step04_write_cleaned_all_channel_rf_period_outputs.py`: writes step-04 ch0-ref RF timing, folded 1x profiles, folded 3x profiles, and folded phase-vs-ToT products.
- `step04_write_rf_period_scan_outputs.py`: writes the dedicated RF-period-scan diagnostics in the old-style merit/sigma/residual format from the step-04 scan text outputs.
- `step04_run_cleaned_all_channel_rf_period.py`: step-04 runner.
- [extract_native_trigger_timestamps.cpp](extract_native_trigger_timestamps.cpp) and [run_trigger_timestamp_envelope.py](run_trigger_timestamp_envelope.py): unbiased native-Prague board1049 trigger-envelope diagnostic, run only on padiwa. The extractor retains every selected-board header, including empty headers; the runner measures dense/sparse groups, absolute-time distributions, blockwise repetition and duty proxies without changing detector-hit selection. Grouping/density sensitivity, rare pauses and hardware-interpretation limits are reported explicitly.
- `step04_alt_trigger_referenced_rf_modulation.py`: parallel RF-modulation test for runs with no channel-0 RF reference; folds trigger-referenced hit times instead.
  - For very long/pathological runs, `--score-stride`, `--single-rf-segment`, and `--skip-cycle-residual-diagnostics` can recover the RF summary and full-statistics folded outputs while keeping optional diagnostic memory bounded.
  - Standalone Step04 defaults are the original Prague real-data settings: ch0 rise -410 to -395 ns (reference ToT 16.5-19.5 ns), both timing axes -6000 to 6000 ns, ToT 0-128 ns in 0.1 ns bins, scoring channel 2 with scoring ToT 0-128 ns, score stride 1, automatic per-segment phase alignment, and folded 2D storage for the scoring channel only. The Zebra 0115 two-generator values are kept as comments beside every Python/C++ default so they can be restored for pulser tests: ch0 rise -392 to -374 ns, timing axes -11000 to 2000 ns in 2600 bins, ToT 0-80 ns in 800 bins, score stride 10, `--single-rf-segment`, and `--histogram-channels 0,1,2,3,6`.
  - The default coarse period range is the Prague 37-41 ns window at 0.05 ns steps, followed by the existing local refinement. The Zebra exploratory range was 5-200 ns at 0.25 ns steps. This is not a measured generator period and does not cover microsecond generators. A coarse grid can miss narrow optima: inspect aliases, scan edges and flat/degenerate scores, then confirm a candidate with a narrower scan and `--score-stride 1`. The default score stride is 10 for search cost only; final histograms still use all eligible pulses.
  - `--single-rf-segment` uses one global fitted phase origin instead of automatic per-segment alignment, and was the Zebra setting. Either way this is peak-recentered folding, not unshifted modulo timing. `--auto-rf-segments` (the default) is the Prague automatic alignment.
  - RF scoring and saved 2D channels are independent: use e.g. `--score-channel 0 --histogram-channels 0,2,6` only when additional ch0 pulses genuinely supply periodic information. Explicit `--histogram-channels` replaces the default list. Selecting five detailed channels and 5 ns unfolded timing bins bounds wide-period phase-time matrix memory; expanding the range or channel list requires checking allocation sizes.
  - When scoring ch0, the guarded reference hit itself is excluded from the search: its self-relative time is always zero and cannot constrain a period. Other accepted ch0 pulses are used. A narrow peak or a scan-edge optimum is not proof of the programmed RF period; verify repeated pulse spacing and scan diagnostics.
  - Folded phase-ToT exports record `stored_channels`. Steps05-07/09-10 reject a channel with no entries instead of drawing a success-shaped empty plot. Confirm board identity and input occupancy before running these steps.
  - `run_all_cleaned_analysis_chain.py` remains the real-data batch entrypoint and still passes the Prague Step04 settings explicitly, consistent with its Step02/03 gates. Standalone real-data Step04 now needs no extra flags; the explicit equivalent is:

    ```bash
    --ch0ref-min-ns -6000 --ch0ref-max-ns 6000 --ch0ref-bins 12000 \
    --trigger-ref-min-ns -6000 --trigger-ref-max-ns 6000 --trigger-ref-bins 2400 \
    --tot-max-ns 128 --tot-bins 1280 --score-channel 2 --score-tot-max-ns 128 \
    --histogram-channels 2 --initial-period-min-ns 37 --initial-period-max-ns 41 \
    --initial-step-ns 0.05 --score-stride 1 \
    --ch0-valid-rise-min-ns -410 --ch0-valid-rise-max-ns -395 --auto-rf-segments
    ```
- `step05_estimate_gamma_neutron_background.py`: first-pass channel-2 phase-vs-ToT sideband estimator for prompt-like ROI, broad neutron-candidate band, and local background control regions.
- `step06_quantify_gamma_neutron_components.py`: TOF-guided prompt-gamma, background-control, and neutron-candidate quantification step. It measures per-run centroids, widths, local/background-control regions, and first-pass background-subtracted rates while explicitly flagging gamma/neutron overlap cases.
- `step07_phase_rate_slice_diagnostics.py`: gamma-relative one-period phase-rate diagnostics by ToT slice, used to inspect delayed excess structure in 1D before committing to 2D population cuts.
- `step08_fit_gamma_tot_peak.py`: narrow Ch02 gamma-ToT peak fitter using a Gaussian plus constant background in a fixed ROI to extract peak position and its statistical fit uncertainty.
- `step09_gamma_window_sensitivity.py`: 5 m-only sensitivity scan for gamma-window definition. It centers fixed-width gamma ToT gates on the fitted peak from step08, rebuilds the gamma phase reference and 1D phase-rate diagnostics, and reports how much the reference shifts under window variations.
- `step10_gamma_reference_cfd.py`: reference comparison and CFD diagnostics; discovers profiles for the requested `--channel` and rejects missing/empty references.

## Runs Without A Channel-0 RF Reference

`step04_alt_trigger_referenced_rf_modulation.py` is a parallel tool for
acquisitions where the channel-0 RF reference pulse was never cabled into the
TDC, so the ch0-referenced step-04 chain cannot produce a phase at all. The
2026-04-27 `NCAL_Calibration_Source_Prescaler1024_*` control run is one such
file: TDC_001049 has zero channel-0 hits.

Because the DOGMA trigger is the cyclotron RF divided by 1024, every trigger
edge is RF-phase-locked and the hit times in the decoded `.dat` are already
measured relative to it. The tool therefore folds trigger-referenced hit times
directly, after dividing out the trigger-window acceptance, and excludes channel
0 so the RF reference cannot inject the signal by construction. It only reads
the decoded `.dat` and leaves all existing step code and outputs untouched.

The validated Prague legacy-format decoder is based on `stream` commit
`6a3b6e1`; the other inspected builds either refuse the non-TDC5 data or report
a zero `GlobalTriggerTime`. Its environment is
`/home/padiwa/compat_bin/apr2026_env.sh`.

The April 27 external-trigger pulser also exposed a first-buffer epoch bug:
the trigger timestamp was converted before the epoch reference was initialized,
producing an incorrect first-block time offset and losing precision in the
float hit timestamps. The isolated build
`/home/padiwa/stream_prague_pulser_build` initializes the reference from the
DOGMA trigger header before converting it. The original build is preserved.
The pulser run's step-00 artifacts retain the patch, environment, runner, and
pilot regression checks. Use that corrected environment for this pulser file;
do not repair the already-rounded float timestamps after export.
The [verified decode status](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step00_decode/decode_status.txt>)
also records that all complete events were exported, but the recording's final
event is missing 24 bytes. Preserve that input warning when using the export.

For this pulser's Step01, use the run-specific
[launcher](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step01_cleaning/run_pulser_step01.sh>).
It preserves the original Zebra reference gates as commented settings and uses
Prague pulser rise -410 to -395 ns and ToT 20.0 to 21.5 ns. Signal and trigger
both use TDC ordinal 1 (TDC_001049). Cleaning retains all valid pulse pairs with
ToT 1-5000 ns; the reference gates only identify the first ch0 reference in the
ROOT export, not a filter deleting the second ch0 pulse. General analysis
defaults are unchanged. Prescaler 4096 is a user assumption recorded for later
Step04, not a Step01 calculation or a measured hardware setting.

Execute only on padiwa. Deploy the launcher and its verifier together under
`Results/<pulser_run>/Step01_Cleaning_20261014/`, with matching copies of
`step01_run_dogma_raw_data_refiner.py`, `step01_dogma_raw_data_refiner.cpp`, and
`step01_write_cleaned_root_tree.py` in its `analysis_source/` subdirectory.
Run `bash run_pulser_step01.sh pilot`, inspect the verification, then run
`bash run_pulser_step01.sh full`. Each stage refuses existing output directories,
records source hashes and input file metadata, and verifies every cleaned ROOT
entry against the cleaner's channel accounting and the configured reference
gates. Empty cleaned windows and missing or ambiguous guarded references are
reported explicitly in `verification.json` and `reference_quality_windows.tsv`;
successful cleaning does not imply that every window is usable by Step04.
Do not silently drop or renumber these windows. Full pulse tables, refined
streams, ROOT parts and binaries stay remote.

The selected local pulser artifacts in [Pulser Runs](<../Results/Pulser Runs/>) use `step00_decode/` and
`step01_cleaning/`: numbered steps with short descriptions. As in the Zebra
results, Step01 separates board products (`TDC_001049/pilot/` and
`TDC_001049/full/`) from execution and verification records (`Reports/`).
Its `COPY_INFO.txt` explains the selection and `SOURCE_SHA256SUMS.txt` verifies
the local paths. This local organization does not rename remote outputs.

The pulser [Step02 launcher](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step02_rates/run_pulser_step02.sh>)
uses the same ch0 cuts, 0.1 s absolute-rate display bins, and a -300 to 1100 ns
ch0-relative axis in 1 ns bins. All accepted pulses contribute to absolute
profiles; ch0-relative profiles require a selected reference, with only that
reference contributing on channel 0. Missing-reference accounting is checked
against Step01. Run only on padiwa with a frozen Step02 source snapshot.
Selected small plots, tables and reports go to local `step02_rates/TDC_001049/`
and `step02_rates/Reports/`; ROOT outputs and cleaned input stay remote.
Rates are recorded accepted counts per acquisition time, not corrected detector
rates or an RF-frequency measurement.

The pulser [Step03 launcher](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step03_time_tot/run_pulser_step03.sh>)
uses trigger-relative time -700 to 700 ns and ch0-relative time -300 to
1100 ns, both in 1 ns bins. Its 2D ToT axis and distribution display are
0-128 ns in 0.1 ns bins; full ToT TSV/ROOT distributions retain the observed
upper tail. The same pulser reference gates select only the first passing
ch0 pulse in file order, not the minimum-time pulse. All accepted pulses,
including both ch0 families, contribute to the full width distribution.
Missing-reference windows still contribute to trigger-relative products.
Independent cleaned-ROOT verification reconciles map counts and separate
time-only, ToT-only, and joint overflows with Steps01-02, and checks every
stored ROOT bin and the generated image set. The frozen source snapshot must
include the Step02 Python runner imported by Step03. Tests in
`tests/test_cleaned_time_tot.py` cover the pulser doublet, both reference
gates, file-order selection, missing references, and full-ToT preservation
outside the 2D display; `tests/test_step03_tot_display.py` checks display-only
cropping. Execute only on padiwa. Small copies go to
`step03_time_tot/TDC_001049/` and `step03_time_tot/Reports/` under the local
pulser acquisition; ROOT and bulk input remain remote. Prescaler 4096 is
recorded but unused, and no RF reconstruction or folding occurs in Step03.

The verified full pulser [Step03 report](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step03_time_tot/Reports/step03_report.txt>)
records 53,226,331 pulses retained in the full width distributions and no
timing-axis overflow. Only 2,231 pulses (0.004192%) exceed the 128 ns display,
all outside NCAL ch2-8. The small local copy contains 123 PNG and 123 PDF.
Individual ToT plots annotate their own channel's display-tail counts;
combined plots annotate the all-channel counts. The final display-label
correction and six passing tests are recorded separately from the preserved
original source snapshot; neither cleaning nor 2D map construction was rerun.

For pulser Step04, use the separate
[ch0-only clock launcher](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/run_pulser_step04.sh>),
not the detector-feature RF-period scanner. The
[native override documentation](zebra_ch0_clock.md#native-prague-external-trigger-pulser)
lists the retained Zebra values beside the Prague settings. Prague uses
TDC_001049, selected ordinal 1, native 5 ns trigger-anchor quantum, and the
explicitly assumed divider 4096. The -410 to -395 ns / 20.0 to 21.5 ns
first-passing ch0 gates and Step03 display axes are unchanged. All 32 channels
have folded profiles; detailed maps cover channels 0-8, including NCAL 2-8.
The local axis is normalized by the measured ch0 period; the 3x product is
actual modulo 3P, not three repeated copies.

Decoded windows are matched using selected-board headers, not board-0 headers;
the export starts with board-0-only events. The 372 reference-gap IDs remain
explicit, with no fabricated periods. Their pulses and nonreference ch0 pulses
are separate folding exclusions, not new cleaning cuts. The
[independent verifier](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/verify_pulser_step04.py>)
reconstructs every profile/detailed-map bin from cleaned ROOT pulses, compares
exact reference/gap IDs with Step01, and reconciles folded eligibility with
Steps02-03. Source/input hashes, original-ID period tables, and the 24-byte raw
tail warning are retained. Small copies go to `step04_ch0_clock/TDC_001049/`
and `step04_ch0_clock/Reports/`; all ROOT and bulk timing tables stay remote.
Frequency remains conditional on the native decoder time scale and assumed
prescaler. Instantaneous RF phase between prescaled edges is not certified.

The preserved [single-clock attempt report](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/attempt_01_single_clock/Reports/full/step04_report.txt>)
records a passing 2,849-window pilot and all 14,374,171 full-run
original/cleaned reference matches, but **full clock calibration and folding
are blocked**, not accepted. In addition to the dominant 163,840.759 ns
reference cadence, the full export contains 57,799 intervals near
81,920.380 ns and 681 near 40,960.192 ns. Adjacent original reference IDs
and sustained faster blocks rule out explaining these populations solely
as missing references. Fixed divider 4096 would imply approximately
40, 20, and 10 ns source periods; a near-40 ns source with changed dividers
is another unverified interpretation. Choose an explicit model before
continuing rather than forcing one clock or silently excluding blocks.
The 372 missing-reference IDs/classes still match Step01 exactly.
All 136 selected first-attempt artifacts are SHA-256 verified and preserved
under `attempt_01_single_clock/`; its 50 PNG and 50 PDF are **pilot only**.
No full RF fold was produced by that attempt.
The archived [cadence diagnostic](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/attempt_01_single_clock/Reports/full/clock_interval_diagnostic.json>)
and [faster-block table](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/attempt_01_single_clock/Reports/full/fast_cadence_blocks.tsv>)
preserve the measured populations and original IDs.

The run-specific launcher explicitly enables the
[fixed-divider cadence-aware model](zebra_ch0_clock.md#explicit-fixed-divider-cadence-model):
divider 4096 throughout, measured cadence frequency factors 1,2,4, separate
fits within contiguous regions, and no averaging across transitions. Regions
with fewer than 32 references and shared transition references have no
assigned period; their IDs and pulse counts are explicit folding exclusions.
All upstream data and original reference IDs remain intact. The generalized
Zebra runner keeps its strict single-clock defaults unless this mode is
explicitly selected. The current full outputs use their local measured
periods and a mixed-regime display normalization, not one imposed 40 ns clock.

The [completed Step04 report](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/Reports/full/step04_report.txt>)
records a passing conditional cadence-aware full run and 20 passing regression
cases. Measured periods/frequencies under the fixed-divider assumption are
40.000186025 ns / 24.999883735 MHz, 20.000093003 ns / 49.999767494 MHz, and
10.000046467 ns / 99.999535332 MHz. All 14,374,171 references were verified;
14,374,136 have supported local periods. The 35 explicitly unresolved reference
windows contain 105 pulses, of which 70 were previously relative-eligible.
Folding retains 38,871,689 pulses; adding those 70 reconciles exactly with
Steps02-03's 38,871,759 eligible pulses. All 53,226,331 cleaned pulses reconcile
with folded, missing-reference, nonreference-ch0, and unresolved-clock categories.

Independent full cleaned-ROOT period/mask/bin, ROOT-axis, image-inventory and
source/input checks pass. Verification-only corrections accept the repository's
qualified cleaning-pass status and distinguish 14,374,728 complete decoded events
from 14,374,712 board-0 headers (16 events lack that header); no analysis products
were recomputed for these metadata-check corrections. Original failed logs and
frozen sources are preserved. Full products include 182 PNG and 182 PDF plus two
PNG/PDF report figures; pilot products are separate. All 555 selected artifacts
(32,771,849 bytes) were copied and source-SHA verified; all ROOT and bulk timing
tables remain remote. No later physics-analysis step was run.

### Pulser diagnostic presentation and interpretation

Numerical Step04 PASS checks are not a physical RF certificate. Prefer the
[measured-delay overview](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/TDC_001049/full/Pulser_Diagnostics/measured_delay_tot_overview.png>)
and [main-pulse zooms](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/TDC_001049/full/Pulser_Diagnostics/>)
for pulser characterization, rather than treating a compact spot on a broad
folded axis as evidence of an RF-correlated population.

[run_cleaned_pulser_diagnostics.py](run_cleaned_pulser_diagnostics.py) streams the
cleaned ROOT parts on padiwa, separates observed ch0 arrival cadences, retains
both ch0 families, counts missing/unresolved references, and writes measured
delay-ToT views without modulo or period rescaling. Its explicitly labeled
main-pulse display gate is **40 <= ch0-relative delay <= 80 ns, with no ToT gate**;
this is not a cleaning/folding cut. Full-population views preserve late/low-ToT
populations, which can dominate a combined width (notably ch4 and ch8).
Display-zoom tails and gate exclusions are counted, not deleted.

The [diagnostic report](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/Reports/pulser_diagnostics/pulser_diagnostics_report.txt>)
and [main-pulse table](<../Results/Pulser Runs/run_NCAL_Calibration_Pulser_ExternalTrig_Prescaler_4_TrigWindow_0x80800080_20260427_0000/step04_ch0_clock/Reports/pulser_diagnostics/main_pulse_display_summary.tsv>)
record the added checks. All 2,590,711 referenced ch2-8 pulses belong to the
dominant ~163.840762 us cadence; neither faster supported cadence contains
referenced ch2-8 pulses. This establishes membership, not the cause of cadence
changes. Ch2's main-display cohort has median delay 60.582 ns and descriptive
central-68% halfwidth 0.179 ns; these are not a Gaussian fit or detector resolution.

Seven targeted tests, including the complete CLI/image inventory, pass.
Independent exact cleaned-pulse quantiles agree with the 0.05 ns histogram
estimates within 0.003559 ns; main-display cohort counts/medians/widths are checked
independently. The additional local copy contains 65 SHA-verified small files
(4,952,829 bytes), including 27 PNG/PDF pairs. The original 555 Step04 and 136
archived artifacts are unchanged. The remote diagnostic stage is
`Step04_PulserDiagnostics_20261015`; no physical RF sign-off or later physics step
is implied by the diagnostic PASS.

## Zebra Mock Tests

For a ch0-only clock estimator that does not require an RF-correlated gamma
feature, use the separate [parallel Zebra calibration](zebra_ch0_clock.md).
It leaves the existing step code and earlier outputs untouched, measures each
board independently, reports local averaging-span sensitivity and timing limits,
and writes to a new `step04_ch0` result directory. Execute it only on `padiwa`.

Do not infer the wired board from a TDC's mere presence in a decoded file. Keep
hardware IDs explicit, compare decoded edge counts with cleaner validation, and
confirm the actual wiring. Never pool boards or relabel a populated channel as
channel 2 to obtain a plot. Sparse counts already present before cleaning are not
recovered by widening downstream cuts. Gamma/neutron labels, distance models,
and NCAL exposure assumptions are not calibrated physics for mock data.

Regression tests for the RF-score/output-channel separation and empty-input
guards are in `tests/test_rf_channel_pipeline.py`. Run them on `padiwa`:

```bash
python -m unittest discover -s tests -p 'test_rf_channel_pipeline.py' -v
```

The batch runner always requests channel 2's 2D products alongside its RF scoring
channel, since its comparison collection consumes channel 2.

`tests/test_mock_chain_integration.py` runs steps01-10 and CleanedHits export on a
synthetic known-count, 40 ns RF fixture. It checks detector occupancy and finite
CFD references, and removes its temporary products. Run it on `padiwa`, not on
the Mac:

```bash
python -m unittest discover -s tests -p 'test_mock_chain_integration.py' -v
```

## Utility Files Kept Nearby

- `README.md`: this file.
- `Prague_Cyclotron_Exp_2.code-workspace`: local VS Code workspace file for opening this analysis folder.

## Archived and Legacy Code

- Replaced raw-analysis code has been moved into `Archive/` with `legacy_` prefixes.
- Old mock-data generators, stress-test helpers, and local ROOT writer sanity checks have also been moved into `Archive/` so the top-level `Analysis/` directory only keeps the active production chain.
- Older exploratory trigger-reference, edge-check, and CSV-digest utilities stay in `Archive/` because they are not part of the current cleaned-data production chain.
- See `Archive/README.md` for the archived-file summary.
