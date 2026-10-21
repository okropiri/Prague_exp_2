# Prague RF-triggered control run — zebra pattern test

> **RETRACTED CONCLUSION (2026-10-13):** This run is not a valid control for
> deciding whether the folded zebra structure originates in the beam, in the
> RF-conditioning chain, or inside DOGMA. The calibration run has no usable
> channel-0 prescaled RF reference, so the same ch0-referenced folding that
> reveals the zebra pattern cannot be performed. The ch0-free,
> trigger-referenced modulation comparison below answers a different question
> and must not be interpreted as proof that the zebra structure is beam
> physics. A valid control must preserve the complete chain
> `25 MHz RF -> fan-in/fan-out -> discriminator -> prescaler -> DOGMA ch0`
> while removing or decorrelating only the detector/beam signal.

**Run key:** `NCAL_Calibration_Source_Prescaler1024_TrigWindow_0x83e803e8_20us_20260427_0000`
**Date analysed:** 2026-10-13
**Remote results:** `/data6/Dogma_analysis_by_Dachi/Results/NCAL_Calibration_Source_Prescaler1024_TrigWindow_0x83e803e8_20us_20260427_0000/`

---

## 1. The question

The Zebra two-generator pulser tests did not reproduce the diagonal "zebra"
stripes seen in the Prague beam data. That left two possibilities:

- the stripes are **real beam physics** (cyclotron RF bunching), or
- the stripes are an **electronics or analysis artifact** that the Zebra setup
  simply failed to reproduce.

To separate these, we needed a Prague acquisition that **keeps the RF-derived
trigger but removes the RF-correlated beam physics**. Exactly one such file
exists in the whole 137-file Prague inventory: the 2026-04-27 radioactive-source
calibration run taken with the `Prescaler1024` RF trigger. Everything else is
either a beam run or an internal-trigger run (no RF at all).

---

## 2. Answer

**The 38.851 ns modulation is present in the selected Prague beam file and is
not detected in the selected radioactive-source calibration file.** This is
strong evidence that it is associated with the beam acquisition rather than an
unavoidable RF-trigger/DAQ effect. Because these are different acquisitions on
different days, this comparison alone does not prove that every possible
electronics or run-condition difference has been excluded.

| | `NCAL_Test_..._Pos_5m_20260425_0000` | `NCAL_Calibration_Source_..._20260427_0000` |
|---|---|---|
| Detector hits analysed | 10,814,070 | 10,503,722 |
| Best-fit period | **38.851 ns** (25.74 MHz) | 48.357 ns (noise) |
| Modulation amplitude | **58.0 %** | 1.1 % |
| Significance over null | **50 σ** | 2.5 σ (scan artifact) |
| Fold depth (peak-to-peak) | **127.7 %** | 3.9 % |
| Amplitude at 38.851 ns | 58.0 % | **0.91 % (1.2 σ — not detected)** |
| 95 % upper limit at 38.851 ns | — | **1.05 %** |

Both files used the RF-derived `Prescaler1024` trigger and were processed with
the same detector-hit selection. The calibration file shows no statistically
significant modulation at the beam period: its fitted amplitude is 0.91 %, with
a 1.05 % upper limit, compared with 58 % in the selected beam sample.

The beam period of 38.851 ns (25.74 MHz) matches the cyclotron RF, and the
folded profile is a clean single sinusoid — the signature of beam bunching.

See `step04_rf_ch0free/zebra_fold_detectors.png`.

---

## 3. Two problems that had to be solved first

### 3.1 No installed software could decode the file

The decoded `.dat` for this run was 0 bytes and the old log ended with
`Cannot use TDC5 data for non TDC5 processor`. Root cause:

- The `stream` library gained a TDC5-only DOGMA path on 2025-09-01. Every
  currently installed build (`dabcgo4`, `dabcgo4_bak20260923`) contains it and
  therefore **refuses to read the Prague data**, which is not TDC5. Forcing
  `CreateTDC5` produced empty subevents and nonsense timestamps — proof the data
  is genuinely not TDC5.
- The one older build, `dabcgo4_bak20260330`, predates the 2025-08-01
  introduction of the per-subevent trigger timestamp, so it decoded hits
  correctly but reported `GlobalTriggerTime 0.000000000`.

**Fix:** rebuild `stream` at commit `6a3b6e1` (2025-09-01, the parent of the
TDC5 rewrite) against the current ROOT/Go4. This is the only revision with both
the non-TDC5 DOGMA path *and* the trigger timestamp.

**Validation:** re-decoding the April beam file reproduced the April reference
`.dat` **byte-for-byte over 174,598 lines**, including every
`GlobalTriggerTime` value. The environment is saved at
`/home/padiwa/compat_bin/apr2026_env.sh`.

### 3.2 The control run has no channel-0 RF reference

The standard `step04` folds hits against the channel-0 RF reference pulse. In
this control run that channel was **not cabled**:

| | Beam run | Control run |
|---|---|---|
| ch0 hits on TDC_001049 | 686,562 (≈100 % occupancy) | **0** |
| ch0 hits on TDC_001048 | 2,906 | 5,203, all inside a 12 s blip of a 416 s run |

So `step02`/`step03` ch0-referenced outputs and `step04` **cannot run as
designed** on this file. This is a cabling fact about the run, not a software
problem.

**Fix:** the DOGMA trigger *is* the cyclotron RF divided by 1024, so every
trigger edge is already RF-phase-locked, and hit times in the `.dat` are
measured relative to it. Folding trigger-referenced hit times therefore tests
for the RF modulation **without needing ch0 at all**. This is implemented as a
new parallel tool, `Analysis/step04_alt_trigger_referenced_rf_modulation.py`.
The existing step code was not modified.

---

## 4. Method (and why it is trustworthy)

1. Histogram trigger-referenced rising-edge times in 0.5 ns bins.
2. **Exclude channel 0.** Including it would inject the RF line by
   construction, since ch0 *is* the RF pulse. This matters: the first pass
   included ch0 and had to be redone.
3. Divide out a 400 ns running mean. The trigger-window acceptance varies on a
   microsecond scale, so this removes the envelope while leaving the ~39 ns
   modulation untouched.
4. Restrict both datasets to the **same** +200…+5000 ns window. This avoids the
   control's large prompt spike at −525 ns (57 % of its hits) and the beam's ch0
   region, and leaves nearly equal statistics (10.8 M vs 10.5 M hits).
5. Scan the trial period; take the first-harmonic amplitude of the residual.
6. Calibrate against the same statistic at trial periods far from the peak. The
   amplitude is positive-definite, so the excess is measured **against the null
   mean, not against zero**.

Two earlier approaches were tried and discarded as invalid, and are documented
here so they are not repeated:

- A **Rayleigh test on raw hit times** gave mean Z ≈ 5×10⁵. The Rayleigh test
  assumes a uniform null, but the hit-time distribution is strongly peaked, so
  it reports enormous "significance" at *every* trial period. Useless.
- A **raw FFT power spectrum** was contaminated by TDC fine-time quantisation,
  which produces comb structure at all frequencies and inflates the noise floor
  by a factor of ~10⁶ above Poisson.

---

## 5. What this does and does not prove

**Shows:** the selected April 25 beam sample has a strong 38.851 ns modulation,
while the April 27 radioactive-source calibration sample does not. This argues
against an RF-trigger/DAQ effect that would necessarily appear in both runs.

**Does not prove:** that every electronics effect or run-condition difference
has been eliminated. The files are different types of acquisition taken on
different days, not a perfectly matched beam-on/beam-off pair. The 1.05 % limit
applies specifically to modulation in this calibration sample under the stated
selection.

**Relation to the Zebra pulser tests:** their lack of stripes is consistent
with this result, but the calibration comparison does not by itself prove why
the Zebra setup lacked modulation.

---

## 6. Run data

- **Beam input:** `/data6/Data/NCAL_Test_Prescaler1024_TrigWindow_0x83e803e8_20us_Pos_5m_20260425_0000.dld.dat`; the histogram used the first 60,000,000 decoded text lines.
- **Calibration input:** `/data6/Data/NCAL_Calibration_Source_Prescaler1024_TrigWindow_0x83e803e8_20us_20260427_0000.dld` (1,912,909,824 bytes) — verified unchanged before and after; its full decoded `.dld.dat` file was histogrammed.
- **Common selection:** board 1 rising edges, channel 0 excluded, and hit times restricted to +200…+5000 ns, leaving 10,814,070 beam hits and 10,503,722 calibration hits.
- **Decode:** 20,868,779 events, 2,079,859,992 bytes processed, 6.1 GB `.dat`, 416.2 s of run time.
- **Structural validation:** passed. Boards paired per event with a constant 13.140 µs epoch offset; trigger timestamps strictly monotonic on both boards.
- **Board occupancy:** TDC_001049 carries the detectors (73.8 M edges); TDC_001048 is nearly empty (1.2 M edges).
- **Calibration files:** unchanged by the run (`config_drift.txt` is empty).

## 7. Local layout

```
step00/              decode logs, validation, channel census, macros, scripts
step04_rf_ch0free/   histograms, period scans, folded profiles, plots, scripts
  verify/            reproduction of the result with the packaged repo tool
```

Every copied file is hash-verified against `SOURCE_SHA256SUMS.txt`.

## 8. Parameter handling

Step 02/03/04 defaults were returned to their Prague values, with the Zebra
two-generator values preserved as inline comments next to each one (for example
`double ch0RefMinNs = -6000.0;  // Zebra: -11000.0`), so both parameter sets
remain available. `Analysis/README.md` documents this.

Note that `step02`, `step03` (ch0-referenced parts) and `step04` were **not**
run on this file, because it has no ch0 reference channel. `step01` cleaning and
the trigger-referenced parts of `step03` would run, but add nothing to the
question above, which the ch0-free test answers directly and with better
statistics.

## 9. Regression-test status

The full `Analysis/tests/` suite was re-run on `padiwa` after the parameter
restoration: **30 tests, all passing**.

Getting there surfaced one issue worth recording. `/data6/.../Analysis` is not a
git checkout, it is kept in sync by hand, and five modules
(`step05`, `step06`, `step07`, `step09`, `step10`) had drifted to stale copies.
That drift, not any code defect, caused the failures and errors first seen in
`test_rf_channel_pipeline.py` — for example the stale `step05` was missing the
guard that raises `No folded phase-ToT entries for channel N`. The stale files
were backed up to `Analysis_bak_20261013_stale_step05_10/` and replaced from the
local repository.

After the sync, only three tests still failed, and all three legitimately
asserted the Zebra defaults. They were updated to assert the restored Prague
defaults, with the Zebra parameter set retained as explicit opt-in coverage:

- `test_python_prague_defaults_and_zebra_segment_override`
- `test_default_saves_score_channel_only_and_zebra_list_is_opt_in`
- `test_cpp_zebra_parameter_set_covers_broad_scan_without_losing_output_counts`
- `test_cpp_prague_defaults_score_only_on_prague_shaped_pulses` (new, guards the
  restored C++ defaults against Prague-shaped mock pulses)

**Before running anything on `padiwa`, verify the remote `Analysis` tree matches
the local repository.**
