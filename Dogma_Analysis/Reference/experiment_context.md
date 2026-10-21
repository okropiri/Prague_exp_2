# DOGMA Experiment Context

This note captures the working detector and timing assumptions for the NCAL DOGMA analysis in this repository.

## Detector and Timing Chain

- `NCAL` is the neutron detector of interest.
- `ch0` is not a detector pulse from NCAL. It is a timing reference derived from the machine RF chain.
- The original source is a `25 MHz` sinusoidal signal.
- That sinusoid is sent through a quad discriminator and converted into a NIM logic signal.
- The discriminated signal is then prescaled by `2^10` before it is used by the DOGMA DAQ.
- The same prescaled signal is fed both to the trigger path and to channel `0`.

## Why ch0 Referencing Matters

- The absolute trigger timestamp is only precise at about `5 ns`.
- Edge timestamps inside the trigger window are more precise than the trigger absolute timestamp.
- Because channel `0` carries the same prescaled timing signal that participates in the trigger, `ch0` to detector-channel timing should provide a more precise relative time measurement than trigger-absolute timing alone.
- In practice, `ch0` is the preferred timing reference for precise relative timing studies against NCAL, Lstilbene, and Sstilbene.

## Trigger Window Model

- Current working assumption from the DAQ setup is that the trigger window is approximately `+-5000 ns`, so about `10 us` wide in total.
- This should correspond to recording roughly one quarter of the underlying RF time structure because the prescaled trigger arrives at about `25 MHz / 2^10`.
- Numerically, the prescaled timing period is about `40.96 us`, so a `10 us` acquisition window covers about `24%` of that period.

## Physical Timescale

- Even though the DAQ trigger is prescaled, the underlying machine-related physics still happens on the original RF timescale of about `40 ns`.
- A follow-up analysis goal is to determine the exact cyclotron RF frequency more precisely.

## Measured Millisecond Trigger Envelope: Canonical 5 m Run

- **Established in the data:** the canonical `NCAL_20us_Pos_5m_0000` recording
  contains a dominant approximately 150 Hz burst envelope with approximately
  20% recorded-trigger occupancy, in addition to the fast RF timescale.
- All **10,516,211** original board1049 trigger-window timestamps were checked
  on padiwa, before pulse cleaning or ch0-reference selection. Dense bursts
  typically contain 33-35 triggers, spaced by about **39.78 us**.
- Measured burst repetition is approximately **6.665 ms**; the median
  first-to-last trigger span is **1.313 ms**, or **1.353 ms** with a
  one-trigger-spacing endpoint correction. The corresponding median
  trigger-duty proxy is **20.29%**. This is not an exact hardware gate duty.
- The pattern agrees with the supplied cyclotron-control-room scope image
  and persists across ten recording segments and multiple grouping thresholds.
  Sparse inter-burst triggers and one **1.2803 s** exceptional gap also exist;
  the inter-burst interval is not literally free of all recorded triggers.
- **Do not automatically interpret ordinary inter-burst gaps as data loss,
  corrupted timestamps, machine restarts or phase resets. Do not discard
  otherwise valid trigger windows/detector hits because a symmetric clock
  estimate lacks neighbours at a natural burst boundary.**
- Timestamp evidence establishes the recorded envelope, not the exact
  RF/discriminator/trigger-gate mechanism, RF carrier continuity during the
  low-occupancy interval, or beam-current duty. Other runs require their own
  checks; this result does not apply automatically to pulser/Zebra data.
- See [measured trigger-envelope context and analysis guardrails](prague_trigger_envelope.md)
  for methods, measurements, figures and reproducible evidence.
- The subsequent [full-population burst-aware ch0 mode](../Analysis/zebra_ch0_clock.md#full-population-burst-aware-mode)
  separates period-calibration support from data eligibility: no passing-reference
  window is removed for lacking symmetric neighbours. Unreferenced/empty windows
  remain in the all-window inventory; unavailable phases are explicit, never fabricated.
- The [completed corrected 5 m run](<../Results/Position Runs/Ch0 Burst Folding/NCAL_20us_Pos_5m_0000/README.md>)
  retains all 10,516,211 windows and all 31,217,712 cleaned NCAL pulses in local
  TTrees. All 10,450,559 passing-reference windows receive the shared ch0-only
  RF period, approximately 38.850071 ns in the nominal DAQ scale. NCAL folding
  includes 31,159,659 pulses; the remaining 58,053 lack a passing ch0 and remain
  explicitly phase-invalid, not deleted. No clock-support exclusion remains.
- The subsequent [canonical matched RF-method comparison](<../Results/Position Runs/NCAL_20us_Pos_5m_0000/README.md>)
  reran shared cleaning/rate/timing Steps01-03, both Step04 alternatives and
  downstream Steps05-10 with frozen settings. Step00's single-copy decode
  was hash-verified and reused. Both methods use the same unshifted,
  first-passing-ch0 phase origin and the same folding implementation;
  detector-derived global/regional recentering is not applied.
  Independently verified counts remain exactly those above, with no
  clock-support exclusions. Ch0 gives 38.8500713060152 ns and ch2 gives
  38.84950000000001 ns: ch2 minus ch0 is -0.000571306 ns, about -14.705 ppm
  in nominal DAQ units, not an absolute clock calibration or a significance.
- Its [common-bin comparison](<../Results/Position Runs/NCAL_20us_Pos_5m_0000/Comparison/NCAL_3P_ToT_comparison.png>)
  uses 0.25 ns phase bins and 0.1 ns ToT bins with a shared logarithmic
  color scale. ToT display limits are 0-60 ns only; stored maps cover
  0-128 ns and complete NCAL pulse TTrees retain all cleaned widths.
  The 77 phase-valid NCAL pulses at/above 60 ns are retained, not rejected.
  Similar visible structures persist in both methods; this agreement alone
  does not establish their physical origin. Missing background-control
  positions and the assumed 10000 ns normalization window remain explicit
  limitations. The comparison has not yet been added to the LaTeX write-up.

## Current Analysis Question

- In the `ch0`-referenced analysis, the observed folded timing width is not exactly the expected `10 us`.
- This mismatch needs investigation.
- Important possibilities include trigger-window configuration details, edge-selection effects, ch0 selection logic, prescaler behavior, or DAQ/event-building effects.

## Current Measured Findings

- For the run used in the current `ch0`-referenced outputs, the input file is `NCAL_Test_Prescaler1024_TrigWindow_0x83e803e8_20us_Pos_3.4m_20260425_0000.dld.dat`.
- The filename carries a `20us` label, but the direct raw timing structure of the DOGMA file does not support a populated detector window anywhere near `20 us`.
- A direct sample of `200000` signal windows from that raw file gives the earliest `ch0` rise at about `-401.4 ns` relative to trigger zero, with a very narrow spread of about `-406.5 ns` to `-396.4 ns` for the `0.1%` to `99.9%` interval.
- In the same raw-file sample, the occupied `Ncal1` trigger-relative window is about `-4995.7 ns` to `+4927.1 ns` for the `0.1%` to `99.9%` interval.
- In the same raw-file sample, the occupied `Lstilbene` trigger-relative window is about `-4976.2 ns` to `+4935.4 ns` for the `0.1%` to `99.9%` interval.
- After shifting by the measured `ch0` timing, the implied `ch0`-referenced occupied window becomes about `-4594.3 ns` to `+5328.5 ns` for `Ncal1`, and about `-4574.8 ns` to `+5336.8 ns` for `Lstilbene`.
- This matches the existing folded plots and explains why the `ch0`-referenced window looks offset and asymmetric around zero even though the total occupied width is still about `10 us`.
- Current best interpretation: the non-centered `ch0`-referenced width is primarily caused by the fixed offset between `ch0` and trigger zero, not by evidence for a truly wider populated detector window.

## Provenance and Confidence

- The `+-5000 ns` trigger-window statement is currently treated as an operational assumption passed on from the colleague who set up the DOGMA DAQ.
- Until independently verified from DAQ configuration or raw timing structure, treat it as a working model rather than a confirmed hardware specification.
- The currently available raw timing structure supports an occupied detector window of about `10 us`, but the exact hardware or firmware meaning of the filename token `TrigWindow_0x83e803e8_20us` remains unresolved in the repository materials.