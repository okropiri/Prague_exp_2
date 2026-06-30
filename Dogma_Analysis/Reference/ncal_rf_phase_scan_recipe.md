# NCAL RF Phase Scan Recipe

## Goal

Recover the fine RF period and phase inside the natural `ch0`-referenced `10 us` window using `NCAL` only.

## Inputs

- Raw DOGMA `.dld.dat` file.
- `ch0` as the per-window timing reference.
- `NCAL` rising edges and their rise/fall-paired `ToT` values.

## Procedure

1. Parse the raw DOGMA stream and keep only TDC-1 signal windows.
2. In each signal window, take the earliest `ch0` rise as the local timing reference.
3. Pair `NCAL` rising and falling edges to build `ToT` and compute raw `ch0`-referenced pulse times.
4. Keep the raw `ch0`-referenced times directly rather than starting from an already binned histogram.
5. Scan trial periods over a coarse range around `40 ns`, then refine the scan around the best region.
6. For each trial period:
   - fold raw times modulo the trial period,
   - estimate the dominant phase peak,
   - measure the residual width of that peak,
   - measure the selected fraction of pulses inside the peak window,
   - fit residual vs cycle index as a drift check.
7. Choose the best period from the scan metrics, then build the final `RF phase vs ToT` map.

## Output Logic

- Keep the main result set first on `padiwa:/data6/Dogma_analysis_by_Dachi/Analysis/<dir_name>/`.
- Copy the same result set into the local workspace mirror under `Dogma_Analysis/Results/<dir_name>/`.
- Mark raw-analysis result directories explicitly, for example `_ch0Ref_rawdata_rfPhaseScan`.

## Current Tooling

- `Analysis/Archive/legacy_dogma_ch0ref_ncal_rf_scan.cpp`: archived raw NCAL-only RF scan and text-output producer.
- `Analysis/step04_write_rf_period_scan_outputs.py`: active plotting helper for RF-period scan diagnostics and final phase/ToT figures.