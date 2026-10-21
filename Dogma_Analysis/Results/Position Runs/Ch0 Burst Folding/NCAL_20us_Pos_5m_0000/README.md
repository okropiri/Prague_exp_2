# Canonical 5 m - full-population ch0 RF folding

## Start here

- [NCAL 3P versus ToT](<RF folding 3P/ch02_rf_phase3_vs_tot.png>)
- [NCAL 1P versus ToT](<RF folding 1P/ch02_rf_phase_vs_tot.png>)
- [All original trigger windows and ROOT replay/pulse trees](<ROOT Trees/>)
- [Independent full-data QC](QC/independent_verification.json)
- [Estimator assessment](QC/assessment.json)

## What changed

One shared RF period, **38.850071306015 ns** in the nominal
DAQ time scale, is estimated from ch0 intervals inside stable burst fragments.
Each fragment has its own intercept. No cycles or carrier phase are inferred
through inactive gaps. The maximum conditional regional phase bound at
6000 ns delay is **0.053006 ns**.
This validates the declared shared-period model, not instantaneous carrier phase.

**Every passing-reference window is folded, including both burst boundaries,
short fragments and isolated references. No symmetric-neighbour mask remains.**
Calibration support is not a detector-quality cut. All ch0 pulses are retained,
not just the reference pulse; the origin remains the guarded first passing ch0
in file order, never the minimum-time ch0 or a detector-derived alignment.

| Accounting | Count |
|---|---:|
| Original board1049 trigger windows in window TTree | 10,516,211 |
| Passing-reference windows | 10,450,559 |
| Clock-support-excluded windows | 0 |
| All cleaned pulses, every channel | 186,894,978 |
| All cleaned NCAL ch2 pulses | 31,217,712 |
| NCAL pulses with justified ch0-referenced phase | 31,159,659 |
| NCAL pulses without passing ch0, retained with invalid phase | 58,053 |

All 10,516,211 windows remain in `TriggerWindows`, including
empty and unreferenced windows. Without a passing ch0 a justified ch0 phase
cannot be manufactured. Those pulse rows have `phase_valid=false` and NaN phase.
ToT display bounds affect displayed bins, not the pulse trees.

The new run reuses the independently verified fresh Steps00-03 inputs and
reruns RF calibration, folding, output writing and independent full-data QC.
Original ch2-estimator results and the earlier support-selected diagnostic
were not overwritten. The legacy ch2 export is duplicated; its raw entry
counts are not directly comparable to a single fresh decode.

## ROOT TTrees: exact reproduction, not mystery ROOT files

The complete NCAL pulse set is copied, including invalid-phase rows; no sampling.

- `Trigger_windows.root`, TTree `TriggerWindows`: every original window,
  native 5 ns timestamp, cleaned-pulse count, passing-reference status and
  direct calibration-support flag.
- `NCAL_ch02_pulses_part*.root`, TTree `Pulses`: original trigger/ch0-relative
  timing, ToT, 1P/3P phases and validity for every cleaned NCAL pulse.
- `All_channel_plot_replay.root`: eight **weighted-bin TTrees**, reproducing
  every 1P/3P channel profile and detailed ToT/time/trigger map. The `count`
  branch is the bin weight, not one event. Axis JSON objects accompany trees.
  These preserve exact plotted bins, not unbinned data or arbitrary finer binning.

For the NCAL 3P/ToT plot in ROOT, using the complete pulse set:

```cpp
TChain pulses("Pulses");
pulses.Add("ROOT Trees/NCAL_ch02_pulses_part*.root");
pulses.Draw("tot_ns:rf_phase3_ns>>h_pulses(468,0,116.550213918054908,1280,0,128)", "phase_valid", "COLZ");
```

Or reproduce the exact binned plot from the small all-channel replay:

```cpp
auto file = TFile::Open("ROOT Trees/All_channel_plot_replay.root");
auto tree = file->Get<TTree>("phase3_tot");
tree->Draw("value_ns:phase_ns>>h_replay(468,0,116.550213918054908,1280,0,128)", "count*(channel==2)", "COLZ");
```

All exported ROOT objects above are real TTrees, not RNTuples or histogram-only
files. Persisted NCAL trees were read back and every NCAL profile/map bin was
checked again; every all-channel replay bin was also checked independently.

## Clean local layout and copy provenance

Only final plots, the TTrees, concise QC and this navigation are copied.
No pilots, rejected trials, build directories, source snapshots, raw tables
or verbose execution logs are brought into this results folder.

ROOT transfer: **921.11 MiB total**, capped at
250 MiB/file and 1000 MiB total. Every copied source artifact is SHA-256 verified
in [COPY_MANIFEST.json](COPY_MANIFEST.json).
Complete inputs and generated products remain on padiwa at:

`/data6/Dogma_analysis_by_Dachi/Results/NCAL_20us_Pos_5m_0000/Ch0_Burst_Folding_20261015/full/Results/NCAL_20us_Pos_5m_0000`
