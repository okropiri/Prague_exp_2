# DOGMA Analysis Workflow Status - 2026-06-30

This note summarizes the current DOGMA/NCAL analysis workflow after mirroring the Dell workspace to the Mac and checking the active padiwa result tree.

## Operating Model

- Mac workspace: edit code, documentation, LaTeX, and small reference artifacts.
- GitHub mirror: `https://github.com/okropiri/Prague_exp_2.git` tracks code, notes, LaTeX source, figures, setup screenshots, and small selected outputs.
- DOGMA execution target: `padiwa`.
- The Mac workspace is for editing, review, git operations, and documentation; DOGMA Python runners, C++ analyzers, compilers, notebooks, and data-processing commands should not be run locally on the Mac.
- Raw data and full generated outputs live under `/data6`, especially:
  - raw inputs: `/data6/Data/*.dld.dat`
  - refined pulse tables: `/data6/Data/Refined_Data/`
  - analysis outputs: `/data6/Dogma_analysis_by_Dachi/Results/`
- Full DOGMA processing should not be moved back to the Mac unless explicitly requested.

## Production Workflow

The active production path is the cleaned-data chain in `Dogma_Analysis/Analysis/`.

1. `step01_run_dogma_raw_data_refiner.py`
   - Compiles/runs `step01_dogma_raw_data_refiner.cpp`.
   - Converts raw DOGMA `.dld.dat` data into a cleaned pulse table.
   - Output layout on padiwa: `/data6/Data/Refined_Data/<run_key>_rawRefined/<run_key>_rawRefined_pulses.tsv` plus summary.

2. `step02_run_cleaned_all_channel_rates.py`
   - Consumes the cleaned pulse table.
   - Produces absolute-rate and guarded ch0-referenced rate outputs.
   - Output folders: `Abs_rates/` and `Ch0_ref_Rates/`.
   - Writes odd-ch0 reference diagnostics.

3. `step03_run_cleaned_all_channel_time_tot.py`
   - Consumes the cleaned pulse table.
   - Produces trigger-referenced and guarded ch0-referenced time-vs-ToT outputs.
   - Output folders: `Trigger_ref_ToT/`, `Ch0_ref_TOT/`, and `TOT_distrib/`.
   - Reuses the same guarded first-valid ch0 selection policy.

4. `step04_run_cleaned_all_channel_rf_period.py`
   - Consumes the cleaned pulse table.
   - Deduces the RF period, folds all channels in 1x and 3x RF phase, and writes RF-period scan diagnostics.
   - Output folders: `Ch0_ref_Rates/`, `Folded_RF/`, `Folded_RF_3x/`, and `RF_period_scan/`.
   - Current scoring channel is `ch02` by default.

5. `step01_write_cleaned_root_tree.py`
   - Final chain step in `run_all_cleaned_analysis_chain.py`.
   - Writes RF-annotated `CleanedHits` ROOT files plus `RunMetadata` and a parts manifest.
   - Uses the step-04 RF summary so the final ROOT export carries `rf_period_ns`, `rf_phase_origin_ns`, and RF segment metadata.

## Current Processing State On Padiwa

Checked on 2026-06-30 from the Mac with direct `ssh padiwa`.

- `/data6/Data` contains `26` top-level `.dld.dat` files. Several calibration/source files in that list are zero-byte placeholders.
- `/data6/Data/Refined_Data` contains `16` refined directories, including real NCAL runs, mock data, the old first test, and one subset directory.
- `/data6/Dogma_analysis_by_Dachi/Results` contains a complete mock-data chain and the main NCAL run outputs.
- Real NCAL runs with final RF-annotated `CleanedHits` summaries: `14`.
- `NCAL_20us_Pos_4m_0000` was recovered on 2026-06-30 after bounding the memory-heavy step-04 diagnostic path.

## Completed Real Runs

These runs have step-02, step-03, step-04, and final RF-annotated `CleanedHits` output on padiwa:

| Run key | RF period ns | Phase origin ns | CleanedHits entries | ROOT parts |
|---|---:|---:|---:|---:|
| `NCAL_20us_Pos_2.8m_0001` | 38.85 | 26.856041 | 128037544 | 7 |
| `NCAL_20us_Pos_2.8m_0002` | 38.8495 | 26.925074 | 127359298 | 7 |
| `NCAL_20us_Pos_2.8m_0003` | 38.8495 | 26.976751 | 128393920 | 7 |
| `NCAL_20us_Pos_2.8m_0004` | 38.85 | 26.990013 | 74149822 | 4 |
| `NCAL_20us_Pos_2.8m_0005` | 38.85 | 27.151211 | 515990510 | 26 |
| `NCAL_20us_Pos_2.8m_0006` | 38.85 | 27.24949 | 54764560 | 3 |
| `NCAL_20us_Pos_3.4m_0000` | 38.85 | 30.562578 | 453246926 | 23 |
| `NCAL_20us_Pos_3.4m_0001` | 38.85 | 30.600537 | 18145014 | 1 |
| `NCAL_20us_Pos_4m_0000` | 38.8545 | 36.55429 | 466856858 | 24 |
| `NCAL_20us_Pos_4m_0001` | 38.85 | 36.753742 | 95861910 | 5 |
| `NCAL_20us_Pos_5m_0000` | 38.8495 | 37.612242 | 375580494 | 19 |
| `NCAL_20us_Pos_5m_0001` | 38.8495 | 37.584363 | 145525832 | 8 |
| `NCAL_20us_Pos_5.8m_0000` | 38.8495 | 1.368496 | 401078638 | 21 |
| `NCAL_20us_Pos_6.6m_0000` | 38.85 | 2.267314 | 255908347 | 13 |

The deduced RF period is stable at about `38.8495-38.8545 ns` across completed runs. Phase origin changes with run/position, as expected for the chosen folding origin.

## Recovered Run

`NCAL_20us_Pos_4m_0000` now has:

- refined pulse table: yes
- step-02 cleaned rates: yes
- step-03 time-vs-ToT: yes
- step-04 RF-period products: yes
- final RF-annotated `CleanedHits` ROOT summary: yes

The old overnight batch report says this run failed at `step04_rf_period`. A direct full rerun on 2026-06-30 also reached the RF analyzer but was killed after the post-scan diagnostic path grew to nearly all available RAM. The successful recovery used `--score-stride 4`, `--single-rf-segment`, and `--skip-cycle-residual-diagnostics` for step 04, then wrote final `CleanedHits` with the recovered RF summary.

### Step-04 Diagnostic On 2026-06-30

Diagnostic commands were run only on `padiwa`.

- Direct `ssh padiwa` works from the Mac; `/data6/Dogma_analysis_by_Dachi/Results` is accessible.
- No DOGMA processing jobs were running during the check.
- `padiwa` had ample space for the result target: `/data6` had about `9.5T` free; `/tmp` had about `16G` free.
- The failed `NCAL_20us_Pos_4m_0000` result directory contained step-02/step-03 products only: `Abs_rates/`, `Ch0_ref_Rates/`, `Trigger_ref_ToT/`, `Ch0_ref_TOT/`, and `TOT_distrib/`.
- There were no step-04 leftovers for this run: no `Folded_RF/`, no `Folded_RF_3x/`, no `RF_period_scan/`, and no `*_cleaned_rf_period_summary.txt`.
- A padiwa-only subset diagnostic was made from the first `20000` cleaned pulse-table windows under `/tmp` and processed with the current `step04_run_cleaned_all_channel_rf_period.py`.
- The subset step-04 diagnostic succeeded, including the C++ analyzer, RF-period scan, plotting, ROOT output, and odd-ch0 diagnostics.
- The subset diagnostic found `deduced_period_ns=38.8465`, `selected_fraction=0.595852`, `sigma_ns=4.371178`, and `score_pulses_stored=23479` for that small subset.
- Temporary diagnostic files under `/tmp/dogma_diag_4m0000_397860` were removed after the check.
- The active padiwa `step04_dogma_cleaned_all_channel_rf_period.cpp` has timestamp `2026-06-15 01:05`, after the incomplete `4m_0000` outputs from `2026-06-12`.

Interpretation: the current padiwa step-04 code path can process real `NCAL_20us_Pos_4m_0000` data on a subset. The most likely next action is a padiwa-only full step-04 rerun for `NCAL_20us_Pos_4m_0000` using the current code, followed by final `CleanedHits` writing if step 04 succeeds.

### Recovery On 2026-06-30

- Full step-04 recovery succeeded on `padiwa` using a bounded diagnostic path:
  - `--score-stride 4`
  - `--single-rf-segment`
  - `--skip-cycle-residual-diagnostics`
- Step-04 RF summary:
  - `deduced_period_ns=38.8545`
  - `phase_origin_ns=36.55429`
  - `score_pulses_before_stride=77767392`
  - `score_pulses_stored=19441848`
  - `selected_pulses=11522836`
  - `selected_fraction=0.592682`
  - `sigma_ns=4.390381`
  - `rf_segment_count=1`
- The step-04 C++ analyzer completed in about `45 min` with maximum resident set size about `528 MB`.
- The normal padiwa system Python lacked `matplotlib`/`uproot`; the step-04 writer scripts were run with `/home/padiwa/dogsoft/.venv/bin/python`.
- Final `CleanedHits` export succeeded:
  - `root_tree_entries=466856858`
  - `root_tree_file_count=24`
  - `root_tree_rf_period_ns=38.8545`
  - `root_tree_rf_phase_origin_ns=36.55429`
  - `root_tree_rf_segment_count=1`
  - manifest: `/data6/Dogma_analysis_by_Dachi/Results/NCAL_20us_Pos_4m_0000/NCAL_20us_Pos_4m_0000_cleaned_hits_parts.txt`
  - runtime: about `1 h 13 min`

## Interpretation

The project is past raw exploration and into a cleaned-data production workflow. All listed real NCAL position/tuning runs have now been processed through the full chain, including RF-period deduction and final RF-annotated ROOT export. The main technical lesson from the `4m_0000` recovery is that long runs can make optional RF residual/segmentation diagnostics much larger than the physics outputs; use the bounded diagnostic flags when recovering similarly pathological runs.

The local Mac/GitHub mirror is suitable for code/documentation work and small artifacts, but the authoritative large data/output state is still padiwa `/data6`.
