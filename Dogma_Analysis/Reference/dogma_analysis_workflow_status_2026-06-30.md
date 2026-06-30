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
- Real NCAL runs with final RF-annotated `CleanedHits` summaries: `13`.
- One real run is incomplete at the cleaned production-chain level: `NCAL_20us_Pos_4m_0000` has step-02/step-03 outputs but no step-04 RF-period summary and no final cleaned-hit ROOT summary.

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
| `NCAL_20us_Pos_4m_0001` | 38.85 | 36.753742 | 95861910 | 5 |
| `NCAL_20us_Pos_5m_0000` | 38.8495 | 37.612242 | 375580494 | 19 |
| `NCAL_20us_Pos_5m_0001` | 38.8495 | 37.584363 | 145525832 | 8 |
| `NCAL_20us_Pos_5.8m_0000` | 38.8495 | 1.368496 | 401078638 | 21 |
| `NCAL_20us_Pos_6.6m_0000` | 38.85 | 2.267314 | 255908347 | 13 |

The deduced RF period is stable at about `38.8495-38.85 ns` across completed runs. Phase origin changes with run/position, as expected for the chosen folding origin.

## Incomplete Run

`NCAL_20us_Pos_4m_0000` has:

- refined pulse table: yes
- step-02 cleaned rates: yes
- step-03 time-vs-ToT: yes
- step-04 RF-period products: missing
- final RF-annotated `CleanedHits` ROOT summary: missing

The old overnight batch report says this run failed at `step04_rf_period`. The current results tree still shows no `*_cleaned_rf_period_summary.txt` or `*_cleaned_hits_summary.txt` for this run.

## Interpretation

The project is past raw exploration and into a cleaned-data production workflow. Most real NCAL position/tuning runs have been processed through the full chain, including RF-period deduction and final RF-annotated ROOT export. The immediate workflow gap is to rerun or debug step 04 for `NCAL_20us_Pos_4m_0000`, then write its final `CleanedHits` output if step 04 succeeds.

The local Mac/GitHub mirror is suitable for code/documentation work and small artifacts, but the authoritative large data/output state is still padiwa `/data6`.