# Repository Instructions

## Workspace layout
- `Dogma_test_Data/` contains raw input data and should be treated as read-only.
- `Analysis/` is where active analysis scripts, utilities, and supporting code should be kept.
- `Analysis/Archive/` is where older or replaced analysis code should be moved, both locally and in mirrored remote analysis trees.
- `Results/` is reserved for generated outputs, plots, tables, and derived result files.

## Processing order and result layout
- Cleaning and validation always happen first.
- Downstream analysis codes should use cleaned data, not the raw DOGMA export, unless the user explicitly asks to work on raw input.
- Derive one parent run directory from the input file, for example `NCAL_20us_Pos_3.4m_0000`.
- Write all outputs under `Results/<run_dir>/` instead of naming result directories from the full raw input filename.
- Cleaner outputs belong under `Results/<run_dir>/cleaned/`, optionally with one tool-specific child directory if needed.
- Downstream cleaned-data analyses belong under `Results/<run_dir>/cleaned/<analysis_name>/`, for example `Results/NCAL_20us_Pos_3.4m_0000/cleaned/ch0Ref_rfPhaseScan_ncalOnly_37to41/`.
- If the user requests a top-level result layout for a downstream cleaned-data analysis, follow that explicit override instead of forcing the generic `cleaned/<analysis_name>/` pattern.
- For cleaned all-channel rate products, write directly to `Results/<run_dir>/Abs_rates/` and `Results/<run_dir>/Ch0_ref_Rates/`.
- For cleaned all-channel time-vs-ToT products, write directly to `Results/<run_dir>/Trigger_ref_ToT/`, `Results/<run_dir>/Ch0_ref_TOT/`, and `Results/<run_dir>/TOT_distrib/`.
- For cleaned all-channel RF-period products, write directly to `Results/<run_dir>/Ch0_ref_Rates/`, `Results/<run_dir>/Folded_RF/`, `Results/<run_dir>/Folded_RF_3x/`, and `Results/<run_dir>/RF_period_scan/`.
- When step-02 and step-04 share the same `run_key`, keep both under `Ch0_ref_Rates/` but use the step-04 `*_Ch0_ref_Rates_scan*` filename stem so the RF-scan products do not overwrite the step-02 ch0-reference rate products.
- Keep RF-period deduction diagnostics together under `RF_period_scan/`; keep folded phase-profile and folded phase-vs-ToT outputs under `Folded_RF/` and `Folded_RF_3x/`.
- Keep the same relative result layout on remote storage.

## Remote and local storage
- Write DOGMA result sets first to the remote tree under `/data6/Dogma_analysis_by_Dachi/Results/<run_dir>/...`.
- Keep archived result sets under `/data6/Dogma_analysis_by_Dachi/Results/Archive/`.
- Do not copy large result files back to the local machine by default.
- Only mirror small artifacts locally unless explicitly requested: PNG, PDF, small text tables/summaries, and small ROOT files.
- Create or edit code locally first, do cheap local smoke checks there, then sync the active code to `padiwa` and let `padiwa` do the full data processing.
- Do not shift the full cleaned-data processing workload back to the local machine when `padiwa` is the intended execution target.
- When a cleaned-data analysis uses a ch0-referenced frame, apply the guarded ch0-reference selection consistently: use the first valid ch0 hit in file order, not simply the minimum-time ch0 hit.
- When a runner executes a cleaned-data ch0-referenced analysis, save odd-ch0 diagnostic outputs in the run directory so the rejected and multi-hit ch0 populations can be audited later.

## Python environment
- Before running any Python command in the terminal, activate the conda environment with `conda activate conda_env`.
- Do not run `python`, `pip`, notebooks, or Python-based tooling until `conda_env` is active.
- Prefer storing reusable analysis code in `Analysis/` rather than mixing code into the raw data or results directories.
- For this workspace, use the local `conda_env` for code creation, smoke tests, and support tooling before syncing to `padiwa`.
- Current cleaned-data entrypoints are `step02_run_cleaned_all_channel_rates.py` for rate products, `step03_run_cleaned_all_channel_time_tot.py` for time-vs-ToT products, and `step04_run_cleaned_all_channel_rf_period.py` for cleaned RF-period products.

## Data handling
- Do not modify files inside `Dogma_test_Data/` unless explicitly requested.
- Write derived files and exports to `Results/` using the run-key-based layout described above.