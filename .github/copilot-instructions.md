# Repository Instructions

## Project Context

- This repository is the Mac-side mirror of the Prague Cyclotron Experiment 2 DOGMA/NCAL workspace formerly organized on Dell.
- The Dell source path was `/home/dachi/Documents/NCAL_FSD/Prague_Cyclotron_Exp_2/`.
- The heavy compute/data target is `padiwa` with result storage under `/data6/Dogma_analysis_by_Dachi/Results/`.
- Do not move full DOGMA processing back to the Mac unless explicitly requested.

## Layout

- `Dogma_Analysis/Analysis/` contains active analysis scripts, C++ analyzers, and Python runners.
- `Dogma_Analysis/Archive/` and `Dogma_Analysis/Analysis/Archive/` contain older or replaced analysis code.
- `Dogma_Analysis/Reference/experiment_context.md` is the first timing/detector context note to read.
- `Dogma_Analysis/Reference/ncal_rf_phase_scan_walkthrough.md` explains the RF phase scan at a human level.
- `Documents/Latex/` contains the write-up.
- `data6/` is a local mount/reference path and must not be committed.

## Analysis Rules

- Treat `Dogma_Analysis/Dogma_test_Data/` as read-only raw/test input.
- Cleaning and validation happen before downstream analysis.
- Downstream analyses should consume cleaned data, not raw DOGMA export, unless explicitly requested.
- Keep result layout run-key based: `Results/<run_dir>/...`.
- Preserve the remote result layout on `/data6`.
- For ch0-referenced analysis, use the guarded first-valid ch0 hit in file order rather than simply the minimum-time ch0 hit.

## Execution Model

- Edit locally, perform cheap local checks where possible, then sync/run heavy jobs on `padiwa`.
- On Linux environments used for this analysis, activate `conda_env` before Python commands.
- Current cleaned-data entrypoints are:
  - `step01_run_dogma_raw_data_refiner.py`
  - `step02_run_cleaned_all_channel_rates.py`
  - `step03_run_cleaned_all_channel_time_tot.py`
  - `step04_run_cleaned_all_channel_rf_period.py`
  - `run_all_cleaned_analysis_chain.py`

## Git/Data Hygiene

- Do not commit `.conda/`, `data6/`, Python caches, LaTeX build products, or generated bulk result trees.
- GitHub should synchronize code, notes, LaTeX sources, reference docs, and small selected artifacts.
- Full raw data and full generated result sets stay on `padiwa`/`/data6` unless the user explicitly asks for a different storage plan.