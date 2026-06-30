# DOGMA Analysis Scripts

## Working Model

- Edit code locally in the mounted workspace.
- Run full processing on `padiwa` against the real `/data6` filesystem.
- Keep the active chain in step order so the cleaned-data workflow is readable from filenames alone.
- Experiment-specific timing context is summarized in `../Reference/experiment_context.md`.
- A human-oriented explanation of the NCAL RF scan is in `../Reference/ncal_rf_phase_scan_walkthrough.md`.
- Activate the existing Python environment before running any helper scripts:

```bash
conda activate conda_env
```

## Active Cleaned-Data Chain

- `step01_dogma_raw_data_refiner.cpp`: raw DOGMA export to cleaned pulse table. This is the gatekeeper step; downstream analyses are expected to consume its cleaned pulse output.
- `step01_run_dogma_raw_data_refiner.py`: step-01 wrapper that compiles and runs the cleaner and can hand the cleaned pulse table to the ROOT writer.
- `step01_write_cleaned_root_tree.py`: converts the cleaned pulse table into the split `CleanedHits` ROOT export plus `RunMetadata` and the parts manifest.
- `step02_dogma_cleaned_all_channel_rates.cpp`: cleaned-data all-channel rate analyzer.
- `step02_write_cleaned_all_channel_rates_outputs.py`: writes the cleaned absolute-rate and ch0-referenced rate plots plus ROOT outputs.
- `step02_run_cleaned_all_channel_rates.py`: step-02 runner and shared helper module used by later cleaned-data steps.
- `step03_dogma_cleaned_all_channel_time_tot.cpp`: cleaned-data all-channel time-vs-ToT analyzer in trigger and guarded ch0 frames.
- `step03_write_cleaned_all_channel_time_tot_outputs.py`: writes step-03 trigger-ref and ch0-ref time-vs-ToT plots plus ROOT outputs.
- `step03_run_cleaned_all_channel_time_tot.py`: step-03 runner.
- `step04_dogma_cleaned_all_channel_rf_period.cpp`: cleaned-data RF-period deduction and all-channel RF folding analyzer.
- `step04_write_cleaned_all_channel_rf_period_outputs.py`: writes step-04 ch0-ref RF timing, folded 1x profiles, folded 3x profiles, and folded phase-vs-ToT products.
- `step04_write_rf_period_scan_outputs.py`: writes the dedicated RF-period-scan diagnostics in the old-style merit/sigma/residual format from the step-04 scan text outputs.
- `step04_run_cleaned_all_channel_rf_period.py`: step-04 runner.
- `step05_estimate_gamma_neutron_background.py`: first-pass channel-2 phase-vs-ToT sideband estimator for prompt-like ROI, broad neutron-candidate band, and local background control regions.

## Utility Files Kept Nearby

- `README.md`: this file.
- `Prague_Cyclotron_Exp_2.code-workspace`: local VS Code workspace file for opening this analysis folder.

## Archived and Legacy Code

- Replaced raw-analysis code has been moved into `Archive/` with `legacy_` prefixes.
- Old mock-data generators, stress-test helpers, and local ROOT writer sanity checks have also been moved into `Archive/` so the top-level `Analysis/` directory only keeps the active production chain.
- Older exploratory trigger-reference, edge-check, and CSV-digest utilities stay in `Archive/` because they are not part of the current cleaned-data production chain.
- See `Archive/README.md` for the archived-file summary.