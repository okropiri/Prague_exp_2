# Archived DOGMA Analysis Code

## Legacy Raw-Analysis Line

- `legacy_dogma_absolute_rate_profile.cpp`: replaced by the cleaned-data rate chain in step 02.
- `legacy_write_absolute_rate_outputs.py`: writer for the legacy absolute-rate raw-analysis outputs.
- `legacy_dogma_ch0ref_full_analysis.cpp`: replaced by the cleaned-data split into step 02 rates, step 03 time-vs-ToT, and step 04 RF analysis.
- `legacy_write_ch0ref_full_analysis_outputs.py`: writer for the legacy ch0-referenced raw-analysis outputs.
- `legacy_dogma_ch0ref_ncal_rf_scan.cpp`: archived raw NCAL-only RF-period scan producer; the active production chain now uses cleaned pulses in step 04.
- `legacy_run_abs_and_ch0ref_rates.py`: legacy combined raw-analysis runner.
- `run_abs_and_ch0ref_rates.py`: superseded pre-step combined runner that still points at retired raw-analysis sources and writers, so it is kept only for reference.

## Mock And Local Test Helpers

- `mock_generate_cleaned_pulse_table.py`: deterministic mock cleaned-pulse generator for end-to-end chain validation.
- `mock_run_step01_cleaner_stress_tests.py`: synthetic raw DOGMA stress-test harness for the step-01 cleaner, with plain-text QC output.
- `mock_write_analysis_chain_qc_report.py`: markdown report writer for the mock cleaned-data chain validation outputs.
- `test_root_writer.py` and `test_root.root`: local uproot TH2 write/read sanity check and its test artifact.

## Older Exploratory Utilities

- `dogma_lstilbene_edge_checks.cpp` and `write_lstilbene_edge_checks_outputs.py`: archived exploratory Lstilbene edge studies.
- `dogma_trgref_artifact_hunting.cpp` and `write_trgref_artifact_hunting_outputs.py`: archived trigger-reference artifact-hunting workflow.
- `dogma_trgref_asymmetry_investigation.cpp` and `write_trgref_asymmetry_investigation_outputs.py`: archived trigger-reference asymmetry study.
- `digest_dld_data.py` and `trigger_hist_summary.py`: archived CSV-digest support utilities, not part of the current cleaned-data production chain.