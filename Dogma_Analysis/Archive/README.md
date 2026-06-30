# Archived DOGMA Scripts

These files are kept for reference and comparison, but they are not the main active analysis line anymore.

## Archived Analysis Programs

- `dogma_channel_scan.cpp`: early scanner for rise/fall pairing anomalies, sentinel rows, and TOT statistics on selected channels.
- `dogma_trigger_time_phase_tot_map.cpp`: early trigger-time phase versus TOT map for a single signal channel.
- `dogma_phase_tot_map.cpp`: parent-ch0 matched phase versus TOT map with modulo phase folding.
- `dogma_ch2_start_time_tot_map.cpp`: ch2 start-time versus TOT grid relative to the trigger block start.
- `dogma_ch2_parent_time_tot_no_modulo.cpp`: no-modulo parent-time versus TOT grid for ch2, with delta profiling and period search.
- `dogma_trgref_colleague_style.cpp`: older trigger-referenced histogramming variant kept as a reference implementation.

## Archived Helper and Plot Scripts

- `write_phase_tot_histograms.py`: helper for turning the older phase/TOT text histograms into ROOT and image outputs.
- `write_time_tot_grid_root.py`: helper for converting older time-versus-TOT grid text files into ROOT TH2 histograms.
- `write_trgref_colleague_outputs.py`: plotting and ROOT conversion helper for the colleague-style trigger-referenced outputs.
- `plot_ch2_start_time_tot.gnuplot`: quick gnuplot view for the archived ch2 start-time versus TOT grid.
- `plot_ch2_parent_time_tot_no_modulo.gnuplot`: quick gnuplot view for the archived no-modulo parent-time versus TOT grid.

## Note

If one of these archived scripts becomes active again, move it back into `Analysis/` together with its helper script and update the README files.