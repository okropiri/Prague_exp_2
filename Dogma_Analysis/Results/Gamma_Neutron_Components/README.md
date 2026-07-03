# Step06 Gamma/Neutron Component Quantification

This directory contains the first quantitative gamma/background/neutron-candidate protocol built from cleaned `Folded_RF_3x` phase-vs-ToT sparse histograms.

Channel: Ch02

Purpose:
- detect the compact prompt-gamma centroid and width independently for each run;
- predict the neutron-candidate phase from the prompt gamma plus the gamma-neutron TOF delay for a 30 MeV reference neutron;
- derive a position-dependent neutron-candidate phase width from separated reference positions when possible;
- define background sidebands by excluding both the gamma ROI and the expected neutron-candidate ROI;
- report raw counts, background estimates, net counts, and rates for gamma and neutron-candidate ROIs.

Important interpretation notes:
- The 3.4 m runs are treated as background-control candidates because the neutron-candidate band is expected to overlap or hide near the prompt-gamma region, not because neutrons are absent.
- The neutron-candidate ROI is not a final neutron identification. It is a TOF-guided candidate region for later template or mixture fitting.
- Neutron-candidate width is allowed to vary with detector distance. Current width model source: `constant_median_separated_runs_nonmonotonic_widths`.
- Overlay plots use linear Z by default. The color scale is clipped at quantile 0.995 to make the gamma and neutron-candidate regions auditable without logarithmic visual inflation of sparse tails.
- Broad neutron widths are used only as exclusion bands for background protection. Direct neutron counts use a compact neutron-core phase-ToT ROI when separation is sufficient.
- Runs flagged as `gamma_neutron_roi_overlap` have no direct neutron rate reported; the expected neutron phase is drawn only as a dashed guide line.
- Runs flagged as `close_to_gamma_tail` report a conservative neutron-core rate, but they should still be treated as template-fit candidates rather than final yields.

Default regions:
- Gamma ToT band: [10.0, 18.0) ns.
- Signal/neutron-candidate ToT band: [8.0, 30.0) ns.
- Neutron core ToT interval: measured from the valid neutron-core phase ROI using quantiles [0.16, 0.84], with minimum width 4.0 ns.
- Background-control positions: 3.4 m.
- Neutron-width reference positions: 5 m, 5.8 m, 6.6 m.
