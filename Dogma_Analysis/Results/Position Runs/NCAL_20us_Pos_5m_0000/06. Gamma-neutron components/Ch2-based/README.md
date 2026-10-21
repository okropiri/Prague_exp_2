# Step06 Gamma/Neutron Component Quantification

This directory contains the first quantitative gamma/background/neutron-candidate protocol built from cleaned `Folded_RF_3x` phase-vs-ToT sparse histograms.

Channel: Ch02

Purpose:
- detect the compact prompt-gamma centroid and width independently for each run;
- draw the nominal neutron-candidate phase expected from the prompt gamma plus the gamma-neutron TOF delay for a 30 MeV reference neutron;
- identify a compact observed fast-neutron core near that TOF guide when it is sufficiently separated from the prompt gamma;
- define background sidebands by excluding prompt-gamma and neutron-candidate phase regions;
- report raw counts, background estimates, net counts, and rates for prompt gamma and for a constrained fast-neutron TOF/ToT slice.

Important interpretation notes:
- The 3.4 m runs are treated as background-control candidates because the neutron-candidate band is expected to overlap or hide near the prompt-gamma region, not because neutrons are absent.
- The broad neutron-candidate distribution is not fully contained by the direct-counting ROI and may be wider than one RF period after folding. The direct neutron number reported here is only a fast-neutron TOF slice, not the total neutron yield.
- The orange dotted line is the nominal 30 MeV neutron TOF guide. The orange box, when present, is centered on the observed fast-neutron core near that guide.
- The orange box is constrained both in TOF/RF phase and in ToT. Its ToT interval is an energy-deposition proxy and must be applied consistently to simulation or control-detector comparisons before converting rates into flux or efficiency.
- Neutron-candidate width is allowed to vary with detector distance. Current width model source: `constant_median_separated_runs`.
- Overlay plots use gamma-relative RF phase by default, so the prompt gamma ROI is centered rather than split at an RF-period boundary.
- Overlay plots use linear Z by default. The color scale is clipped at quantile 0.995 to make the gamma and neutron-candidate regions auditable without logarithmic visual inflation of sparse tails.
- Broad neutron widths are used only as exclusion bands for background protection. Direct neutron counts use a compact observed fast-neutron phase-ToT ROI when separation is sufficient.
- Runs flagged as `gamma_neutron_roi_overlap` have no direct neutron rate reported; the expected neutron phase is drawn only as a dashed guide line.
- Runs flagged as `close_to_gamma_tail` report a conservative neutron-core rate, but they should still be treated as template-fit candidates rather than final yields.

Default regions:
- Gamma ToT band: [10.0, 18.0) ns.
- Signal/neutron-candidate ToT band: [8.0, 30.0) ns.
- Fast-neutron TOF window: 5.0 ns full width.
- Neutron core ToT interval: measured from the valid neutron-core phase ROI using quantiles [0.16, 0.84], with minimum width 4.0 ns.
- Background-control positions: 3.4 m.
- Neutron-width reference positions: 5 m, 5.8 m, 6.6 m.
