# Step09 Efficiency Finding

This directory captures the current efficiency-finding view built on the step09 gamma-reference work.

Method:
- use the fitted gamma ToT peak position `mu_gamma` from step08;
- compare the original broad gamma band [10.0, 18.0] ns against centered cut gates `[mu_gamma - half_width, mu_gamma + half_width]`;
- scan half widths: 1 ns;
- rebuild the prompt-gamma phase reference and the 1D phase-rate slice diagnostics for each width;
- use a local quadratic peak fit around the cut-band profile maximum as the final phase anchor.
- keep the narrow centered cut only for anchor finding;
- define both efficiency ROIs with the shared ToT cut [8.0, inf] ns and separate them only by the dip-bounded phase regions.

Interpretation:
- the spread in `phase_center_ns` across the width scan is a direct measure of the selected local-fit gamma-reference sensitivity to the ToT gate choice;
- the `*_gamma_reference.png` plots keep the selected gamma anchor line on the 1D panel and show the current efficiency ROIs only on the 2D panel and the integral summary panel;
- the right-hand 1D profiles are displayed as mean `counts / accepted RF cycle / ns`;
- an additional high-ToT band [8.0, inf] ns is shown in green to track the candidate neutron population in the same accepted-cycle normalization;
- the integration ROIs are shown on the 2D panel, and the third panel reports their integrated counts per accepted RF cycle;
- the left-hand 2D panels are built from the canonical one-period fold and tiled 3x only for visibility;
- the stored TSV outputs retain `Hz/ns` rate-density columns for wall-clock normalization and add the accepted-cycle-normalized profile columns and validation-efficiency metadata;
- accepted-cycle normalization uses the explicit working-model trigger-window width assumption `assumed_trigger_window_width_ns=10000`;
- when `--apply-zebra-timewalk-correction` is enabled, the correction slope is derived from the phase-peak shift between the `2-3 ns` and `6-7 ns` low-ToT zebra bands before rebuilding the same efficiency-style plots;
- this is a systematic-selection study, not a statistical fit error.
