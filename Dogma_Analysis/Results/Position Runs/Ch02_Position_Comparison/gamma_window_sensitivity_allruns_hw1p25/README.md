# Step09 Gamma Window Sensitivity

This directory tests the 5 m prompt-gamma reference against centered fixed-width gamma ToT windows.

Method:
- use the fitted gamma ToT peak position `mu_gamma` from step08;
- compare the original broad gamma band [10.0, 18.0] ns against centered cut gates `[mu_gamma - half_width, mu_gamma + half_width]`;
- scan half widths: 1.25 ns;
- rebuild the prompt-gamma phase reference and the 1D phase-rate slice diagnostics for each width;
- use a local quadratic peak fit around the cut-band profile maximum as the final phase anchor.

Interpretation:
- the spread in `phase_center_ns` across the width scan is a direct measure of the selected local-fit gamma-reference sensitivity to the ToT gate choice;
- the `*_gamma_reference.png` plots keep the selected gamma anchor line on the 1D panel and use the broader ROI overlays only on the 2D panel;
- the right-hand 1D profiles are displayed as mean `counts / accepted RF cycle / ns`;
- an additional high-ToT band [8.0, inf] ns is shown in green to track the candidate neutron population in the same accepted-cycle normalization;
- the integration ROIs are shown on the 2D panel, and the third panel reports their integrated counts per accepted RF cycle;
- the left-hand 2D panels are built from the canonical one-period fold and tiled 3x only for visibility;
- the stored TSV outputs retain `Hz/ns` rate-density columns for wall-clock normalization and add the accepted-cycle-normalized profile columns and validation-efficiency metadata;
- accepted-cycle normalization uses the explicit working-model trigger-window width assumption `assumed_trigger_window_width_ns=10000`;
- this is a systematic-selection study, not a statistical fit error.
