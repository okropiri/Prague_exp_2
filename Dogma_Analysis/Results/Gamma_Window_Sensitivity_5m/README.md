# Step09 Gamma Window Sensitivity

This directory tests the 5 m prompt-gamma reference against centered fixed-width gamma ToT windows.

Method:
- use the fitted gamma ToT peak position `mu_gamma` from step08;
- define gamma ToT gates as `[mu_gamma - half_width, mu_gamma + half_width]`;
- scan half widths: 1 ns, 1.25 ns, 1.5 ns;
- rebuild the prompt-gamma phase reference and the 1D phase-rate slice diagnostics for each width.

Interpretation:
- the spread in `phase_center_ns` across the width scan is a direct measure of gamma-reference sensitivity to the ToT gate choice;
- this is a systematic-selection study, not a statistical fit error.
