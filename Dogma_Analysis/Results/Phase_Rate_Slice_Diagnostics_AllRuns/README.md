# Step07 Phase-Rate Slice Diagnostics

This directory contains gamma-relative one-period rate-vs-phase diagnostics by ToT slice.

Goal:
- normalize each phase profile to rate per unit time;
- inspect abrupt excess regions in 1D without committing yet to a full 2D ROI;
- derive mathematically auditable candidate intervals from excess-over-baseline significance.

Method:
- Prompt gamma is centered at phase 0 in the displayed one-period coordinate.
- For each ToT slice, a folded one-period profile is converted to rate (Hz per phase bin).
- A robust baseline is estimated from the low quantile of bins outside the prompt-gamma exclusion region.
- Excess significance is approximated by z = (counts - baseline_counts) / sqrt(max(baseline_counts, 1)).
- Candidate intervals are contiguous bins with z >= 5 for at least 3 bins.

Important note:
- This is a selection-design diagnostic for the 5 m separated case. It is more rigorous than choosing boxes by eye, but it still needs 2D confirmation before defining final efficiency cuts.

Gamma reference definition:
- The prompt-gamma reference phase is derived from the compact gamma ToT band [10.0, 18.0) ns.
- Counts in that ToT band are folded to one RF period.
- A small circular smoothing is applied.
- The selected gamma reference is the phase bin at the maximum of the smoothed one-period gamma-band profile.
- All phase-rate slice profiles are then displayed relative to that selected prompt-gamma phase.
