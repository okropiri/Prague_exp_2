# Step10 Gamma Reference CFD Comparison

This directory compares prompt-gamma reference definitions for the processed runs in the supplied gamma-window scan input.

References compared:
- broad-band maximum,
- broad-band left-edge constant-fraction crossing at 0.50 of the net peak height,
- centered-cut maximum,
- centered-cut local quadratic peak fit around the maximum,
- centered-cut left-edge constant-fraction crossing at 0.50 of the net peak height.

Interpretation:
- the orange dash-dot line is the selected local-fit prompt-gamma reference,
- the cyan line is the actual raw cut-band maximum from the unsmoothed cut-band profile,
- the green line is the CFD-style left-edge reference on the cut-band profile,
- the purple line is the CFD-style left-edge reference on the broad-band profile,
- the blue dotted line is the broad-band maximum reference.
- `cut_local_peak_shift_from_max_ns` is the signed local-fit minus raw-cut-maximum shift.
- The right-hand 1D profiles are displayed as mean `counts / accepted RF cycle / ns`.
- The left-hand 2D panels are built from the canonical one-period fold and tiled 3x only for visibility.
- The TSV summary retains `Hz/ns` reference-shape columns and also records the accepted-window validation metadata carried from step09.
