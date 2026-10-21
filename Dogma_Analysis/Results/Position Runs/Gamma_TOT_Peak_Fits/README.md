# Step08 Sharp Gamma ToT Peak Fits

This directory contains Gaussian fits to the sharp Ch02 gamma ToT peak.

Method:
- use the Ch02 whole-run ToT distribution already produced by step03;
- convert counts to rate per ToT bin by dividing by `covered_duration_seconds`;
- fit only the narrow ROI [12.5, 15.0] ns;
- model = constant background + Gaussian peak.

Interpretation:
- `peak_mu_ns` is the fitted ToT peak position.
- `peak_mu_err_ns` is the statistical fit uncertainty on that mean under the chosen model and fit ROI.
- it is not the full uncertainty budget; ROI choice, model mismatch, binning, and residual time-walk effects are systematic terms.
