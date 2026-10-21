# Gamma/Neutron Phase-ToT Background Estimate

This directory contains a first-pass sideband estimate from existing Folded_RF_3x phase-vs-ToT sparse histograms.

Last refreshed: 2026-07-01 from padiwa `/data6` outputs.

Included runs: 14 completed NCAL position-scan runs, including the recovered `NCAL_20us_Pos_4m_0000` run.

Channel: Ch02

Regions used:
- Signal ToT band: [8.0, 30.0) ns.
- Prompt-like gamma ROI: detected automatically in [10.0, 18.0) ns and repeated over the three RF periods with half-width 2.5 ns.
- Empty-phase sideband: lowest 0.20 quantile of the one-period signal-band phase projection after excluding +/- 7.5 ns around the prompt-like phase.
- High-ToT sideband: [45.0, 80.0) ns over all RF phases.

Interpretation:
- The sideband estimate is a local background/pedestal estimate, not a particle-identification label.
- In overlap regions, gamma and neutron yields should be extracted with template or mixture fits after background treatment.
- The broad neutron-candidate band in the summary table is the signal ToT band after removing the prompt-like gamma ROI; it is not yet a pure neutron sample.
