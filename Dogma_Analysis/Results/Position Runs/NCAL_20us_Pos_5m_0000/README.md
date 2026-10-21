# Canonical 5 m: matched RF-method analysis

## Start here

- [Side-by-side 3P phase versus ToT](Comparison/NCAL_3P_ToT_comparison.png)
- [Side-by-side 1P phase versus ToT](Comparison/NCAL_1P_ToT_comparison.png)
- [Equal-bin 3P profile overlay](Comparison/NCAL_3P_equal_bin_profiles.png)
- [Ch2 period-search procedure](Comparison/Ch2_period_search_procedure.png)
- [Frozen run settings and executed source identities](Run_settings.json)
- [Equal-bin scientific verification](Comparison/equal_bin_verification.json)
- [Local copy verification](LOCAL_VERIFICATION.json)

## What is controlled

Both methods use exactly the same newly cleaned, single-copy recording on
TDC_001049 and the same guarded first passing ch0 reference in file order.
Steps01-03 are shared. The two Step04 methods differ in period evidence only:
ch0 burst timing versus ch2 detector timing. Both then use the same unshifted
folding implementation. No fitted detector offset or regional alignment is
applied to final RF phases. Step08's unfolded ToT fit is also shared.

All primary ToT displays use **0-60 ns**. Stored maps retain **0-128 ns**,
and the full NCAL pulse trees retain the complete cleaned ToT range.
Background/fit selections are not silently narrowed to the display.

| Accounting | Count |
|---|---:|
| All cleaned NCAL pulses in the shared pulse trees | 31,217,712 |
| Phase-valid NCAL pulses, identical for both methods | 31,159,659 |
| Missing-reference NCAL pulses retained phase-invalid | 58,053 |
| Phase-valid NCAL pulses at/above 60 ns, hidden only by the ToT zoom | 77 |
| Ch0-derived period (nominal DAQ ns) | 38.850071306015 |
| Ch2-derived period (nominal DAQ ns) | 38.849500000000 |

Ch2 minus ch0 is -0.000571306 ns, or approximately -14.705 ppm
(a magnitude of 0.571306 ps). These reported digits describe the computation,
not absolute hardware timing accuracy or an estimator-disagreement significance.
The selected detector-period optimum is interior to every search round;
see the [scan assessment](Comparison/period_scan_assessment.json).
Similar visible structures remain in both matched maps. Their agreement
does not by itself prove the physical origin of those structures.

## Numbered steps and reports

Steps01-03 contain shared cleaning, rate and timing/ToT results.
The two labelled Step04 directories contain calibration/folding results.
Steps05-07 and09-10 contain `Ch0-based` and `Ch2-based` outputs.
Step08 is the shared unfolded ToT fit.
Every step retains a `Reports` area with resolved parameters (including
defaults), exact commands, logs, exit status, validation and output hashes.
Reports reference immutable remote processing paths; do not execute their
commands locally on the Mac. Bulk decoded and all-channel cleaned inputs
remain on padiwa. Step00 was hash-verified and reused, not decoded again.

Some downstream plots intentionally use a fitted gamma-relative coordinate.
They are labelled diagnostics and do not replace the unshifted RF products.
The 10000 ns acquisition-window width used for accepted-cycle normalization
is an explicit assumption, not a measurement from plotting axes. No 20%
duty-factor correction is applied. A 5m-only analysis cannot validate a
cross-position neutron-width model or absent background-control run.

## ROOT reproduction

The complete NCAL pulse set is stored **once** in
[the ch0 ROOT directory](<04. Ch0-based RF folding/ROOT Trees/>).
It includes invalid-phase rows. Its stored `rf_phase_ns`/`rf_phase3_ns`
branches use the ch0-derived period. For either method, reconstruct phase
from `rise_time_ch0_ref_ns` using that method's period in the table above.
The following is an event-level exploratory view, not exact histogram replay:

```cpp
TChain pulses("Pulses");
pulses.Add("04. Ch0-based RF folding/ROOT Trees/NCAL_ch02_pulses_part*.root");
// Choose the ch0 or ch2 measured period above; no peak-centering subtraction.
double T = 38.849500000000013;
pulses.Draw(Form("tot_ns:((fmod(fmod(rise_time_ch0_ref_ns/%.17g,3)+3,3))*%.17g)",T,T),
            "phase_valid", "COLZ");
```

Each Step04 branch has exact all-channel weighted-bin replay **TTrees**.
The comparison ROOT file has equal-0.25ns-bin NCAL replay TTrees for both
periods, with 0.1ns ToT bins over the full stored 0-128ns range.
These are weighted-bin trees, not unbinned events. Use `count` as a weight
and the stored axis metadata when reproducing maps. Set the display y range
to 0-60ns without deleting the upper bins.
For exact comparison-map replay in ROOT (the common 3P axis has
467 bins):

```cpp
TFile replay("Comparison/Equal_bin_NCAL_replay.root");
auto bins = replay.Get<TTree>("ch2_phase3_tot");
auto map = new TH2D("ch2_replay", "Ch2-derived period;Unshifted phase (ns);ToT (ns)",
                   467,
                   0, 116.75,
                   1280, 0, 128);
bins->Draw("tot_ns:phase_ns>>ch2_replay", "count", "goff");
map->GetYaxis()->SetRangeUser(0, 60);
map->Draw("COLZ");
gPad->SetLogz();
```

Use `ch0_phase3_tot` with the same axes for the ch0 map. The one-period
trees are named `ch0_phase1_tot` and `ch2_phase1_tot`; their own stored
axis metadata gives the shorter x range. Count weights are integers,
and every persisted comparison bin was read back and checked.
The two methods use the same hits, so their histograms are statistically
correlated; no independent-sample chi-square significance is claimed.

## Preservation and copy provenance

Previous canonical contents are retained unchanged under
[Archive](Archive/2026-10-19_Previous_5m_results/).
The [remote archival verification](Archive/remote_archival_verification.json)
and [per-file original hashes](Archive/remote_original_hashes.json) additionally
record preservation of all 4,048 original remote files. Historical remote
paths remain compatibility symlinks to their archived originals.
Older separately named diagnostic runs remain historical; this folder is
the canonical matched comparison.
All 1,341 copied source artifacts were SHA-256 checked.
The [copy manifest](COPY_MANIFEST.json) maps every local artifact to its
remote source. No pilots, build trees, bulk inputs or source snapshots are
copied here. Complete execution products and frozen sources remain remote.
