import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ANALYSIS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS))
from step05_estimate_gamma_neutron_background import read_channel_histogram


class MockChainIntegrationTest(unittest.TestCase):
    def test_numbered_chain_with_known_rf_and_detector_counts(self):
        # Synthetic software fixture only: not a replacement for the Zebra input.
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "synthetic.dld.dat"
            run = "Synthetic_20us_Pos_5m_0000"
            windows = 300
            rng = np.random.default_rng(20261008)
            with raw.open("w") as handle:
                for window in range(windows):
                    handle.write("TDC 0 of total 1 TDCs sub 0x1 finding sub-events 0x1 TDC TDC_MOCK size 28\n")
                    handle.write(f"GlobalTriggerTime {100 + window * .001:.9f}\n")
                    hits = [(0, -401 + cycle * 40, 17) for cycle in range(4)]
                    for cycle in range(3):
                        hits.extend([
                            (2, -391 + cycle * 40 + rng.normal(0, .8), rng.normal(13.6, .35)),
                            (2, -381 + cycle * 40 + rng.normal(0, 1), rng.normal(26, 1)),
                            (2, -401 + cycle * 40 + rng.uniform(0, 40), rng.uniform(45, 55)),
                        ])
                    hits.append((6, -380, 16))
                    for channel, rise, tot in hits:
                        handle.write(f"1 {channel} 1 {rise:.6f}\n1 {channel} 0 {rise + tot:.6f}\n")

            def execute(script, *args):
                result = subprocess.run(
                    [sys.executable, str(ANALYSIS / script), *map(str, args)],
                    text=True, capture_output=True,
                )
                self.assertEqual(result.returncode, 0, f"{script}\n{result.stdout}\n{result.stderr}")

            results = root / "Results"
            refined = root / "Refined"
            build = root / "Build"
            execute("step01_run_dogma_raw_data_refiner.py",
                    "--input", raw, "--output-root", refined, "--run-key", run,
                    "--signal-tdc-ordinal", 0, "--trigger-tdc-ordinal", 0,
                    "--binary-path", build / "cleaner", "--skip-root-tree")
            pulses = refined / f"{run}_rawRefined" / f"{run}_rawRefined_pulses.tsv"
            common = ["--input", pulses, "--output-root", results, "--run-key", run, "--build-dir", build]
            execute("step02_run_cleaned_all_channel_rates.py", *common)
            execute("step03_run_cleaned_all_channel_time_tot.py", *common)
            execute("step04_run_cleaned_all_channel_rf_period.py", *common,
                    "--score-channel", 0, "--histogram-channels", "0,2,6",
                    "--initial-period-min-ns", 39.5, "--initial-period-max-ns", 40.5,
                    "--initial-step-ns", .1, "--refine-rounds", 1,
                    "--ch0-valid-rise-min-ns", -410, "--ch0-valid-rise-max-ns", -395,
                    "--score-stride", 1,
                    "--min-selected-pulses", 100, "--single-rf-segment",
                    "--skip-cycle-residual-diagnostics")
            histogram = read_channel_histogram(
                results / run / "Folded_RF_3x" / f"{run}_Folded_RF_3x_phase_tot_sparse.tsv", 2
            )
            self.assertEqual(histogram.counts.sum(), windows * 9)
            self.assertAlmostEqual(float(histogram.metadata["deduced_period_ns"]), 40, delta=.15)

            root_file = results / run / f"{run}_cleaned_hits.root"
            execute("step01_write_cleaned_root_tree.py", "--pulse-table", pulses,
                    "--output-root-file", root_file,
                    "--summary-file", results / run / "cleaned_hits_summary.txt")
            import uproot
            with uproot.open(root_file) as handle:
                self.assertEqual(handle["CleanedHits"].num_entries, windows * 14)

            downstream = ["--results-root", results, "--run-glob", run, "--channel", 2]
            for step, script in (
                (5, "step05_estimate_gamma_neutron_background.py"),
                (6, "step06_quantify_gamma_neutron_components.py"),
                (7, "step07_phase_rate_slice_diagnostics.py"),
                (8, "step08_fit_gamma_tot_peak.py"),
                (9, "step09_gamma_window_sensitivity.py"),
            ):
                args = [*downstream, "--output-dir", root / f"Step{step:02d}"]
                if step == 9:
                    args.extend(["--gamma-fit-summary", root / "Step08/gamma_tot_peak_fit_summary.tsv"])
                execute(script, *args)
                self.assertTrue(list((root / f"Step{step:02d}").rglob("*.png")), script)
            execute("step10_gamma_reference_cfd.py", "--results-root", results,
                    "--input-dir", root / "Step09", "--output-dir", root / "Step10", "--channel", 2)
            with (root / "Step10/gamma_reference_cfd50_summary.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 3)
            for row in rows:
                self.assertTrue(np.isfinite(float(row["cut_cfd50_ref_ns"])), row)
            print("Verified synthetic steps01-10: RF=40ns, ch2=2700 folded entries, "
                  "ROOT=4200 hits, three finite CFD references. Temporary outputs removed.")


if __name__ == "__main__":
    unittest.main()
