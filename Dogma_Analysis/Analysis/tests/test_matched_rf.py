from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import uproot

ANALYSIS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS))
from analysis_provenance import file_identity, record_arguments
from step04_write_cleaned_all_channel_rf_period_outputs import read_matrix, write_sparse_2d_outputs
from step03_write_cleaned_all_channel_time_tot_outputs import write_mode_outputs


class MatchedRFTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        for name in ("zebra_ch0_clock", "step01_dogma_raw_data_refiner"):
            subprocess.run(["g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-pedantic",
                            str(ANALYSIS / f"{name}.cpp"), "-o", str(cls.root / name)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_equal_periods_give_identical_unshifted_bins_and_retain_all_ch0(self):
        source = self.root / "matched.tsv"
        source.write_text(
            "window_index\ttdc_ordinal\tglobal_trigger_seconds\tchannel\t"
            "rise_tdc_index_one_based\tfall_tdc_index_one_based\trise_ns\tfall_ns\ttot_ns\n"
            "1\t1\t100\t0\t2\t2\t-400\t-383\t17\n"
            "1\t1\t100\t0\t2\t2\t-395\t-378\t17\n"
            "1\t1\t100\t2\t2\t2\t-300\t-285\t15\n"
            "1\t1\t100\t2\t2\t2\t-410\t-340\t70\n"
            "2\t1\t101\t2\t2\t2\t-300\t-285\t15\n")
        periods = self.root / "matched_periods.tsv"
        periods.write_text("1\t40\n")
        for method, model in (("ch0", "shared_burst_rf_period"), ("ch2", "shared_ch2_detector_period")):
            command = [str(self.root / "zebra_ch0_clock"), "fold", "--input", str(source),
                       "--periods", str(periods), "--output", str(self.root / method),
                       "--run-key", "matched", "--ordinal", "1", "--allow-reference-gaps", "1",
                       "--expected-windows", "2", "--include-all-ch0-hits", "1",
                       "--period-model", model, "--histogram-channels", "0,2",
                       "--ch0-valid-rise-min-ns", "-410", "--ch0-valid-rise-max-ns", "-395",
                       "--ch0ref-min-ns", "-200", "--ch0ref-max-ns", "200", "--ch0ref-bins", "80",
                       "--trigger-ref-min-ns", "-600", "--trigger-ref-max-ns", "200", "--trigger-ref-bins", "160",
                       "--tot-max-ns", "128", "--tot-bins", "1280"]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        for name in ("folded_matrix", "folded3_matrix", "phase_tot", "phase3_tot",
                     "phase_time", "phase3_time", "phase_trigger", "phase3_trigger"):
            a, b = [(self.root / method / f"{name}.tsv").read_text().splitlines()
                    for method in ("ch0", "ch2")]
            self.assertEqual([row for row in a if not row.startswith("#")],
                             [row for row in b if not row.startswith("#")])
        _, centers, counts, _ = read_matrix(self.root / "ch0/folded_matrix.tsv")
        self.assertEqual(counts[:, 0].sum(), 2)
        self.assertEqual(counts[:, 2].sum(), 2)
        self.assertEqual(centers[np.flatnonzero(counts[:, 2])[0]], 20.125)

    def test_duplicate_recording_restart_is_rejected_not_stitched(self):
        source = self.root / "repeated.dld.dat"
        source.write_text("".join(
            "TDC 0 of total 1 TDCs sub 0x1 finding sub-events 0x1 TDC TDC_MOCK size 28\n"
            f"GlobalTriggerTime {time:.9f}\n1 0 1 -400\n1 0 0 -383\n"
            for time in (100, 102, 100, 102)))
        result = subprocess.run([str(self.root / "step01_dogma_raw_data_refiner"),
            "--input", str(source), "--output-root", str(self.root / "rejected"),
            "--run-key", "duplicate", "--signal-tdc-ordinal", "0", "--trigger-tdc-ordinal", "0",
            "--reject-trigger-restarts"], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("possible duplicated recording", result.stderr)

    def test_resolved_settings_include_defaults_and_nonfinite_values_are_explicit(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.dict(os.environ, {"DOGMA_STEP_REPORT_DIR": directory}), patch.object(
                    sys, "argv", [str(ANALYSIS / "run_matched_five_m.py")]):
                record_arguments(argparse.Namespace(display=60, enabled=False, upper=float("inf")))
            files = list(Path(directory).glob("*_settings.json"))
            self.assertEqual(len(files), 1)
            settings = json.loads(files[0].read_text())
            self.assertEqual(settings["effective_arguments"], {"display": 60, "enabled": False, "upper": "inf"})
            self.assertEqual(settings["source"]["sha256"], file_identity(ANALYSIS / "run_matched_five_m.py")["sha256"])

    def test_both_tot_map_zooms_preserve_high_tot_root_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            source = folder / "sparse.tsv"
            source.write_text("# run_key=zoom\n# x_min_ns=0\n# x_max_ns=40\n# x_bins=4\n"
                              "# y_min_ns=0\n# y_max_ns=128\n# y_bins=1280\n"
                              "# tot_min_ns=0\n# tot_max_ns=128\n# tot_bins=1280\n"
                              "# channel_count=32\n# stored_channels=2\n2 1 150 10\n2 2 700 3\n")
            saved = []
            def check(figure, *args, **kwargs):
                for axis in figure.axes:
                    if axis.get_ylabel() == "ToT (ns)":
                        self.assertEqual(tuple(axis.get_ylim()), (0, 60))
                        saved.append(True)
            with patch("matplotlib.figure.Figure.savefig", new=check):
                write_sparse_2d_outputs(source, folder, combined_filename_suffix="map",
                    per_channel_filename_suffix="map", root_name_suffix="rf", hist_group="RF",
                    x_label="phase", y_label="ToT (ns)", title_prefix="RF", display_y_max_ns=60)
                write_mode_outputs(source, folder, mode_name="time", filename_suffix="time",
                    hist_group="Time", x_label="time", title_prefix="time", display_max_ns=60)
            self.assertTrue(saved)
            for filename, group in (("zoom_rf.root", "RF"), ("zoom_time.root", "Time")):
                with uproot.open(folder / filename) as root:
                    hist = root[f"{group}/hCh02"]
                    self.assertEqual(hist.values().sum(), 13)
                    self.assertEqual(hist.axis(1).edges()[-1], 128)
                    self.assertEqual(hist.values()[2, 700], 3)


if __name__ == "__main__":
    unittest.main()
