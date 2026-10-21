from __future__ import annotations

from decimal import Decimal
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ANALYSIS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS))
from run_zebra_ch0_clock import calibrate_clock, tick_numbers
from step02_run_cleaned_all_channel_rates import compile_cpp
from step04_write_cleaned_all_channel_rf_period_outputs import read_matrix


class ClockFitTests(unittest.TestCase):
    def test_off_grid_period_without_detector_feature(self) -> None:
        n = np.arange(4000)
        times = n * 1024 * 39.125371 + np.random.default_rng(72).normal(0, .08, n.size)
        model, assessment = calibrate_clock(times, 1024, 11000, .1)
        self.assertLess(abs(model.periods.mean() - 39.125371), 1e-7)
        self.assertLess(assessment["max_3sigma_conditional_phase_error_ns"], .1)

    def test_quadratic_clock_drift(self) -> None:
        n = np.arange(4000)
        times = 1024 * (39.12 * n + .5 * 8e-7 * n * n)
        times += np.random.default_rng(7).normal(0, .08, n.size)
        model, _ = calibrate_clock(times, 1024, 11000, .1)
        midpoint = (np.maximum(0, n - 4) + np.minimum(n.size - 1, n + 4)) / 2
        self.assertLess(float(np.max(abs(model.periods - (39.12 + 8e-7 * midpoint)))), 1e-4)

    def test_smooth_modulation_and_phase_precision_limit_are_explicit(self) -> None:
        n = np.arange(4000)
        times = 1024 * (39.12 * n - .5 * np.cos(n / 100))
        model, assessment = calibrate_clock(times, 1024, 11000, .1)
        phase_error = 11000 * abs(model.periods - (39.12 + .005 * np.sin(n / 100))) / 39.12
        self.assertLess(float(phase_error.max()), .1)
        self.assertEqual(assessment["half_span"], 4)
        self.assertIn("NOT_CERTIFIED", assessment["phase_precision_status"])
        self.assertEqual(len(assessment["interval_span_sensitivity"]), 3)

    def test_missing_prescaled_ticks_are_counted(self) -> None:
        full = np.arange(4000)
        kept = np.delete(full, np.r_[100:110, 300])
        times = kept * 1024 * 39.125371
        ticks = tick_numbers(times, 1024)
        np.testing.assert_array_equal(ticks, kept)
        _, assessment = calibrate_clock(times, 1024, 11000, .1)
        self.assertEqual(assessment["missed_prescaled_ticks"], 11)

    def test_invalid_clock_is_explicitly_rejected(self) -> None:
        for times, prescaler in ((np.arange(64) * 40000, 512),
                                 (np.r_[np.arange(63), 61] * 40000, 1024),
                                 (np.arange(64) * np.nan, 1024)):
            with self.assertRaises(RuntimeError):
                calibrate_clock(times, prescaler, 11000, .1)


class BackendTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.workspace = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.workspace.name) / "zebra_ch0_clock"
        compile_cpp(ANALYSIS / "zebra_ch0_clock.cpp", cls.binary, "g++")

    @classmethod
    def tearDownClass(cls) -> None:
        cls.workspace.cleanup()

    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.folder = Path(self.temp.name)
        self.cleaned = self.folder / "cleaned.tsv"
        self.decoded = self.folder / "decoded.dat"
        self.periods = self.folder / "periods.tsv"
        self.quantum = Decimal(40) / 3
        self.true_period = 39.125371
        self.rises = []
        with self.cleaned.open("w") as clean, self.decoded.open("w") as raw, self.periods.open("w") as periods:
            clean.write("# run_key=synthetic_rawRefined\n")
            for i in range(128):
                time = i * 1024 * self.true_period
                tick = round((time + 383) / float(self.quantum))
                rise = time - tick * float(self.quantum)
                self.rises.append(rise)
                stamp = (Decimal(65000000000000 + tick) * self.quantum / Decimal(10**9))
                raw.write(f"TDC 0 of total 2 TDCs TDC AD size 0\nGlobalTriggerTime {stamp:.9f}\n")
                raw.write(f"TDC 1 of total 2 TDCs TDC AE size 6\nGlobalTriggerTime {stamp:.9f}\n")
                for ch, edge, value in ((0, 1, rise), (0, 0, rise + 17),
                                        (0, 1, rise + 42.9), (0, 0, rise + 95.9),
                                        (2, 1, rise - 3000), (2, 0, rise - 2983)):
                    raw.write(f"2 {ch} {edge} {value:.3f}\n")
                for ch, r, f, tot in ((0, rise, rise + 17, 17),
                                     (0, rise + 42.9, rise + 95.9, 53),
                                     (2, rise - 3000, rise - 2983, 17)):
                    clean.write(f"{i + 1}\t1\t{stamp:.6f}\t{ch}\t2\t2\t{r:.3f}\t{f:.3f}\t{tot}\n")
                periods.write(f"{i + 1}\t{self.true_period}\n")

    def tearDown(self) -> None:
        self.temp.cleanup()

    def extract(self, output: Path, quantum: str = "13.333333333333333333") -> subprocess.CompletedProcess[str]:
        return subprocess.run([str(self.binary), "extract", "--input", str(self.cleaned),
                               "--decoded", str(self.decoded), "--ordinal", "1",
                               "--coarse-quantum-ns", quantum, "--output", str(output)],
                              capture_output=True, text=True)

    def test_fine_ch0_reconstruction_undoes_anchor_rounding(self) -> None:
        output = self.folder / "references.tsv"
        result = self.extract(output)
        self.assertEqual(result.returncode, 0, result.stderr)
        data = np.loadtxt(output, comments="#")
        expected = np.arange(128) * 1024 * self.true_period
        self.assertLess(float(np.max(abs(data[:, 4] - expected))), .002)
        np.testing.assert_array_equal(data[:, 5], np.ones(128))

    def test_wrong_clock_lattice_is_rejected(self) -> None:
        result = self.extract(self.folder / "bad.tsv", "14")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("coarse lattice", result.stderr)

    def test_low_precision_trigger_metadata_is_rejected(self) -> None:
        text = self.decoded.read_text()
        lines = [f"GlobalTriggerTime {Decimal(line.split()[1]):.6f}"
                 if line.startswith("GlobalTriggerTime ") else line for line in text.splitlines()]
        self.decoded.write_text("\n".join(lines) + "\n")
        result = self.extract(self.folder / "bad.tsv")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("9-decimal", result.stderr)

    def test_fold_counts_and_actual_modulo_three_periods(self) -> None:
        output = self.folder / "folded"
        result = subprocess.run([str(self.binary), "fold", "--input", str(self.cleaned),
                                 "--periods", str(self.periods), "--output", str(output),
                                 "--run-key", "synthetic"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        for factor, filename in ((1, "folded_matrix.tsv"), (3, "folded3_matrix.tsv")):
            metadata, _, counts, _ = read_matrix(output / filename)
            self.assertEqual(counts[:, 0].sum(), 128)
            self.assertEqual(counts[:, 2].sum(), 128)
            self.assertEqual(counts[0, 0], 128)
            self.assertAlmostEqual(float(metadata["x_max_ns"]), factor * self.true_period, places=8)
            expected_phase = ((-3000 / self.true_period) % factor) * self.true_period
            expected_bin = int(expected_phase / float(metadata["x_max_ns"]) * counts.shape[0])
            self.assertEqual(counts[expected_bin, 2], 128)
        second = subprocess.run([str(self.binary), "fold", "--input", str(self.cleaned),
                                 "--periods", str(self.periods), "--output", str(output),
                                 "--run-key", "synthetic"], capture_output=True, text=True)
        self.assertNotEqual(second.returncode, 0)
        self.assertIn("replace existing", second.stderr)

    def test_variable_period_normalization_and_long_tot_accounting(self) -> None:
        local_periods = self.true_period + np.linspace(-.01, .01, 128)
        self.periods.write_text("".join(f"{i + 1}\t{value:.12f}\n" for i, value in enumerate(local_periods)))
        lines = self.cleaned.read_text().splitlines()
        fields = lines[1].split("\t")
        for channel, tot in ((6, 130), (13, 81)):
            lines.insert(4, f"1\t1\t{fields[2]}\t{channel}\t2\t2\t-1000\t{-1000 + tot}\t{tot}")
        self.cleaned.write_text("\n".join(lines) + "\n")
        output = self.folder / "variable"
        result = subprocess.run([str(self.binary), "fold", "--input", str(self.cleaned),
                                 "--periods", str(self.periods), "--output", str(output),
                                 "--run-key", "synthetic"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        for factor, filename in ((1, "folded_matrix.tsv"), (3, "folded3_matrix.tsv")):
            metadata, _, counts, _ = read_matrix(output / filename)
            expected_phase = ((-3000 / local_periods) % factor) * local_periods.mean()
            expected_bins = (expected_phase / float(metadata["x_max_ns"]) * counts.shape[0]).astype(int)
            np.testing.assert_array_equal(counts[:, 2], np.bincount(expected_bins, minlength=counts.shape[0]))
            self.assertEqual(counts[:, 6].sum(), 1)
            self.assertEqual(counts[:, 13].sum(), 1)
        stats = np.loadtxt(output / "channel_counts.tsv", skiprows=1)
        self.assertEqual(stats[6, 3], 1)
        self.assertEqual(stats[13, 3], 1)


if __name__ == "__main__":
    unittest.main()
