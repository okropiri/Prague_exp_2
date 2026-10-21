from decimal import Decimal
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ANALYSIS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS))
from run_trigger_timestamp_envelope import burst_groups, fit_repetition, measure


class TriggerEnvelopeTests(unittest.TestCase):
    def make_ticks(self, period_ns=6666665, bursts=300, count=33):
        starts = np.arange(bursts, dtype=np.int64) * period_ns
        return ((starts[:, None] + np.arange(count) * 39780) // 5).ravel().astype(np.uint64)

    def test_grouping_and_occupancy_do_not_discard_any_trigger(self):
        ticks = self.make_ticks()
        result, starts, stops = measure(ticks, 500000)
        self.assertEqual(result["trigger_windows"], 9900)
        self.assertEqual(result["bursts"], 300)
        self.assertEqual(int(np.sum(stops - starts)), len(ticks))
        self.assertAlmostEqual(result["repetition_fit"]["fitted_period_ns"], 6666665, places=6)
        self.assertAlmostEqual(result["regular_burst_span_plus_one_trigger_interval_duty_proxy"]["median"],
                               33 * 39780 / 6666665, places=10)
        for threshold in (200000, 1000000):
            np.testing.assert_array_equal(burst_groups(ticks, threshold)[0], starts)

    def test_frequency_is_measured_not_forced_to_150_hz(self):
        for period in (8000000, 5555555):
            result = measure(self.make_ticks(period_ns=period), 500000)[0]
            self.assertAlmostEqual(result["repetition_fit"]["frequency_hz"], 1e9 / period, places=8)

    def test_long_pause_is_reported_and_not_bridged_by_fit(self):
        ticks = self.make_ticks()
        ticks[4950:] += np.uint64(200000000)
        result = measure(ticks, 500000)[0]
        self.assertEqual(result["gaps_above_100ms"], 1)
        self.assertEqual(result["repetition_fit"]["fit_regions"], 2)
        self.assertEqual(result["repetition_fit"]["excluded_start_intervals"], 1)
        self.assertAlmostEqual(result["repetition_fit"]["fitted_period_ns"], 6666665, places=6)

    def test_missing_recorded_trigger_does_not_remove_other_windows(self):
        ticks = np.delete(self.make_ticks(), 10)
        result = measure(ticks, 500000)[0]
        self.assertEqual(result["trigger_windows"], 9899)
        self.assertEqual(result["bursts"], 300)
        self.assertEqual(result["dense_burst_trigger_count_inventory"]["32"], 1)

    def test_no_burst_structure_fails_instead_of_claiming_modulation(self):
        ticks = np.arange(300, dtype=np.uint64) * 7956
        with self.assertRaisesRegex(RuntimeError, "Insufficient burst starts"):
            measure(ticks, 500000)

    def test_burst_fit_uses_no_cross_pause_cycle_count(self):
        starts = np.r_[np.arange(100) * 7000000., 900123456 + np.arange(100) * 7000000.]
        fitted = fit_repetition(starts)
        self.assertEqual(fitted["fit_regions"], 2)
        self.assertAlmostEqual(fitted["fitted_period_ns"], 7000000, places=6)

    def test_sparse_between_burst_triggers_remain_accounted(self):
        ticks = self.make_ticks()
        sparse = ((np.arange(299) * 6666665 + 3500000) // 5).astype(np.uint64)
        ticks = np.sort(np.r_[ticks, sparse])
        result, starts, stops = measure(ticks, 200000)
        self.assertEqual(result["bursts"], 300)
        self.assertEqual(result["sparse_group_trigger_windows"], 299)
        self.assertEqual(result["dense_group_trigger_windows"], 9900)
        self.assertTrue(result["all_trigger_accounting_verified"])
        self.assertEqual(int(np.sum(stops - starts)) + 299, len(ticks))

    def test_full_cli_retains_last_trigger_in_segment_accounting(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            ticks = self.make_ticks()
            timestamps = folder / "timestamps"
            timestamps.mkdir()
            ticks.astype("<u8").tofile(timestamps / "trigger_ticks_5ns.bin")
            (timestamps / "extraction_summary.txt").write_text(
                f"status=PASS\ntrigger_timestamps={len(ticks)}\nfirst_tick={ticks[0]}\nlast_tick={ticks[-1]}\n")
            with (timestamps / "timestamp_samples.tsv").open("w") as handle:
                handle.write("window_index\tdecoded_line\toriginal_seconds\tcoarse_tick_5ns\n")
                for index in (0, len(ticks) - 1):
                    seconds = Decimal(int(ticks[index])) * Decimal("0.000000005")
                    handle.write(f"{index + 1}\t1\t{seconds:.9f}\t{ticks[index]}\n")
            result = subprocess.run([sys.executable, str(ANALYSIS / "run_trigger_timestamp_envelope.py"),
                                     "--timestamps", str(timestamps), "--output", str(folder / "output")],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            import json
            report = json.loads((folder / "output/measurement.json").read_text())
            self.assertEqual(sum(row["trigger_windows"] for row in report["equal_duration_segments"]), len(ticks))
            self.assertEqual(report["absolute_histogram_integral"], len(ticks))
            self.assertEqual(len(list((folder / "output").glob("*.png"))), 3)


class NativeTriggerExtractionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.workspace.name) / "extract"
        subprocess.run(["g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-pedantic",
                        str(ANALYSIS / "extract_native_trigger_timestamps.cpp"),
                        "-o", str(cls.binary)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.workspace.cleanup()

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)

    def run_extractor(self, stamps, expected=None, hardware="TDC_001049"):
        raw = self.folder / "decoded.dat"
        lines = []
        for stamp in stamps:
            lines.extend(["TDC 0 of total 2 TDCs TDC TDC_001048 size 0",
                          "GlobalTriggerTime 1.000000000",
                          f"TDC 1 of total 2 TDCs TDC {hardware} size 0",
                          f"GlobalTriggerTime {stamp}"])
        raw.write_text("\n".join(lines) + "\n")
        return subprocess.run([str(self.binary), str(raw), str(self.folder / "output"),
                               str(len(stamps) if expected is None else expected)],
                              capture_output=True, text=True)

    def test_empty_headers_and_noninteger_burst_gaps_are_retained(self):
        stamps = ["902.178410820", "902.178450600", "902.183850735"]
        result = self.run_extractor(stamps)
        self.assertEqual(result.returncode, 0, result.stderr)
        expected = [int(Decimal(value) * Decimal(1000000000) / 5) for value in stamps]
        np.testing.assert_array_equal(np.fromfile(self.folder / "output/trigger_ticks_5ns.bin",
                                                 dtype="<u8"), expected)

    def test_duplicate_or_backward_timestamps_are_rejected(self):
        result = self.run_extractor(["1.000000000", "1.000000000"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Duplicate/backwards", result.stderr)

    def test_wrong_board_is_rejected(self):
        result = self.run_extractor(["1.000000000", "1.000040000"], hardware="TDC_001048")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not TDC_001049", result.stderr)

    def test_truncated_precision_is_rejected(self):
        result = self.run_extractor(["1.000000", "1.000040"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("nine timestamp decimals", result.stderr)

    def test_count_mismatch_is_rejected(self):
        result = self.run_extractor(["1.000000000", "1.000040000"], expected=3)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("count differs", result.stderr)


if __name__ == "__main__":
    unittest.main()
