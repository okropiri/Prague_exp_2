from __future__ import annotations

import csv
from decimal import Decimal
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import uproot

ANALYSIS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS))
from run_zebra_ch0_clock import calibrate_burst_clock, calibrate_cadence_clock, calibrate_clock, interval_periods, tick_numbers, verify_folded
from step02_run_cleaned_all_channel_rates import compile_cpp, parse_metadata
from step04_write_cleaned_all_channel_rf_period_outputs import read_matrix


class BurstClockTests(unittest.TestCase):
    def fit(self, times, ids=None):
        if ids is None:
            ids = np.arange(len(times)) + 1
        return calibrate_burst_clock(times, 1024, 6000, .1, window_ids=ids)

    def test_all_33_references_including_both_boundaries_are_calibrated(self):
        model, result = self.fit(np.arange(33) * 40960.)
        np.testing.assert_allclose(model.periods, 40, rtol=0, atol=1e-12)
        self.assertEqual(result["calibrated_reference_windows"], 33)
        self.assertEqual(result["unresolved_reference_windows"], 0)

    def test_arbitrary_inactive_gap_cannot_change_measured_period(self):
        original = np.arange(34) * 40960.
        a, _ = self.fit(np.r_[original, original + 5300000.])
        b, result = self.fit(np.r_[original, original + 1234567890.123])
        np.testing.assert_allclose(a.periods, b.periods, rtol=0, atol=1e-10)
        self.assertEqual(result["unbridged_clock_gaps"], 1)
        self.assertIn("unknown", result["clock_tick_scope"])

    def test_short_and_isolated_regions_keep_shared_period_without_fake_local_fit(self):
        times = np.r_[np.arange(34) * 40960., 5300000.,
                      7000000. + np.arange(4) * 40960.,
                      10000000. + np.arange(34) * 40960.]
        model, result = self.fit(times)
        self.assertTrue(np.isfinite(model.periods).all())
        self.assertEqual(result["fit_reference_windows"], 68)
        self.assertEqual(result["calibrated_reference_windows"], len(times))

    def test_missing_reference_never_removes_other_window_periods(self):
        ids = np.arange(100) + 1
        ids[50:] += 1
        model, _ = self.fit(np.arange(100) * 40960., ids)
        np.testing.assert_allclose(model.periods, 40, atol=1e-12)

    def test_real_frequency_change_rejects_shared_model_not_selected_windows(self):
        times = np.r_[np.arange(80) * 1024 * 40.,
                      5300000 + np.arange(80) * 1024 * 40.005]
        with self.assertRaisesRegex(RuntimeError, "stability budget"):
            self.fit(times)

    def test_short_only_and_reordered_timelines_fail_explicitly(self):
        with self.assertRaisesRegex(RuntimeError, "Insufficient"):
            self.fit(np.r_[np.arange(20) * 40960., 5300000 + np.arange(20) * 40960.])
        with self.assertRaisesRegex(RuntimeError, "increase"):
            self.fit(np.r_[np.arange(33), 30] * 40960.)


class GapAwareClockTests(unittest.TestCase):
    def fit(self, times, ids=None):
        from run_zebra_ch0_clock import calibrate_gap_clock
        if ids is None:
            ids = np.arange(len(times)) + 1
        return calibrate_gap_clock(times, 1024, 6000, .1, 5, window_ids=ids)

    def test_noninteger_acquisition_gap_is_never_bridged(self):
        times = np.r_[np.arange(100) * 40960., 5300000 + np.arange(100) * 40960.]
        model, result = self.fit(times)
        expected = np.r_[np.arange(8, 92), np.arange(108, 192)]
        np.testing.assert_array_equal(np.flatnonzero(np.isfinite(model.periods)), expected)
        np.testing.assert_allclose(model.periods[expected], 40, rtol=0, atol=1e-12)
        self.assertEqual(result["unbridged_clock_gaps"], 1)
        self.assertEqual(result["unresolved_reference_windows"], 32)
        self.assertEqual(result["endpoint_period_ns"], 40)
        self.assertIn("unknown", result["clock_tick_scope"])

    def test_missing_reference_id_splits_an_otherwise_regular_timeline(self):
        times = np.arange(200) * 40960.
        ids = np.arange(200) + 1
        ids[100:] += 1
        model, result = self.fit(times, ids)
        self.assertTrue(np.isnan(model.periods[92:108]).all())
        self.assertEqual(result["unbridged_clock_gaps"], 1)
        self.assertEqual(result["clock_support_reason_counts"]["resolved"], 168)

    def test_phase_disagreement_excludes_only_explicitly_unsupported_windows(self):
        times = np.arange(200) * 40960.
        times[50] += 2
        model, result = self.fit(times)
        self.assertTrue(np.isnan(model.periods[50]))
        self.assertTrue(np.isfinite(model.periods[70]))
        self.assertGreater(result["clock_support_reason_counts"]["phase_consistency"], 0)
        self.assertLessEqual(result["max_half_window_phase_disagreement_ns"], .1)
        self.assertTrue(all(row["max_phase_change_ns"] <= .1 for row in result["interval_span_sensitivity"]))

    def test_short_regions_are_not_assigned_a_period(self):
        times = np.r_[np.arange(100) * 40960., 5300000 + np.arange(20) * 40960.,
                      10000000 + np.arange(100) * 40960.]
        model, result = self.fit(times)
        self.assertTrue(np.isnan(model.periods[100:120]).all())
        self.assertEqual(result["unbridged_clock_gaps"], 2)
        self.assertEqual(result["clock_regions"][1]["sufficient_support"], 0)

    def test_insufficient_or_invalid_support_fails_explicitly(self):
        with self.assertRaisesRegex(RuntimeError, "Insufficient supported"):
            self.fit(np.arange(32) * 40960.)
        with self.assertRaisesRegex(RuntimeError, "increase"):
            self.fit(np.r_[np.arange(99), 90] * 40960.)


class PragueClockTests(unittest.TestCase):
    def test_mixed_frequency_cadences_keep_divider_and_exclude_boundaries(self) -> None:
        periods = np.r_[np.full(400, 40.00018), np.full(80, 20.00009),
                        np.full(80, 10.000045), np.full(400, 40.00018)]
        times = np.r_[0, np.cumsum(periods * 4096)]
        ids = np.arange(len(times)) + 1
        with self.assertRaisesRegex(RuntimeError, "short/duplicate"):
            calibrate_clock(times, 4096, 1100, .1, window_ids=ids)
        model, result = calibrate_cadence_clock(times, 4096, 1100, .1, [1, 2, 4], window_ids=ids)
        np.testing.assert_array_equal(np.flatnonzero(~np.isfinite(model.periods)), [400, 480, 560])
        for start, stop, expected in ((0, 400, 40.00018), (401, 480, 20.00009),
                                      (481, 560, 10.000045), (561, 961, 40.00018)):
            np.testing.assert_allclose(model.periods[start:stop], expected, rtol=0, atol=1e-9)
        self.assertEqual(result["prescaler"], 4096)
        self.assertEqual(result["missed_prescaled_ticks"], 0)
        self.assertEqual(result["unresolved_reference_windows"], 3)

    def test_missing_tick_inside_fast_region_is_not_a_slow_rate_change(self) -> None:
        periods = np.r_[np.full(400, 40.00018), np.full(80, 20.00009), np.full(400, 40.00018)]
        times = np.r_[0, np.cumsum(periods * 4096)]
        ids = np.arange(len(times)) + 1
        times[441:] += 20.00009 * 4096
        ids[441:] += 1
        model, result = calibrate_cadence_clock(times, 4096, 1100, .1, [1, 2, 4], window_ids=ids)
        self.assertEqual(result["missed_prescaled_ticks"], 1)
        self.assertEqual(result["unresolved_reference_windows"], 2)
        np.testing.assert_allclose(model.periods[401:480], 20.00009, rtol=0, atol=1e-9)

    def test_short_rate_region_is_explicitly_unresolved(self) -> None:
        periods = np.r_[np.full(100, 40.00018), np.full(5, 10.000045), np.full(100, 40.00018)]
        times = np.r_[0, np.cumsum(periods * 4096)]
        ids = np.arange(len(times)) + 1
        model, result = calibrate_cadence_clock(times, 4096, 1100, .1, [1, 2, 4], window_ids=ids)
        np.testing.assert_array_equal(np.flatnonzero(~np.isfinite(model.periods)), np.arange(100, 106))
        self.assertEqual(result["unresolved_reference_windows"], 6)

    def test_chunk_boundaries_preserve_clock_formula_and_original_window_ids(self) -> None:
        ids = np.arange(65553) + 1
        ids[100:] += 3
        times = (ids - 1) * 4096 * 40.000183 + np.sin(ids / 100) * .05
        ticks = tick_numbers(times, 4096)
        expected, _, _ = interval_periods(times, ticks, 4096, 4)
        model, result = calibrate_clock(times, 4096, 1100, .1, window_ids=ids)
        np.testing.assert_allclose(model.periods, expected, rtol=0, atol=1e-12)
        self.assertEqual(result["missed_prescaled_ticks"], 3)
        self.assertEqual(model.blocks[-1]["first_window"], int(ids[65536]))
        self.assertEqual(model.blocks[-1]["last_window"], int(ids[-1]))
        self.assertIn("NOT_CERTIFIED", result["phase_precision_status"])


class PragueBackendTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.workspace = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.workspace.name) / "clock"
        compile_cpp(ANALYSIS / "zebra_ch0_clock.cpp", cls.binary, "g++")

    @classmethod
    def tearDownClass(cls) -> None:
        cls.workspace.cleanup()

    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.cleaned = self.folder / "cleaned.tsv"
        self.decoded = self.folder / "decoded.dat"
        self.periods = self.folder / "periods.tsv"
        self.gates = [
            "--ch0-valid-rise-min-ns", "-410", "--ch0-valid-rise-max-ns", "-395",
            "--ch0-valid-tot-min-ns", "20", "--ch0-valid-tot-max-ns", "21.5",
        ]
        self.make_fixture()

    def make_fixture(self, first_passing: bool = False) -> None:
        with self.cleaned.open("w") as clean, self.decoded.open("w") as raw, self.periods.open("w") as periods:
            clean.write("# run_key=Prague_rawRefined\n")
            for event in range(5):
                raw.write(f"TDC 0 of total 2 TDCs TDC TDC_001048 size 0\n"
                          f"GlobalTriggerTime {Decimal('1.0') + Decimal(event) * Decimal('.000163840'):.9f}\n")
            for window in range(1, 129):
                stamp = Decimal("1.1") + Decimal(window - 1) * Decimal(".000163840")
                hits = [(0, -401.0, 20.66), (0, -316.5, 14.57)]
                hits += [(ch, -401 + 59.5, 17.9) for ch in range(2, 9)]
                if window == 10:
                    hits = []
                elif window == 20:
                    hits = [hit for hit in hits if hit[0] != 0]
                elif window == 30:
                    hits[0] = (0, -401.0, 14.57)
                elif window == 1 and first_passing:
                    hits[:2] = [(0, -402, 14.57), (0, -399, 20.66), (0, -405, 20.9)]
                raw.write(f"TDC 0 of total 2 TDCs TDC TDC_001048 size 0\nGlobalTriggerTime {stamp:.9f}\n")
                raw.write(f"TDC 1 of total 2 TDCs TDC TDC_001049 size {2 * len(hits)}\nGlobalTriggerTime {stamp:.9f}\n")
                for channel, rise, tot in hits:
                    raw.write(f"2 {channel} 1 {rise:.3f}\n2 {channel} 0 {rise + tot:.3f}\n")
                    clean.write(f"{window}\t1\t{stamp:.6f}\t{channel}\t2\t2\t"
                                f"{rise:.3f}\t{rise + tot:.3f}\t{tot:.3f}\n")
                if window not in (10, 20, 30):
                    periods.write(f"{window}\t40\n")

    def extract(self, *, gaps: bool = True, hardware: str = "TDC_001049") -> subprocess.CompletedProcess[str]:
        return subprocess.run([
            str(self.binary), "extract", "--input", str(self.cleaned),
            "--decoded", str(self.decoded), "--output", str(self.folder / "references.tsv"),
            "--ordinal", "1", "--hardware-id", hardware, "--coarse-quantum-ns", "5",
            "--allow-reference-gaps", "1" if gaps else "0", "--expected-windows", "128",
            *self.gates,
        ], capture_output=True, text=True)

    def fold(self, name: str, extra: list[str] | None = None) -> subprocess.CompletedProcess[str]:
        return subprocess.run([
            str(self.binary), "fold", "--input", str(self.cleaned),
            "--periods", str(self.periods), "--output", str(self.folder / name),
            "--run-key", "Prague", "--ordinal", "1",
            "--allow-reference-gaps", "1", "--expected-windows", "128", *self.gates,
            "--trigger-ref-min-ns", "-700", "--trigger-ref-max-ns", "700",
            "--trigger-ref-bins", "1400", "--ch0ref-min-ns", "-300",
            "--ch0ref-max-ns", "1100", "--ch0ref-bins", "1400",
            "--tot-max-ns", "128", "--tot-bins", "1280",
            "--histogram-channels", "0,2,3,4,5,6,7,8", *(extra or []),
        ], capture_output=True, text=True)

    def test_explicit_unresolved_clock_and_native_period_guard(self) -> None:
        self.periods.write_text("".join(f"{index}\t10\n" for index in range(1, 129)
                                       if index not in (10, 20, 30, 35)))
        unresolved = self.folder / "unresolved.tsv"
        unresolved.write_text("# transition\n35\n")
        failed = self.fold("strict_guard", ["--unresolved-clock-windows", str(unresolved)])
        self.assertNotEqual(failed.returncode, 0)
        result = self.fold("cadence", ["--unresolved-clock-windows", str(unresolved), "--period-min-ns", "8.75"])
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.folder / "cadence/channel_counts.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(int(rows[0]["without_clock_pulses"]), 2)
        self.assertEqual(int(rows[2]["without_clock_pulses"]), 1)
        self.assertEqual(int(rows[2]["folded_pulses"]), 124)
        self.assertEqual(int(rows[2]["input_pulses"]), 127)

    def test_native_scale_selected_board_matching_and_all_gap_types(self) -> None:
        result = self.extract()
        self.assertEqual(result.returncode, 0, result.stderr)
        path = self.folder / "references.tsv"
        data = np.loadtxt(path, comments="#")
        ids = np.delete(np.arange(1, 129), [9, 19, 29])
        np.testing.assert_array_equal(data[:, 0], ids)
        np.testing.assert_allclose(data[:, 4], (ids - 1) * 163840, rtol=0, atol=.001)
        metadata = parse_metadata(path)
        self.assertEqual(int(metadata["decoded_selected_board_windows"]), 128)
        self.assertEqual(int(metadata["decoded_board0_headers"]), 133)
        self.assertEqual(int(metadata["reference_gap_windows"]), 3)
        self.assertEqual(int(metadata["empty_cleaned_windows"]), 1)
        with (self.folder / "reference_gap_windows.tsv").open() as handle:
            gaps = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["reason"] for row in gaps],
                         ["empty_cleaned_window", "no_accepted_ch0", "no_passing_ch0"])

    def test_gap_handling_requires_explicit_opt_in_and_hardware_match(self) -> None:
        strict = self.extract(gaps=False)
        self.assertNotEqual(strict.returncode, 0)
        self.assertIn("Non-contiguous", strict.stderr)
        hardware = self.extract(hardware="TDC_001048")
        self.assertNotEqual(hardware.returncode, 0)
        self.assertIn("hardware identity mismatch", hardware.stderr)

    def test_first_passing_reference_is_not_minimum_time(self) -> None:
        self.make_fixture(first_passing=True)
        result = self.extract()
        self.assertEqual(result.returncode, 0, result.stderr)
        data = np.loadtxt(self.folder / "references.tsv", comments="#")
        self.assertEqual(data[0, 2], -399)
        self.assertEqual(data[0, 5], 2)

    def test_fold_accounts_for_every_pulse_and_true_modulo_three(self) -> None:
        result = self.extract()
        self.assertEqual(result.returncode, 0, result.stderr)
        metadata = parse_metadata(self.folder / "references.tsv")
        expected = np.asarray([int(value) for value in metadata["input_channel_counts"].split(",")])
        folded = self.fold("folded")
        self.assertEqual(folded.returncode, 0, folded.stderr)
        rows = verify_folded(self.folder / "folded", expected, 125)
        self.assertEqual(rows[0]["folded_pulses"], 125)
        self.assertEqual(rows[0]["without_reference_pulses"], 2)
        self.assertEqual(rows[0]["nonreference_ch0_pulses"], 125)
        for ch in range(2, 9):
            self.assertEqual(rows[ch]["input_pulses"], 127)
            self.assertEqual(rows[ch]["folded_pulses"], 125)
            self.assertEqual(rows[ch]["without_reference_pulses"], 2)
        for factor, name in ((1, "folded_matrix.tsv"), (3, "folded3_matrix.tsv")):
            meta, _, counts, _ = read_matrix(self.folder / "folded" / name)
            phase = ((59.5 / 40) % factor) * 40
            index = int(phase / float(meta["x_max_ns"]) * len(counts))
            self.assertEqual(counts[index, 2], 125)
            self.assertEqual(counts[:, 2].sum(), 125)

    def test_reference_window_without_period_is_rejected(self) -> None:
        self.periods.write_text("".join(line for line in self.periods.read_text().splitlines(True)
                                       if not line.startswith("5\t")))
        result = self.fold("bad")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no calibrated period", result.stderr)

    def test_full_population_fold_retains_nonreference_ch0_pulses(self) -> None:
        self.assertEqual(self.extract().returncode, 0)
        metadata = parse_metadata(self.folder / "references.tsv")
        expected = np.asarray([int(value) for value in metadata["input_channel_counts"].split(",")])
        result = self.fold("all_ch0", ["--include-all-ch0-hits", "1"])
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = verify_folded(self.folder / "all_ch0", expected, 125, include_all_ch0=True)
        self.assertEqual(rows[0]["folded_pulses"], 250)
        self.assertEqual(rows[0]["nonreference_ch0_pulses"], 0)
        self.assertEqual(rows[0]["without_reference_pulses"], 2)

    def test_weighted_root_ttrees_reproduce_all_plot_bins(self) -> None:
        from verify_ch0_burst_fold import write_replay
        self.assertEqual(self.extract().returncode, 0)
        result = self.fold("replay", ["--include-all-ch0-hits", "1",
                                      "--period-model", "shared_burst_rf_period"])
        self.assertEqual(result.returncode, 0, result.stderr)
        path = self.folder / "replay.root"
        write_replay(self.folder / "replay", path)
        with uproot.open(path) as root:
            self.assertEqual(sum(value == "TTree" for value in root.classnames().values()), 8)
            data = root["phase3_tot"].arrays(library="np")
            self.assertEqual(data["count"][data["channel"] == 2].sum(), 125)
        self.assertEqual(parse_metadata(self.folder / "replay/phase3_tot.tsv")["period_model"],
                         "shared_burst_rf_period")

    def test_pulse_ttree_retains_invalid_phase_rows(self) -> None:
        from verify_ch0_burst_fold import PulseTrees
        writer = PulseTrees(self.folder)
        writer.extend({"trigger_window_id": np.array([1, 2, 3], dtype=np.uint32),
                       "phase_valid": np.array([True, False, True]),
                       "rf_phase_ns": np.array([1., np.nan, 2.])})
        writer.close()
        with uproot.open(writer.parts[0]) as root:
            self.assertEqual(root["Pulses"].classname, "TTree")
            self.assertEqual(root["Pulses"].num_entries, 3)
            self.assertTrue(np.isnan(root["Pulses"]["rf_phase_ns"].array(library="np")[1]))


if __name__ == "__main__":
    unittest.main()
