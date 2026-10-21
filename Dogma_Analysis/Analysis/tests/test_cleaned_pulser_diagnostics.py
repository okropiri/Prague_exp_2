from __future__ import annotations

from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import uproot

ANALYSIS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS))
from run_cleaned_pulser_diagnostics import (
    cadence_labels, crop_histogram, eligibility, histogram_quantile, zoom_edges,
)


class PulserDiagnosticsTests(unittest.TestCase):
    def test_cadence_membership_keeps_transition_and_short_regions_unresolved(self) -> None:
        regions = [
            {"first_reference_window": "1", "last_reference_window": "40",
             "cadence_factor": "1", "calibrated_reference_windows": "39"},
            {"first_reference_window": "40", "last_reference_window": "80",
             "cadence_factor": "2", "calibrated_reference_windows": "39"},
            {"first_reference_window": "80", "last_reference_window": "84",
             "cadence_factor": "4", "calibrated_reference_windows": "0"},
        ]
        with self.assertRaisesRegex(RuntimeError, "boundary"):
            cadence_labels(regions, np.arange(80, 85), 84)
        labels = cadence_labels(regions, np.r_[40, np.arange(80, 85)], 84)
        np.testing.assert_array_equal(labels[1:40], np.ones(39))
        np.testing.assert_array_equal(labels[41:80], np.full(39, 2))
        self.assertTrue(np.all(labels[np.r_[0, 40, np.arange(80, 85)]] == 0))

    def test_invalid_regions_and_duplicate_uncertainty_ids_fail(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "unresolved IDs"):
            cadence_labels([], np.array([1, 1]), 4)
        with self.assertRaisesRegex(RuntimeError, "outside"):
            cadence_labels([], np.array([5]), 4)
        with self.assertRaisesRegex(RuntimeError, "overlapping"):
            cadence_labels([{"first_reference_window": "1", "last_reference_window": "5",
                             "cadence_factor": "3", "calibrated_reference_windows": "4"}],
                           np.array([], dtype=np.int64), 5)

    def test_fold_eligibility_does_not_remove_delayed_ch0_from_diagnostic_inputs(self) -> None:
        channel = np.array([0, 0, 2, 2])
        relative = np.array([0., 84., 61., np.nan])
        reference = np.array([True, False, False, False])
        np.testing.assert_array_equal(eligibility(channel, relative, reference),
                                      [True, False, True, False])
        self.assertEqual(int(np.isfinite(relative).sum()), 3)

    def test_quantiles_are_weighted_and_empty_histograms_fail(self) -> None:
        counts, edges = np.array([1, 8, 1]), np.arange(4.)
        self.assertAlmostEqual(histogram_quantile(counts, edges, .5), 1.5)
        self.assertAlmostEqual(histogram_quantile(counts, edges, .16), 1.075)
        with self.assertRaisesRegex(RuntimeError, "empty histogram"):
            histogram_quantile(np.zeros(3), edges, .5)

    def test_zoom_uses_core_not_single_tail_and_bounds_memory(self) -> None:
        edges = np.linspace(0, 128, 2561)
        counts = np.zeros(2560, dtype=np.int64)
        counts[360], counts[2400] = 10000, 1
        zoom = zoom_edges(counts, edges)
        self.assertLess(zoom[-1], 21)
        self.assertGreater(zoom[0], 15)
        self.assertLessEqual(len(zoom), 385)
        counts[:] = 1
        self.assertLessEqual(len(zoom_edges(counts, edges)), 385)

    def test_fold_zoom_preserves_original_bins_and_explicit_tail_counts(self) -> None:
        counts = np.arange(12).reshape(3, 4)
        cropped, x, y = crop_histogram(counts, np.arange(4.), np.arange(5.),
                                      np.array([1.1, 2.8]), np.array([.2, 2.]))
        np.testing.assert_array_equal(cropped, counts[1:3, :2])
        np.testing.assert_array_equal(x, [1., 2., 3.])
        np.testing.assert_array_equal(y, [0., 1., 2.])
        self.assertEqual(int(cropped.sum()), 26)
        self.assertEqual(int(counts.sum() - cropped.sum()), 40)

    def test_full_cli_preserves_counts_and_writes_all_twenty_seven_figure_pairs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            cleaned, clock, folded = folder / "cleaned", folder / "Clock_Calibration", folder / "Folded_RF_3x"
            for directory in (cleaned, clock, folded):
                directory.mkdir()
            window, channel, delay, tot, reference = [], [], [], [], []
            for wid in range(1, 51):
                for ch, time, width, selected in (
                    [(0, 0., 20., True), (0, 84. + wid * .003, 14., False)]
                    + [(ch, 60. + ch + wid * .005, 18., False) for ch in range(1, 9)]
                    + [(17, -26., 15., False)]
                ):
                    window.append(wid)
                    channel.append(ch)
                    delay.append(time)
                    tot.append(width)
                    reference.append(selected)
            data = {"trigger_window_id": np.array(window, dtype=np.uint64),
                    "channel": np.array(channel, dtype=np.int32),
                    "rise_time_ch0_ref_ns": np.array(delay), "tot_ns": np.array(tot),
                    "is_ch0_reference_hit": np.array(reference)}
            with uproot.recreate(cleaned / "fixture_pulses_part0001.root") as root:
                root.mktree("CleanedPulses", data)
            total = np.bincount(data["channel"], minlength=32)
            eligible = eligibility(data["channel"], data["rise_time_ch0_ref_ns"],
                                   data["is_ch0_reference_hit"])
            supported = data["trigger_window_id"] < 50
            folded_counts = np.bincount(data["channel"][eligible & supported], minlength=32)
            unresolved_counts = np.bincount(data["channel"][~supported], minlength=32)
            excluded = np.bincount(data["channel"][eligible & ~supported], minlength=32)
            assessment = {"prescaler": 4096, "unreferenced_channel_counts": [0] * 32,
                          "unresolved_clock_channel_counts": unresolved_counts.tolist()}
            verified = {"status": "PASS", "selected_board_windows": 50,
                        "cleaned_pulses": int(total.sum()), "folded_channel_counts": folded_counts.tolist(),
                        "eligible_clock_uncertainty_channel_counts": excluded.tolist(),
                        "unresolved_clock_reference_windows": 1}
            (folder / "assessment.json").write_text(json.dumps(assessment))
            (folder / "verification.json").write_text(json.dumps(verified))
            (clock / "unresolved_clock_windows.tsv").write_text("# window_index\n50\n")
            (clock / "cadence_regions.tsv").write_text(
                "first_reference_window\tlast_reference_window\tcadence_factor\tcalibrated_reference_windows\n"
                "1\t50\t1\t49\n")
            (clock / "cadence_regimes.tsv").write_text(
                "cadence_factor\tperiod_ns\tcalibrated_reference_windows\n"
                "1\t40.0\t49\n2\t20.0\t0\n4\t10.0\t0\n")
            with uproot.recreate(folded / "fixture_rf_phase3_vs_tot.root") as root:
                for ch in range(9):
                    mask = eligible & supported & (data["channel"] == ch)
                    root[f"RFPhase3xToT/hCh{ch:02d}"] = np.histogram2d(
                        data["rise_time_ch0_ref_ns"][mask] % 120, data["tot_ns"][mask],
                        (np.linspace(0, 120, 481), np.linspace(0, 128, 1281)))
            subprocess.run(
                [sys.executable, str(ANALYSIS / "run_cleaned_pulser_diagnostics.py"),
                 "--cleaned-root-dir", str(cleaned), "--clock-result-dir", str(folder),
                 "--verification", str(folder / "verification.json"), "--output", str(folder / "output")],
                check=True, stdout=subprocess.DEVNULL,
            )
            result = json.loads((folder / "output/Reports/verification.json").read_text())
            self.assertEqual(result["cleaned_pulses"], 550)
            self.assertEqual(result["png_files"], 27)
            self.assertEqual(result["pdf_files"], 27)
            self.assertEqual(result["references_by_factor"], [1, 49, 0, 0, 0])
            self.assertEqual(result["main_pulse_display_gate_ns"], [40, 80])
            self.assertEqual([row["core_display_pulses"] for row in result["main_pulse_display_summary"]],
                             [50] * 7)


if __name__ == "__main__":
    unittest.main()
