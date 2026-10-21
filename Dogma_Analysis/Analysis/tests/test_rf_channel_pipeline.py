import argparse
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ANALYSIS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ANALYSIS))

import step05_estimate_gamma_neutron_background as background
import step04_run_cleaned_all_channel_rf_period as rf_runner
import run_all_cleaned_analysis_chain as batch_runner
import step06_quantify_gamma_neutron_components as components
import step07_phase_rate_slice_diagnostics as slices
import step09_gamma_window_sensitivity as sensitivity
import step10_gamma_reference_cfd as cfd


class RFChannelPipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.TemporaryDirectory()
        cls.root = Path(cls.work.name)
        cls.binary = cls.root / "rf_analyzer"
        subprocess.run([
            "g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-pedantic",
            str(ANALYSIS / "step04_dogma_cleaned_all_channel_rf_period.cpp"),
            "-o", str(cls.binary),
        ], check=True)
        cls.pulses = cls.root / "pulses.tsv"
        with cls.pulses.open("w") as handle:
            handle.write("window_index\ttdc_ordinal\tglobal_trigger_seconds\tchannel\t"
                         "rise_tdc_index_one_based\tfall_tdc_index_one_based\trise_ns\tfall_ns\ttot_ns\n")
            for window in range(40):
                for channel, times in ((0, (0, 40, 80, 120)), (2, (5, 45, 85)), (4, (10, 50, 90)), (6, (15,))):
                    for rise in times:
                        handle.write(f"{window}\t0\t{100 + window * .001:.6f}\t{channel}\t1\t1\t{rise}\t{rise + 17}\t17\n")

    @classmethod
    def tearDownClass(cls):
        cls.work.cleanup()

    def analyze(self, run_key, score_channel, extra=()):
        command = [
            str(self.binary), "--input", str(self.pulses), "--output-root", str(self.root),
            "--run-key", run_key, "--score-channel", str(score_channel),
            "--initial-period-min-ns", "39.5", "--initial-period-max-ns", "40.5",
            "--initial-step-ns", ".1", "--refine-rounds", "1", "--min-selected-pulses", "20",
            "--score-stride", "1",
            "--ch0-valid-rise-min-ns", "-1", "--ch0-valid-rise-max-ns", "1",
            "--ch0-valid-tot-min-ns", "16", "--ch0-valid-tot-max-ns", "18",
            "--ch0ref-min-ns", "-20", "--ch0ref-max-ns", "160", "--ch0ref-bins", "180",
            "--trigger-ref-min-ns", "-20", "--trigger-ref-max-ns", "160", "--trigger-ref-bins", "180",
            "--tot-max-ns", "30", "--tot-bins", "300",
            "--single-rf-segment", "--skip-cycle-residual-diagnostics", *extra,
        ]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        return self.root / run_key

    def histogram(self, folder, channel):
        return background.read_channel_histogram(
            folder / "Folded_RF_3x" / f"{folder.name}_Folded_RF_3x_phase_tot_sparse.tsv", channel
        )

    def test_explicit_score_channel_only_preserves_prague_output(self):
        folder = self.analyze("PragueOutput", 4, ("--histogram-channels", "4"))
        self.assertEqual(self.histogram(folder, 4).counts.sum(), 120)
        with self.assertRaisesRegex(ValueError, "No folded phase-ToT entries for channel 2"):
            self.histogram(folder, 2)

    def test_estimate_only_preserves_period_search_without_aligned_fold(self):
        historical = self.analyze("HistoricalPeriod", 2)
        estimate = self.analyze("EstimateOnly", 2, ("--estimate-only",))
        original = rf_runner.read_key_value_file(historical / "HistoricalPeriod_cleaned_rf_period_summary.txt")
        measured = rf_runner.read_key_value_file(estimate / "period_estimate.txt")
        self.assertAlmostEqual(float(measured["deduced_period_ns"]), float(original["deduced_period_ns"]), places=9)
        self.assertEqual(measured["phase_offset_applied_to_final_folding"], "false")
        self.assertFalse((estimate / "Folded_RF").exists())

    def test_default_saves_score_channel_only_and_zebra_list_is_opt_in(self):
        folder = self.analyze("PragueDefaultChannels", 2)
        histogram = self.histogram(folder, 2)
        self.assertEqual(histogram.counts.sum(), 120)
        self.assertEqual(histogram.metadata["stored_channels"], "2")
        for absent in (4, 6):
            with self.assertRaisesRegex(ValueError, f"No folded phase-ToT entries for channel {absent}"):
                self.histogram(folder, absent)

        folder = self.analyze("ZebraChannels", 2, ("--histogram-channels", "0,1,2,3,6"))
        histogram = self.histogram(folder, 2)
        self.assertEqual(histogram.counts.sum(), 120)
        self.assertEqual(histogram.metadata["stored_channels"], "0,1,2,3,6")
        self.assertEqual(self.histogram(folder, 6).counts.sum(), 40)
        with self.assertRaisesRegex(ValueError, "No folded phase-ToT entries for channel 4"):
            self.histogram(folder, 4)

    def test_python_prague_defaults_and_zebra_segment_override(self):
        with patch.object(sys, "argv", ["step04", "--input", "input.tsv", "--output-root", "output"]):
            args = rf_runner.parse_args()
        self.assertEqual((args.ch0_valid_rise_min_ns, args.ch0_valid_rise_max_ns), (-410, -395))
        self.assertEqual((args.ch0ref_min_ns, args.ch0ref_max_ns, args.ch0ref_bins), (-6000, 6000, 12000))
        self.assertEqual((args.trigger_ref_min_ns, args.trigger_ref_max_ns, args.trigger_ref_bins), (-6000, 6000, 2400))
        self.assertEqual((args.tot_max_ns, args.tot_bins, args.score_tot_max_ns), (128, 1280, 128))
        self.assertEqual((args.initial_period_min_ns, args.initial_period_max_ns, args.initial_step_ns), (37, 41, .05))
        self.assertIsNone(args.histogram_channels)
        self.assertEqual(args.score_stride, 1)
        self.assertEqual(args.score_channel, 2)
        self.assertFalse(args.single_rf_segment)
        with patch.object(sys, "argv", [
            "step04", "--input", "input.tsv", "--output-root", "output", "--single-rf-segment",
        ]):
            self.assertTrue(rf_runner.parse_args().single_rf_segment)

    def test_cpp_auto_segments_restores_prague_setting(self):
        folder = self.analyze("PragueSegments", 2, ("--auto-rf-segments", "--histogram-channels", "2"))
        summary = rf_runner.read_key_value_file(folder / "PragueSegments_cleaned_rf_period_summary.txt")
        self.assertEqual(summary["single_rf_segment"], "false")
        self.assertEqual(self.histogram(folder, 2).counts.sum(), 120)

    def shifted_pulses(self, name, offset_ns):
        source = self.root / name
        lines = self.pulses.read_text().splitlines()
        shifted = [lines[0]]
        for line in lines[1:]:
            fields = line.split("\t")
            fields[6] = str(float(fields[6]) + offset_ns)
            fields[7] = str(float(fields[7]) + offset_ns)
            shifted.append("\t".join(fields))
        source.write_text("\n".join(shifted) + "\n")
        return source

    def test_cpp_prague_defaults_score_only_on_prague_shaped_pulses(self):
        source = self.shifted_pulses("prague.tsv", -400)
        result = subprocess.run([
            str(self.binary), "--input", str(source), "--output-root", str(self.root),
            "--run-key", "PragueDefaults", "--min-selected-pulses", "10",
            "--skip-cycle-residual-diagnostics",
        ], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        folder = self.root / "PragueDefaults"
        summary = rf_runner.read_key_value_file(folder / "PragueDefaults_cleaned_rf_period_summary.txt")
        self.assertEqual(summary["valid_ch0_windows"], "40")
        self.assertEqual(summary["score_pulses_before_stride"], "120")
        self.assertEqual(summary["score_pulses_stored"], "120")
        self.assertEqual(summary["single_rf_segment"], "false")
        histogram = self.histogram(folder, 2)
        self.assertEqual(histogram.counts.sum(), 120)
        self.assertEqual(histogram.metadata["stored_channels"], "2")
        self.assertEqual((histogram.metadata["y_min_ns"], histogram.metadata["y_max_ns"]), ("0", "128"))

    def test_cpp_zebra_parameter_set_covers_broad_scan_without_losing_output_counts(self):
        source = self.shifted_pulses("zebra.tsv", -381)
        result = subprocess.run([
            str(self.binary), "--input", str(source), "--output-root", str(self.root),
            "--run-key", "ZebraDefaults", "--min-selected-pulses", "10",
            "--skip-cycle-residual-diagnostics",
            # Zebra parameter set, now opt-in because the shipped defaults are the Prague ones.
            "--ch0-valid-rise-min-ns", "-392", "--ch0-valid-rise-max-ns", "-374",
            "--ch0ref-min-ns", "-11000", "--ch0ref-max-ns", "2000", "--ch0ref-bins", "2600",
            "--trigger-ref-min-ns", "-11000", "--trigger-ref-max-ns", "2000", "--trigger-ref-bins", "2600",
            "--tot-max-ns", "80", "--tot-bins", "800", "--score-tot-max-ns", "80",
            "--initial-period-min-ns", "5", "--initial-period-max-ns", "200", "--initial-step-ns", "0.25",
            "--score-stride", "10", "--single-rf-segment", "--histogram-channels", "0,1,2,3,6",
        ], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        folder = self.root / "ZebraDefaults"
        summary = rf_runner.read_key_value_file(folder / "ZebraDefaults_cleaned_rf_period_summary.txt")
        self.assertEqual(summary["valid_ch0_windows"], "40")
        self.assertEqual(summary["score_pulses_before_stride"], "120")
        self.assertEqual(summary["score_pulses_stored"], "12")
        self.assertEqual(summary["single_rf_segment"], "true")
        histogram = self.histogram(folder, 2)
        self.assertEqual(histogram.counts.sum(), 120)
        self.assertEqual((histogram.metadata["y_min_ns"], histogram.metadata["y_max_ns"]), ("0", "80"))
        tables = list((folder / "RF_period_scan").glob("*_scan.txt"))
        self.assertEqual(len(tables), 1)
        periods = [float(line.split()[1]) for line in tables[0].read_text().splitlines()
                   if not line.startswith("#") and line.split()[0] == "0"]
        self.assertEqual(len(periods), 781)
        np.testing.assert_allclose(periods, np.arange(5, 200.25, .25))

    def test_ch0_score_and_detector_outputs_are_independent(self):
        folder = self.analyze("Independent", 0, ("--histogram-channels", "0,2,6"))
        self.assertEqual(self.histogram(folder, 2).counts.sum(), 120)
        self.assertEqual(self.histogram(folder, 6).counts.sum(), 40)
        histogram = self.histogram(folder, 2)
        self.assertEqual(histogram.metadata["stored_channels"], "0,2,6")
        self.assertAlmostEqual(float(histogram.metadata["deduced_period_ns"]), 40, delta=.15)
        summary = (folder / f"{folder.name}_cleaned_rf_period_summary.txt").read_text()
        self.assertIn("score_pulses_before_stride=120\n", summary)
        for axis in ("ch0_ref_time", "trigger_time"):
            table = folder / "Folded_RF_3x" / f"{folder.name}_Folded_RF_3x_phase_vs_{axis}_sparse.tsv"
            counts = {}
            for line in table.read_text().splitlines():
                if line.startswith("#"):
                    continue
                channel, _, _, count = map(int, line.split())
                counts[channel] = counts.get(channel, 0) + count
            self.assertEqual(counts[2], 120)
            self.assertEqual(counts[6], 40)

    def test_no_other_ch0_pulses_cannot_score_self_reference(self):
        single = self.root / "single.tsv"
        header, *rows = self.pulses.read_text().splitlines()
        single.write_text(header + "\n" + "\n".join(
            row for row in rows if row.split("\t")[3] != "0" or row.split("\t")[6] == "0"
        ) + "\n")
        with patch.object(self, "pulses", single):
            with self.assertRaises(AssertionError):
                self.analyze("OnlyReference", 0)

    def test_invalid_channels_rejected(self):
        with self.assertRaisesRegex(AssertionError, "histogram-channels"):
            self.analyze("InvalidOutput", 0, ("--histogram-channels", "32"))
        with self.assertRaisesRegex(AssertionError, "score-channel"):
            self.analyze("InvalidScore", 32)

    def test_empty_reference_window_rejected(self):
        folder = self.analyze("EmptyWindow", 0, ("--histogram-channels", "2"))
        histogram = self.histogram(folder, 2)
        with self.assertRaisesRegex(ValueError, "No reference entries"):
            sensitivity.detect_gamma_center_for_window(histogram, 20, 25, 2)
        args = argparse.Namespace(gamma_tot_min_ns=20, gamma_tot_max_ns=25)
        with self.assertRaisesRegex(ValueError, "No reference entries"):
            slices.detect_gamma_center(histogram, args)
        args = components.parse_args
        with patch.object(sys, "argv", ["step06"]):
            settings = args()
        settings.gamma_tot_min_ns = 20
        settings.gamma_tot_max_ns = 25
        with self.assertRaisesRegex(ValueError, "No measured prompt reference"):
            components.measure_run(histogram, settings)

    def test_peak_outside_fit_roi_rejected(self):
        path = self.root / "invalid_peak.tsv"
        fields = ["run_key", "fit_success", "peak_mu_ns", "peak_mu_err_ns", "peak_sigma_ns", "fit_tot_min_ns", "fit_tot_max_ns"]
        with path.open("w") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t")
            writer.writeheader()
            writer.writerow(dict(zip(fields, ["Mock", "True", "11", "240", "5", "12.5", "15"])))
        with self.assertRaisesRegex(ValueError, "outside its fitted ToT ROI"):
            sensitivity.read_peak_summary(path)

    def test_cfd_requested_channel_discovered_and_empty_profile_rejected(self):
        folder = self.root / "cfd_input"
        folder.mkdir()
        (folder / "gamma_window_sensitivity_summary.tsv").write_text("run_key\thalf_width_ns\nMock\t1.25\n")
        (folder / "Mock_ch06_hw1p25_gamma_reference_profile.tsv").touch()
        args = argparse.Namespace(input_dir=str(folder), results_root=str(self.root),
                                  output_dir=str(self.root / "cfd_output"), channel=6)
        zeros = np.zeros(5)
        with patch.object(cfd, "parse_args", return_value=args), patch.object(
            cfd, "parse_profile", return_value=(zeros, zeros, zeros, zeros, zeros, zeros, zeros, {})
        ):
            with self.assertRaisesRegex(ValueError, "Empty broad or selected reference"):
                cfd.main()

    def test_wrapper_passes_histogram_channels(self):
        class AnalyzerReached(Exception):
            pass

        def command_handler(command):
            if "--score-channel" in command:
                self.assertIn("--histogram-channels", command)
                self.assertEqual(command[command.index("--histogram-channels") + 1], "0,2,6")
                raise AnalyzerReached

        with patch.object(sys, "argv", [
            "step04", "--input", str(self.pulses), "--output-root", str(self.root / "wrapper"),
            "--build-dir", str(self.root / "wrapper_build"), "--score-channel", "0",
            "--histogram-channels", "0,2,6",
        ]), patch.object(rf_runner, "compile_cpp"), patch.object(rf_runner, "run_command", side_effect=command_handler):
            with self.assertRaises(AnalyzerReached):
                rf_runner.main()

    def test_batch_keeps_ch2_output_when_scoring_ch0(self):
        commands = []

        def command_handler(command, step_name):
            commands.append(command)
            if step_name == "step04_rf_period":
                raise RuntimeError("Stop after capturing step04 command")

        with patch.object(batch_runner, "run_command", side_effect=command_handler), patch.object(batch_runner, "ensure_exists"):
            batch_runner.run_chain_for_file(
                self.pulses, ANALYSIS, self.root / "batch_refined", self.root / "batch_results",
                self.root / "batch_comparison", sys.executable, "g++", self.root / "batch_build", 0, False,
            )
        rf_commands = [command for command in commands if "--score-channel" in command]
        self.assertEqual(len(rf_commands), 1)
        command = rf_commands[0]
        self.assertEqual(command[command.index("--histogram-channels") + 1], "0,2")
        for flag, value in (
            ("--ch0ref-min-ns", "-6000"), ("--ch0ref-max-ns", "6000"), ("--ch0ref-bins", "12000"),
            ("--trigger-ref-min-ns", "-6000"), ("--trigger-ref-max-ns", "6000"), ("--trigger-ref-bins", "2400"),
            ("--tot-max-ns", "128"), ("--tot-bins", "1280"), ("--score-tot-max-ns", "128"),
            ("--initial-period-min-ns", "37"), ("--initial-period-max-ns", "41"), ("--initial-step-ns", "0.05"),
            ("--score-stride", "1"), ("--ch0-valid-rise-min-ns", "-410"), ("--ch0-valid-rise-max-ns", "-395"),
        ):
            self.assertEqual(command[command.index(flag) + 1], value)
        self.assertIn("--auto-rf-segments", command)


if __name__ == "__main__":
    unittest.main()
