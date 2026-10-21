from pathlib import Path
import subprocess
import tempfile
import unittest


ANALYSIS = Path(__file__).resolve().parents[1]


class CleanedTimeTotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.TemporaryDirectory()
        cls.root = Path(cls.work.name)
        cls.binary = cls.root / "time_tot"
        subprocess.run([
            "g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-pedantic",
            str(ANALYSIS / "step03_dogma_cleaned_all_channel_time_tot.cpp"),
            "-o", str(cls.binary),
        ], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.work.cleanup()

    def analyze(self, name, hits):
        pulses = self.root / f"{name}.tsv"
        with pulses.open("w") as handle:
            handle.write("window_index\ttdc_ordinal\tglobal_trigger_seconds\tchannel\t"
                         "rise_tdc_index_one_based\tfall_tdc_index_one_based\trise_ns\tfall_ns\ttot_ns\n")
            for window, channel, rise, tot in hits:
                handle.write(f"{window}\t1\t{100 + window * .001:.6f}\t{channel}\t"
                             f"2\t2\t{rise}\t{rise + tot}\t{tot}\n")
        subprocess.run([
            str(self.binary), "--input", str(pulses), "--output-root", str(self.root),
            "--run-key", name,
            "--trigger-ref-min-ns", "-700", "--trigger-ref-max-ns", "700",
            "--trigger-ref-bins", "1400",
            "--ch0ref-min-ns", "-300", "--ch0ref-max-ns", "1100",
            "--ch0ref-bins", "1400",
            "--tot-min-ns", "0", "--tot-max-ns", "128", "--tot-bins", "1280",
            "--ch0-valid-rise-min-ns", "-410", "--ch0-valid-rise-max-ns", "-395",
            "--ch0-valid-tot-min-ns", "20", "--ch0-valid-tot-max-ns", "21.5",
        ], check=True, capture_output=True, text=True)
        folder = self.root / name
        summary = dict(line.split("=", 1) for line in
                       (folder / f"{name}_cleaned_time_tot_summary.txt").read_text().splitlines()
                       if "=" in line)
        return folder, summary

    def read_rows(self, path):
        return [tuple(map(int, line.split())) for line in path.read_text().splitlines()
                if line and not line.startswith("#")]

    def test_doublet_reference_gaps_and_display_overflows_preserve_full_tot(self):
        folder, summary = self.analyze("Pulser", [
            (1, 0, -401, 20.66), (1, 0, -317, 14.57), (1, 2, -342, 10),
            (1, 3, -350, 4999.5),
            (2, 0, -401, 14.57), (2, 2, -340, 10),
            (3, 2, -339, 10),
            (4, 0, -401, 20.66), (4, 2, 700, 128), (4, 4, 650, 128),
        ])
        for key, value in {
            "parsed_pulse_rows": 10, "invalid_rows": 0, "windows_seen": 4,
            "valid_ch0_windows": 2, "windows_without_ch0": 1,
            "windows_without_valid_ch0_candidate": 1,
            "trigger_ref_accepted_pulses": 7,
            "trigger_ref_rejected_outside_time_range": 0,
            "trigger_ref_rejected_outside_tot_range": 2,
            "trigger_ref_rejected_outside_time_and_tot_range": 1,
            "ch0_ref_accepted_pulses": 3,
            "ch0_ref_rejected_outside_time_range": 0,
            "ch0_ref_rejected_outside_tot_range": 2,
            "ch0_ref_rejected_outside_time_and_tot_range": 1,
            "whole_run_tot_accepted_pulses": 10,
            "whole_run_tot_rejected_below_min": 0,
        }.items():
            self.assertEqual(int(summary[key]), value, key)
        trigger = self.read_rows(folder / "Trigger_ref_ToT/Pulser_Trigger_ref_ToT_sparse.tsv")
        relative = self.read_rows(folder / "Ch0_ref_TOT/Pulser_Ch0_ref_TOT_sparse.tsv")
        widths = self.read_rows(folder / "TOT_distrib/Pulser_TOT_distrib.tsv")
        self.assertEqual(sum(row[-1] for row in trigger if row[0] == 0), 4)
        self.assertEqual(sum(row[-1] for row in relative if row[0] == 0), 2)
        self.assertEqual(sum(row[-1] for row in widths), 10)
        self.assertIn((3, 49995, 1), widths)

    def test_reference_uses_both_inclusive_gates_and_file_order_not_minimum_time(self):
        folder, summary = self.analyze("FileOrder", [
            (1, 0, -400, 14.57), (1, 0, -395, 21.5),
            (1, 0, -410, 20), (1, 2, -350, 10),
        ])
        self.assertEqual(int(summary["total_valid_ch0_candidates"]), 2)
        self.assertEqual(int(summary["windows_with_multiple_valid_ch0_candidates"]), 1)
        self.assertEqual(int(summary["rejected_ch0_outside_tot_range"]), 1)
        relative = self.read_rows(folder / "Ch0_ref_TOT/FileOrder_Ch0_ref_TOT_sparse.tsv")
        self.assertIn((2, 345, 100, 1), relative)
        self.assertIn((0, 300, 215, 1), relative)
        self.assertEqual(sum(row[-1] for row in relative), 2)


if __name__ == "__main__":
    unittest.main()
