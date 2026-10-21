from pathlib import Path
import subprocess
import tempfile
import unittest


ANALYSIS = Path(__file__).resolve().parents[1]


class CleanedRatesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.work = tempfile.TemporaryDirectory()
        cls.root = Path(cls.work.name)
        cls.binary = cls.root / "rates"
        subprocess.run([
            "g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-pedantic",
            str(ANALYSIS / "step02_dogma_cleaned_all_channel_rates.cpp"),
            "-o", str(cls.binary),
        ], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.work.cleanup()

    def analyze(self, name, times, *, rise_gate=(-392, -374), tot_gate=(16.5, 19.5)):
        pulses = self.root / f"{name}.tsv"
        with pulses.open("w") as handle:
            handle.write("window_index\ttdc_ordinal\tglobal_trigger_seconds\tchannel\t"
                         "rise_tdc_index_one_based\tfall_tdc_index_one_based\trise_ns\tfall_ns\ttot_ns\n")
            for hit in times:
                window, channel, rise = hit[:3]
                tot = hit[3] if len(hit) == 4 else 17
                handle.write(f"{window}\t1\t{100 + window * .001:.6f}\t{channel}\t"
                             f"2\t2\t{rise}\t{rise + tot}\t{tot}\n")
        subprocess.run([
            str(self.binary), "--input", str(pulses), "--output-root", str(self.root),
            "--run-key", name, "--ch0-valid-rise-min-ns", str(rise_gate[0]),
            "--ch0-valid-rise-max-ns", str(rise_gate[1]),
            "--ch0-valid-tot-min-ns", str(tot_gate[0]),
            "--ch0-valid-tot-max-ns", str(tot_gate[1]),
        ], check=True, capture_output=True, text=True)
        folder = self.root / name
        summary = dict(line.split("=", 1) for line in
                       (folder / f"{name}_cleaned_rates_summary.txt").read_text().splitlines()
                       if "=" in line)
        totals = []
        for subdir, suffix in (("Abs_rates", "Abs_rates"), ("Ch0_ref_Rates", "Ch0_ref_Rates")):
            sums = [0] * 32
            path = folder / subdir / f"{name}_{suffix}_matrix.tsv"
            for line in path.read_text().splitlines():
                if line.startswith("#"):
                    continue
                values = line.split()
                for channel, value in enumerate(values[2:-1]):
                    sums[channel] += int(value)
            totals.append(sums)
        return summary, totals

    def test_earlier_detector_in_first_window_is_retained(self):
        summary, (absolute, relative) = self.analyze("Unsorted", [
            (1, 0, -383), (1, 0, -340), (1, 1, -6083), (1, 2, -6000),
            (2, 0, -384), (2, 1, -6084),
        ])
        self.assertEqual(int(summary["negative_absolute_time_rows"]), 0)
        self.assertEqual(absolute[:3], [3, 2, 1])
        self.assertEqual(relative[:3], [2, 2, 1])
        self.assertEqual(int(summary["valid_ch0_windows"]), 2)
        self.assertEqual(int(summary["parsed_pulse_rows"]), 6)

    def test_sorted_input_counts_are_unchanged(self):
        summary, (absolute, relative) = self.analyze("Sorted", [
            (1, 1, -6083), (1, 2, -6000), (1, 0, -383), (1, 0, -340),
            (2, 1, -6084), (2, 0, -384),
        ])
        self.assertEqual(int(summary["negative_absolute_time_rows"]), 0)
        self.assertEqual(absolute[:3], [3, 2, 1])
        self.assertEqual(relative[:3], [2, 2, 1])
        self.assertEqual(int(summary["valid_ch0_windows"]), 2)

    def test_pulser_doublet_and_missing_reference_preserve_absolute_counts(self):
        summary, (absolute, relative) = self.analyze("PraguePulser", [
            (1, 0, -401, 20.66), (1, 0, -317, 14.57), (1, 2, -200, 10),
            (2, 0, -401, 14.57), (2, 2, -190, 10),
            (3, 2, -180, 10),
        ], rise_gate=(-410, -395), tot_gate=(20, 21.5))
        self.assertEqual(absolute[:3], [3, 0, 3])
        self.assertEqual(relative[:3], [1, 0, 1])
        self.assertEqual(int(summary["parsed_pulse_rows"]), 6)
        self.assertEqual(int(summary["windows_seen"]), 3)
        self.assertEqual(int(summary["valid_ch0_windows"]), 1)
        self.assertEqual(int(summary["windows_without_ch0"]), 1)
        self.assertEqual(int(summary["windows_without_valid_ch0_candidate"]), 1)
        self.assertEqual(int(summary["windows_with_multiple_valid_ch0_candidates"]), 0)


if __name__ == "__main__":
    unittest.main()
