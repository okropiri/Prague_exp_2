from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from step03_write_cleaned_all_channel_time_tot_outputs import write_tot_distribution_outputs


class TotDisplayTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.source = self.directory / "tot.tsv"
        self.source.write_text(
            "# run_key=test\n# tot_min_ns=0\n# tot_max_ns=5000\n"
            "# tot_bins=50000\n# tot_bin_width_ns=0.1\n# channel_count=32\n"
            "0\t170\t100\n6\t150\t1000\n6\t900\t2\n6\t49835\t223\n"
        )

    def test_display_limit_preserves_full_root_counts(self):
        saved_figures = []

        def check_figure(figure, *args, **kwargs):
            for axis in figure.axes:
                if axis.get_xlabel() == "ToT (ns)":
                    np.testing.assert_allclose(axis.get_xlim(), [0, 80])
                for line in axis.lines:
                    self.assertLess(float(np.max(line.get_xdata())), 80)
                for image in axis.images:
                    self.assertEqual(image.get_array().shape[1], 800)
            if figure._suptitle is not None:
                self.assertIn("225 (16.981132%)", figure._suptitle.get_text())
            saved_figures.append(figure)

        with patch("matplotlib.figure.Figure.savefig", new=check_figure):
            write_tot_distribution_outputs(self.source, self.directory, 80)
        self.assertEqual(len(saved_figures), 14)
        with uproot.open(self.directory / "test_TOT_distrib.root") as root:
            self.assertEqual(root["ToTDistrib/hTotal"].values().sum(), 1325)
            self.assertEqual(root["ToTDistrib/hCh06"].values().sum(), 1225)
            self.assertEqual(root["ToTDistrib/hCh06"].axis().edges()[-1], 5000)
            self.assertEqual(root["ToTDistrib/hCh06"].values()[49835], 223)

    def test_default_display_retains_full_extent(self):
        def check_figure(figure, *args, **kwargs):
            for axis in figure.axes:
                for image in axis.images:
                    self.assertEqual(image.get_extent()[1], 5000)
                    self.assertEqual(image.get_array().shape[1], 50000)
        with patch("matplotlib.figure.Figure.savefig", new=check_figure):
            write_tot_distribution_outputs(self.source, self.directory)

    def test_channel_display_exclusions_use_that_channels_counts(self):
        checked = set()

        def check_figure(figure, *args, **kwargs):
            title = figure.axes[0].get_title()
            if "| Ch00 |" in title:
                self.assertIn("outside displayed bins: 0 (0.000000%)", title)
                checked.add(0)
            elif "| Ch06 |" in title:
                self.assertIn("outside displayed bins: 225 (18.367347%)", title)
                checked.add(6)
            elif figure._suptitle is not None:
                self.assertIn("outside displayed bins: 225 (16.981132%)",
                              figure._suptitle.get_text())

        with patch("matplotlib.figure.Figure.savefig", new=check_figure):
            write_tot_distribution_outputs(self.source, self.directory, 80)
        self.assertEqual(checked, {0, 6})

    def test_invalid_display_limit_is_an_error(self):
        for maximum in (0, -1, float("inf"), float("nan")):
            with self.subTest(maximum=maximum), self.assertRaises(ValueError):
                write_tot_distribution_outputs(self.source, self.directory, maximum)


if __name__ == "__main__":
    unittest.main()
