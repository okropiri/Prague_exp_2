from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import uproot
from uproot.writing.identify import to_TAxis, to_TH2x


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a ROOT TH2 from a DOGMA no-modulo time-vs-TOT grid file."
    )
    parser.add_argument("--histogram-input", required=True)
    parser.add_argument("--root-output", required=True)
    parser.add_argument("--hist-name", default="h_time_tot")
    parser.add_argument("--title", default="DOGMA time vs TOT")
    parser.add_argument("--x-title", default="time since parent ch0 [us]")
    parser.add_argument("--y-title", default="TOT [ns]")
    parser.add_argument("--x-scale", type=float, default=1.0)
    parser.add_argument("--y-scale", type=float, default=1.0)
    return parser.parse_args()


def read_grid(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] != 3:
        raise ValueError(f"Expected three columns in {path}, found shape {data.shape}")

    x_values = np.unique(data[:, 0])
    y_values = np.unique(data[:, 1])
    x_bins = x_values.size
    y_bins = y_values.size

    expected_rows = x_bins * y_bins
    if data.shape[0] != expected_rows:
        raise ValueError(
            f"Grid is incomplete: expected {expected_rows} rows from unique axes, found {data.shape[0]}"
        )

    counts = data[:, 2].reshape(y_bins, x_bins)
    return x_values, y_values, counts


def make_edges(centers: np.ndarray) -> np.ndarray:
    if centers.size < 2:
        raise ValueError("Need at least two bin centers to infer edges")
    deltas = np.diff(centers)
    if not np.allclose(deltas, deltas[0], rtol=1e-6, atol=1e-9):
        raise ValueError("Bin centers are not evenly spaced")
    half_width = 0.5 * deltas[0]
    edges = np.empty(centers.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - half_width
    edges[-1] = centers[-1] + half_width
    return edges


def write_root_histogram(
    x_centers: np.ndarray,
    y_centers: np.ndarray,
    counts: np.ndarray,
    output_path: Path,
    hist_name: str,
    title: str,
    x_title: str,
    y_title: str,
    x_scale: float,
    y_scale: float,
) -> None:
    scaled_x_centers = x_centers * x_scale
    scaled_y_centers = y_centers * y_scale
    x_edges = make_edges(scaled_x_centers)
    y_edges = make_edges(scaled_y_centers)

    x_bins = scaled_x_centers.size
    y_bins = scaled_y_centers.size

    storage = np.zeros((x_bins + 2, y_bins + 2), dtype=np.float64)
    storage[1:-1, 1:-1] = counts.T

    x_centers_2d = scaled_x_centers[np.newaxis, :]
    y_centers_2d = scaled_y_centers[:, np.newaxis]
    entries = float(counts.sum())

    histogram = to_TH2x(
        fName=None,
        fTitle=title,
        data=storage.ravel(),
        fEntries=entries,
        fTsumw=entries,
        fTsumw2=entries,
        fTsumwx=float((counts * x_centers_2d).sum()),
        fTsumwx2=float((counts * x_centers_2d * x_centers_2d).sum()),
        fTsumwy=float((counts * y_centers_2d).sum()),
        fTsumwy2=float((counts * y_centers_2d * y_centers_2d).sum()),
        fTsumwxy=float((counts * x_centers_2d * y_centers_2d).sum()),
        fSumw2=None,
        fXaxis=to_TAxis("xaxis", x_title, x_bins, float(x_edges[0]), float(x_edges[-1])),
        fYaxis=to_TAxis("yaxis", y_title, y_bins, float(y_edges[0]), float(y_edges[-1])),
    )

    with uproot.recreate(output_path) as root_file:
        root_file[hist_name] = histogram


def main() -> None:
    args = parse_args()
    histogram_input = Path(args.histogram_input)
    root_output = Path(args.root_output)

    x_centers, y_centers, counts = read_grid(histogram_input)
    write_root_histogram(
        x_centers,
        y_centers,
        counts,
        root_output,
        args.hist_name,
        args.title,
        args.x_title,
        args.y_title,
        args.x_scale,
        args.y_scale,
    )

    print(f"Histogram input: {histogram_input}")
    print(f"ROOT output: {root_output}")
    print(f"Histogram name: {args.hist_name}")
    print(f"X bins: {x_centers.size} [{x_centers[0]:.6f}, {x_centers[-1]:.6f}]")
    print(f"Y bins: {y_centers.size} [{y_centers[0]:.6f}, {y_centers[-1]:.6f}]")
    print(f"Total entries: {int(counts.sum())}")


if __name__ == "__main__":
    main()