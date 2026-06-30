from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.colors import LogNorm
from uproot.writing.identify import to_TAxis, to_TH2x


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write ROOT and PNG outputs from a DOGMA phase-TOT histogram TSV."
    )
    parser.add_argument("--histogram-input", required=True)
    parser.add_argument("--root-output", required=True)
    parser.add_argument("--png-output", required=True)
    parser.add_argument("--title", default="DOGMA phase vs TOT")
    return parser.parse_args()


def read_histogram(path: Path) -> tuple[np.ndarray, dict[str, float]]:
    metadata: dict[str, float] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            key, value = line[1:].strip().split("=", 1)
            metadata[key.strip()] = float(value.strip())
    counts = np.loadtxt(path, comments="#", delimiter="\t")
    if counts.ndim == 1:
        counts = counts[np.newaxis, :]
    return counts, metadata


def write_root_histogram(
    counts: np.ndarray,
    metadata: dict[str, float],
    output_path: Path,
    title: str,
) -> None:
    phase_bins = int(metadata["phase_bins"])
    tot_bins = int(metadata["tot_bins"])
    phase_min = metadata["phase_min_ns"]
    phase_max = metadata["phase_max_ns"]
    tot_min = metadata["tot_min_ns"]
    tot_max = metadata["tot_max_ns"]

    phase_edges = np.linspace(phase_min, phase_max, phase_bins + 1)
    tot_edges = np.linspace(tot_min, tot_max, tot_bins + 1)

    storage = np.zeros((phase_bins + 2, tot_bins + 2), dtype=np.float64)
    storage[1:-1, 1:-1] = counts.T

    phase_centers = ((phase_edges[:-1] + phase_edges[1:]) / 2.0)[np.newaxis, :]
    tot_centers = ((tot_edges[:-1] + tot_edges[1:]) / 2.0)[:, np.newaxis]
    entries = float(counts.sum())

    histogram = to_TH2x(
        fName=None,
        fTitle=title,
        data=storage.ravel(),
        fEntries=entries,
        fTsumw=entries,
        fTsumw2=entries,
        fTsumwx=float((counts * phase_centers).sum()),
        fTsumwx2=float((counts * phase_centers * phase_centers).sum()),
        fTsumwy=float((counts * tot_centers).sum()),
        fTsumwy2=float((counts * tot_centers * tot_centers).sum()),
        fTsumwxy=float((counts * phase_centers * tot_centers).sum()),
        fSumw2=None,
        fXaxis=to_TAxis("xaxis", "Phase [ns]", phase_bins, phase_min, phase_max),
        fYaxis=to_TAxis("yaxis", "TOT [ns]", tot_bins, tot_min, tot_max),
    )

    with uproot.recreate(output_path) as root_file:
        root_file["h_phase_tot"] = histogram


def write_png(
    counts: np.ndarray,
    metadata: dict[str, float],
    output_path: Path,
    title: str,
) -> None:
    phase_min = metadata["phase_min_ns"]
    phase_max = metadata["phase_max_ns"]
    tot_min = metadata["tot_min_ns"]
    tot_max = metadata["tot_max_ns"]

    fig, ax = plt.subplots(figsize=(10, 7))
    nonzero = counts[counts > 0]
    norm = None
    if nonzero.size > 0:
        vmax = float(np.quantile(nonzero, 0.999))
        vmax = max(vmax, 1.0)
        norm = LogNorm(vmin=1.0, vmax=vmax)

    image = ax.imshow(
        counts,
        origin="lower",
        aspect="auto",
        extent=[phase_min, phase_max, tot_min, tot_max],
        cmap="magma",
        norm=norm,
    )
    ax.set_xlabel("Phase relative to parent ch0 [ns]")
    ax.set_ylabel("ch2 TOT [ns]")
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label="counts")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    histogram_input = Path(args.histogram_input)
    root_output = Path(args.root_output)
    png_output = Path(args.png_output)

    counts, metadata = read_histogram(histogram_input)
    write_root_histogram(counts, metadata, root_output, args.title)
    write_png(counts, metadata, png_output, args.title)

    print(f"Histogram input: {histogram_input}")
    print(f"ROOT output: {root_output}")
    print(f"PNG output: {png_output}")
    print(f"Total histogrammed entries: {int(counts.sum())}")


if __name__ == "__main__":
    main()