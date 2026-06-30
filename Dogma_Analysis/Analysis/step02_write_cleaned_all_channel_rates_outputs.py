from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["agg.path.chunksize"] = 10000

import matplotlib.pyplot as plt
import numpy as np
import uproot
from uproot.writing.identify import to_TAxis, to_TH1x


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write PDF, PNG, and ROOT outputs for cleaned all-channel absolute and ch0-referenced rates."
    )
    parser.add_argument("--abs-matrix-file", required=True)
    parser.add_argument("--ch0ref-matrix-file", required=True)
    parser.add_argument("--abs-output-dir")
    parser.add_argument("--ch0ref-output-dir")
    return parser.parse_args()


def parse_metadata_line(line: str) -> tuple[str, str] | None:
    if not line.startswith("#") or "=" not in line:
        return None
    key, value = line[1:].split("=", 1)
    return key.strip(), value.strip()


def read_matrix(path: Path) -> tuple[dict[str, str], np.ndarray, np.ndarray, np.ndarray]:
    metadata: dict[str, str] = {}
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            parsed = parse_metadata_line(line)
            if parsed is not None:
                key, value = parsed
                metadata[key] = value
                continue
            if not line or line.startswith("#"):
                continue
            rows.append([float(value) for value in line.split()])

    if not rows:
        return metadata, np.asarray([], dtype=np.float64), np.zeros((0, 0), dtype=np.float64), np.asarray([], dtype=np.float64)

    data = np.asarray(rows, dtype=np.float64)
    time_axis = data[:, 1]
    counts = data[:, 2:-1]
    total = data[:, -1]
    return metadata, time_axis, counts, total


def channel_labels(channel_count: int) -> list[str]:
    return [f"Ch{index:02d}" for index in range(channel_count)]


def choose_time_unit(max_seconds: float) -> tuple[float, str]:
    if max_seconds >= 3600.0:
        return 3600.0, "h"
    if max_seconds >= 120.0:
        return 60.0, "min"
    return 1.0, "s"


def make_uniform_histogram(title: str, x_title: str, values: np.ndarray, x_min: float, x_max: float):
    bins = int(values.size)
    storage = np.zeros(bins + 2, dtype=np.float64)
    if bins > 0:
        storage[1:-1] = values.astype(np.float64)

    if bins > 0 and x_max > x_min:
        bin_width = (x_max - x_min) / bins
        centers = x_min + (np.arange(bins, dtype=np.float64) + 0.5) * bin_width
    else:
        centers = np.asarray([], dtype=np.float64)

    sumw = float(values.sum(dtype=np.float64))
    sumw2 = float(np.square(values, dtype=np.float64).sum(dtype=np.float64))
    sumwx = float((values.astype(np.float64) * centers).sum(dtype=np.float64)) if centers.size else 0.0
    sumwx2 = float((values.astype(np.float64) * np.square(centers)).sum(dtype=np.float64)) if centers.size else 0.0

    return to_TH1x(
        fName=None,
        fTitle=title,
        data=storage,
        fEntries=float(sumw),
        fTsumw=sumw,
        fTsumw2=sumw2,
        fTsumwx=sumwx,
        fTsumwx2=sumwx2,
        fSumw2=None,
        fXaxis=to_TAxis("xaxis", x_title, bins if bins > 0 else 1, float(x_min), float(x_max) if x_max > x_min else float(x_min + 1.0)),
    )


def write_abs_outputs(matrix_path: Path, output_dir: Path) -> None:
    metadata, time_seconds, raw_counts, raw_total = read_matrix(matrix_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_key = metadata.get("run_key", matrix_path.stem.replace("_Abs_rates_matrix", ""))
    display_bin_width_seconds = float(metadata["display_bin_width_seconds"])
    display_bin_label = f"{display_bin_width_seconds:g}"
    rate_counts_per_second = raw_counts / display_bin_width_seconds if raw_counts.size else raw_counts
    total_counts_per_second = raw_total / display_bin_width_seconds if raw_total.size else raw_total
    rate_khz = rate_counts_per_second / 1e3 if raw_counts.size else raw_counts
    total_rate_khz = total_counts_per_second / 1e3 if raw_total.size else raw_total
    labels = channel_labels(raw_counts.shape[1] if raw_counts.ndim == 2 else 0)
    nonzero = np.where(raw_counts.sum(axis=0) > 0)[0] if raw_counts.size else np.asarray([], dtype=np.int64)

    total_run_seconds = float(time_seconds[-1] + 0.5 * display_bin_width_seconds) if time_seconds.size else 0.0
    unit_scale, unit_label = choose_time_unit(total_run_seconds)

    combined_figure, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=False, height_ratios=(1, 2))
    axes[0].plot(time_seconds / unit_scale, total_rate_khz, color="#1f77b4", linewidth=0.8)
    axes[0].set_ylabel("Rate (kHz)")
    axes[0].set_title(f"Absolute rate vs run time | all channels | {display_bin_label} s display bins")
    axes[0].grid(True, linewidth=0.3, alpha=0.4)
    if nonzero.size > 0:
        heatmap = rate_khz[:, nonzero].T
        image = axes[1].imshow(
            heatmap,
            aspect="auto",
            origin="lower",
            extent=(0.0, total_run_seconds / unit_scale if total_run_seconds > 0.0 else 1.0, -0.5, nonzero.size - 0.5),
            interpolation="nearest",
            cmap="viridis",
        )
        axes[1].set_yticks(np.arange(nonzero.size, dtype=np.float64))
        axes[1].set_yticklabels([labels[index] for index in nonzero])
        combined_figure.colorbar(image, ax=axes[1], label="Rate (kHz)")
    else:
        axes[1].text(0.5, 0.5, "No nonzero channels", transform=axes[1].transAxes, ha="center", va="center")
    axes[1].set_xlabel(f"Time since run start ({unit_label})")
    axes[1].set_ylabel("Channel")
    combined_figure.tight_layout()
    combined_png = output_dir / "combined_abs_rates.png"
    combined_pdf = output_dir / "combined_abs_rates.pdf"
    combined_figure.savefig(combined_png, dpi=180)
    combined_figure.savefig(combined_pdf)
    plt.close(combined_figure)

    for channel_index in nonzero:
        figure, axis = plt.subplots(figsize=(13, 5))
        axis.plot(time_seconds / unit_scale, rate_khz[:, channel_index], color="#1f77b4", linewidth=0.8)
        axis.set_title(f"Absolute rate vs run time | {labels[channel_index]} | {display_bin_label} s display bins")
        axis.set_xlabel(f"Time since run start ({unit_label})")
        axis.set_ylabel("Rate (kHz)")
        axis.grid(True, linewidth=0.3, alpha=0.4)
        figure.tight_layout()
        png_path = output_dir / f"ch{channel_index:02d}_abs_rate.png"
        pdf_path = output_dir / f"ch{channel_index:02d}_abs_rate.pdf"
        figure.savefig(png_path, dpi=180)
        figure.savefig(pdf_path)
        plt.close(figure)

    root_path = output_dir / f"{run_key}_Abs_rates.root"
    with uproot.recreate(root_path) as root_file:
        x_min = 0.0
        x_max = total_run_seconds if total_run_seconds > 0.0 else max(display_bin_width_seconds, 1.0)
        root_file["AbsRates/hTotal_display"] = make_uniform_histogram(
            f"Absolute total rate ({display_bin_label} s display bins, kHz)",
            "Time since run start (s)",
            total_rate_khz,
            x_min,
            x_max,
        )
        for channel_index in nonzero:
            root_file[f"AbsRates/h{labels[channel_index]}_display"] = make_uniform_histogram(
                f"Absolute rate {labels[channel_index]} ({display_bin_label} s display bins, kHz)",
                "Time since run start (s)",
                rate_khz[:, channel_index],
                x_min,
                x_max,
            )


def write_ch0ref_outputs(matrix_path: Path, output_dir: Path) -> None:
    metadata, time_ns, counts, total = read_matrix(matrix_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_key = metadata.get("run_key", matrix_path.stem.replace("_Ch0_ref_Rates_matrix", ""))
    x_min = float(metadata["x_min_ns"])
    x_max = float(metadata["x_max_ns"])
    reference_label = "validated first ch0 reference hit in window"
    labels = channel_labels(counts.shape[1] if counts.ndim == 2 else 0)
    nonzero = np.where(counts.sum(axis=0) > 0)[0] if counts.size else np.asarray([], dtype=np.int64)

    combined_figure, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=False, height_ratios=(1, 2))
    axes[0].plot(time_ns, total, color="#d62728", linewidth=0.8)
    axes[0].set_ylabel("counts")
    axes[0].set_title("ch0-referenced pulse distributions | all channels | validated ch0 reference")
    axes[0].grid(True, linewidth=0.3, alpha=0.4)
    if nonzero.size > 0:
        heatmap = counts[:, nonzero].T
        image = axes[1].imshow(
            heatmap,
            aspect="auto",
            origin="lower",
            extent=(x_min, x_max, -0.5, nonzero.size - 0.5),
            interpolation="nearest",
            cmap="magma",
        )
        axes[1].set_yticks(np.arange(nonzero.size, dtype=np.float64))
        axes[1].set_yticklabels([labels[index] for index in nonzero])
        combined_figure.colorbar(image, ax=axes[1], label="counts")
    else:
        axes[1].text(0.5, 0.5, "No nonzero channels", transform=axes[1].transAxes, ha="center", va="center")
    axes[1].set_xlabel(f"Time relative to {reference_label} (ns)")
    axes[1].set_ylabel("Channel")
    combined_figure.tight_layout()
    combined_png = output_dir / "combined_ch0_ref_rates.png"
    combined_pdf = output_dir / "combined_ch0_ref_rates.pdf"
    combined_figure.savefig(combined_png, dpi=180)
    combined_figure.savefig(combined_pdf)
    plt.close(combined_figure)

    for channel_index in nonzero:
        figure, axis = plt.subplots(figsize=(13, 5))
        axis.plot(time_ns, counts[:, channel_index], color="#d62728", linewidth=0.8)
        axis.set_title(f"ch0-referenced distribution | {labels[channel_index]} | validated ch0 reference")
        axis.set_xlabel(f"Time relative to {reference_label} (ns)")
        axis.set_ylabel("counts")
        axis.grid(True, linewidth=0.3, alpha=0.4)
        figure.tight_layout()
        png_path = output_dir / f"ch{channel_index:02d}_ch0_ref_rate.png"
        pdf_path = output_dir / f"ch{channel_index:02d}_ch0_ref_rate.pdf"
        figure.savefig(png_path, dpi=180)
        figure.savefig(pdf_path)
        plt.close(figure)

    root_path = output_dir / f"{run_key}_Ch0_ref_Rates.root"
    with uproot.recreate(root_path) as root_file:
        root_file["Ch0Ref/hTotal"] = make_uniform_histogram(
            "ch0-referenced total pulse distribution",
            f"Time relative to {reference_label} (ns)",
            total,
            x_min,
            x_max,
        )
        for channel_index in nonzero:
            root_file[f"Ch0Ref/h{labels[channel_index]}"] = make_uniform_histogram(
                f"ch0-referenced pulse distribution {labels[channel_index]}",
                f"Time relative to {reference_label} (ns)",
                counts[:, channel_index],
                x_min,
                x_max,
            )


def main() -> None:
    args = parse_args()
    abs_matrix_path = Path(args.abs_matrix_file)
    ch0ref_matrix_path = Path(args.ch0ref_matrix_file)
    abs_output_dir = Path(args.abs_output_dir) if args.abs_output_dir else abs_matrix_path.parent
    ch0ref_output_dir = Path(args.ch0ref_output_dir) if args.ch0ref_output_dir else ch0ref_matrix_path.parent

    write_abs_outputs(abs_matrix_path, abs_output_dir)
    write_ch0ref_outputs(ch0ref_matrix_path, ch0ref_output_dir)


if __name__ == "__main__":
    main()