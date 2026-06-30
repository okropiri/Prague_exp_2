from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["agg.path.chunksize"] = 10000

from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import uproot
from uproot.writing.identify import to_TAxis, to_TH1x, to_TH2x


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write PNG, PDF, and ROOT outputs for cleaned all-channel RF period profiles and folded RF-phase 2D maps."
    )
    parser.add_argument("--ch0ref-matrix-file", required=True)
    parser.add_argument("--folded-matrix-file", required=True)
    parser.add_argument("--folded3x-matrix-file", required=True)
    parser.add_argument("--folded-phase-tot-file", required=True)
    parser.add_argument("--folded3x-phase-tot-file", required=True)
    parser.add_argument("--folded-phase-ch0-time-file", required=True)
    parser.add_argument("--folded-phase-trigger-time-file", required=True)
    parser.add_argument("--folded3x-phase-ch0-time-file", required=True)
    parser.add_argument("--folded3x-phase-trigger-time-file", required=True)
    parser.add_argument("--ch0ref-output-dir")
    parser.add_argument("--folded-output-dir")
    parser.add_argument("--folded3x-output-dir")
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
    x_axis = data[:, 1]
    counts = data[:, 2:-1]
    total = data[:, -1]
    return metadata, x_axis, counts, total


def read_sparse_histogram(path: Path) -> tuple[dict[str, str], dict[int, list[tuple[int, int, int]]], dict[tuple[int, int], int], int]:
    metadata: dict[str, str] = {}
    entries_by_channel: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    aggregate_entries: dict[tuple[int, int], int] = defaultdict(int)
    max_cell_count = 0
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
            parts = line.split()
            if len(parts) != 4:
                raise ValueError(f"Expected 4 columns in {path}, found: {line}")
            channel = int(parts[0])
            x_bin = int(parts[1])
            y_bin = int(parts[2])
            count = int(parts[3])
            entries_by_channel[channel].append((x_bin, y_bin, count))
            aggregate_entries[(x_bin, y_bin)] += count
            max_cell_count = max(max_cell_count, count)
    return metadata, entries_by_channel, aggregate_entries, max_cell_count


def build_dense(entries: list[tuple[int, int, int]], x_bins: int, y_bins: int) -> np.ndarray:
    counts = np.zeros((y_bins, x_bins), dtype=np.float64)
    for x_bin, y_bin, count in entries:
        counts[y_bin, x_bin] = count
    return counts


def channel_labels(channel_count: int) -> list[str]:
    return [f"Ch{index:02d}" for index in range(channel_count)]


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
        fXaxis=to_TAxis(
            "xaxis",
            x_title,
            bins if bins > 0 else 1,
            float(x_min),
            float(x_max) if x_max > x_min else float(x_min + 1.0),
        ),
    )


def make_uniform_histogram2d(
    counts: np.ndarray,
    title: str,
    x_title: str,
    y_title: str,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
):
    x_bins = counts.shape[1]
    y_bins = counts.shape[0]
    # ROOT TH2 expects x bins to vary fastest in the flattened storage.
    storage = np.zeros((y_bins + 2, x_bins + 2), dtype=np.float64)
    storage[1:-1, 1:-1] = counts

    x_centers = x_min + (np.arange(x_bins, dtype=np.float64) + 0.5) * ((x_max - x_min) / x_bins)
    y_centers = y_min + (np.arange(y_bins, dtype=np.float64) + 0.5) * ((y_max - y_min) / y_bins)
    x_centers_2d = x_centers[np.newaxis, :]
    y_centers_2d = y_centers[:, np.newaxis]
    entries = float(counts.sum(dtype=np.float64))

    return to_TH2x(
        fName=None,
        fTitle=title,
        data=storage.ravel(),
        fEntries=entries,
        fTsumw=entries,
        fTsumw2=float(np.square(counts, dtype=np.float64).sum(dtype=np.float64)),
        fTsumwx=float((counts * x_centers_2d).sum(dtype=np.float64)),
        fTsumwx2=float((counts * np.square(x_centers_2d)).sum(dtype=np.float64)),
        fTsumwy=float((counts * y_centers[:, np.newaxis]).sum(dtype=np.float64)),
        fTsumwy2=float((counts * np.square(y_centers[:, np.newaxis])).sum(dtype=np.float64)),
        fTsumwxy=float((counts * x_centers_2d * y_centers[:, np.newaxis]).sum(dtype=np.float64)),
        fSumw2=None,
        fXaxis=to_TAxis("xaxis", x_title, x_bins, float(x_min), float(x_max)),
        fYaxis=to_TAxis("yaxis", y_title, y_bins, float(y_min), float(y_max)),
    )


def nonzero_channel_indices(counts: np.ndarray) -> np.ndarray:
    if counts.size == 0:
        return np.asarray([], dtype=np.int64)
    return np.where(counts.sum(axis=0) > 0)[0]


def build_heatmap_cmap() -> colors.Colormap:
    base = mpl.colormaps["magma"]
    samples = base(np.linspace(0.08, 1.0, 256))
    cmap = colors.ListedColormap(samples, name="magma_visible")
    cmap.set_bad("#f7f7f7")
    return cmap


def build_norm(max_cell_count: float, *, logarithmic: bool):
    if max_cell_count <= 0.0:
        return colors.Normalize(vmin=0.0, vmax=1.0)
    if logarithmic:
        if max_cell_count > 1.0:
            return colors.LogNorm(vmin=1.0, vmax=max_cell_count)
        return colors.Normalize(vmin=0.0, vmax=1.0)
    return colors.Normalize(vmin=0.0, vmax=max(max_cell_count, 1.0))


def draw_heatmap(
    axis: plt.Axes,
    counts: np.ndarray,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    title: str,
    norm,
    cmap: colors.Colormap,
) -> None:
    if counts.sum() <= 0:
        axis.text(0.5, 0.5, "No entries", transform=axis.transAxes, ha="center", va="center")
        axis.set_title(title)
        return
    masked_counts = np.ma.masked_less_equal(counts, 0.0)
    axis.imshow(
        masked_counts,
        aspect="auto",
        origin="lower",
        extent=(x_min, x_max, y_min, y_max),
        interpolation="nearest",
        cmap=cmap,
        norm=norm,
    )
    axis.set_title(title)


def write_profile_outputs(
    matrix_path: Path,
    output_dir: Path,
    *,
    combined_png_name: str,
    combined_pdf_name: str,
    per_channel_suffix: str,
    root_name_suffix: str,
    color: str,
    cmap: str,
    title_prefix: str,
    x_label: str,
    root_dir_name: str,
) -> None:
    metadata, x_axis, counts, total = read_matrix(matrix_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_key = metadata.get("run_key", matrix_path.stem)
    x_min = float(metadata["x_min_ns"])
    x_max = float(metadata["x_max_ns"])
    labels = channel_labels(counts.shape[1] if counts.ndim == 2 else 0)
    nonzero = nonzero_channel_indices(counts)
    deduced_period_ns = float(metadata.get("deduced_period_ns", "0"))
    phase_origin_ns = float(metadata.get("phase_origin_ns", "0"))

    figure, axes = plt.subplots(
        2,
        1,
        figsize=(14, 8),
        sharex=False,
        height_ratios=(1, 2),
        constrained_layout=True,
    )
    axes[0].plot(x_axis, total, color=color, linewidth=0.8)
    axes[0].set_ylabel("counts")
    axes[0].set_title(
        f"{title_prefix} | all channels | period {deduced_period_ns:.6f} ns | phase origin {phase_origin_ns:.6f} ns"
    )
    axes[0].grid(True, linewidth=0.3, alpha=0.4)
    if nonzero.size > 0:
        heatmap = counts[:, nonzero].T
        image = axes[1].imshow(
            heatmap,
            aspect="auto",
            origin="lower",
            extent=(x_min, x_max, -0.5, nonzero.size - 0.5),
            interpolation="nearest",
            cmap=cmap,
        )
        axes[1].set_yticks(np.arange(nonzero.size, dtype=np.float64))
        axes[1].set_yticklabels([labels[index] for index in nonzero])
        figure.colorbar(image, ax=axes[1], label="counts")
    else:
        axes[1].text(0.5, 0.5, "No nonzero channels", transform=axes[1].transAxes, ha="center", va="center")
    axes[1].set_xlabel(x_label)
    axes[1].set_ylabel("Channel")
    figure.savefig(output_dir / combined_png_name, dpi=180)
    figure.savefig(output_dir / combined_pdf_name)
    plt.close(figure)

    for channel_index in nonzero:
        channel_figure, channel_axis = plt.subplots(figsize=(13, 5), constrained_layout=True)
        channel_axis.plot(x_axis, counts[:, channel_index], color=color, linewidth=0.8)
        channel_axis.set_title(f"{title_prefix} | {labels[channel_index]}")
        channel_axis.set_xlabel(x_label)
        channel_axis.set_ylabel("counts")
        channel_axis.grid(True, linewidth=0.3, alpha=0.4)
        channel_figure.savefig(output_dir / f"ch{channel_index:02d}_{per_channel_suffix}.png", dpi=180)
        channel_figure.savefig(output_dir / f"ch{channel_index:02d}_{per_channel_suffix}.pdf")
        plt.close(channel_figure)

    root_path = output_dir / f"{run_key}_{root_name_suffix}.root"
    with uproot.recreate(root_path) as root_file:
        root_file[f"{root_dir_name}/hTotal"] = make_uniform_histogram(
            f"{title_prefix} total",
            x_label,
            total,
            x_min,
            x_max,
        )
        for channel_index in nonzero:
            root_file[f"{root_dir_name}/h{labels[channel_index]}"] = make_uniform_histogram(
                f"{title_prefix} {labels[channel_index]}",
                x_label,
                counts[:, channel_index],
                x_min,
                x_max,
            )


def write_sparse_2d_outputs(
    sparse_path: Path,
    output_dir: Path,
    *,
    combined_filename_suffix: str,
    per_channel_filename_suffix: str,
    root_name_suffix: str,
    hist_group: str,
    x_label: str,
    y_label: str,
    title_prefix: str,
) -> None:
    metadata, channel_entries, aggregate_entries, max_cell_count = read_sparse_histogram(sparse_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_key = metadata.get("run_key", sparse_path.stem)
    x_min = float(metadata["x_min_ns"])
    x_max = float(metadata["x_max_ns"])
    x_bins = int(metadata["x_bins"])
    y_min = float(metadata.get("y_min_ns", metadata.get("tot_min_ns", "0")))
    y_max = float(metadata.get("y_max_ns", metadata.get("tot_max_ns", "1")))
    y_bins = int(metadata.get("y_bins", metadata.get("tot_bins", "1")))
    channel_count = int(metadata.get("channel_count", "32"))

    heatmap_cmap = build_heatmap_cmap()

    log_norm = build_norm(float(max_cell_count), logarithmic=True)
    linear_norm = build_norm(float(max_cell_count), logarithmic=False)

    aggregate_counts = build_dense(
        [(x_bin, y_bin, count) for (x_bin, y_bin), count in aggregate_entries.items()],
        x_bins,
        y_bins,
    )

    for suffix, norm, scale_label in (
        (combined_filename_suffix, log_norm, "log Z"),
        (combined_filename_suffix + "_linearZ", linear_norm, "linear Z"),
    ):
        combined_figure, axes = plt.subplots(4, 8, figsize=(26, 13), sharex=True, sharey=True, constrained_layout=True)
        image = None
        for channel in range(channel_count):
            axis = axes.flat[channel]
            counts = build_dense(channel_entries.get(channel, []), x_bins, y_bins)
            draw_heatmap(axis, counts, x_min, x_max, y_min, y_max, f"Ch{channel:02d}", norm, heatmap_cmap)
            if counts.sum() > 0 and image is None and axis.images:
                image = axis.images[0]
            if channel // 8 == 3:
                axis.set_xlabel(x_label)
            if channel % 8 == 0:
                axis.set_ylabel(y_label)
        combined_figure.suptitle(f"{title_prefix} | all channels | {scale_label}")
        if image is not None:
            combined_figure.colorbar(image, ax=axes.ravel().tolist(), label="counts")
        combined_figure.savefig(output_dir / f"combined_{suffix}.png", dpi=180)
        combined_figure.savefig(output_dir / f"combined_{suffix}.pdf")
        plt.close(combined_figure)

        for channel in range(channel_count):
            counts = build_dense(channel_entries.get(channel, []), x_bins, y_bins)
            if counts.sum() <= 0:
                continue
            panel_norm = build_norm(float(np.max(counts)), logarithmic=(scale_label == "log Z"))
            figure, axis = plt.subplots(figsize=(12, 6), constrained_layout=True)
            draw_heatmap(
                axis,
                counts,
                x_min,
                x_max,
                y_min,
                y_max,
                f"{title_prefix} | Ch{channel:02d} | {scale_label}",
                panel_norm,
                heatmap_cmap,
            )
            axis.set_xlabel(x_label)
            axis.set_ylabel(y_label)
            if axis.images:
                figure.colorbar(axis.images[0], ax=axis, label="counts")
            figure.savefig(output_dir / f"ch{channel:02d}_{per_channel_filename_suffix}{'_linearZ' if scale_label == 'linear Z' else ''}.png", dpi=180)
            figure.savefig(output_dir / f"ch{channel:02d}_{per_channel_filename_suffix}{'_linearZ' if scale_label == 'linear Z' else ''}.pdf")
            plt.close(figure)

    root_path = output_dir / f"{run_key}_{root_name_suffix}.root"
    with uproot.recreate(root_path) as root_file:
        root_file[f"{hist_group}/hAllChannels"] = make_uniform_histogram2d(
            aggregate_counts,
            f"{title_prefix} all channels",
            x_label,
            y_label,
            x_min,
            x_max,
            y_min,
            y_max,
        )
        for channel in range(channel_count):
            counts = build_dense(channel_entries.get(channel, []), x_bins, y_bins)
            if counts.sum() <= 0:
                continue
            root_file[f"{hist_group}/hCh{channel:02d}"] = make_uniform_histogram2d(
                counts,
                f"{title_prefix} Ch{channel:02d}",
                x_label,
                y_label,
                x_min,
                x_max,
                y_min,
                y_max,
            )


def main() -> None:
    args = parse_args()
    ch0ref_matrix_path = Path(args.ch0ref_matrix_file)
    folded_matrix_path = Path(args.folded_matrix_file)
    folded3x_matrix_path = Path(args.folded3x_matrix_file)
    folded_phase_tot_path = Path(args.folded_phase_tot_file)
    folded3x_phase_tot_path = Path(args.folded3x_phase_tot_file)
    folded_phase_ch0_time_path = Path(args.folded_phase_ch0_time_file)
    folded_phase_trigger_time_path = Path(args.folded_phase_trigger_time_file)
    folded3x_phase_ch0_time_path = Path(args.folded3x_phase_ch0_time_file)
    folded3x_phase_trigger_time_path = Path(args.folded3x_phase_trigger_time_file)

    ch0ref_output_dir = Path(args.ch0ref_output_dir) if args.ch0ref_output_dir else ch0ref_matrix_path.parent
    folded_output_dir = Path(args.folded_output_dir) if args.folded_output_dir else folded_matrix_path.parent
    folded3x_output_dir = Path(args.folded3x_output_dir) if args.folded3x_output_dir else folded3x_matrix_path.parent

    write_profile_outputs(
        ch0ref_matrix_path,
        ch0ref_output_dir,
        combined_png_name="combined_ch0_ref_rates_scan.png",
        combined_pdf_name="combined_ch0_ref_rates_scan.pdf",
        per_channel_suffix="ch0_ref_rates_scan",
        root_name_suffix="Ch0_ref_Rates_scan",
        color="#d62728",
        cmap="magma",
        title_prefix="Counts vs ch0-referenced time",
        x_label="Time relative to validated first ch0 reference hit (ns)",
        root_dir_name="Ch0RefRatesScan",
    )
    write_profile_outputs(
        folded_matrix_path,
        folded_output_dir,
        combined_png_name="combined_folded_rf.png",
        combined_pdf_name="combined_folded_rf.pdf",
        per_channel_suffix="folded_rf",
        root_name_suffix="Folded_RF",
        color="#1f77b4",
        cmap="viridis",
        title_prefix="Counts vs RF phase | 1x period",
        x_label="RF phase (ns)",
        root_dir_name="FoldedRF",
    )
    write_sparse_2d_outputs(
        folded_phase_tot_path,
        folded_output_dir,
        combined_filename_suffix="rf_phase_vs_tot",
        per_channel_filename_suffix="rf_phase_vs_tot",
        root_name_suffix="RF_phase_vs_tot",
        hist_group="RFPhaseToT",
        x_label="RF phase (ns)",
        y_label="ToT (ns)",
        title_prefix="ToT vs RF phase | 1x period",
    )
    write_sparse_2d_outputs(
        folded_phase_ch0_time_path,
        folded_output_dir,
        combined_filename_suffix="rf_phase_vs_ch0_ref_time",
        per_channel_filename_suffix="rf_phase_vs_ch0_ref_time",
        root_name_suffix="RF_phase_vs_ch0_ref_time",
        hist_group="RFPhaseCh0RefTime",
        x_label="RF phase (ns)",
        y_label="Time relative to validated first ch0 reference hit (ns)",
        title_prefix="ch0-referenced time vs RF phase | 1x period",
    )
    write_sparse_2d_outputs(
        folded_phase_trigger_time_path,
        folded_output_dir,
        combined_filename_suffix="rf_phase_vs_trigger_time",
        per_channel_filename_suffix="rf_phase_vs_trigger_time",
        root_name_suffix="RF_phase_vs_trigger_time",
        hist_group="RFPhaseTriggerTime",
        x_label="RF phase (ns)",
        y_label="Time relative to trigger window (ns)",
        title_prefix="Trigger-referenced time vs RF phase | 1x period",
    )
    write_profile_outputs(
        folded3x_matrix_path,
        folded3x_output_dir,
        combined_png_name="combined_folded_rf_3x.png",
        combined_pdf_name="combined_folded_rf_3x.pdf",
        per_channel_suffix="folded_rf_3x",
        root_name_suffix="Folded_RF_3x",
        color="#2ca02c",
        cmap="cividis",
        title_prefix="Counts vs RF phase | 3x period",
        x_label="RF phase (ns)",
        root_dir_name="FoldedRF3x",
    )
    write_sparse_2d_outputs(
        folded3x_phase_tot_path,
        folded3x_output_dir,
        combined_filename_suffix="rf_phase_3x_vs_tot",
        per_channel_filename_suffix="rf_phase_3x_vs_tot",
        root_name_suffix="RF_phase_3x_vs_tot",
        hist_group="RFPhase3xToT",
        x_label="RF phase (ns)",
        y_label="ToT (ns)",
        title_prefix="ToT vs RF phase | 3x period",
    )
    write_sparse_2d_outputs(
        folded3x_phase_ch0_time_path,
        folded3x_output_dir,
        combined_filename_suffix="rf_phase_3x_vs_ch0_ref_time",
        per_channel_filename_suffix="rf_phase_3x_vs_ch0_ref_time",
        root_name_suffix="RF_phase_3x_vs_ch0_ref_time",
        hist_group="RFPhase3xCh0RefTime",
        x_label="RF phase (ns)",
        y_label="Time relative to validated first ch0 reference hit (ns)",
        title_prefix="ch0-referenced time vs RF phase | 3x period",
    )
    write_sparse_2d_outputs(
        folded3x_phase_trigger_time_path,
        folded3x_output_dir,
        combined_filename_suffix="rf_phase_3x_vs_trigger_time",
        per_channel_filename_suffix="rf_phase_3x_vs_trigger_time",
        root_name_suffix="RF_phase_3x_vs_trigger_time",
        hist_group="RFPhase3xTriggerTime",
        x_label="RF phase (ns)",
        y_label="Time relative to trigger window (ns)",
        title_prefix="Trigger-referenced time vs RF phase | 3x period",
    )


if __name__ == "__main__":
    main()