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
        description="Write PNG, PDF, and ROOT outputs for cleaned all-channel time-vs-ToT histograms."
    )
    parser.add_argument("--trigger-ref-file", required=True)
    parser.add_argument("--ch0ref-file", required=True)
    parser.add_argument("--tot-distrib-file")
    parser.add_argument("--trigger-ref-output-dir")
    parser.add_argument("--ch0ref-output-dir")
    parser.add_argument("--tot-distrib-output-dir")
    return parser.parse_args()


def parse_metadata_line(line: str) -> tuple[str, str] | None:
    if not line.startswith("#") or "=" not in line:
        return None
    key, value = line[1:].split("=", 1)
    return key.strip(), value.strip()


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


def read_tot_distribution(path: Path) -> tuple[dict[str, str], np.ndarray]:
    metadata: dict[str, str] = {}
    entries_by_channel: dict[int, list[tuple[int, int]]] = defaultdict(list)
    max_tot_bin = -1
    max_channel = -1

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
            if len(parts) != 3:
                raise ValueError(f"Expected 3 columns in {path}, found: {line}")
            channel = int(parts[0])
            tot_bin = int(parts[1])
            count = int(parts[2])
            entries_by_channel[channel].append((tot_bin, count))
            max_tot_bin = max(max_tot_bin, tot_bin)
            max_channel = max(max_channel, channel)

    metadata_tot_bins = int(metadata.get("tot_bins", "0"))
    metadata_channel_count = int(metadata.get("channel_count", "32"))
    tot_bins = max(metadata_tot_bins, max_tot_bin + 1, 1)
    channel_count = max(metadata_channel_count, max_channel + 1, 1)
    counts_by_channel = np.zeros((tot_bins, channel_count), dtype=np.float64)
    for channel, entries in entries_by_channel.items():
        for tot_bin, count in entries:
            counts_by_channel[tot_bin, channel] = count
    return metadata, counts_by_channel


def build_dense(entries: list[tuple[int, int, int]], x_bins: int, y_bins: int) -> np.ndarray:
    counts = np.zeros((y_bins, x_bins), dtype=np.float64)
    for x_bin, y_bin, count in entries:
        counts[y_bin, x_bin] = count
    return counts


def channel_labels(channel_count: int) -> list[str]:
    return [f"Ch{index:02d}" for index in range(channel_count)]


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
        fTsumwy=float((counts * y_centers_2d).sum(dtype=np.float64)),
        fTsumwy2=float((counts * np.square(y_centers_2d)).sum(dtype=np.float64)),
        fTsumwxy=float((counts * x_centers_2d * y_centers_2d).sum(dtype=np.float64)),
        fSumw2=None,
        fXaxis=to_TAxis("xaxis", x_title, x_bins, float(x_min), float(x_max)),
        fYaxis=to_TAxis("yaxis", y_title, y_bins, float(y_min), float(y_max)),
    )


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


def write_mode_outputs(
    sparse_path: Path,
    output_dir: Path,
    mode_name: str,
    filename_suffix: str,
    hist_group: str,
    x_label: str,
    title_prefix: str,
) -> None:
    metadata, channel_entries, aggregate_entries, max_cell_count = read_sparse_histogram(sparse_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_key = metadata.get("run_key", sparse_path.stem)
    x_min = float(metadata["x_min_ns"])
    x_max = float(metadata["x_max_ns"])
    x_bins = int(metadata["x_bins"])
    tot_min = float(metadata["tot_min_ns"])
    tot_max = float(metadata["tot_max_ns"])
    tot_bins = int(metadata["tot_bins"])
    channel_count = int(metadata.get("channel_count", "32"))

    heatmap_cmap = build_heatmap_cmap()
    norm = build_norm(float(max_cell_count), logarithmic=True)

    aggregate_counts = build_dense(
        [(x_bin, y_bin, count) for (x_bin, y_bin), count in aggregate_entries.items()],
        x_bins,
        tot_bins,
    )

    combined_figure, axes = plt.subplots(4, 8, figsize=(26, 13), sharex=True, sharey=True, constrained_layout=True)
    image = None
    for channel in range(channel_count):
        axis = axes.flat[channel]
        counts = build_dense(channel_entries.get(channel, []), x_bins, tot_bins)
        draw_heatmap(axis, counts, x_min, x_max, tot_min, tot_max, f"Ch{channel:02d}", norm, heatmap_cmap)
        if counts.sum() > 0 and image is None:
            image = axis.images[0]
        if channel // 8 == 3:
            axis.set_xlabel("Time (ns)")
        if channel % 8 == 0:
            axis.set_ylabel("ToT (ns)")
    combined_figure.suptitle(f"{title_prefix} | all channels")
    if image is not None:
        combined_figure.colorbar(image, ax=axes.ravel().tolist(), label="counts")
    combined_png = output_dir / f"combined_{filename_suffix}.png"
    combined_pdf = output_dir / f"combined_{filename_suffix}.pdf"
    combined_figure.savefig(combined_png, dpi=180)
    combined_figure.savefig(combined_pdf)
    plt.close(combined_figure)

    for channel in range(channel_count):
        counts = build_dense(channel_entries.get(channel, []), x_bins, tot_bins)
        if counts.sum() <= 0:
            continue
        figure, axis = plt.subplots(figsize=(12, 6))
        panel_norm = build_norm(float(np.max(counts)), logarithmic=True)
        draw_heatmap(axis, counts, x_min, x_max, tot_min, tot_max, f"{title_prefix} | Ch{channel:02d}", panel_norm, heatmap_cmap)
        axis.set_xlabel(x_label)
        axis.set_ylabel("ToT (ns)")
        if axis.images:
            figure.colorbar(axis.images[0], ax=axis, label="counts")
        figure.tight_layout()
        png_path = output_dir / f"ch{channel:02d}_{filename_suffix}.png"
        pdf_path = output_dir / f"ch{channel:02d}_{filename_suffix}.pdf"
        figure.savefig(png_path, dpi=180)
        figure.savefig(pdf_path)
        plt.close(figure)

    root_path = output_dir / f"{run_key}_{mode_name}.root"
    with uproot.recreate(root_path) as root_file:
        root_file[f"{hist_group}/hAllChannels"] = make_uniform_histogram2d(
            aggregate_counts,
            f"{title_prefix} all channels",
            x_label,
            "ToT (ns)",
            x_min,
            x_max,
            tot_min,
            tot_max,
        )
        for channel in range(channel_count):
            counts = build_dense(channel_entries.get(channel, []), x_bins, tot_bins)
            if counts.sum() <= 0:
                continue
            root_file[f"{hist_group}/hCh{channel:02d}"] = make_uniform_histogram2d(
                counts,
                f"{title_prefix} Ch{channel:02d}",
                x_label,
                "ToT (ns)",
                x_min,
                x_max,
                tot_min,
                tot_max,
            )


def write_tot_distribution_outputs(tot_distribution_path: Path, output_dir: Path) -> None:
    metadata, counts_by_channel = read_tot_distribution(tot_distribution_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_key = metadata.get("run_key", tot_distribution_path.stem)
    tot_min = float(metadata["tot_min_ns"])
    tot_max = float(metadata["tot_max_ns"])
    tot_bins = int(metadata["tot_bins"])
    channel_count = int(metadata.get("channel_count", "32"))
    labels = channel_labels(channel_count)

    bin_width = float(metadata.get("tot_bin_width_ns", "0"))
    if bin_width <= 0.0:
        bin_width = (tot_max - tot_min) / tot_bins if tot_bins > 0 else 1.0
    tot_axis = tot_min + (np.arange(tot_bins, dtype=np.float64) + 0.5) * bin_width

    total = counts_by_channel.sum(axis=1)
    nonzero = nonzero_channel_indices(counts_by_channel)
    source_label = f"all parsed pulses with ToT >= {tot_min:g} ns"
    heatmap_cmap = build_heatmap_cmap()

    def save_combined_figure(suffix: str, *, logarithmic_total: bool) -> None:
        figure, axes = plt.subplots(
            2,
            1,
            figsize=(14, 8),
            sharex=True,
            height_ratios=(1, 2),
            constrained_layout=True,
        )

        line_values = np.where(total > 0.0, total, np.nan) if logarithmic_total else total
        axes[0].plot(tot_axis, line_values, color="#ff7f0e", linewidth=0.9)
        axes[0].set_ylabel("counts")
        axes[0].set_title("Top: total counts summed over all channels")
        axes[0].grid(True, linewidth=0.3, alpha=0.4)
        if logarithmic_total:
            axes[0].set_yscale("log")

        if nonzero.size > 0:
            heatmap = counts_by_channel[:, nonzero].T
            masked_heatmap = np.ma.masked_less_equal(heatmap, 0.0)
            image = axes[1].imshow(
                masked_heatmap,
                aspect="auto",
                origin="lower",
                extent=(tot_min, tot_max, -0.5, nonzero.size - 0.5),
                interpolation="nearest",
                cmap=heatmap_cmap,
                norm=build_norm(float(np.max(heatmap)), logarithmic=True),
            )
            axes[1].set_yticks(np.arange(nonzero.size, dtype=np.float64))
            axes[1].set_yticklabels([labels[index] for index in nonzero])
            figure.colorbar(image, ax=axes[1], label="counts")
        else:
            axes[1].text(0.5, 0.5, "No nonzero channels", transform=axes[1].transAxes, ha="center", va="center")

        figure.suptitle(
            "Whole-run ToT distribution"
            f" | {source_label}"
            f" | {'log-y total counts' if logarithmic_total else 'linear total counts'}"
            " | log-z channel heatmap"
        )
        axes[1].set_title("Bottom: channel-resolved counts")
        axes[1].set_xlabel("ToT (ns)")
        axes[1].set_ylabel("Channel")
        figure.savefig(output_dir / f"combined_tot_distrib{suffix}.png", dpi=180)
        figure.savefig(output_dir / f"combined_tot_distrib{suffix}.pdf")
        plt.close(figure)

    save_combined_figure("", logarithmic_total=False)
    save_combined_figure("_logz", logarithmic_total=False)
    save_combined_figure("_logy", logarithmic_total=True)

    for channel in nonzero:
        channel_counts = counts_by_channel[:, channel]
        for suffix, logarithmic_total in (("", False), ("_logy", True)):
            figure, axis = plt.subplots(figsize=(13, 5), constrained_layout=True)
            line_values = np.where(channel_counts > 0.0, channel_counts, np.nan) if logarithmic_total else channel_counts
            axis.plot(tot_axis, line_values, color="#ff7f0e", linewidth=0.9)
            axis.set_title(
                "Whole-run ToT counts"
                f" | {source_label}"
                f" | {labels[channel]}"
                f" | {'log-y' if logarithmic_total else 'linear'}"
            )
            axis.set_xlabel("ToT (ns)")
            axis.set_ylabel("counts")
            axis.grid(True, linewidth=0.3, alpha=0.4)
            if logarithmic_total:
                axis.set_yscale("log")
            figure.savefig(output_dir / f"ch{channel:02d}_tot_distrib{suffix}.png", dpi=180)
            figure.savefig(output_dir / f"ch{channel:02d}_tot_distrib{suffix}.pdf")
            plt.close(figure)

    root_path = output_dir / f"{run_key}_TOT_distrib.root"
    with uproot.recreate(root_path) as root_file:
        root_file["ToTDistrib/hTotal"] = make_uniform_histogram(
            "Whole-run ToT counts total",
            "ToT (ns)",
            total,
            tot_min,
            tot_max,
        )
        for channel in nonzero:
            root_file[f"ToTDistrib/h{labels[channel]}"] = make_uniform_histogram(
                f"Whole-run ToT counts {labels[channel]}",
                "ToT (ns)",
                counts_by_channel[:, channel],
                tot_min,
                tot_max,
            )


def main() -> None:
    args = parse_args()
    trigger_ref_path = Path(args.trigger_ref_file)
    ch0ref_path = Path(args.ch0ref_file)
    tot_distrib_path = Path(args.tot_distrib_file) if args.tot_distrib_file else ch0ref_path
    trigger_ref_output_dir = Path(args.trigger_ref_output_dir) if args.trigger_ref_output_dir else trigger_ref_path.parent
    ch0ref_output_dir = Path(args.ch0ref_output_dir) if args.ch0ref_output_dir else ch0ref_path.parent
    tot_distrib_output_dir = Path(args.tot_distrib_output_dir) if args.tot_distrib_output_dir else tot_distrib_path.parent

    write_mode_outputs(
        trigger_ref_path,
        trigger_ref_output_dir,
        mode_name="Trigger_ref_ToT",
        filename_suffix="trigger_ref_tot",
        hist_group="TriggerRefToT",
        x_label="Rising time relative to trigger window (ns)",
        title_prefix="Trigger-referenced ToT vs rising time",
    )
    write_mode_outputs(
        ch0ref_path,
        ch0ref_output_dir,
        mode_name="Ch0_ref_TOT",
        filename_suffix="ch0_ref_tot",
        hist_group="Ch0RefToT",
        x_label="Rising time relative to validated first ch0 reference hit (ns)",
        title_prefix="ch0-referenced ToT vs rising time",
    )
    write_tot_distribution_outputs(tot_distrib_path, tot_distrib_output_dir)


if __name__ == "__main__":
    main()