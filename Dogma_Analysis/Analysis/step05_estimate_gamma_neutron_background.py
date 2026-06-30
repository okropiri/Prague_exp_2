from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

from matplotlib import colors, patches
import matplotlib.pyplot as plt
import numpy as np


POSITION_PATTERN = re.compile(r"Pos_([0-9.]+)m")


@dataclass(frozen=True)
class SparseHistogram:
    path: Path
    metadata: dict[str, str]
    counts: np.ndarray
    x_centers_ns: np.ndarray
    y_centers_ns: np.ndarray
    x_bin_width_ns: float
    y_bin_width_ns: float


@dataclass(frozen=True)
class RunResult:
    run_key: str
    position_m: float
    input_file: Path
    total_counts: float
    signal_band_counts: float
    gamma_phase_center_ns: float
    gamma_roi_counts: float
    gamma_roi_bins: int
    neutron_candidate_counts: float
    neutron_candidate_bins: int
    empty_sideband_counts: float
    empty_sideband_bins: int
    empty_sideband_density_per_bin: float
    high_tot_sideband_counts: float
    high_tot_sideband_bins: int
    high_tot_sideband_density_per_bin: float
    gamma_expected_background: float
    gamma_background_fraction: float
    gamma_background_subtracted: float
    neutron_expected_background: float
    neutron_background_fraction: float
    neutron_background_subtracted: float
    empty_phase_fraction: float
    empty_phase_density_cv: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate local background sidebands and first-pass gamma/neutron ROI counts "
            "from cleaned Folded_RF_3x phase-vs-ToT sparse histograms."
        )
    )
    parser.add_argument("--results-root", default="/data6/Dogma_analysis_by_Dachi/Results")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--channel", type=int, default=2)
    parser.add_argument("--run-glob", default="NCAL_20us_Pos_*")
    parser.add_argument("--tot-signal-min-ns", type=float, default=8.0)
    parser.add_argument("--tot-signal-max-ns", type=float, default=30.0)
    parser.add_argument("--tot-gamma-min-ns", type=float, default=10.0)
    parser.add_argument("--tot-gamma-max-ns", type=float, default=18.0)
    parser.add_argument("--tot-background-min-ns", type=float, default=45.0)
    parser.add_argument("--tot-background-max-ns", type=float, default=80.0)
    parser.add_argument("--gamma-phase-half-width-ns", type=float, default=2.5)
    parser.add_argument("--gamma-phase-exclusion-half-width-ns", type=float, default=7.5)
    parser.add_argument("--empty-phase-quantile", type=float, default=0.20)
    parser.add_argument("--display-tot-max-ns", type=float, default=55.0)
    return parser.parse_args()


def parse_metadata_line(line: str) -> tuple[str, str] | None:
    if not line.startswith("#") or "=" not in line:
        return None
    key, value = line[1:].split("=", 1)
    return key.strip(), value.strip()


def parse_run_key(path: Path) -> str:
    name = path.name
    suffix = "_Folded_RF_3x_phase_tot_sparse.tsv"
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return path.stem


def parse_position_m(run_key: str) -> float:
    match = POSITION_PATTERN.search(run_key)
    if match is None:
        return float("nan")
    return float(match.group(1))


def find_sparse_files(results_root: Path, run_glob: str) -> list[Path]:
    paths = sorted(
        results_root.glob(f"{run_glob}/Folded_RF_3x/*_Folded_RF_3x_phase_tot_sparse.tsv"),
        key=lambda path: (parse_position_m(parse_run_key(path)), parse_run_key(path)),
    )
    return [path for path in paths if path.is_file()]


def read_channel_histogram(path: Path, channel: int) -> SparseHistogram:
    metadata: dict[str, str] = {}
    entries: list[tuple[int, int, int]] = []
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
            entry_channel = int(parts[0])
            if entry_channel != channel:
                continue
            x_bin = int(parts[1])
            y_bin = int(parts[2])
            count = int(parts[3])
            entries.append((x_bin, y_bin, count))

    x_bins = int(metadata["x_bins"])
    y_bins = int(metadata["y_bins"])
    x_min_ns = float(metadata["x_min_ns"])
    x_max_ns = float(metadata["x_max_ns"])
    y_min_ns = float(metadata["y_min_ns"])
    y_max_ns = float(metadata["y_max_ns"])
    x_bin_width_ns = (x_max_ns - x_min_ns) / x_bins
    y_bin_width_ns = (y_max_ns - y_min_ns) / y_bins

    counts = np.zeros((y_bins, x_bins), dtype=np.float64)
    for x_bin, y_bin, count in entries:
        counts[y_bin, x_bin] = float(count)

    x_centers_ns = x_min_ns + (np.arange(x_bins, dtype=np.float64) + 0.5) * x_bin_width_ns
    y_centers_ns = y_min_ns + (np.arange(y_bins, dtype=np.float64) + 0.5) * y_bin_width_ns
    return SparseHistogram(
        path=path,
        metadata=metadata,
        counts=counts,
        x_centers_ns=x_centers_ns,
        y_centers_ns=y_centers_ns,
        x_bin_width_ns=x_bin_width_ns,
        y_bin_width_ns=y_bin_width_ns,
    )


def circular_delta_ns(values_ns: np.ndarray, center_ns: float, period_ns: float) -> np.ndarray:
    return np.abs((values_ns - center_ns + 0.5 * period_ns) % period_ns - 0.5 * period_ns)


def circular_smooth(values: np.ndarray, half_width_bins: int = 2) -> np.ndarray:
    smoothed = np.zeros_like(values, dtype=np.float64)
    weight_count = 0
    for offset in range(-half_width_bins, half_width_bins + 1):
        smoothed += np.roll(values, offset)
        weight_count += 1
    return smoothed / max(1, weight_count)


def fold_projection_to_one_period(projection: np.ndarray, period_bins: int) -> np.ndarray:
    folded = np.zeros(period_bins, dtype=np.float64)
    for x_bin, count in enumerate(projection):
        folded[x_bin % period_bins] += count
    return folded


def mask_between(values: np.ndarray, lower: float, upper: float) -> np.ndarray:
    return (values >= lower) & (values < upper)


def count_in_mask(counts: np.ndarray, y_mask: np.ndarray, x_mask: np.ndarray) -> float:
    return float(counts[np.ix_(y_mask, x_mask)].sum())


def bin_count(y_mask: np.ndarray, x_mask: np.ndarray) -> int:
    return int(np.count_nonzero(y_mask) * np.count_nonzero(x_mask))


def expected_background(density_per_bin: float, target_bins: int) -> float:
    return density_per_bin * target_bins


def safe_fraction(numerator: float, denominator: float) -> float:
    if denominator <= 0.0:
        return float("nan")
    return numerator / denominator


def contiguous_phase_ranges(mask: np.ndarray, bin_width_ns: float) -> list[tuple[float, float]]:
    ranges: list[tuple[float, float]] = []
    start_index: int | None = None
    for index, selected in enumerate(mask):
        if selected and start_index is None:
            start_index = index
        if start_index is not None and (not selected or index == len(mask) - 1):
            end_index = index + 1 if selected and index == len(mask) - 1 else index
            ranges.append((start_index * bin_width_ns, end_index * bin_width_ns))
            start_index = None
    return ranges


def periodic_ranges(center_ns: float, half_width_ns: float, period_ns: float, fold_factor: int) -> list[tuple[float, float]]:
    phase_start = center_ns - half_width_ns
    phase_end = center_ns + half_width_ns
    base_ranges: list[tuple[float, float]]
    if phase_start < 0.0:
        base_ranges = [(0.0, phase_end), (period_ns + phase_start, period_ns)]
    elif phase_end > period_ns:
        base_ranges = [(phase_start, period_ns), (0.0, phase_end - period_ns)]
    else:
        base_ranges = [(phase_start, phase_end)]

    ranges: list[tuple[float, float]] = []
    for fold_index in range(fold_factor):
        offset_ns = fold_index * period_ns
        for start_ns, end_ns in base_ranges:
            ranges.append((offset_ns + start_ns, offset_ns + end_ns))
    return ranges


def analyze_histogram(histogram: SparseHistogram, channel: int, args: argparse.Namespace) -> tuple[RunResult, dict[str, np.ndarray | float | int]]:
    run_key = histogram.metadata.get("run_key", parse_run_key(histogram.path))
    position_m = parse_position_m(run_key)
    period_ns = float(histogram.metadata["deduced_period_ns"])
    fold_factor = int(histogram.metadata.get("fold_factor", "3"))
    x_bins = histogram.counts.shape[1]
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    x_phase_ns = histogram.x_centers_ns % period_ns

    signal_tot_mask = mask_between(histogram.y_centers_ns, args.tot_signal_min_ns, args.tot_signal_max_ns)
    gamma_tot_mask = mask_between(histogram.y_centers_ns, args.tot_gamma_min_ns, args.tot_gamma_max_ns)
    high_tot_background_mask = mask_between(histogram.y_centers_ns, args.tot_background_min_ns, args.tot_background_max_ns)

    gamma_projection_3x = histogram.counts[gamma_tot_mask, :].sum(axis=0)
    gamma_projection_1x = fold_projection_to_one_period(gamma_projection_3x, period_bins)
    gamma_projection_smoothed = circular_smooth(gamma_projection_1x, half_width_bins=2)
    gamma_phase_center_ns = float(phase_centers_ns[int(np.argmax(gamma_projection_smoothed))])

    gamma_phase_mask = circular_delta_ns(x_phase_ns, gamma_phase_center_ns, period_ns) <= args.gamma_phase_half_width_ns
    signal_projection_3x = histogram.counts[signal_tot_mask, :].sum(axis=0)
    signal_projection_1x = fold_projection_to_one_period(signal_projection_3x, period_bins)
    candidate_empty_phase_mask = circular_delta_ns(phase_centers_ns, gamma_phase_center_ns, period_ns) > args.gamma_phase_exclusion_half_width_ns
    candidate_values = signal_projection_1x[candidate_empty_phase_mask]
    if candidate_values.size == 0:
        raise ValueError(f"No candidate empty phase bins for {run_key}")
    threshold = float(np.quantile(candidate_values, args.empty_phase_quantile))
    empty_phase_mask_1x = candidate_empty_phase_mask & (signal_projection_1x <= threshold)
    minimum_empty_bins = max(5, int(round(args.empty_phase_quantile * np.count_nonzero(candidate_empty_phase_mask))))
    if np.count_nonzero(empty_phase_mask_1x) < minimum_empty_bins:
        candidate_indices = np.flatnonzero(candidate_empty_phase_mask)
        ordered_indices = candidate_indices[np.argsort(signal_projection_1x[candidate_indices])]
        empty_phase_mask_1x = np.zeros(period_bins, dtype=bool)
        empty_phase_mask_1x[ordered_indices[:minimum_empty_bins]] = True
    empty_phase_mask_3x = empty_phase_mask_1x[np.arange(x_bins) % period_bins]

    all_x_mask = np.ones(x_bins, dtype=bool)
    gamma_roi_counts = count_in_mask(histogram.counts, gamma_tot_mask, gamma_phase_mask)
    gamma_roi_bins = bin_count(gamma_tot_mask, gamma_phase_mask)
    signal_band_counts = count_in_mask(histogram.counts, signal_tot_mask, all_x_mask)
    signal_band_bins = bin_count(signal_tot_mask, all_x_mask)
    neutron_candidate_counts = signal_band_counts - gamma_roi_counts
    neutron_candidate_bins = signal_band_bins - gamma_roi_bins
    empty_sideband_counts = count_in_mask(histogram.counts, signal_tot_mask, empty_phase_mask_3x)
    empty_sideband_bins = bin_count(signal_tot_mask, empty_phase_mask_3x)
    high_tot_sideband_counts = count_in_mask(histogram.counts, high_tot_background_mask, all_x_mask)
    high_tot_sideband_bins = bin_count(high_tot_background_mask, all_x_mask)

    empty_density_per_bin = safe_fraction(empty_sideband_counts, empty_sideband_bins)
    high_tot_density_per_bin = safe_fraction(high_tot_sideband_counts, high_tot_sideband_bins)
    gamma_expected_bg = expected_background(empty_density_per_bin, gamma_roi_bins)
    neutron_expected_bg = expected_background(empty_density_per_bin, neutron_candidate_bins)

    empty_phase_counts_1x = signal_projection_1x[empty_phase_mask_1x]
    empty_phase_density_cv = float(np.std(empty_phase_counts_1x) / np.mean(empty_phase_counts_1x)) if np.mean(empty_phase_counts_1x) > 0 else float("nan")
    total_counts = float(histogram.counts.sum())

    result = RunResult(
        run_key=run_key,
        position_m=position_m,
        input_file=histogram.path,
        total_counts=total_counts,
        signal_band_counts=signal_band_counts,
        gamma_phase_center_ns=gamma_phase_center_ns,
        gamma_roi_counts=gamma_roi_counts,
        gamma_roi_bins=gamma_roi_bins,
        neutron_candidate_counts=neutron_candidate_counts,
        neutron_candidate_bins=neutron_candidate_bins,
        empty_sideband_counts=empty_sideband_counts,
        empty_sideband_bins=empty_sideband_bins,
        empty_sideband_density_per_bin=empty_density_per_bin,
        high_tot_sideband_counts=high_tot_sideband_counts,
        high_tot_sideband_bins=high_tot_sideband_bins,
        high_tot_sideband_density_per_bin=high_tot_density_per_bin,
        gamma_expected_background=gamma_expected_bg,
        gamma_background_fraction=safe_fraction(gamma_expected_bg, gamma_roi_counts),
        gamma_background_subtracted=gamma_roi_counts - gamma_expected_bg,
        neutron_expected_background=neutron_expected_bg,
        neutron_background_fraction=safe_fraction(neutron_expected_bg, neutron_candidate_counts),
        neutron_background_subtracted=neutron_candidate_counts - neutron_expected_bg,
        empty_phase_fraction=safe_fraction(np.count_nonzero(empty_phase_mask_1x), period_bins),
        empty_phase_density_cv=empty_phase_density_cv,
    )

    diagnostic = {
        "period_ns": period_ns,
        "fold_factor": fold_factor,
        "gamma_phase_mask": gamma_phase_mask,
        "empty_phase_mask_3x": empty_phase_mask_3x,
        "empty_phase_mask_1x": empty_phase_mask_1x,
        "gamma_tot_mask": gamma_tot_mask,
        "signal_tot_mask": signal_tot_mask,
        "high_tot_background_mask": high_tot_background_mask,
        "signal_projection_1x": signal_projection_1x,
        "gamma_projection_1x": gamma_projection_1x,
    }
    return result, diagnostic


def draw_overlay(histogram: SparseHistogram, result: RunResult, diagnostic: dict[str, np.ndarray | float | int], output_path: Path, args: argparse.Namespace) -> None:
    period_ns = float(diagnostic["period_ns"])
    fold_factor = int(diagnostic["fold_factor"])
    display_tot_mask = histogram.y_centers_ns <= args.display_tot_max_ns
    display_counts = histogram.counts[display_tot_mask, :]
    masked_counts = np.ma.masked_less_equal(display_counts, 0.0)
    positive_counts = display_counts[display_counts > 0]
    vmax = float(np.quantile(positive_counts, 0.995)) if positive_counts.size else 1.0

    figure, axis = plt.subplots(figsize=(13.5, 5.2), constrained_layout=True)
    axis.imshow(
        masked_counts,
        aspect="auto",
        origin="lower",
        extent=(float(histogram.x_centers_ns[0] - 0.5 * histogram.x_bin_width_ns), float(histogram.x_centers_ns[-1] + 0.5 * histogram.x_bin_width_ns), 0.0, args.display_tot_max_ns),
        interpolation="nearest",
        cmap="magma",
        norm=colors.LogNorm(vmin=1.0, vmax=max(1.0, vmax)),
    )

    for start_ns, end_ns in periodic_ranges(result.gamma_phase_center_ns, args.gamma_phase_half_width_ns, period_ns, fold_factor):
        axis.add_patch(
            patches.Rectangle(
                (start_ns, args.tot_gamma_min_ns),
                end_ns - start_ns,
                args.tot_gamma_max_ns - args.tot_gamma_min_ns,
                linewidth=1.4,
                edgecolor="#00e5ff",
                facecolor="none",
                label="gamma ROI" if start_ns == periodic_ranges(result.gamma_phase_center_ns, args.gamma_phase_half_width_ns, period_ns, fold_factor)[0][0] else None,
            )
        )

    empty_phase_mask_1x = diagnostic["empty_phase_mask_1x"]
    assert isinstance(empty_phase_mask_1x, np.ndarray)
    for phase_start_ns, phase_end_ns in contiguous_phase_ranges(empty_phase_mask_1x, histogram.x_bin_width_ns):
        for fold_index in range(fold_factor):
            start_ns = fold_index * period_ns + phase_start_ns
            end_ns = fold_index * period_ns + phase_end_ns
            axis.add_patch(
                patches.Rectangle(
                    (start_ns, args.tot_signal_min_ns),
                    end_ns - start_ns,
                    args.tot_signal_max_ns - args.tot_signal_min_ns,
                    linewidth=0.0,
                    facecolor="#5ec2ff",
                    alpha=0.18,
                )
            )

    axis.axhspan(args.tot_background_min_ns, min(args.tot_background_max_ns, args.display_tot_max_ns), color="#9cff57", alpha=0.10)
    axis.axhline(args.tot_signal_min_ns, color="white", linewidth=0.7, alpha=0.7)
    axis.axhline(args.tot_signal_max_ns, color="white", linewidth=0.7, alpha=0.7)
    axis.set_xlim(0.0, float(histogram.metadata["x_max_ns"]))
    axis.set_ylim(0.0, args.display_tot_max_ns)
    axis.set_xlabel("RF phase, three periods (ns)")
    axis.set_ylabel("ToT (ns)")
    axis.set_title(
        f"{result.run_key} | detector face {result.position_m:g} m | Ch{args.channel:02d} background sidebands"
    )
    axis.text(
        0.01,
        0.98,
        "cyan: prompt-like ROI\nblue: empty-phase sideband\ngreen: high-ToT sideband",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none", "pad": 4},
    )
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_summary_table(results: list[RunResult], output_path: Path) -> None:
    fieldnames = [field for field in RunResult.__dataclass_fields__]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for result in results:
            row = {field: getattr(result, field) for field in fieldnames}
            row["input_file"] = str(result.input_file)
            writer.writerow(row)


def write_position_summary_table(results: list[RunResult], output_path: Path) -> None:
    fieldnames = [
        "position_m",
        "run_count",
        "gamma_roi_counts",
        "gamma_expected_background",
        "gamma_background_fraction",
        "gamma_background_subtracted",
        "neutron_candidate_counts",
        "neutron_expected_background",
        "neutron_background_fraction",
        "neutron_background_subtracted",
        "empty_sideband_density_per_bin",
        "high_tot_sideband_density_per_bin",
        "run_keys",
    ]
    positions = sorted({result.position_m for result in results})
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for position_m in positions:
            group = [result for result in results if result.position_m == position_m]
            gamma_roi_counts = sum(result.gamma_roi_counts for result in group)
            gamma_expected_background = sum(result.gamma_expected_background for result in group)
            neutron_candidate_counts = sum(result.neutron_candidate_counts for result in group)
            neutron_expected_background = sum(result.neutron_expected_background for result in group)
            empty_sideband_counts = sum(result.empty_sideband_counts for result in group)
            empty_sideband_bins = sum(result.empty_sideband_bins for result in group)
            high_tot_sideband_counts = sum(result.high_tot_sideband_counts for result in group)
            high_tot_sideband_bins = sum(result.high_tot_sideband_bins for result in group)
            writer.writerow(
                {
                    "position_m": position_m,
                    "run_count": len(group),
                    "gamma_roi_counts": gamma_roi_counts,
                    "gamma_expected_background": gamma_expected_background,
                    "gamma_background_fraction": safe_fraction(gamma_expected_background, gamma_roi_counts),
                    "gamma_background_subtracted": gamma_roi_counts - gamma_expected_background,
                    "neutron_candidate_counts": neutron_candidate_counts,
                    "neutron_expected_background": neutron_expected_background,
                    "neutron_background_fraction": safe_fraction(neutron_expected_background, neutron_candidate_counts),
                    "neutron_background_subtracted": neutron_candidate_counts - neutron_expected_background,
                    "empty_sideband_density_per_bin": safe_fraction(empty_sideband_counts, empty_sideband_bins),
                    "high_tot_sideband_density_per_bin": safe_fraction(high_tot_sideband_counts, high_tot_sideband_bins),
                    "run_keys": ",".join(result.run_key for result in group),
                }
            )


def draw_summary_plot(results: list[RunResult], output_path: Path) -> None:
    if not results:
        return
    positions = np.asarray([result.position_m for result in results], dtype=np.float64)
    gamma_background_fraction = np.asarray([result.gamma_background_fraction for result in results], dtype=np.float64)
    neutron_background_fraction = np.asarray([result.neutron_background_fraction for result in results], dtype=np.float64)
    empty_density = np.asarray([result.empty_sideband_density_per_bin for result in results], dtype=np.float64)
    high_density = np.asarray([result.high_tot_sideband_density_per_bin for result in results], dtype=np.float64)

    figure, axes = plt.subplots(2, 1, figsize=(9.2, 7.0), sharex=True, constrained_layout=True)
    axes[0].scatter(positions, gamma_background_fraction, label="gamma ROI sideband fraction", color="#0096c7")
    axes[0].scatter(positions, neutron_background_fraction, label="broad band sideband fraction", color="#d97706")
    unique_positions = np.asarray(sorted({result.position_m for result in results}), dtype=np.float64)
    gamma_means = np.asarray([np.nanmean(gamma_background_fraction[positions == position]) for position in unique_positions])
    neutron_means = np.asarray([np.nanmean(neutron_background_fraction[positions == position]) for position in unique_positions])
    axes[0].plot(unique_positions, gamma_means, color="#0096c7", linewidth=1.2, alpha=0.7)
    axes[0].plot(unique_positions, neutron_means, color="#d97706", linewidth=1.2, alpha=0.7)
    axes[0].set_ylabel("expected background / raw counts")
    axes[0].grid(True, linewidth=0.3, alpha=0.5)
    axes[0].legend(fontsize=8)

    axes[1].scatter(positions, empty_density, label="empty phase sideband", color="#0096c7")
    axes[1].scatter(positions, high_density, label="high-ToT sideband", color="#65a30d")
    empty_means = np.asarray([np.nanmean(empty_density[positions == position]) for position in unique_positions])
    high_means = np.asarray([np.nanmean(high_density[positions == position]) for position in unique_positions])
    axes[1].plot(unique_positions, empty_means, color="#0096c7", linewidth=1.2, alpha=0.7)
    axes[1].plot(unique_positions, high_means, color="#65a30d", linewidth=1.2, alpha=0.7)
    axes[1].set_xlabel("Detector-face distance (m)")
    axes[1].set_ylabel("counts per histogram bin")
    axes[1].grid(True, linewidth=0.3, alpha=0.5)
    axes[1].legend(fontsize=8)
    figure.suptitle("First-pass phase-ToT background sideband estimates | Ch02")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_readme(output_dir: Path, args: argparse.Namespace) -> None:
    text = f"""# Gamma/Neutron Phase-ToT Background Estimate

This directory contains a first-pass sideband estimate from existing Folded_RF_3x phase-vs-ToT sparse histograms.

Channel: Ch{args.channel:02d}

Regions used:
- Signal ToT band: [{args.tot_signal_min_ns}, {args.tot_signal_max_ns}) ns.
- Prompt-like gamma ROI: detected automatically in [{args.tot_gamma_min_ns}, {args.tot_gamma_max_ns}) ns and repeated over the three RF periods with half-width {args.gamma_phase_half_width_ns} ns.
- Empty-phase sideband: lowest {args.empty_phase_quantile:.2f} quantile of the one-period signal-band phase projection after excluding +/- {args.gamma_phase_exclusion_half_width_ns} ns around the prompt-like phase.
- High-ToT sideband: [{args.tot_background_min_ns}, {args.tot_background_max_ns}) ns over all RF phases.

Interpretation:
- The sideband estimate is a local background/pedestal estimate, not a particle-identification label.
- In overlap regions, gamma and neutron yields should be extracted with template or mixture fits after background treatment.
- The broad neutron-candidate band in the summary table is the signal ToT band after removing the prompt-like gamma ROI; it is not yet a pure neutron sample.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    results_root = Path(args.results_root).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else results_root / "Ch02_Position_Comparison" / "gamma_neutron_background"
    output_dir.mkdir(parents=True, exist_ok=True)
    overlay_dir = output_dir / "overlays"
    overlay_dir.mkdir(parents=True, exist_ok=True)

    sparse_files = find_sparse_files(results_root, args.run_glob)
    if not sparse_files:
        raise FileNotFoundError(f"No Folded_RF_3x phase-ToT sparse files found under {results_root}")

    results: list[RunResult] = []
    for sparse_path in sparse_files:
        histogram = read_channel_histogram(sparse_path, args.channel)
        result, diagnostic = analyze_histogram(histogram, args.channel, args)
        results.append(result)
        draw_overlay(histogram, result, diagnostic, overlay_dir / f"{result.run_key}_ch{args.channel:02d}_background_rois.png", args)
        print(
            f"{result.run_key}: position={result.position_m:g} m, "
            f"gamma_bg_fraction={result.gamma_background_fraction:.4f}, "
            f"neutron_candidate_bg_fraction={result.neutron_background_fraction:.4f}"
        )

    results.sort(key=lambda result: (result.position_m, result.run_key))
    write_summary_table(results, output_dir / "phase_tot_background_summary.tsv")
    write_position_summary_table(results, output_dir / "phase_tot_background_by_position.tsv")
    draw_summary_plot(results, output_dir / "phase_tot_background_summary.png")
    write_readme(output_dir, args)
    print(f"Wrote {output_dir}")


if __name__ == "__main__":
    main()