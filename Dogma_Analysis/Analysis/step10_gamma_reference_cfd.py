from __future__ import annotations
from analysis_provenance import record_arguments

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from step05_estimate_gamma_neutron_background import read_channel_histogram
from step07_phase_rate_slice_diagnostics import parse_position_m, run_duration_seconds, signed_delta_ns


@dataclass(frozen=True)
class ReferenceResult:
    run_key: str
    position_m: float
    width_tag: str
    accepted_window_count: int
    total_window_count: int
    validation_efficiency_window: float
    assumed_trigger_window_width_ns: float
    accepted_rf_cycles: float
    gamma_tot_min_ns: float
    gamma_tot_max_ns: float
    broad_max_ref_ns: float
    broad_cfd50_ref_ns: float
    cut_max_ref_ns: float
    cut_local_peak_ref_ns: float
    cut_local_peak_shift_from_max_ns: float
    cut_cfd50_ref_ns: float
    broad_baseline_hz_per_ns: float
    broad_peak_hz_per_ns: float
    broad_threshold50_hz_per_ns: float
    cut_baseline_hz_per_ns: float
    cut_peak_hz_per_ns: float
    cut_threshold50_hz_per_ns: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare max-bin and left-edge constant-fraction gamma references for gamma-window scans."
    )
    parser.add_argument("--results-root", default="/data6/Dogma_analysis_by_Dachi/Results")
    parser.add_argument(
        "--input-dir",
        default="/data6/Dogma_analysis_by_Dachi/Results/Ch02_Position_Comparison/gamma_window_sensitivity_5m",
    )
    parser.add_argument("--output-dir")
    parser.add_argument("--channel", type=int, default=2)
    parser.add_argument("--broad-gamma-tot-min-ns", type=float, default=10.0)
    parser.add_argument("--broad-gamma-tot-max-ns", type=float, default=18.0)
    parser.add_argument("--baseline-exclusion-ns", type=float, default=6.0)
    parser.add_argument("--fraction", type=float, default=0.5)
    parser.add_argument("--gamma-display-tot-max-ns", type=float, default=30.0)
    parser.add_argument("--local-fit-half-window-bins", type=int, default=2)
    return record_arguments(parser.parse_args())


def parse_profile(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    phase = []
    broad_counts = []
    cut_counts = []
    broad_rate = []
    cut_rate = []
    broad_cycle_density = []
    cut_cycle_density = []
    metadata: dict[str, float] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            phase.append(float(row["phase_relative_ns"]))
            broad_counts.append(float(row["broad_gamma_counts"]))
            cut_counts.append(float(row["cut_gamma_counts"]))
            broad_rate.append(float(row["broad_gamma_rate_density_hz_per_ns"]))
            cut_rate.append(float(row["cut_gamma_rate_density_hz_per_ns"]))
            broad_cycle_density.append(float(row["broad_gamma_mean_counts_per_accepted_cycle_per_ns"]))
            cut_cycle_density.append(float(row["cut_gamma_mean_counts_per_accepted_cycle_per_ns"]))
            if not metadata:
                metadata = {
                    "accepted_window_count": float(row["accepted_window_count"]),
                    "total_window_count": float(row["total_window_count"]),
                    "validation_efficiency_window": float(row["validation_efficiency_window"]),
                    "assumed_trigger_window_width_ns": float(row["assumed_trigger_window_width_ns"]),
                    "accepted_rf_cycles": float(row["accepted_rf_cycles"]),
                }
    return (
        np.asarray(phase),
        np.asarray(broad_counts),
        np.asarray(cut_counts),
        np.asarray(broad_rate),
        np.asarray(cut_rate),
        np.asarray(broad_cycle_density),
        np.asarray(cut_cycle_density),
        metadata,
    )


def fold_histogram_to_one_period(histogram, period_ns: float) -> np.ndarray:
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    folded = np.zeros((histogram.counts.shape[0], period_bins), dtype=np.float64)
    for x_bin in range(histogram.counts.shape[1]):
        folded[:, x_bin % period_bins] += histogram.counts[:, x_bin]
    return folded


def left_fraction_crossing(phase: np.ndarray, rate: np.ndarray, exclusion_ns: float, fraction: float) -> tuple[float, float, float, float]:
    peak_idx = int(np.argmax(rate))
    side_vals = [value for value, phase_ns in zip(rate, phase) if abs(phase_ns) > exclusion_ns]
    baseline = float(np.quantile(side_vals, 0.20)) if side_vals else float(np.min(rate))
    peak = float(rate[peak_idx])
    threshold = baseline + fraction * (peak - baseline)
    crossing = float("nan")
    for index in range(1, peak_idx + 1):
        if rate[index - 1] < threshold <= rate[index]:
            x0, x1 = phase[index - 1], phase[index]
            y0, y1 = rate[index - 1], rate[index]
            if y1 != y0:
                crossing = float(x0 + (threshold - y0) * (x1 - x0) / (y1 - y0))
            else:
                crossing = float(x1)
            break
    return crossing, baseline, peak, threshold


def local_quadratic_peak_fit(
    phase: np.ndarray,
    rate: np.ndarray,
    half_window_bins: int,
) -> tuple[float, np.ndarray | None, np.ndarray | None]:
    peak_idx = int(np.argmax(rate))
    left = max(0, peak_idx - half_window_bins)
    right = min(len(phase), peak_idx + half_window_bins + 1)
    window_phase = phase[left:right]
    window_rate = rate[left:right]
    if window_phase.size < 3:
        return float(phase[peak_idx]), None, None
    coeffs = np.polyfit(window_phase, window_rate, deg=2)
    a, b, c = coeffs
    if a >= 0.0:
        return float(phase[peak_idx]), window_phase, np.polyval(coeffs, window_phase)
    vertex = float(-b / (2.0 * a))
    if not (float(window_phase[0]) <= vertex <= float(window_phase[-1])):
        vertex = float(phase[peak_idx])
    dense_phase = np.linspace(float(window_phase[0]), float(window_phase[-1]), 200)
    dense_rate = np.polyval(coeffs, dense_phase)
    return vertex, dense_phase, dense_rate


def draw_comparison(
    run_key: str,
    position_m: float,
    histogram,
    phase: np.ndarray,
    broad_cycle_density: np.ndarray,
    cut_cycle_density: np.ndarray,
    accepted_rf_cycles: float,
    broad_rate: np.ndarray,
    cut_rate: np.ndarray,
    broad_max_ref_raw: float,
    broad_cfd50_ref_raw: float,
    cut_max_ref_raw: float,
    cut_local_peak_ref_raw: float,
    cut_cfd50_ref_raw: float,
    broad_max_ref_relative: float,
    broad_cfd50_ref_relative: float,
    cut_max_ref_relative: float,
    cut_local_peak_ref_relative: float,
    cut_cfd50_ref_relative: float,
    cut_local_fit_phase: np.ndarray | None,
    cut_local_fit_count_density: np.ndarray | None,
    gamma_tot_min_ns: float,
    gamma_tot_max_ns: float,
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    display_tot_mask = histogram.y_centers_ns <= args.gamma_display_tot_max_ns
    period_ns = float(histogram.metadata["deduced_period_ns"])
    one_period_counts = fold_histogram_to_one_period(histogram, period_ns)
    display_counts = (np.tile(one_period_counts, (1, 3)) / max(accepted_rf_cycles, 1.0))[display_tot_mask, :]
    vmax = float(np.quantile(display_counts[display_counts > 0.0], 0.995)) if np.any(display_counts > 0.0) else 1.0
    figure, axes = plt.subplots(1, 2, figsize=(13.0, 4.8), constrained_layout=True)
    x_min_relative = float(-cut_local_peak_ref_raw)
    x_max_relative = float(3.0 * period_ns - cut_local_peak_ref_raw)
    image = axes[0].imshow(
        np.ma.masked_less_equal(display_counts, 0.0),
        aspect="auto",
        origin="lower",
        extent=(
            x_min_relative,
            x_max_relative,
            0.0,
            args.gamma_display_tot_max_ns,
        ),
        interpolation="nearest",
        cmap="magma",
        vmin=0.0,
        vmax=max(1.0e-12, vmax),
    )
    figure.colorbar(image, ax=axes[0], label="mean counts / accepted RF cycle / 2D bin, one-period fold tiled 3x")
    axes[0].axhspan(args.broad_gamma_tot_min_ns, args.broad_gamma_tot_max_ns, color="#2563eb", alpha=0.10)
    axes[0].axhspan(gamma_tot_min_ns, gamma_tot_max_ns, color="#f59e0b", alpha=0.20)
    for fold_index in range(3):
        base_relative = fold_index * period_ns
        axes[0].axvline(base_relative + broad_max_ref_relative, color="#2563eb", linestyle=":", linewidth=1.1)
        axes[0].axvline(base_relative + broad_cfd50_ref_relative, color="#7c3aed", linestyle="--", linewidth=1.1)
        axes[0].axvline(base_relative + cut_max_ref_relative, color="#00e5ff", linestyle="--", linewidth=1.2)
        axes[0].axvline(base_relative + cut_local_peak_ref_relative, color="#f59e0b", linestyle="-.", linewidth=1.2)
        axes[0].axvline(base_relative + cut_cfd50_ref_relative, color="#10b981", linestyle="--", linewidth=1.2)
    axes[0].set_xlabel("RF phase relative to cut-band local-fit reference, three periods (ns)")
    axes[0].set_ylabel("ToT (ns)")
    axes[0].set_title(f"{run_key} | one-period fold tiled 3x, local-fit centered")

    axes[1].plot(phase, broad_cycle_density, color="#2563eb", linewidth=1.2, label="broad gamma band mean density")
    axes[1].plot(phase, cut_cycle_density, color="#dc2626", linewidth=1.2, label="centered cut gamma mean density")
    if cut_local_fit_phase is not None and cut_local_fit_count_density is not None:
        axes[1].plot(cut_local_fit_phase, cut_local_fit_count_density, color="#f59e0b", linewidth=1.6, linestyle="-.", label="cut-band local quadratic fit")
    axes[1].axvline(broad_max_ref_relative, color="#2563eb", linestyle=":", linewidth=1.1, label="broad-band max")
    axes[1].axvline(broad_cfd50_ref_relative, color="#7c3aed", linestyle="--", linewidth=1.1, label="broad-band CFD50")
    axes[1].axvline(cut_max_ref_relative, color="#00e5ff", linestyle="--", linewidth=1.2, label="cut-band max")
    axes[1].axvline(0.0, color="#f59e0b", linestyle="-.", linewidth=1.2, label="cut-band local fit peak")
    axes[1].axvline(cut_cfd50_ref_relative, color="#10b981", linestyle="--", linewidth=1.2, label="cut-band CFD50")
    axes[1].set_xlabel("RF phase relative to cut-band local-fit reference (ns)")
    axes[1].set_ylabel("mean counts / accepted RF cycle / ns")
    axes[1].set_title(f"{run_key} | max vs local fit vs CFD50")
    axes[1].grid(True, linewidth=0.3, alpha=0.5)
    axes[1].legend(fontsize=8)
    annotation = (
        f"broad max: {broad_max_ref_relative:+.3f} ns\n"
        f"broad CFD50: {broad_cfd50_ref_relative:+.3f} ns\n"
        f"cut max: {cut_max_ref_relative:+.3f} ns\n"
        f"cut local fit: {0.0:+.3f} ns\n"
        f"cut CFD50: {cut_cfd50_ref_relative:+.3f} ns"
    )
    axes[0].text(
        0.02,
        0.98,
        annotation,
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=8,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.50, "edgecolor": "none", "pad": 4},
    )

    figure.suptitle(f"Prompt-gamma reference comparison | detector face {position_m:g} m")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_results(results: list[ReferenceResult], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[field.name for field in ReferenceResult.__dataclass_fields__.values()],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for result in results:
            writer.writerow({field: getattr(result, field) for field in ReferenceResult.__dataclass_fields__})


def write_readme(output_dir: Path, args: argparse.Namespace) -> None:
    text = f"""# Step10 Gamma Reference CFD Comparison

This directory compares prompt-gamma reference definitions for the processed runs in the supplied gamma-window scan input.

References compared:
- broad-band maximum,
- broad-band left-edge constant-fraction crossing at {args.fraction:.2f} of the net peak height,
- centered-cut maximum,
- centered-cut local quadratic peak fit around the maximum,
- centered-cut left-edge constant-fraction crossing at {args.fraction:.2f} of the net peak height.

Interpretation:
- the orange dash-dot line is the selected local-fit prompt-gamma reference,
- the cyan line is the actual raw cut-band maximum from the unsmoothed cut-band profile,
- the green line is the CFD-style left-edge reference on the cut-band profile,
- the purple line is the CFD-style left-edge reference on the broad-band profile,
- the blue dotted line is the broad-band maximum reference.
- `cut_local_peak_shift_from_max_ns` is the signed local-fit minus raw-cut-maximum shift.
- The right-hand 1D profiles are displayed as mean `counts / accepted RF cycle / ns`.
- The left-hand 2D panels are built from the canonical one-period fold and tiled 3x only for visibility.
- The TSV summary retains `Hz/ns` reference-shape columns and also records the accepted-window validation metadata carried from step09.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir).expanduser().resolve()
    results_root = Path(args.results_root).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else input_dir.parent / "gamma_reference_cfd_5m"
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_path = input_dir / "gamma_window_sensitivity_summary.tsv"
    summary_rows = {}
    with summary_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            summary_rows[(row["run_key"], float(row["half_width_ns"]))] = row

    results: list[ReferenceResult] = []
    channel_tag = f"_ch{args.channel:02d}_"
    profile_paths = sorted(input_dir.glob(f"*{channel_tag}hw*_gamma_reference_profile.tsv"))
    if not profile_paths:
        raise FileNotFoundError(f"No channel {args.channel} reference profiles found in {input_dir}")
    for profile_path in profile_paths:
        name = profile_path.name
        run_key, width_tag = name.rsplit(channel_tag, 1)
        width_tag = width_tag.replace("_gamma_reference_profile.tsv", "")
        try:
            half_width_ns = float(width_tag.replace("hw", "").replace("p", "."))
        except ValueError as error:
            raise ValueError(f"Invalid reference profile filename: {profile_path}") from error
        row = summary_rows.get((run_key, half_width_ns))
        if row is None:
            raise ValueError(f"No sensitivity summary for {run_key}, half width {half_width_ns}")
        position_m = parse_position_m(run_key)
        phase, broad_counts, cut_counts, broad_rate, cut_rate, broad_cycle_density, cut_cycle_density, profile_meta = parse_profile(profile_path)
        if broad_counts.sum() <= 0 or cut_counts.sum() <= 0:
            raise ValueError(f"Empty broad or selected reference in {profile_path}; CFD cannot be measured.")
        accepted_rf_cycles = float(row["accepted_rf_cycles"])
        accepted_window_count = int(float(row["accepted_window_count"]))
        total_window_count = int(float(row["total_window_count"]))
        validation_efficiency_window = float(row["validation_efficiency_window"])
        assumed_trigger_window_width_ns = float(row["assumed_trigger_window_width_ns"])
        broad_max_ref_relative_raw = float(phase[int(np.argmax(broad_rate))])
        broad_cfd50_ref_relative_raw, broad_baseline, broad_peak, broad_threshold = left_fraction_crossing(phase, broad_rate, args.baseline_exclusion_ns, args.fraction)
        cut_max_ref_relative_raw = float(phase[int(np.argmax(cut_rate))])
        cut_local_peak_ref_relative_raw, cut_local_fit_phase, cut_local_fit_rate = local_quadratic_peak_fit(phase, cut_rate, args.local_fit_half_window_bins)
        cut_cfd50_ref_relative_raw, baseline, peak, threshold = left_fraction_crossing(phase, cut_rate, args.baseline_exclusion_ns, args.fraction)
        broad_max_ref_relative = broad_max_ref_relative_raw - cut_local_peak_ref_relative_raw
        broad_cfd50_ref_relative = broad_cfd50_ref_relative_raw - cut_local_peak_ref_relative_raw
        cut_max_ref_relative = cut_max_ref_relative_raw - cut_local_peak_ref_relative_raw
        cut_cfd50_ref_relative = cut_cfd50_ref_relative_raw - cut_local_peak_ref_relative_raw
        cut_local_peak_ref_relative = 0.0
        if cut_local_fit_phase is not None:
            cut_local_fit_phase = cut_local_fit_phase - cut_local_peak_ref_relative_raw

        histogram = read_channel_histogram(results_root / run_key / "Folded_RF_3x" / f"{run_key}_Folded_RF_3x_phase_tot_sparse.tsv", args.channel)
        duration_seconds = run_duration_seconds(histogram, run_key)
        cut_local_fit_count_density = cut_local_fit_rate * duration_seconds / max(accepted_rf_cycles, 1.0) if cut_local_fit_rate is not None else None
        input_center_raw = float(row["phase_center_ns"])
        broad_max_ref_raw = input_center_raw + broad_max_ref_relative_raw
        broad_cfd50_ref_raw = input_center_raw + broad_cfd50_ref_relative_raw
        cut_max_ref_raw = input_center_raw + cut_max_ref_relative_raw
        cut_local_peak_ref_raw = input_center_raw + cut_local_peak_ref_relative_raw
        cut_cfd50_ref_raw = input_center_raw + cut_cfd50_ref_relative_raw
        gamma_tot_min_ns = float(row["gamma_tot_min_ns"])
        gamma_tot_max_ns = float(row["gamma_tot_max_ns"])

        draw_comparison(
            run_key,
            position_m,
            histogram,
            phase,
            broad_cycle_density,
            cut_cycle_density,
            accepted_rf_cycles,
            broad_rate,
            cut_rate,
            broad_max_ref_raw,
            broad_cfd50_ref_raw,
            cut_max_ref_raw,
            cut_local_peak_ref_raw,
            cut_cfd50_ref_raw,
            broad_max_ref_relative,
            broad_cfd50_ref_relative,
            cut_max_ref_relative,
            cut_local_peak_ref_relative,
            cut_cfd50_ref_relative,
            cut_local_fit_phase,
            cut_local_fit_count_density,
            gamma_tot_min_ns,
            gamma_tot_max_ns,
            output_dir / f"{run_key}_ch02_{width_tag}_gamma_reference_cfd50.png",
            args,
        )
        results.append(
            ReferenceResult(
                run_key=run_key,
                position_m=position_m,
                width_tag=width_tag,
                accepted_window_count=accepted_window_count,
                total_window_count=total_window_count,
                validation_efficiency_window=validation_efficiency_window,
                assumed_trigger_window_width_ns=assumed_trigger_window_width_ns,
                accepted_rf_cycles=accepted_rf_cycles,
                gamma_tot_min_ns=gamma_tot_min_ns,
                gamma_tot_max_ns=gamma_tot_max_ns,
                broad_max_ref_ns=broad_max_ref_relative,
                broad_cfd50_ref_ns=broad_cfd50_ref_relative,
                cut_max_ref_ns=cut_max_ref_relative,
                cut_local_peak_ref_ns=cut_local_peak_ref_relative,
                cut_local_peak_shift_from_max_ns=cut_local_peak_ref_relative_raw - cut_max_ref_relative_raw,
                cut_cfd50_ref_ns=cut_cfd50_ref_relative,
                broad_baseline_hz_per_ns=broad_baseline,
                broad_peak_hz_per_ns=broad_peak,
                broad_threshold50_hz_per_ns=broad_threshold,
                cut_baseline_hz_per_ns=baseline,
                cut_peak_hz_per_ns=peak,
                cut_threshold50_hz_per_ns=threshold,
            )
        )
        print(
            f"{run_key} {width_tag}: broad_max={broad_max_ref_relative:.4f} ns, "
            f"cut_max={cut_max_ref_relative:.4f} ns, cut_local_fit={cut_local_peak_ref_relative:.4f} ns, "
            f"broad_cfd50={broad_cfd50_ref_relative:.4f} ns, cut_cfd50={cut_cfd50_ref_relative:.4f} ns"
        )

    write_results(results, output_dir / "gamma_reference_cfd50_summary.tsv")
    write_readme(output_dir, args)
    print(f"Wrote {output_dir}")


if __name__ == "__main__":
    main()
