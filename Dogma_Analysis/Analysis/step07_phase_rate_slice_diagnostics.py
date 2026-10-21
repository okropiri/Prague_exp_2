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

from step05_estimate_gamma_neutron_background import (
    SparseHistogram,
    circular_smooth,
    find_sparse_files,
    fold_projection_to_one_period,
    mask_between,
    parse_run_key,
    read_channel_histogram,
)


@dataclass(frozen=True)
class SliceProfile:
    run_key: str
    position_m: float
    tot_min_ns: float
    tot_max_ns: float
    phase_center_ns: float
    gamma_relative_phase_ns: np.ndarray
    rate_hz: np.ndarray
    counts: np.ndarray
    baseline_rate_hz: float
    baseline_counts_per_bin: float
    z_score: np.ndarray


@dataclass(frozen=True)
class ExcessInterval:
    run_key: str
    position_m: float
    tot_min_ns: float
    tot_max_ns: float
    start_phase_ns: float
    end_phase_ns: float
    peak_phase_ns: float
    peak_rate_hz: float
    peak_z: float
    integrated_excess_counts: float
    integrated_excess_rate_hz: float


@dataclass(frozen=True)
class AcceptedExposure:
    run_key: str
    covered_duration_seconds: float
    accepted_window_count: int
    total_window_count: int
    validation_efficiency_window: float
    rejected_window_fraction: float
    assumed_trigger_window_width_ns: float
    accepted_exposure_seconds: float
    deduced_period_ns: float
    accepted_rf_cycles: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build normalized one-period rate-vs-phase diagnostics by ToT slice, "
            "for robust fast-neutron/gamma selection design."
        )
    )
    parser.add_argument("--results-root", default="/data6/Dogma_analysis_by_Dachi/Results")
    parser.add_argument("--output-dir")
    parser.add_argument("--channel", type=int, default=2)
    parser.add_argument("--run-glob", default="NCAL_20us_Pos_5m_*")
    parser.add_argument("--gamma-tot-min-ns", type=float, default=10.0)
    parser.add_argument("--gamma-tot-max-ns", type=float, default=18.0)
    parser.add_argument("--tot-slice-edges-ns", type=float, nargs="*", default=[8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 22.0, 26.0, 30.0])
    parser.add_argument("--gamma-exclusion-half-width-ns", type=float, default=4.0)
    parser.add_argument("--gamma-display-tot-max-ns", type=float, default=30.0)
    parser.add_argument("--baseline-quantile", type=float, default=0.20)
    parser.add_argument("--z-threshold", type=float, default=5.0)
    parser.add_argument("--min-contiguous-bins", type=int, default=3)
    return record_arguments(parser.parse_args())


def parse_position_m(run_key: str) -> float:
    if "Pos_" not in run_key:
        return float("nan")
    return float(run_key.split("Pos_", 1)[1].split("m_", 1)[0])


def read_key_value_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip().strip('"')
    return values


def run_duration_seconds(histogram: SparseHistogram, run_key: str) -> float:
    run_dir = histogram.path.parents[1]
    summary = read_key_value_file(run_dir / f"{run_key}_cleaned_rates_summary.txt")
    return float(summary["covered_duration_seconds"])


def accepted_exposure(histogram: SparseHistogram, run_key: str, assumed_trigger_window_width_ns: float) -> AcceptedExposure:
    run_dir = histogram.path.parents[1]
    rf_summary = read_key_value_file(run_dir / f"{run_key}_cleaned_rf_period_summary.txt")
    rates_summary = read_key_value_file(run_dir / f"{run_key}_cleaned_rates_summary.txt")
    covered_duration_seconds = float(rates_summary["covered_duration_seconds"])
    accepted_window_count = int(rf_summary["valid_ch0_windows"])
    total_window_count = int(rf_summary["windows_seen"])
    validation_efficiency_window = accepted_window_count / total_window_count if total_window_count > 0 else float("nan")
    rejected_window_fraction = 1.0 - validation_efficiency_window if total_window_count > 0 else float("nan")
    accepted_exposure_seconds = accepted_window_count * assumed_trigger_window_width_ns * 1e-9
    deduced_period_ns = float(rf_summary["deduced_period_ns"])
    accepted_rf_cycles = accepted_window_count * (assumed_trigger_window_width_ns / deduced_period_ns)
    return AcceptedExposure(
        run_key=run_key,
        covered_duration_seconds=covered_duration_seconds,
        accepted_window_count=accepted_window_count,
        total_window_count=total_window_count,
        validation_efficiency_window=validation_efficiency_window,
        rejected_window_fraction=rejected_window_fraction,
        assumed_trigger_window_width_ns=assumed_trigger_window_width_ns,
        accepted_exposure_seconds=accepted_exposure_seconds,
        deduced_period_ns=deduced_period_ns,
        accepted_rf_cycles=accepted_rf_cycles,
    )


def signed_delta_ns(values_ns: np.ndarray, center_ns: float, period_ns: float) -> np.ndarray:
    return (values_ns - center_ns + 0.5 * period_ns) % period_ns - 0.5 * period_ns


def detect_gamma_center(histogram: SparseHistogram, args: argparse.Namespace) -> tuple[float, np.ndarray, np.ndarray, float]:
    period_ns = float(histogram.metadata["deduced_period_ns"])
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    gamma_tot_mask = mask_between(histogram.y_centers_ns, args.gamma_tot_min_ns, args.gamma_tot_max_ns)
    projection_3x = histogram.counts[gamma_tot_mask, :].sum(axis=0)
    if projection_3x.sum() <= 0:
        raise ValueError(f"No reference entries in the gamma ToT window for {histogram.path}")
    projection_1x = fold_projection_to_one_period(projection_3x, period_bins)
    smoothed = circular_smooth(projection_1x, half_width_bins=2)
    gamma_center_ns = float(phase_centers_ns[int(np.argmax(smoothed))])
    gamma_relative_phase_ns = signed_delta_ns(phase_centers_ns, gamma_center_ns, period_ns)
    order = np.argsort(gamma_relative_phase_ns)
    return gamma_center_ns, gamma_relative_phase_ns[order], projection_1x[order], period_ns


def draw_gamma_reference(
    histogram: SparseHistogram,
    run_key: str,
    position_m: float,
    duration_seconds: float,
    gamma_center_ns: float,
    sorted_relative_phase_ns: np.ndarray,
    sorted_gamma_counts: np.ndarray,
    period_ns: float,
    output_path: Path,
    profile_path: Path,
    args: argparse.Namespace,
) -> None:
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    smoothed_counts = circular_smooth(sorted_gamma_counts, half_width_bins=2)
    sorted_gamma_rate = sorted_gamma_counts / duration_seconds
    smoothed_rate = smoothed_counts / duration_seconds

    with profile_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["run_key", "phase_relative_ns", "gamma_band_counts", "gamma_band_rate_hz", "smoothed_gamma_band_rate_hz"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for phase_ns, counts, rate_hz, smooth_rate_hz in zip(sorted_relative_phase_ns, sorted_gamma_counts, sorted_gamma_rate, smoothed_rate):
            writer.writerow(
                {
                    "run_key": run_key,
                    "phase_relative_ns": phase_ns,
                    "gamma_band_counts": counts,
                    "gamma_band_rate_hz": rate_hz,
                    "smoothed_gamma_band_rate_hz": smooth_rate_hz,
                }
            )

    display_tot_mask = histogram.y_centers_ns <= args.gamma_display_tot_max_ns
    display_counts = histogram.counts[display_tot_mask, :]
    vmax = float(np.quantile(display_counts[display_counts > 0.0], 0.995)) if np.any(display_counts > 0.0) else 1.0

    figure, axes = plt.subplots(1, 2, figsize=(13.0, 4.8), constrained_layout=True)
    image = axes[0].imshow(
        np.ma.masked_less_equal(display_counts, 0.0),
        aspect="auto",
        origin="lower",
        extent=(
            float(histogram.x_centers_ns[0] - 0.5 * histogram.x_bin_width_ns),
            float(histogram.x_centers_ns[-1] + 0.5 * histogram.x_bin_width_ns),
            0.0,
            args.gamma_display_tot_max_ns,
        ),
        interpolation="nearest",
        cmap="magma",
        vmin=0.0,
        vmax=max(1.0, vmax),
    )
    figure.colorbar(image, ax=axes[0], label="counts/bin, linear scale, clipped at q0.995")
    axes[0].axhspan(args.gamma_tot_min_ns, args.gamma_tot_max_ns, color="#00bcd4", alpha=0.14)
    for fold_index in range(int(histogram.metadata.get("fold_factor", "3"))):
        axes[0].axvline(fold_index * period_ns + gamma_center_ns, color="#00e5ff", linestyle="--", linewidth=1.2)
    axes[0].set_xlabel("RF phase, three periods (ns)")
    axes[0].set_ylabel("ToT (ns)")
    axes[0].set_title(f"{run_key} | gamma ToT band on 2D")

    axes[1].plot(sorted_relative_phase_ns, sorted_gamma_rate, color="#0891b2", linewidth=1.0, label="gamma-band rate")
    axes[1].plot(sorted_relative_phase_ns, smoothed_rate, color="#ef4444", linewidth=1.5, label="smoothed gamma-band rate")
    axes[1].axvline(0.0, color="#00e5ff", linestyle="--", linewidth=1.2, label="selected gamma reference")
    axes[1].set_xlabel("RF phase relative to selected gamma reference (ns)")
    axes[1].set_ylabel("rate per phase bin (Hz)")
    axes[1].set_title(f"{run_key} | one-period gamma-band profile")
    axes[1].grid(True, linewidth=0.3, alpha=0.5)
    axes[1].legend(fontsize=8)

    figure.suptitle(f"Prompt-gamma reference definition | detector face {position_m:g} m")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def build_slice_profile(
    histogram: SparseHistogram,
    run_key: str,
    position_m: float,
    duration_seconds: float,
    period_ns: float,
    gamma_center_ns: float,
    sorted_relative_phase_ns: np.ndarray,
    tot_min_ns: float,
    tot_max_ns: float,
    args: argparse.Namespace,
) -> SliceProfile:
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    tot_mask = mask_between(histogram.y_centers_ns, tot_min_ns, tot_max_ns)
    projection_3x = histogram.counts[tot_mask, :].sum(axis=0)
    projection_1x = fold_projection_to_one_period(projection_3x, period_bins)
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    relative_phase_ns = signed_delta_ns(phase_centers_ns, gamma_center_ns, period_ns)
    order = np.argsort(relative_phase_ns)
    sorted_counts = projection_1x[order]
    sorted_rate = sorted_counts / duration_seconds

    baseline_candidates = sorted_counts[np.abs(sorted_relative_phase_ns) > args.gamma_exclusion_half_width_ns]
    if baseline_candidates.size == 0:
        baseline_counts = float(np.quantile(sorted_counts, args.baseline_quantile))
    else:
        baseline_counts = float(np.quantile(baseline_candidates, args.baseline_quantile))
    z_score = (sorted_counts - baseline_counts) / np.sqrt(np.maximum(baseline_counts, 1.0))

    return SliceProfile(
        run_key=run_key,
        position_m=position_m,
        tot_min_ns=tot_min_ns,
        tot_max_ns=tot_max_ns,
        phase_center_ns=gamma_center_ns,
        gamma_relative_phase_ns=sorted_relative_phase_ns,
        rate_hz=sorted_rate,
        counts=sorted_counts,
        baseline_rate_hz=baseline_counts / duration_seconds,
        baseline_counts_per_bin=baseline_counts,
        z_score=z_score,
    )


def find_intervals(profile: SliceProfile, args: argparse.Namespace) -> list[ExcessInterval]:
    selected = (np.abs(profile.gamma_relative_phase_ns) > args.gamma_exclusion_half_width_ns) & (profile.z_score >= args.z_threshold)
    intervals: list[ExcessInterval] = []
    start_index: int | None = None
    for index, is_selected in enumerate(selected):
        if is_selected and start_index is None:
            start_index = index
        if start_index is not None and (not is_selected or index == len(selected) - 1):
            end_index = index if not is_selected else index + 1
            if end_index - start_index >= args.min_contiguous_bins:
                interval_slice = slice(start_index, end_index)
                rates = profile.rate_hz[interval_slice]
                counts = profile.counts[interval_slice]
                phases = profile.gamma_relative_phase_ns[interval_slice]
                z = profile.z_score[interval_slice]
                peak_local = int(np.argmax(rates))
                integrated_excess_counts = float(np.sum(np.maximum(0.0, counts - profile.baseline_counts_per_bin)))
                intervals.append(
                    ExcessInterval(
                        run_key=profile.run_key,
                        position_m=profile.position_m,
                        tot_min_ns=profile.tot_min_ns,
                        tot_max_ns=profile.tot_max_ns,
                        start_phase_ns=float(phases[0]),
                        end_phase_ns=float(phases[-1]),
                        peak_phase_ns=float(phases[peak_local]),
                        peak_rate_hz=float(rates[peak_local]),
                        peak_z=float(z[peak_local]),
                        integrated_excess_counts=integrated_excess_counts,
                        integrated_excess_rate_hz=integrated_excess_counts / max(1.0, np.sum(profile.counts) / np.sum(profile.rate_hz)) if np.sum(profile.rate_hz) > 0 else float("nan"),
                    )
                )
            start_index = None
    return intervals


def write_phase_profiles(profiles: list[SliceProfile], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "run_key",
                "position_m",
                "tot_min_ns",
                "tot_max_ns",
                "phase_relative_ns",
                "counts",
                "rate_hz",
                "baseline_rate_hz",
                "z_score",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for profile in profiles:
            for phase_ns, counts, rate_hz, z_score in zip(profile.gamma_relative_phase_ns, profile.counts, profile.rate_hz, profile.z_score):
                writer.writerow(
                    {
                        "run_key": profile.run_key,
                        "position_m": profile.position_m,
                        "tot_min_ns": profile.tot_min_ns,
                        "tot_max_ns": profile.tot_max_ns,
                        "phase_relative_ns": phase_ns,
                        "counts": counts,
                        "rate_hz": rate_hz,
                        "baseline_rate_hz": profile.baseline_rate_hz,
                        "z_score": z_score,
                    }
                )


def write_intervals(intervals: list[ExcessInterval], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[field.name for field in ExcessInterval.__dataclass_fields__.values()], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for interval in intervals:
            writer.writerow({field: getattr(interval, field) for field in ExcessInterval.__dataclass_fields__})


def draw_profiles(run_key: str, position_m: float, profiles: list[SliceProfile], output_path: Path, gamma_exclusion_half_width_ns: float) -> None:
    figure, axes = plt.subplots(2, 1, figsize=(11.0, 8.0), sharex=True, constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(profiles)))
    for color, profile in zip(colors, profiles):
        label = f"{profile.tot_min_ns:g}-{profile.tot_max_ns:g} ns"
        axes[0].plot(profile.gamma_relative_phase_ns, profile.rate_hz, color=color, linewidth=1.2, label=label)
        axes[0].hlines(profile.baseline_rate_hz, profile.gamma_relative_phase_ns[0], profile.gamma_relative_phase_ns[-1], color=color, linewidth=0.8, linestyle="--", alpha=0.6)
        axes[1].plot(profile.gamma_relative_phase_ns, profile.z_score, color=color, linewidth=1.2, label=label)

    axes[0].axvspan(-gamma_exclusion_half_width_ns, gamma_exclusion_half_width_ns, color="#00bcd4", alpha=0.08)
    axes[1].axvspan(-gamma_exclusion_half_width_ns, gamma_exclusion_half_width_ns, color="#00bcd4", alpha=0.08)
    axes[1].axhline(5.0, color="#d97706", linestyle="--", linewidth=1.0, label="z = 5")

    axes[0].set_ylabel("rate per phase bin (Hz)")
    axes[0].set_title(f"{run_key} | detector face {position_m:g} m | Ch02 phase-rate slices")
    axes[0].grid(True, linewidth=0.3, alpha=0.5)
    axes[0].legend(fontsize=8, ncol=2)

    axes[1].set_xlabel("RF phase relative to prompt gamma (ns)")
    axes[1].set_ylabel("excess significance z")
    axes[1].grid(True, linewidth=0.3, alpha=0.5)

    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_readme(output_dir: Path, args: argparse.Namespace) -> None:
    text = f"""# Step07 Phase-Rate Slice Diagnostics

This directory contains gamma-relative one-period rate-vs-phase diagnostics by ToT slice.

Goal:
- normalize each phase profile to rate per unit time;
- inspect abrupt excess regions in 1D without committing yet to a full 2D ROI;
- derive mathematically auditable candidate intervals from excess-over-baseline significance.

Method:
- Prompt gamma is centered at phase 0 in the displayed one-period coordinate.
- For each ToT slice, a folded one-period profile is converted to rate (Hz per phase bin).
- A robust baseline is estimated from the low quantile of bins outside the prompt-gamma exclusion region.
- Excess significance is approximated by z = (counts - baseline_counts) / sqrt(max(baseline_counts, 1)).
- Candidate intervals are contiguous bins with z >= {args.z_threshold:g} for at least {args.min_contiguous_bins} bins.

Important note:
- This is a selection-design diagnostic for the 5 m separated case. It is more rigorous than choosing boxes by eye, but it still needs 2D confirmation before defining final efficiency cuts.

Gamma reference definition:
- The prompt-gamma reference phase is derived from the compact gamma ToT band [{args.gamma_tot_min_ns}, {args.gamma_tot_max_ns}) ns.
- Counts in that ToT band are folded to one RF period.
- A small circular smoothing is applied.
- The selected gamma reference is the phase bin at the maximum of the smoothed one-period gamma-band profile.
- All phase-rate slice profiles are then displayed relative to that selected prompt-gamma phase.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    results_root = Path(args.results_root).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else results_root / "Ch02_Position_Comparison" / "phase_rate_slice_diagnostics"
    output_dir.mkdir(parents=True, exist_ok=True)

    sparse_files = find_sparse_files(results_root, args.run_glob)
    if not sparse_files:
        raise FileNotFoundError(f"No Folded_RF_3x phase-ToT sparse files found under {results_root}")

    edges = list(args.tot_slice_edges_ns)
    if len(edges) < 2:
        raise ValueError("Need at least two ToT slice edges")

    all_profiles: list[SliceProfile] = []
    all_intervals: list[ExcessInterval] = []
    for sparse_path in sparse_files:
        histogram = read_channel_histogram(sparse_path, args.channel)
        run_key = histogram.metadata.get("run_key", parse_run_key(histogram.path))
        position_m = parse_position_m(run_key)
        duration_seconds = run_duration_seconds(histogram, run_key)
        gamma_center_ns, sorted_relative_phase_ns, _, period_ns = detect_gamma_center(histogram, args)
        gamma_profile_3x = histogram.counts[mask_between(histogram.y_centers_ns, args.gamma_tot_min_ns, args.gamma_tot_max_ns), :].sum(axis=0)
        gamma_profile_1x = fold_projection_to_one_period(gamma_profile_3x, max(1, round(period_ns / histogram.x_bin_width_ns)))
        phase_centers_ns = (np.arange(max(1, round(period_ns / histogram.x_bin_width_ns)), dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
        gamma_relative_phase_ns = signed_delta_ns(phase_centers_ns, gamma_center_ns, period_ns)
        gamma_order = np.argsort(gamma_relative_phase_ns)
        draw_gamma_reference(
            histogram,
            run_key,
            position_m,
            duration_seconds,
            gamma_center_ns,
            gamma_relative_phase_ns[gamma_order],
            gamma_profile_1x[gamma_order],
            period_ns,
            output_dir / f"{run_key}_ch{args.channel:02d}_gamma_reference.png",
            output_dir / f"{run_key}_ch{args.channel:02d}_gamma_reference_profile.tsv",
            args,
        )

        run_profiles: list[SliceProfile] = []
        for tot_min_ns, tot_max_ns in zip(edges[:-1], edges[1:]):
            profile = build_slice_profile(
                histogram,
                run_key,
                position_m,
                duration_seconds,
                period_ns,
                gamma_center_ns,
                sorted_relative_phase_ns,
                tot_min_ns,
                tot_max_ns,
                args,
            )
            run_profiles.append(profile)
            all_profiles.append(profile)
            all_intervals.extend(find_intervals(profile, args))

        draw_profiles(run_key, position_m, run_profiles, output_dir / f"{run_key}_ch{args.channel:02d}_phase_rate_slices.png", args.gamma_exclusion_half_width_ns)
        print(f"{run_key}: wrote normalized phase-rate slice diagnostics")

    write_phase_profiles(all_profiles, output_dir / "phase_rate_slices.tsv")
    write_intervals(all_intervals, output_dir / "phase_rate_excess_intervals.tsv")
    write_readme(output_dir, args)
    print(f"Wrote {output_dir}")


if __name__ == "__main__":
    main()