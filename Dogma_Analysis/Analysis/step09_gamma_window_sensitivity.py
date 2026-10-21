from __future__ import annotations
from analysis_provenance import record_arguments

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

from matplotlib import patches
import matplotlib.pyplot as plt
import numpy as np

from step05_estimate_gamma_neutron_background import SparseHistogram, find_sparse_files, mask_between, parse_run_key, read_channel_histogram
from step07_phase_rate_slice_diagnostics import (
    SliceProfile,
    accepted_exposure,
    draw_profiles,
    find_intervals,
    parse_position_m,
    run_duration_seconds,
    signed_delta_ns,
    write_intervals,
    write_phase_profiles,
)


@dataclass(frozen=True)
class GammaWindowResult:
    run_key: str
    position_m: float
    duration_seconds: float
    accepted_exposure_seconds: float
    accepted_window_count: int
    total_window_count: int
    validation_efficiency_window: float
    rejected_window_fraction: float
    assumed_trigger_window_width_ns: float
    accepted_rf_cycles: float
    peak_mu_ns: float
    half_width_ns: float
    gamma_tot_min_ns: float
    gamma_tot_max_ns: float
    phase_center_ns: float
    phase_center_shift_from_1p25_ns: float
    gamma_total_counts: float
    gamma_total_rate_hz: float
    neutron_tot_min_ns: float
    neutron_tot_max_ns: float
    gamma_integral_phase_start_ns: float
    gamma_integral_phase_end_ns: float
    gamma_integral_counts_per_accepted_cycle: float
    neutron_integral_phase_start_ns: float
    neutron_integral_phase_end_ns: float
    neutron_integral_counts_per_accepted_cycle: float


@dataclass(frozen=True)
class ZebraTimewalkFitResult:
    run_key: str
    period_ns: float
    low_tot_min_ns: float
    low_tot_max_ns: float
    high_tot_min_ns: float
    high_tot_max_ns: float
    fit_phase_min_ns: float
    fit_phase_max_ns: float
    fit_phase_half_width_ns: float
    fit_min_band_counts: float
    refine_tot_min_ns: float
    refine_tot_max_ns: float
    seed_slope_ns_per_ns: float
    slope_ns_per_ns: float
    corrected_artifact_phase_ns: float
    tot_reference_ns: float
    low_peak_tot_ns: float
    high_peak_tot_ns: float
    low_peak_phase_ns: float
    high_peak_phase_ns: float
    low_peak_counts: float
    high_peak_counts: float
    tracked_row_count: int
    tracked_phase_rms_ns: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Test the sensitivity of the prompt-gamma phase reference to fixed-width gamma ToT windows "
            "centered on the fitted gamma ToT peak."
        )
    )
    parser.add_argument("--results-root", default="/data6/Dogma_analysis_by_Dachi/Results")
    parser.add_argument("--gamma-fit-summary", default=None)
    parser.add_argument("--output-dir")
    parser.add_argument("--channel", type=int, default=2)
    parser.add_argument("--run-glob", default="NCAL_20us_Pos_5m_*")
    parser.add_argument("--broad-gamma-tot-min-ns", type=float, default=10.0)
    parser.add_argument("--broad-gamma-tot-max-ns", type=float, default=18.0)
    parser.add_argument("--half-widths-ns", type=float, nargs="*", default=[1.0, 1.25, 1.5])
    parser.add_argument("--tot-slice-edges-ns", type=float, nargs="*", default=[8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 22.0, 26.0, 30.0])
    parser.add_argument("--gamma-exclusion-half-width-ns", type=float, default=4.0)
    parser.add_argument("--gamma-display-tot-max-ns", type=float, default=30.0)
    parser.add_argument("--neutron-tot-min-ns", type=float, default=8.0)
    parser.add_argument("--neutron-tot-max-ns", type=float, default=float("inf"))
    parser.add_argument("--neutron-phase-search-min-ns", type=float, default=6.0)
    parser.add_argument("--integration-smoothing-bins", type=int, default=5)
    parser.add_argument("--baseline-quantile", type=float, default=0.20)
    parser.add_argument("--z-threshold", type=float, default=5.0)
    parser.add_argument("--min-contiguous-bins", type=int, default=3)
    parser.add_argument("--local-fit-half-window-bins", type=int, default=2)
    parser.add_argument("--apply-zebra-timewalk-correction", action="store_true")
    parser.add_argument("--zebra-low-tot-min-ns", type=float, default=1.0)
    parser.add_argument("--zebra-low-tot-max-ns", type=float, default=2.0)
    parser.add_argument("--zebra-high-tot-min-ns", type=float, default=7.0)
    parser.add_argument("--zebra-high-tot-max-ns", type=float, default=8.0)
    parser.add_argument("--zebra-fit-phase-min-ns", type=float, default=0.0)
    parser.add_argument("--zebra-fit-phase-max-ns", type=float, default=float("inf"))
    parser.add_argument("--zebra-fit-phase-half-width-ns", type=float, default=2.0)
    parser.add_argument("--zebra-min-band-counts", type=float, default=300.0)
    parser.add_argument("--zebra-refine-tot-min-ns", type=float, default=1.0)
    parser.add_argument("--zebra-refine-tot-max-ns", type=float, default=8.0)
    parser.add_argument("--zebra-refine-slope-span-ns-per-ns", type=float, default=0.75)
    parser.add_argument("--zebra-refine-slope-step-ns-per-ns", type=float, default=0.01)
    parser.add_argument(
        "--assumed-trigger-window-width-ns",
        type=float,
        default=10000.0,
        help="Working-model trigger window width used to convert accepted windows into accepted RF cycles.",
    )
    return record_arguments(parser.parse_args())


def read_peak_summary(path: Path) -> dict[str, dict[str, float]]:
    summary: dict[str, dict[str, float]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("fit_success", "").lower() != "true":
                raise ValueError(f"Unsuccessful step08 fit for {row['run_key']} in {path}")
            peak_mu = float(row["peak_mu_ns"])
            peak_sigma = float(row["peak_sigma_ns"])
            if not np.isfinite(peak_mu) or not np.isfinite(peak_sigma) or peak_sigma <= 0:
                raise ValueError(f"Invalid step08 peak parameters for {row['run_key']} in {path}")
            if not float(row["fit_tot_min_ns"]) <= peak_mu <= float(row["fit_tot_max_ns"]):
                raise ValueError(f"Step08 peak outside its fitted ToT ROI for {row['run_key']} in {path}")
            summary[row["run_key"]] = {
                "peak_mu_ns": peak_mu,
                "peak_mu_err_ns": float(row["peak_mu_err_ns"]),
                "peak_sigma_ns": peak_sigma,
            }
    return summary


def local_quadratic_peak_fit(phase: np.ndarray, values: np.ndarray, half_window_bins: int) -> float:
    peak_idx = int(np.argmax(values))
    left = max(0, peak_idx - half_window_bins)
    right = min(len(phase), peak_idx + half_window_bins + 1)
    window_phase = phase[left:right]
    window_values = values[left:right]
    if window_phase.size < 3:
        return float(phase[peak_idx])
    coeffs = np.polyfit(window_phase, window_values, deg=2)
    a, b, _c = coeffs
    if a >= 0.0:
        return float(phase[peak_idx])
    vertex = float(-b / (2.0 * a))
    if not (float(window_phase[0]) <= vertex <= float(window_phase[-1])):
        return float(phase[peak_idx])
    return vertex


def weighted_phase_center_absolute(
    phase_centers_ns: np.ndarray,
    counts: np.ndarray,
    center_ns: float,
    period_ns: float,
    half_width_ns: float,
) -> tuple[float, float]:
    relative_phase_ns = signed_delta_ns(phase_centers_ns, center_ns, period_ns)
    mask = np.abs(relative_phase_ns) <= half_width_ns
    selected_counts = counts[mask]
    total_counts = float(np.sum(selected_counts))
    if total_counts <= 0.0:
        return float("nan"), 0.0
    mean_relative_ns = float(np.sum(relative_phase_ns[mask] * selected_counts) / total_counts)
    return float((center_ns + mean_relative_ns) % period_ns), total_counts


def local_peak_phase_absolute(
    phase_centers_ns: np.ndarray,
    counts: np.ndarray,
    center_ns: float,
    period_ns: float,
    half_width_ns: float,
) -> tuple[float, float]:
    relative_phase_ns = signed_delta_ns(phase_centers_ns, center_ns, period_ns)
    mask = np.abs(relative_phase_ns) <= half_width_ns
    if not np.any(mask):
        return float("nan"), 0.0
    selected_phase_ns = relative_phase_ns[mask]
    selected_counts = counts[mask]
    total_counts = float(np.sum(selected_counts))
    if total_counts <= 0.0:
        return float("nan"), 0.0

    order = np.argsort(selected_phase_ns)
    ordered_phase_ns = selected_phase_ns[order]
    ordered_counts = selected_counts[order]
    smoothed_counts = circular_smooth(ordered_counts, kernel_bins=5)
    refined_relative_ns = local_quadratic_peak_fit(ordered_phase_ns, smoothed_counts, half_window_bins=2)
    return float((center_ns + refined_relative_ns) % period_ns), total_counts


def band_peak_phase(
    one_period_counts: np.ndarray,
    phase_centers_ns: np.ndarray,
    tot_centers_ns: np.ndarray,
    period_ns: float,
    tot_min_ns: float,
    tot_max_ns: float,
    phase_min_ns: float,
    phase_max_ns: float,
    phase_half_width_ns: float,
) -> tuple[float, float, float]:
    tot_mask = mask_between(tot_centers_ns, tot_min_ns, tot_max_ns)
    projection = np.sum(one_period_counts[tot_mask, :], axis=0)
    if np.isfinite(phase_max_ns) and phase_max_ns > phase_min_ns:
        phase_mask = mask_between(phase_centers_ns, phase_min_ns, phase_max_ns)
        phase_indices = np.where(phase_mask)[0]
    else:
        phase_indices = np.arange(len(phase_centers_ns), dtype=int)
    if phase_indices.size == 0:
        return float("nan"), 0.0, float("nan")
    selected_projection = projection[phase_indices]
    selected_total = float(np.sum(selected_projection))
    if selected_total <= 0.0:
        return float("nan"), 0.0, float("nan")
    peak_index = int(phase_indices[int(np.argmax(selected_projection))])
    peak_center_ns = float(phase_centers_ns[peak_index])
    peak_phase_ns, peak_counts = weighted_phase_center_absolute(
        phase_centers_ns,
        projection,
        peak_center_ns,
        period_ns,
        phase_half_width_ns,
    )
    band_rows = one_period_counts[tot_mask, :]
    band_tot_centers_ns = tot_centers_ns[tot_mask]
    band_relative_phase_ns = signed_delta_ns(phase_centers_ns, peak_phase_ns, period_ns)
    band_phase_mask = np.abs(band_relative_phase_ns) <= phase_half_width_ns
    band_phase_counts = np.sum(band_rows[:, band_phase_mask], axis=1)
    if np.sum(band_phase_counts) > 0.0:
        peak_tot_ns = float(np.average(band_tot_centers_ns, weights=band_phase_counts))
    else:
        peak_tot_ns = float(0.5 * (tot_min_ns + tot_max_ns))
    return peak_phase_ns, peak_counts, peak_tot_ns


def shift_one_period_histogram(
    one_period_counts: np.ndarray,
    tot_centers_ns: np.ndarray,
    phase_bin_width_ns: float,
    period_ns: float,
    slope_ns_per_ns: float,
    tot_reference_ns: float,
) -> np.ndarray:
    phase_bins = one_period_counts.shape[1]
    source_indices = np.arange(phase_bins, dtype=np.float64)
    corrected = np.zeros_like(one_period_counts, dtype=np.float64)
    for tot_index, tot_ns in enumerate(tot_centers_ns):
        row = one_period_counts[tot_index, :]
        if not np.any(row > 0.0):
            continue
        shift_bins = slope_ns_per_ns * (tot_ns - tot_reference_ns) / phase_bin_width_ns
        target = (source_indices - shift_bins) % phase_bins
        left = np.floor(target).astype(int)
        frac = target - left
        right = (left + 1) % phase_bins
        np.add.at(corrected[tot_index], left, row * (1.0 - frac))
        np.add.at(corrected[tot_index], right, row * frac)
    return corrected


def phase_model_value(
    tot_ns: float,
    tracked_tot_ns: np.ndarray,
    tracked_phase_ns: np.ndarray,
    slope_ns_per_ns: float,
) -> float:
    if tracked_tot_ns.size == 0:
        return float("nan")
    if tracked_tot_ns.size == 1:
        return float(tracked_phase_ns[0] + slope_ns_per_ns * (tot_ns - tracked_tot_ns[0]))
    if tot_ns <= tracked_tot_ns[0]:
        return float(tracked_phase_ns[0] + slope_ns_per_ns * (tot_ns - tracked_tot_ns[0]))
    if tot_ns >= tracked_tot_ns[-1]:
        return float(tracked_phase_ns[-1] + slope_ns_per_ns * (tot_ns - tracked_tot_ns[-1]))
    return float(np.interp(tot_ns, tracked_tot_ns, tracked_phase_ns))


def shift_one_period_histogram_with_phase_model(
    one_period_counts: np.ndarray,
    tot_centers_ns: np.ndarray,
    phase_bin_width_ns: float,
    target_phase_ns: float,
    tracked_tot_ns: np.ndarray,
    tracked_phase_ns: np.ndarray,
    slope_ns_per_ns: float,
) -> np.ndarray:
    phase_bins = one_period_counts.shape[1]
    source_indices = np.arange(phase_bins, dtype=np.float64)
    corrected = np.zeros_like(one_period_counts, dtype=np.float64)
    for tot_index, tot_ns in enumerate(tot_centers_ns):
        row = one_period_counts[tot_index, :]
        if not np.any(row > 0.0):
            continue
        source_phase_ns = phase_model_value(tot_ns, tracked_tot_ns, tracked_phase_ns, slope_ns_per_ns)
        shift_bins = (source_phase_ns - target_phase_ns) / phase_bin_width_ns
        target = (source_indices - shift_bins) % phase_bins
        left = np.floor(target).astype(int)
        frac = target - left
        right = (left + 1) % phase_bins
        np.add.at(corrected[tot_index], left, row * (1.0 - frac))
        np.add.at(corrected[tot_index], right, row * frac)
    return corrected


def rotate_tiled_histogram_for_display(
    one_period_counts: np.ndarray,
    tot_centers_ns: np.ndarray,
    phase_bin_width_ns: float,
    tot_bin_width_ns: float,
    period_ns: float,
    anchor_phase_ns: float,
    tot_reference_ns: float,
    slope_ns_per_ns: float,
) -> np.ndarray:
    tiled_counts = np.tile(one_period_counts, (1, 3))
    rotated = np.zeros_like(tiled_counts, dtype=np.float64)
    x_centers_ns = (np.arange(tiled_counts.shape[1], dtype=np.float64) + 0.5) * phase_bin_width_ns
    x_min_center_ns = float(x_centers_ns[0])
    y_min_center_ns = float(tot_centers_ns[0])
    theta = float(np.arctan(slope_ns_per_ns))
    cos_theta = float(np.cos(theta))
    sin_theta = float(np.sin(theta))

    for tot_index, tot_ns in enumerate(tot_centers_ns):
        row = tiled_counts[tot_index, :]
        nonzero_indices = np.where(row > 0.0)[0]
        if nonzero_indices.size == 0:
            continue
        y_rel_ns = tot_ns - tot_reference_ns
        for phase_index in nonzero_indices:
            count = float(row[phase_index])
            phase_ns = float(x_centers_ns[phase_index])
            fold_index = int(np.floor(phase_ns / period_ns))
            anchor_phase_abs_ns = anchor_phase_ns + fold_index * period_ns
            x_rel_ns = phase_ns - anchor_phase_abs_ns
            rotated_phase_ns = anchor_phase_abs_ns + cos_theta * x_rel_ns - sin_theta * y_rel_ns
            rotated_tot_ns = tot_reference_ns + sin_theta * x_rel_ns + cos_theta * y_rel_ns
            x_float = (rotated_phase_ns - x_min_center_ns) / phase_bin_width_ns
            y_float = (rotated_tot_ns - y_min_center_ns) / tot_bin_width_ns
            if not (-1.0 < x_float < tiled_counts.shape[1] and -1.0 < y_float < tiled_counts.shape[0]):
                continue
            x_left = int(np.floor(x_float))
            y_low = int(np.floor(y_float))
            x_frac = x_float - x_left
            y_frac = y_float - y_low
            for x_target, x_weight in ((x_left, 1.0 - x_frac), (x_left + 1, x_frac)):
                if x_target < 0 or x_target >= tiled_counts.shape[1]:
                    continue
                for y_target, y_weight in ((y_low, 1.0 - y_frac), (y_low + 1, y_frac)):
                    if y_target < 0 or y_target >= tiled_counts.shape[0]:
                        continue
                    rotated[y_target, x_target] += count * x_weight * y_weight
    return rotated


def shift_one_period_phase_counts(one_period_counts: np.ndarray, shift_bins: float) -> np.ndarray:
    phase_bins = one_period_counts.shape[1]
    source_indices = np.arange(phase_bins, dtype=np.float64)
    shifted = np.zeros_like(one_period_counts, dtype=np.float64)
    for tot_index in range(one_period_counts.shape[0]):
        row = one_period_counts[tot_index, :]
        if not np.any(row > 0.0):
            continue
        target = (source_indices + shift_bins) % phase_bins
        left = np.floor(target).astype(int)
        frac = target - left
        right = (left + 1) % phase_bins
        np.add.at(shifted[tot_index], left, row * (1.0 - frac))
        np.add.at(shifted[tot_index], right, row * frac)
    return shifted


def shift_phase_value_for_display(phase_ns: float, phase_offset_ns: float, period_ns: float) -> float:
    return float((phase_ns + phase_offset_ns) % period_ns)


def shift_segments_for_display(segments: list[tuple[float, float]], phase_offset_ns: float, period_ns: float) -> list[tuple[float, float]]:
    shifted_segments: list[tuple[float, float]] = []
    for start_ns, end_ns in segments:
        shifted_start_ns = shift_phase_value_for_display(start_ns, phase_offset_ns, period_ns)
        shifted_end_ns = shift_phase_value_for_display(end_ns, phase_offset_ns, period_ns)
        if shifted_end_ns > shifted_start_ns:
            shifted_segments.append((shifted_start_ns, shifted_end_ns))
        else:
            shifted_segments.append((0.0, shifted_end_ns))
            shifted_segments.append((shifted_start_ns, period_ns))
    return shifted_segments


def refine_zebra_slope(
    one_period_counts: np.ndarray,
    phase_centers_ns: np.ndarray,
    tot_centers_ns: np.ndarray,
    period_ns: float,
    seed_intercept_ns: float,
    seed_slope_ns_per_ns: float,
    tot_reference_ns: float,
    args: argparse.Namespace,
) -> tuple[float, float, int, float, np.ndarray, np.ndarray, np.ndarray]:
    refine_tot_mask = mask_between(tot_centers_ns, args.zebra_refine_tot_min_ns, args.zebra_refine_tot_max_ns)
    refine_indices = np.where(refine_tot_mask)[0]
    if refine_indices.size < 2:
        empty = np.asarray([], dtype=np.float64)
        return seed_intercept_ns, seed_slope_ns_per_ns, 0, float("nan"), empty, empty, empty

    window_half_width_ns = max(args.zebra_fit_phase_half_width_ns, 1.5)
    min_row_counts = max(20.0, 0.02 * args.zebra_min_band_counts)
    current_intercept_ns = seed_intercept_ns
    current_slope = seed_slope_ns_per_ns
    tracked_tot_ns = np.asarray([], dtype=np.float64)
    tracked_phase_ns = np.asarray([], dtype=np.float64)
    tracked_weights = np.asarray([], dtype=np.float64)
    coeffs = np.asarray([current_slope, current_intercept_ns], dtype=np.float64)

    for _ in range(3):
        tot_samples: list[float] = []
        phase_samples: list[float] = []
        weight_samples: list[float] = []
        for tot_index in refine_indices:
            row_counts = one_period_counts[tot_index, :]
            if not np.any(row_counts > 0.0):
                continue
            tot_ns = float(tot_centers_ns[tot_index])
            predicted_phase_unwrapped_ns = float(current_intercept_ns + current_slope * (tot_ns - tot_reference_ns))
            predicted_phase_ns = float(predicted_phase_unwrapped_ns % period_ns)
            row_phase_ns, row_window_counts = local_peak_phase_absolute(
                phase_centers_ns,
                row_counts,
                predicted_phase_ns,
                period_ns,
                window_half_width_ns,
            )
            if not np.isfinite(row_phase_ns) or row_window_counts < min_row_counts:
                continue
            unwrapped_phase_ns = predicted_phase_unwrapped_ns + float(signed_delta_ns(np.asarray([row_phase_ns]), predicted_phase_ns, period_ns)[0])
            tot_samples.append(tot_ns)
            phase_samples.append(unwrapped_phase_ns)
            weight_samples.append(float(row_window_counts))
        if len(tot_samples) < 2:
            empty = np.asarray([], dtype=np.float64)
            return current_intercept_ns, current_slope, len(tot_samples), float("nan"), empty, empty, empty
        tracked_tot_ns = np.asarray(tot_samples, dtype=np.float64)
        tracked_phase_ns = np.asarray(phase_samples, dtype=np.float64)
        tracked_weights = np.asarray(weight_samples, dtype=np.float64)
        coeffs = np.polyfit(
            tracked_tot_ns - tot_reference_ns,
            tracked_phase_ns,
            deg=1,
            w=np.sqrt(np.maximum(tracked_weights, 1.0)),
        )
        current_slope = float(coeffs[0])
        current_intercept_ns = float(coeffs[1])

    fitted_phase_ns = np.polyval(coeffs, tracked_tot_ns - tot_reference_ns)
    residual_rms_ns = float(np.sqrt(np.average((tracked_phase_ns - fitted_phase_ns) ** 2, weights=np.maximum(tracked_weights, 1.0))))
    return current_intercept_ns, current_slope, int(tracked_tot_ns.size), residual_rms_ns, tracked_tot_ns, tracked_phase_ns, tracked_weights


def fit_zebra_timewalk_model(histogram: SparseHistogram, run_key: str, args: argparse.Namespace) -> tuple[ZebraTimewalkFitResult, np.ndarray, np.ndarray, np.ndarray]:
    period_ns = float(histogram.metadata["deduced_period_ns"])
    one_period_counts = fold_histogram_to_one_period(histogram, period_ns)
    period_bins = one_period_counts.shape[1]
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    low_peak_phase_ns, low_peak_counts, low_peak_tot_ns = band_peak_phase(
        one_period_counts,
        phase_centers_ns,
        histogram.y_centers_ns,
        period_ns,
        args.zebra_low_tot_min_ns,
        args.zebra_low_tot_max_ns,
        args.zebra_fit_phase_min_ns,
        args.zebra_fit_phase_max_ns,
        args.zebra_fit_phase_half_width_ns,
    )
    high_peak_phase_ns, high_peak_counts, high_peak_tot_ns = band_peak_phase(
        one_period_counts,
        phase_centers_ns,
        histogram.y_centers_ns,
        period_ns,
        args.zebra_high_tot_min_ns,
        args.zebra_high_tot_max_ns,
        args.zebra_fit_phase_min_ns,
        args.zebra_fit_phase_max_ns,
        args.zebra_fit_phase_half_width_ns,
    )
    if not np.isfinite(low_peak_phase_ns) or not np.isfinite(high_peak_phase_ns):
        raise ValueError(f"{run_key}: unable to locate zebra-band peaks for timewalk correction")
    if low_peak_counts < args.zebra_min_band_counts or high_peak_counts < args.zebra_min_band_counts:
        raise ValueError(f"{run_key}: insufficient zebra-band counts for timewalk correction")

    delta_phase_short_ns = float(signed_delta_ns(np.asarray([high_peak_phase_ns]), low_peak_phase_ns, period_ns)[0])
    delta_tot_ns = high_peak_tot_ns - low_peak_tot_ns
    if abs(delta_tot_ns) <= 1.0e-12:
        raise ValueError(f"{run_key}: zebra-band ToT separation is zero")
    tot_reference_ns = 0.5 * (low_peak_tot_ns + high_peak_tot_ns)

    # Scan slopes from -3 to +3 ns/ns in steps of 0.25.  For each candidate slope,
    # track the zebra ridge with a FROZEN slope (only the intercept is updated per
    # iteration) so the scan explores the full range without the polyfit pulling every
    # seed to the same attractor.  The winner (most tracked rows; residual RMS as
    # tiebreak) becomes the seed for one final free-slope refinement.
    seed_slope_ns_per_ns = delta_phase_short_ns / delta_tot_ns  # stored for reference only
    scan_step = 0.10
    scan_half_range = 3.0
    scan_slopes = np.arange(-scan_half_range, scan_half_range + 0.5 * scan_step, scan_step)

    # Frozen-intercept-only tracking helper (avoids modifying refine_zebra_slope signature)
    def _track_frozen_slope(
        try_slope: float,
    ) -> tuple[float, int, float]:
        """Return (intercept, tracked_rows, residual_rms) for a fixed slope."""
        refine_tot_mask = mask_between(histogram.y_centers_ns, args.zebra_refine_tot_min_ns, args.zebra_refine_tot_max_ns)
        refine_indices = np.where(refine_tot_mask)[0]
        if refine_indices.size < 2:
            return float("nan"), 0, float("inf")
        window_half_width_ns = 3.5
        min_row_counts = max(20.0, 0.02 * args.zebra_min_band_counts)
        try_intercept = low_peak_phase_ns - try_slope * (low_peak_tot_ns - tot_reference_ns)
        current_intercept = try_intercept
        tot_s: list[float] = []
        phase_s: list[float] = []
        wt_s: list[float] = []
        for _iter in range(3):
            tot_s, phase_s, wt_s = [], [], []
            for tot_index in refine_indices:
                row_counts = one_period_counts[tot_index, :]
                if not np.any(row_counts > 0.0):
                    continue
                tot_ns_local = float(histogram.y_centers_ns[tot_index])
                pred_unwrap = float(current_intercept + try_slope * (tot_ns_local - tot_reference_ns))
                pred_wrap = float(pred_unwrap % period_ns)
                row_phase, row_wt = local_peak_phase_absolute(
                    phase_centers_ns, row_counts, pred_wrap, period_ns, window_half_width_ns,
                )
                if not np.isfinite(row_phase) or row_wt < min_row_counts:
                    continue
                unwrap = pred_unwrap + float(signed_delta_ns(np.asarray([row_phase]), pred_wrap, period_ns)[0])
                tot_s.append(tot_ns_local)
                phase_s.append(unwrap)
                wt_s.append(float(row_wt))
            if len(tot_s) < 2:
                return current_intercept, len(tot_s), float("inf")
            arr_tot = np.asarray(tot_s, dtype=np.float64)
            arr_phase = np.asarray(phase_s, dtype=np.float64)
            arr_wt = np.asarray(wt_s, dtype=np.float64)
            # Only update intercept; slope stays frozen
            intercept_samples = arr_phase - try_slope * (arr_tot - tot_reference_ns)
            current_intercept = float(np.average(intercept_samples, weights=np.maximum(arr_wt, 1.0)))
        if len(tot_s) < 2:
            return current_intercept, 0, float("inf")
        arr_tot = np.asarray(tot_s, dtype=np.float64)
        arr_phase = np.asarray(phase_s, dtype=np.float64)
        arr_wt = np.asarray(wt_s, dtype=np.float64)
        fitted = current_intercept + try_slope * (arr_tot - tot_reference_ns)
        rms = float(np.sqrt(np.average((arr_phase - fitted) ** 2, weights=np.maximum(arr_wt, 1.0))))
        return current_intercept, len(tot_s), rms

    best_scan_slope = seed_slope_ns_per_ns
    best_scan_rows = -1
    best_scan_rms = float("inf")
    for try_slope_val in scan_slopes:
        _ic, _rows, _rms = _track_frozen_slope(float(try_slope_val))
        rms_ok = np.isfinite(_rms)
        better = _rows > best_scan_rows or (
            _rows == best_scan_rows and rms_ok and _rms < best_scan_rms
        )
        if better:
            best_scan_slope = float(try_slope_val)
            best_scan_rows = _rows
            best_scan_rms = _rms if rms_ok else float("inf")

    # One final free-slope refinement starting from the scan winner
    best_intercept_scan = low_peak_phase_ns - best_scan_slope * (low_peak_tot_ns - tot_reference_ns)
    corrected_intercept_ns, slope_ns_per_ns, tracked_row_count, tracked_phase_rms_ns, tracked_tot_ns, tracked_phase_ns, tracked_weights = refine_zebra_slope(
        one_period_counts,
        phase_centers_ns,
        histogram.y_centers_ns,
        period_ns,
        best_intercept_scan,
        best_scan_slope,
        tot_reference_ns,
        args,
    )
    corrected_artifact_phase_ns = float(corrected_intercept_ns % period_ns)
    fit = ZebraTimewalkFitResult(
        run_key=run_key,
        period_ns=period_ns,
        low_tot_min_ns=args.zebra_low_tot_min_ns,
        low_tot_max_ns=args.zebra_low_tot_max_ns,
        high_tot_min_ns=args.zebra_high_tot_min_ns,
        high_tot_max_ns=args.zebra_high_tot_max_ns,
        fit_phase_min_ns=args.zebra_fit_phase_min_ns,
        fit_phase_max_ns=args.zebra_fit_phase_max_ns,
        fit_phase_half_width_ns=args.zebra_fit_phase_half_width_ns,
        fit_min_band_counts=args.zebra_min_band_counts,
        refine_tot_min_ns=args.zebra_refine_tot_min_ns,
        refine_tot_max_ns=args.zebra_refine_tot_max_ns,
        seed_slope_ns_per_ns=seed_slope_ns_per_ns,
        slope_ns_per_ns=slope_ns_per_ns,
        corrected_artifact_phase_ns=corrected_artifact_phase_ns,
        tot_reference_ns=tot_reference_ns,
        low_peak_tot_ns=low_peak_tot_ns,
        high_peak_tot_ns=high_peak_tot_ns,
        low_peak_phase_ns=low_peak_phase_ns,
        high_peak_phase_ns=high_peak_phase_ns,
        low_peak_counts=low_peak_counts,
        high_peak_counts=high_peak_counts,
        tracked_row_count=tracked_row_count,
        tracked_phase_rms_ns=tracked_phase_rms_ns,
    )
    return fit, tracked_tot_ns, tracked_phase_ns, tracked_weights


def apply_zebra_timewalk_correction(
    histogram: SparseHistogram,
    fit: ZebraTimewalkFitResult,
    tracked_tot_ns: np.ndarray,
    tracked_phase_ns: np.ndarray,
) -> SparseHistogram:
    one_period_counts = fold_histogram_to_one_period(histogram, fit.period_ns)
    corrected_one_period_counts = shift_one_period_histogram(
        one_period_counts,
        histogram.y_centers_ns,
        histogram.x_bin_width_ns,
        fit.period_ns,
        fit.slope_ns_per_ns,
        fit.tot_reference_ns,
    )

    # Keep the corrected synthetic histogram in one-period form. The existing
    # efficiency plotting path already folds to one period and tiles 3x only for display.
    corrected_counts = corrected_one_period_counts
    corrected_metadata = dict(histogram.metadata)
    corrected_metadata["timewalk_correction_mode"] = "zebra_linear"
    corrected_metadata["timewalk_zebra_slope_ns_per_ns"] = f"{fit.slope_ns_per_ns:.12f}"
    corrected_metadata["timewalk_zebra_corrected_artifact_phase_ns"] = f"{fit.corrected_artifact_phase_ns:.12f}"
    corrected_metadata["timewalk_zebra_tot_reference_ns"] = f"{fit.tot_reference_ns:.12f}"
    corrected_metadata["timewalk_zebra_low_tot_min_ns"] = f"{fit.low_tot_min_ns:.6f}"
    corrected_metadata["timewalk_zebra_low_tot_max_ns"] = f"{fit.low_tot_max_ns:.6f}"
    corrected_metadata["timewalk_zebra_high_tot_min_ns"] = f"{fit.high_tot_min_ns:.6f}"
    corrected_metadata["timewalk_zebra_high_tot_max_ns"] = f"{fit.high_tot_max_ns:.6f}"
    x_centers_ns = (np.arange(corrected_counts.shape[1], dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    return SparseHistogram(
        path=histogram.path,
        metadata=corrected_metadata,
        counts=corrected_counts,
        x_centers_ns=x_centers_ns,
        y_centers_ns=histogram.y_centers_ns.copy(),
        x_bin_width_ns=histogram.x_bin_width_ns,
        y_bin_width_ns=histogram.y_bin_width_ns,
    )


def write_timewalk_fit_summary(fits: list[ZebraTimewalkFitResult], output_path: Path) -> None:
    if not fits:
        return
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[field.name for field in ZebraTimewalkFitResult.__dataclass_fields__.values()], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for fit in fits:
            writer.writerow({field: getattr(fit, field) for field in ZebraTimewalkFitResult.__dataclass_fields__})


def detect_gamma_center_for_window(histogram, gamma_tot_min_ns: float, gamma_tot_max_ns: float, local_fit_half_window_bins: int):
    period_ns = float(histogram.metadata["deduced_period_ns"])
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    gamma_tot_mask = mask_between(histogram.y_centers_ns, gamma_tot_min_ns, gamma_tot_max_ns)
    projection_3x = histogram.counts[gamma_tot_mask, :].sum(axis=0)
    if projection_3x.sum() <= 0:
        raise ValueError(
            f"No reference entries in ToT window [{gamma_tot_min_ns}, {gamma_tot_max_ns}) ns "
            f"for {histogram.path}; cannot determine a phase reference."
        )
    projection_1x = np.zeros(period_bins, dtype=np.float64)
    for index, count in enumerate(projection_3x):
        projection_1x[index % period_bins] += count
    smoothed = np.convolve(np.r_[projection_1x[-2:], projection_1x, projection_1x[:2]], np.ones(5) / 5.0, mode="same")[2:-2]
    seed_center_ns = float(phase_centers_ns[int(np.argmax(smoothed))])
    seed_relative_phase_ns = signed_delta_ns(phase_centers_ns, seed_center_ns, period_ns)
    seed_order = np.argsort(seed_relative_phase_ns)
    local_offset_ns = local_quadratic_peak_fit(seed_relative_phase_ns[seed_order], projection_1x[seed_order], local_fit_half_window_bins)
    gamma_center_ns = float((seed_center_ns + local_offset_ns) % period_ns)
    gamma_relative_phase_ns = signed_delta_ns(phase_centers_ns, gamma_center_ns, period_ns)
    order = np.argsort(gamma_relative_phase_ns)
    return gamma_center_ns, gamma_relative_phase_ns[order], projection_1x[order], period_ns


def build_slice_profile(histogram, run_key, position_m, duration_seconds, period_ns, gamma_center_ns, sorted_relative_phase_ns, tot_min_ns, tot_max_ns, args) -> SliceProfile:
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    phase_bin_width_ns = period_ns / period_bins
    tot_mask = mask_between(histogram.y_centers_ns, tot_min_ns, tot_max_ns)
    projection_3x = histogram.counts[tot_mask, :].sum(axis=0)
    projection_1x = np.zeros(period_bins, dtype=np.float64)
    for index, count in enumerate(projection_3x):
        projection_1x[index % period_bins] += count
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    relative_phase_ns = signed_delta_ns(phase_centers_ns, gamma_center_ns, period_ns)
    order = np.argsort(relative_phase_ns)
    sorted_counts = projection_1x[order]
    sorted_rate = sorted_counts / duration_seconds / phase_bin_width_ns
    baseline_candidates = sorted_counts[np.abs(sorted_relative_phase_ns) > args.gamma_exclusion_half_width_ns]
    baseline_counts = float(np.quantile(baseline_candidates if baseline_candidates.size else sorted_counts, args.baseline_quantile))
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
        baseline_rate_hz=baseline_counts / duration_seconds / phase_bin_width_ns,
        baseline_counts_per_bin=baseline_counts,
        z_score=z_score,
    )


def build_phase_profile(histogram, duration_seconds, period_ns, gamma_center_ns, tot_min_ns, tot_max_ns):
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    phase_bin_width_ns = period_ns / period_bins
    tot_mask = mask_between(histogram.y_centers_ns, tot_min_ns, tot_max_ns)
    projection_3x = histogram.counts[tot_mask, :].sum(axis=0)
    projection_1x = np.zeros(period_bins, dtype=np.float64)
    for index, count in enumerate(projection_3x):
        projection_1x[index % period_bins] += count
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    relative_phase_ns = signed_delta_ns(phase_centers_ns, gamma_center_ns, period_ns)
    order = np.argsort(relative_phase_ns)
    return relative_phase_ns[order], projection_1x[order], projection_1x[order] / duration_seconds / phase_bin_width_ns, phase_bin_width_ns


def fold_histogram_to_one_period(histogram, period_ns: float) -> np.ndarray:
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    folded = np.zeros((histogram.counts.shape[0], period_bins), dtype=np.float64)
    for x_bin in range(histogram.counts.shape[1]):
        folded[:, x_bin % period_bins] += histogram.counts[:, x_bin]
    return folded


def circular_smooth(values: np.ndarray, kernel_bins: int) -> np.ndarray:
    if values.size == 0 or kernel_bins <= 1:
        return values.astype(np.float64, copy=True)
    kernel_bins = max(1, int(kernel_bins))
    if kernel_bins % 2 == 0:
        kernel_bins += 1
    pad = kernel_bins // 2
    kernel = np.ones(kernel_bins, dtype=np.float64) / kernel_bins
    padded = np.r_[values[-pad:], values, values[:pad]]
    return np.convolve(padded, kernel, mode="valid")


def local_minima_indices(values: np.ndarray) -> np.ndarray:
    if values.size < 3:
        return np.asarray([], dtype=int)
    return np.where((values[1:-1] <= values[:-2]) & (values[1:-1] <= values[2:]))[0] + 1


def wrap_phase_ns(phase_ns: float, phase_min_ns: float, period_ns: float) -> float:
    wrapped = phase_ns
    while wrapped < phase_min_ns:
        wrapped += period_ns
    while wrapped >= phase_min_ns + period_ns:
        wrapped -= period_ns
    return float(wrapped)


def mask_to_spans(phase: np.ndarray, mask: np.ndarray, phase_bin_width_ns: float) -> list[tuple[float, float]]:
    spans: list[tuple[float, float]] = []
    start_idx: int | None = None
    for idx, selected in enumerate(mask):
        if selected and start_idx is None:
            start_idx = idx
        elif not selected and start_idx is not None:
            spans.append((float(phase[start_idx] - 0.5 * phase_bin_width_ns), float(phase[idx - 1] + 0.5 * phase_bin_width_ns)))
            start_idx = None
    if start_idx is not None:
        spans.append((float(phase[start_idx] - 0.5 * phase_bin_width_ns), float(phase[-1] + 0.5 * phase_bin_width_ns)))
    return spans


def select_peak_interval_between_dips(
    phase: np.ndarray,
    density: np.ndarray,
    period_ns: float,
    phase_bin_width_ns: float,
    smoothing_bins: int,
    peak_selector,
) -> tuple[np.ndarray, float, float, float]:
    sample_count = len(phase)
    if sample_count == 0:
        empty_mask = np.zeros(0, dtype=bool)
        return empty_mask, float("nan"), float("nan"), float("nan")
    smoothed = circular_smooth(density, smoothing_bins)
    peak_index = int(peak_selector(smoothed))
    extended_phase = np.r_[phase - period_ns, phase, phase + period_ns]
    extended_smoothed = np.r_[smoothed, smoothed, smoothed]
    minima = local_minima_indices(extended_smoothed)
    center_peak_index = peak_index + sample_count
    left_candidates = minima[minima < center_peak_index]
    right_candidates = minima[minima > center_peak_index]
    left_index = int(left_candidates[-1]) if left_candidates.size else center_peak_index - max(1, sample_count // 4)
    right_index = int(right_candidates[0]) if right_candidates.size else center_peak_index + max(1, sample_count // 4)
    if right_index <= left_index:
        left_index = center_peak_index - max(1, sample_count // 6)
        right_index = center_peak_index + max(1, sample_count // 6)
    selection_mask = np.zeros(sample_count, dtype=bool)
    for extended_index in range(left_index, right_index + 1):
        selection_mask[extended_index % sample_count] = True
    selected_density = density[selection_mask]
    integral_counts_per_cycle = float(np.sum(selected_density) * phase_bin_width_ns)
    phase_min_ns = float(phase[0] - 0.5 * phase_bin_width_ns)
    start_phase_ns = wrap_phase_ns(float(extended_phase[left_index]), phase_min_ns, period_ns)
    end_phase_ns = wrap_phase_ns(float(extended_phase[right_index]), phase_min_ns, period_ns)
    return selection_mask, start_phase_ns, end_phase_ns, integral_counts_per_cycle


def select_prompt_interval_between_outer_dips(
    phase: np.ndarray,
    density: np.ndarray,
    period_ns: float,
    phase_bin_width_ns: float,
    smoothing_bins: int,
    outer_dip_exclusion_half_width_ns: float,
) -> tuple[np.ndarray, float, float, float]:
    sample_count = len(phase)
    if sample_count == 0:
        empty_mask = np.zeros(0, dtype=bool)
        return empty_mask, float("nan"), float("nan"), float("nan")
    smoothed = circular_smooth(density, smoothing_bins)
    minima = local_minima_indices(smoothed)
    left_candidates = minima[phase[minima] <= -outer_dip_exclusion_half_width_ns]
    right_candidates = minima[phase[minima] >= outer_dip_exclusion_half_width_ns]
    if left_candidates.size == 0 or right_candidates.size == 0:
        return select_peak_interval_between_dips(
            phase,
            density,
            period_ns,
            phase_bin_width_ns,
            smoothing_bins,
            lambda values: int(np.argmax(values)),
        )
    left_index = int(left_candidates[-1])
    right_index = int(right_candidates[0])
    if right_index <= left_index:
        return select_peak_interval_between_dips(
            phase,
            density,
            period_ns,
            phase_bin_width_ns,
            smoothing_bins,
            lambda values: int(np.argmax(values)),
        )
    selection_mask = np.zeros(sample_count, dtype=bool)
    selection_mask[left_index : right_index + 1] = True
    integral_counts_per_cycle = float(np.sum(density[selection_mask]) * phase_bin_width_ns)
    start_phase_ns = float(phase[left_index] - 0.5 * phase_bin_width_ns)
    end_phase_ns = float(phase[right_index] + 0.5 * phase_bin_width_ns)
    return selection_mask, start_phase_ns, end_phase_ns, integral_counts_per_cycle


def phase_region_description(start_phase_ns: float, end_phase_ns: float) -> str:
    if not np.isfinite(start_phase_ns) or not np.isfinite(end_phase_ns):
        return "n/a"
    if end_phase_ns < start_phase_ns:
        return f"[{start_phase_ns:+.2f}, {end_phase_ns:+.2f}] ns (wrap)"
    return f"[{start_phase_ns:+.2f}, {end_phase_ns:+.2f}] ns"


def region_from_mask(phase: np.ndarray, mask: np.ndarray, phase_bin_width_ns: float) -> tuple[list[tuple[float, float]], float, float]:
    spans = mask_to_spans(phase, mask, phase_bin_width_ns)
    if not spans:
        return spans, float("nan"), float("nan")
    if len(spans) == 1:
        return spans, spans[0][0], spans[0][1]
    return spans, spans[-1][0], spans[0][1]


def relative_spans_to_absolute_segments(
    spans: list[tuple[float, float]],
    period_ns: float,
    reference_center_ns: float,
) -> list[tuple[float, float]]:
    segments: list[tuple[float, float]] = []
    for start_rel_ns, end_rel_ns in spans:
        absolute_start_ns = (reference_center_ns + start_rel_ns) % period_ns
        absolute_end_ns = (reference_center_ns + end_rel_ns) % period_ns
        if absolute_end_ns > absolute_start_ns:
            segments.append((absolute_start_ns, absolute_end_ns))
        else:
            segments.append((0.0, absolute_end_ns))
            segments.append((absolute_start_ns, period_ns))
    return segments


def tot_range_label(min_ns: float, max_ns: float) -> str:
    if np.isfinite(max_ns):
        return f"{min_ns:g}-{max_ns:g} ns"
    return f">= {min_ns:g} ns"


def draw_gamma_cut_comparison(
    histogram,
    display_histogram,
    run_key: str,
    position_m: float,
    duration_seconds: float,
    accepted_rf_cycles: float,
    accepted_window_count: int,
    total_window_count: int,
    validation_efficiency_window: float,
    assumed_trigger_window_width_ns: float,
    broad_gamma_center_ns: float,
    cut_gamma_center_ns: float,
    period_ns: float,
    broad_tot_min_ns: float,
    broad_tot_max_ns: float,
    cut_tot_min_ns: float,
    cut_tot_max_ns: float,
    neutron_tot_min_ns: float,
    neutron_tot_max_ns: float,
    output_path: Path,
    profile_path: Path,
    args: argparse.Namespace,
) -> tuple[float, float, float, float, float, float]:
    broad_ref_relative_ns = float(signed_delta_ns(np.asarray([broad_gamma_center_ns]), cut_gamma_center_ns, period_ns)[0])
    broad_phase_ns, broad_counts, broad_rate, broad_phase_bin_width_ns = build_phase_profile(
        histogram, duration_seconds, period_ns, cut_gamma_center_ns, broad_tot_min_ns, broad_tot_max_ns
    )
    cut_phase_ns, cut_counts, cut_rate, cut_phase_bin_width_ns = build_phase_profile(
        histogram, duration_seconds, period_ns, cut_gamma_center_ns, cut_tot_min_ns, cut_tot_max_ns
    )
    neutron_phase_ns, neutron_counts, neutron_rate, neutron_phase_bin_width_ns = build_phase_profile(
        histogram, duration_seconds, period_ns, cut_gamma_center_ns, neutron_tot_min_ns, neutron_tot_max_ns
    )
    broad_cycle_density = broad_counts / max(accepted_rf_cycles, 1.0) / broad_phase_bin_width_ns
    cut_cycle_density = cut_counts / max(accepted_rf_cycles, 1.0) / cut_phase_bin_width_ns
    neutron_cycle_density = neutron_counts / max(accepted_rf_cycles, 1.0) / neutron_phase_bin_width_ns

    prompt_phase_candidates = np.where(np.abs(broad_phase_ns) <= args.neutron_phase_search_min_ns)[0]
    prompt_gamma_selection_mask = np.zeros_like(broad_phase_ns, dtype=bool)
    if prompt_phase_candidates.size > 0 and np.any(broad_cycle_density[prompt_phase_candidates] > 0.0):
        prompt_gamma_selection_mask, _, _, _ = select_prompt_interval_between_outer_dips(
            broad_phase_ns,
            broad_cycle_density,
            period_ns,
            broad_phase_bin_width_ns,
            args.integration_smoothing_bins,
            args.gamma_exclusion_half_width_ns,
        )
    else:
        prompt_gamma_selection_mask, _, _, _ = select_peak_interval_between_dips(
            cut_phase_ns,
            cut_cycle_density,
            period_ns,
            cut_phase_bin_width_ns,
            args.integration_smoothing_bins,
            lambda values: int(np.argmax(values)),
        )

    gamma_selection_mask = prompt_gamma_selection_mask.copy()
    neutron_selection_mask = ~gamma_selection_mask

    gamma_selection_spans, gamma_start_phase_ns, gamma_end_phase_ns = region_from_mask(neutron_phase_ns, gamma_selection_mask, neutron_phase_bin_width_ns)
    neutron_selection_spans, neutron_start_phase_ns, neutron_end_phase_ns = region_from_mask(neutron_phase_ns, neutron_selection_mask, neutron_phase_bin_width_ns)
    gamma_integral_counts_per_cycle = float(np.sum(neutron_cycle_density[gamma_selection_mask]) * neutron_phase_bin_width_ns)
    neutron_integral_counts_per_cycle = float(np.sum(neutron_cycle_density[neutron_selection_mask]) * neutron_phase_bin_width_ns) if np.any(neutron_selection_mask) else float("nan")

    gamma_absolute_segments = relative_spans_to_absolute_segments(gamma_selection_spans, period_ns, cut_gamma_center_ns)
    neutron_absolute_segments = relative_spans_to_absolute_segments(neutron_selection_spans, period_ns, cut_gamma_center_ns)

    with profile_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "run_key",
                "phase_relative_ns",
                "broad_gamma_counts",
                "broad_gamma_rate_density_hz_per_ns",
                "broad_gamma_mean_counts_per_accepted_cycle_per_ns",
                "cut_gamma_counts",
                "cut_gamma_rate_density_hz_per_ns",
                "cut_gamma_mean_counts_per_accepted_cycle_per_ns",
                "neutron_high_tot_counts",
                "neutron_high_tot_rate_density_hz_per_ns",
                "neutron_high_tot_mean_counts_per_accepted_cycle_per_ns",
                "in_gamma_integral_region",
                "in_neutron_integral_region",
                "accepted_window_count",
                "total_window_count",
                "validation_efficiency_window",
                "assumed_trigger_window_width_ns",
                "accepted_rf_cycles",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for phase_ns, broad_c, broad_r, broad_cd, cut_c, cut_r, cut_cd, neutron_c, neutron_r, neutron_cd, in_gamma_region, in_neutron_region in zip(
            broad_phase_ns,
            broad_counts,
            broad_rate,
            broad_cycle_density,
            cut_counts,
            cut_rate,
            cut_cycle_density,
            neutron_counts,
            neutron_rate,
            neutron_cycle_density,
            gamma_selection_mask,
            neutron_selection_mask,
        ):
            writer.writerow(
                {
                    "run_key": run_key,
                    "phase_relative_ns": phase_ns,
                    "broad_gamma_counts": broad_c,
                    "broad_gamma_rate_density_hz_per_ns": broad_r,
                    "broad_gamma_mean_counts_per_accepted_cycle_per_ns": broad_cd,
                    "cut_gamma_counts": cut_c,
                    "cut_gamma_rate_density_hz_per_ns": cut_r,
                    "cut_gamma_mean_counts_per_accepted_cycle_per_ns": cut_cd,
                    "neutron_high_tot_counts": neutron_c,
                    "neutron_high_tot_rate_density_hz_per_ns": neutron_r,
                    "neutron_high_tot_mean_counts_per_accepted_cycle_per_ns": neutron_cd,
                    "in_gamma_integral_region": int(in_gamma_region),
                    "in_neutron_integral_region": int(in_neutron_region),
                    "accepted_window_count": accepted_window_count,
                    "total_window_count": total_window_count,
                    "validation_efficiency_window": validation_efficiency_window,
                    "assumed_trigger_window_width_ns": assumed_trigger_window_width_ns,
                    "accepted_rf_cycles": accepted_rf_cycles,
                }
            )

    timewalk_applied = histogram.metadata.get("timewalk_correction_mode", "") == "zebra_linear"
    display_source_histogram = histogram if timewalk_applied else display_histogram
    display_tot_mask = histogram.y_centers_ns <= args.gamma_display_tot_max_ns
    one_period_counts = fold_histogram_to_one_period(display_source_histogram, period_ns)
    tiled_counts = np.tile(one_period_counts, (1, 3))
    display_counts = tiled_counts[display_tot_mask, :]
    vmax = float(np.quantile(display_counts[display_counts > 0.0], 0.995)) if np.any(display_counts > 0.0) else 1.0
    neutron_display_max_ns = min(args.gamma_display_tot_max_ns, neutron_tot_max_ns) if np.isfinite(neutron_tot_max_ns) else args.gamma_display_tot_max_ns
    neutron_label = tot_range_label(neutron_tot_min_ns, neutron_tot_max_ns)
    zebra_corrected_artifact_phase_ns = float(histogram.metadata.get("timewalk_zebra_corrected_artifact_phase_ns", "nan"))
    zebra_slope_ns_per_ns = float(histogram.metadata.get("timewalk_zebra_slope_ns_per_ns", "nan"))
    zebra_tot_reference_ns = float(histogram.metadata.get("timewalk_zebra_tot_reference_ns", "nan"))
    zebra_low_tot_min_ns = float(histogram.metadata.get("timewalk_zebra_low_tot_min_ns", "nan"))
    zebra_high_tot_max_ns = float(histogram.metadata.get("timewalk_zebra_high_tot_max_ns", "nan"))
    display_phase_offset_ns = 0.0
    display_broad_gamma_center_ns = broad_gamma_center_ns
    display_cut_gamma_center_ns = cut_gamma_center_ns
    display_zebra_phase_ns = zebra_corrected_artifact_phase_ns
    display_gamma_segments = gamma_absolute_segments
    display_neutron_segments = neutron_absolute_segments
    if timewalk_applied and np.isfinite(zebra_corrected_artifact_phase_ns):
        display_phase_offset_ns = 0.5 * period_ns - zebra_corrected_artifact_phase_ns
        shifted_one_period_counts = shift_one_period_phase_counts(one_period_counts, display_phase_offset_ns / histogram.x_bin_width_ns)
        display_counts = np.tile(shifted_one_period_counts, (1, 3))[display_tot_mask, :]
        vmax = float(np.quantile(display_counts[display_counts > 0.0], 0.995)) if np.any(display_counts > 0.0) else 1.0
        display_broad_gamma_center_ns = shift_phase_value_for_display(broad_gamma_center_ns, display_phase_offset_ns, period_ns)
        display_cut_gamma_center_ns = shift_phase_value_for_display(cut_gamma_center_ns, display_phase_offset_ns, period_ns)
        display_zebra_phase_ns = shift_phase_value_for_display(zebra_corrected_artifact_phase_ns, display_phase_offset_ns, period_ns)
        display_gamma_segments = shift_segments_for_display(gamma_absolute_segments, display_phase_offset_ns, period_ns)
        display_neutron_segments = shift_segments_for_display(neutron_absolute_segments, display_phase_offset_ns, period_ns)

    figure, axes = plt.subplots(1, 3, figsize=(17.5, 4.8), constrained_layout=True)
    image = axes[0].imshow(
        np.ma.masked_less_equal(display_counts, 0.0),
        aspect="auto",
        origin="lower",
        extent=(
            0.0,
            3.0 * period_ns,
            0.0,
            args.gamma_display_tot_max_ns,
        ),
        interpolation="nearest",
        cmap="magma",
        vmin=0.0,
        vmax=max(1.0, vmax),
    )
    figure.colorbar(image, ax=axes[0], label="counts / 2D bin, one-period fold tiled 3x")
    anchor_band_legend_added = False
    gamma_roi_legend_added = False
    neutron_roi_legend_added = False
    for fold_index in range(3):
        axes[0].axvline(fold_index * period_ns + display_broad_gamma_center_ns, color="#2563eb", linestyle="--", linewidth=1.1)
        axes[0].axvline(fold_index * period_ns + display_cut_gamma_center_ns, color="#00e5ff", linestyle="--", linewidth=1.2)
        axes[0].add_patch(
            patches.Rectangle(
                (fold_index * period_ns, cut_tot_min_ns),
                period_ns,
                cut_tot_max_ns - cut_tot_min_ns,
                linewidth=1.0,
                edgecolor="#f59e0b",
                facecolor="#f59e0b",
                alpha=0.06,
                label=f"anchor ToT cut ({cut_tot_min_ns:.2f}-{cut_tot_max_ns:.2f} ns)" if not anchor_band_legend_added else None,
            )
        )
        anchor_band_legend_added = True
        for start_ns, end_ns in display_gamma_segments:
            axes[0].add_patch(
                patches.Rectangle(
                    (fold_index * period_ns + start_ns, neutron_tot_min_ns),
                    end_ns - start_ns,
                    neutron_display_max_ns - neutron_tot_min_ns,
                    linewidth=1.3,
                    edgecolor="#dc2626",
                    facecolor="#dc2626",
                    alpha=0.18,
                    label=f"gamma efficiency ROI ({neutron_label})" if not gamma_roi_legend_added else None,
                )
            )
            gamma_roi_legend_added = True
        for start_ns, end_ns in display_neutron_segments:
            axes[0].add_patch(
                patches.Rectangle(
                    (fold_index * period_ns + start_ns, neutron_tot_min_ns),
                    end_ns - start_ns,
                    neutron_display_max_ns - neutron_tot_min_ns,
                    linewidth=1.3,
                    edgecolor="#10b981",
                    facecolor="#10b981",
                    alpha=0.16,
                    label=f"neutron efficiency ROI ({neutron_label})" if not neutron_roi_legend_added else None,
                )
            )
            neutron_roi_legend_added = True
    axes[0].set_xlabel("RF phase, three periods (ns)")
    axes[0].set_ylabel("ToT (ns)")
    axes[0].set_title(
        f"{run_key} | {'artifact-centered corrected 2D view' if timewalk_applied else 'raw-count'}"
    )
    axes[0].legend(loc="upper right", fontsize=8, framealpha=0.9)

    axes[1].plot(broad_phase_ns, broad_cycle_density, color="#0891b2", linewidth=1.2, label="broad gamma band mean density")
    axes[1].plot(cut_phase_ns, cut_cycle_density, color="#dc2626", linewidth=1.2, label=f"centered cut gamma mean density ({cut_tot_min_ns:.2f}-{cut_tot_max_ns:.2f} ns ToT)")
    axes[1].plot(neutron_phase_ns, neutron_cycle_density, color="#10b981", linewidth=1.2, label=f"high-ToT {neutron_label} mean density")
    axes[1].axvline(0.0, color="#00e5ff", linestyle="--", linewidth=1.2, label="selected gamma anchor")
    for span_index, (start_ns, end_ns) in enumerate(gamma_selection_spans):
        axes[1].axvspan(start_ns, end_ns, color="#dc2626", alpha=0.08, label="gamma phase ROI" if span_index == 0 else None)
    for span_index, (start_ns, end_ns) in enumerate(neutron_selection_spans):
        axes[1].axvspan(start_ns, end_ns, color="#10b981", alpha=0.08, label="neutron phase ROI" if span_index == 0 else None)
    axes[1].set_xlabel("RF phase relative to selected gamma reference (ns)")
    axes[1].set_ylabel("mean counts / accepted RF cycle / ns")
    axes[1].set_title(f"{run_key} | {'corrected ' if timewalk_applied else ''}cut vs no-cut mean profile")
    axes[1].grid(True, linewidth=0.3, alpha=0.5)
    finite_profile_densities = [
        density[np.isfinite(density)]
        for density in (broad_cycle_density, cut_cycle_density, neutron_cycle_density)
        if np.any(np.isfinite(density))
    ]
    if finite_profile_densities:
        profile_ymax = max(float(np.max(density)) for density in finite_profile_densities)
        axes[1].set_ylim(0.0, profile_ymax * 1.35)
    axes[1].legend(fontsize=8)

    integral_labels = [f"gamma ROI\n{neutron_label}", f"neutron ROI\n{neutron_label}"]
    integral_values = [gamma_integral_counts_per_cycle, neutron_integral_counts_per_cycle]
    integral_colors = ["#dc2626", "#10b981"]
    axes[2].bar(integral_labels, integral_values, color=integral_colors, alpha=0.85)
    axes[2].set_ylabel("integrated counts / accepted RF cycle")
    axes[2].set_title(f"{run_key} | {'corrected ' if timewalk_applied else ''}efficiency ROIs between dips")
    axes[2].grid(True, axis="y", linewidth=0.3, alpha=0.5)
    finite_integrals = [value for value in integral_values if np.isfinite(value)]
    if finite_integrals:
        axes[2].set_ylim(0.0, max(finite_integrals) * 1.22)
    axes[2].text(
        0.02,
        0.90,
        f"gamma: {phase_region_description(gamma_start_phase_ns, gamma_end_phase_ns)}\nneutron: {phase_region_description(neutron_start_phase_ns, neutron_end_phase_ns)}",
        transform=axes[2].transAxes,
        ha="left",
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
    )
    for bar, value in zip(axes[2].patches, integral_values):
        if np.isfinite(value):
            axes[2].text(bar.get_x() + bar.get_width() / 2.0, value, f"{value:.5g}", ha="center", va="bottom", fontsize=8)

    figure.suptitle(
        f"Efficiency finding{' | zebra-fit timewalk correction' if timewalk_applied else ''} | detector face {position_m:g} m"
    )
    figure.savefig(output_path, dpi=180)
    plt.close(figure)
    return gamma_start_phase_ns, gamma_end_phase_ns, gamma_integral_counts_per_cycle, neutron_start_phase_ns, neutron_end_phase_ns, neutron_integral_counts_per_cycle


def write_results(results: list[GammaWindowResult], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[field.name for field in GammaWindowResult.__dataclass_fields__.values()], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerow({field: getattr(result, field) for field in GammaWindowResult.__dataclass_fields__})


def draw_profiles_density(run_key: str, position_m: float, profiles: list[SliceProfile], output_path: Path, gamma_exclusion_half_width_ns: float) -> None:
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

    axes[0].set_ylabel("Hz/ns")
    axes[0].set_title(f"{run_key} | detector face {position_m:g} m | Ch02 phase-rate slices")
    axes[0].grid(True, linewidth=0.3, alpha=0.5)
    axes[0].legend(fontsize=8, ncol=2)

    axes[1].set_xlabel("RF phase relative to prompt gamma (ns)")
    axes[1].set_ylabel("excess significance z")
    axes[1].grid(True, linewidth=0.3, alpha=0.5)

    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def draw_phase_center_summary(results: list[GammaWindowResult], output_path: Path) -> None:
    if not results:
        return
    figure, axis = plt.subplots(figsize=(8.8, 5.2), constrained_layout=True)
    by_run: dict[str, list[GammaWindowResult]] = {}
    for result in results:
        by_run.setdefault(result.run_key, []).append(result)
    colors = ["#2563eb", "#dc2626", "#059669", "#d97706"]
    for color, (run_key, group) in zip(colors, sorted(by_run.items())):
        group = sorted(group, key=lambda item: item.half_width_ns)
        axis.plot([item.half_width_ns for item in group], [item.phase_center_shift_from_1p25_ns for item in group], "o-", color=color, label=run_key)
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_xlabel("Gamma ToT half-width around fitted mu (ns)")
    axis.set_ylabel("phase-center shift relative to ±1.25 ns (ns)")
    axis.set_title("5 m gamma-reference sensitivity to centered fixed-width ToT gate")
    axis.grid(True, linewidth=0.3, alpha=0.5)
    axis.legend(fontsize=8)
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_readme(output_dir: Path, args: argparse.Namespace) -> None:
    text = f"""# Step09 Efficiency Finding

This directory captures the current efficiency-finding view built on the step09 gamma-reference work.

Method:
- use the fitted gamma ToT peak position `mu_gamma` from step08;
- compare the original broad gamma band [{args.broad_gamma_tot_min_ns}, {args.broad_gamma_tot_max_ns}] ns against centered cut gates `[mu_gamma - half_width, mu_gamma + half_width]`;
- scan half widths: {', '.join(f'{width:g} ns' for width in args.half_widths_ns)};
- rebuild the prompt-gamma phase reference and the 1D phase-rate slice diagnostics for each width;
- use a local quadratic peak fit around the cut-band profile maximum as the final phase anchor.
- keep the narrow centered cut only for anchor finding;
- define both efficiency ROIs with the shared ToT cut [{args.neutron_tot_min_ns}, {args.neutron_tot_max_ns}] ns and separate them only by the dip-bounded phase regions.

Interpretation:
- the spread in `phase_center_ns` across the width scan is a direct measure of the selected local-fit gamma-reference sensitivity to the ToT gate choice;
- the `*_gamma_reference.png` plots keep the selected gamma anchor line on the 1D panel and show the current efficiency ROIs only on the 2D panel and the integral summary panel;
- the right-hand 1D profiles are displayed as mean `counts / accepted RF cycle / ns`;
- an additional high-ToT band [{args.neutron_tot_min_ns}, {args.neutron_tot_max_ns}] ns is shown in green to track the candidate neutron population in the same accepted-cycle normalization;
- the integration ROIs are shown on the 2D panel, and the third panel reports their integrated counts per accepted RF cycle;
- the left-hand 2D panels are built from the canonical one-period fold and tiled 3x only for visibility;
- the stored TSV outputs retain `Hz/ns` rate-density columns for wall-clock normalization and add the accepted-cycle-normalized profile columns and validation-efficiency metadata;
- accepted-cycle normalization uses the explicit working-model trigger-window width assumption `assumed_trigger_window_width_ns={args.assumed_trigger_window_width_ns:g}`;
- when `--apply-zebra-timewalk-correction` is enabled, the correction slope is derived from the phase-peak shift between the `2-3 ns` and `6-7 ns` low-ToT zebra bands before rebuilding the same efficiency-style plots;
- this is a systematic-selection study, not a statistical fit error.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    results_root = Path(args.results_root).expanduser().resolve()
    gamma_fit_summary = (
        Path(args.gamma_fit_summary).expanduser().resolve()
        if args.gamma_fit_summary
        else results_root / "Ch02_Position_Comparison" / "gamma_tot_peak_fits" / "gamma_tot_peak_fit_summary.tsv"
    )
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else results_root / "Ch02_Position_Comparison" / "Efficiency_finding"
    output_dir.mkdir(parents=True, exist_ok=True)

    peak_summary = read_peak_summary(gamma_fit_summary)
    sparse_files = find_sparse_files(results_root, args.run_glob)
    if not sparse_files:
        raise FileNotFoundError(f"No Folded_RF_3x phase-ToT sparse files found under {results_root}")

    if len(args.half_widths_ns) == 0:
        raise ValueError("Need at least one half-width")
    reference_half_width = min(args.half_widths_ns, key=lambda value: abs(value - 1.25))

    results: list[GammaWindowResult] = []
    intervals: list[object] = []
    profiles: list[SliceProfile] = []
    timewalk_fits: list[ZebraTimewalkFitResult] = []

    for sparse_path in sparse_files:
        raw_histogram = read_channel_histogram(sparse_path, args.channel)
        histogram = raw_histogram
        run_key = histogram.metadata.get("run_key", parse_run_key(histogram.path))
        if run_key not in peak_summary:
            raise ValueError(f"Missing step08 peak fit for {run_key} in {gamma_fit_summary}")
        if args.apply_zebra_timewalk_correction:
            timewalk_fit, tracked_tot_ns, tracked_phase_ns, _tracked_weights = fit_zebra_timewalk_model(histogram, run_key, args)
            histogram = apply_zebra_timewalk_correction(histogram, timewalk_fit, tracked_tot_ns, tracked_phase_ns)
            timewalk_fits.append(timewalk_fit)
        position_m = parse_position_m(run_key)
        duration_seconds = run_duration_seconds(histogram, run_key)
        exposure = accepted_exposure(histogram, run_key, args.assumed_trigger_window_width_ns)
        mu_gamma = peak_summary[run_key]["peak_mu_ns"]
        broad_gamma_center_ns, _, _, _ = detect_gamma_center_for_window(
            histogram,
            args.broad_gamma_tot_min_ns,
            args.broad_gamma_tot_max_ns,
            args.local_fit_half_window_bins,
        )
        phase_centers_by_width: dict[float, float] = {}

        for half_width_ns in args.half_widths_ns:
            gamma_tot_min_ns = mu_gamma - half_width_ns
            gamma_tot_max_ns = mu_gamma + half_width_ns
            gamma_center_ns, sorted_relative_phase_ns, sorted_counts, period_ns = detect_gamma_center_for_window(
                histogram,
                gamma_tot_min_ns,
                gamma_tot_max_ns,
                args.local_fit_half_window_bins,
            )
            phase_centers_by_width[half_width_ns] = gamma_center_ns
            args.gamma_tot_min_ns = gamma_tot_min_ns
            args.gamma_tot_max_ns = gamma_tot_max_ns

            width_tag = f"hw{half_width_ns:.2f}".replace(".", "p")
            gamma_start_phase_ns, gamma_end_phase_ns, gamma_integral_counts_per_cycle, neutron_start_phase_ns, neutron_end_phase_ns, neutron_integral_counts_per_cycle = draw_gamma_cut_comparison(
                histogram,
                raw_histogram,
                run_key,
                position_m,
                duration_seconds,
                exposure.accepted_rf_cycles,
                exposure.accepted_window_count,
                exposure.total_window_count,
                exposure.validation_efficiency_window,
                exposure.assumed_trigger_window_width_ns,
                broad_gamma_center_ns,
                gamma_center_ns,
                period_ns,
                args.broad_gamma_tot_min_ns,
                args.broad_gamma_tot_max_ns,
                gamma_tot_min_ns,
                gamma_tot_max_ns,
                args.neutron_tot_min_ns,
                args.neutron_tot_max_ns,
                output_dir / f"{run_key}_ch{args.channel:02d}_{width_tag}_gamma_reference.png",
                output_dir / f"{run_key}_ch{args.channel:02d}_{width_tag}_gamma_reference_profile.tsv",
                args,
            )

            run_profiles: list[SliceProfile] = []
            for tot_min_ns, tot_max_ns in zip(args.tot_slice_edges_ns[:-1], args.tot_slice_edges_ns[1:]):
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
                profiles.append(profile)
                intervals.extend(find_intervals(profile, args))
            draw_profiles_density(run_key, position_m, run_profiles, output_dir / f"{run_key}_ch{args.channel:02d}_{width_tag}_phase_rate_slices.png", args.gamma_exclusion_half_width_ns)

            results.append(
                GammaWindowResult(
                    run_key=run_key,
                    position_m=position_m,
                    duration_seconds=duration_seconds,
                    accepted_exposure_seconds=exposure.accepted_exposure_seconds,
                    accepted_window_count=exposure.accepted_window_count,
                    total_window_count=exposure.total_window_count,
                    validation_efficiency_window=exposure.validation_efficiency_window,
                    rejected_window_fraction=exposure.rejected_window_fraction,
                    assumed_trigger_window_width_ns=exposure.assumed_trigger_window_width_ns,
                    accepted_rf_cycles=exposure.accepted_rf_cycles,
                    peak_mu_ns=mu_gamma,
                    half_width_ns=half_width_ns,
                    gamma_tot_min_ns=gamma_tot_min_ns,
                    gamma_tot_max_ns=gamma_tot_max_ns,
                    phase_center_ns=gamma_center_ns,
                    phase_center_shift_from_1p25_ns=0.0,
                    gamma_total_counts=float(np.sum(sorted_counts)),
                    gamma_total_rate_hz=float(np.sum(sorted_counts) / duration_seconds),
                    neutron_tot_min_ns=args.neutron_tot_min_ns,
                    neutron_tot_max_ns=args.neutron_tot_max_ns,
                    gamma_integral_phase_start_ns=gamma_start_phase_ns,
                    gamma_integral_phase_end_ns=gamma_end_phase_ns,
                    gamma_integral_counts_per_accepted_cycle=gamma_integral_counts_per_cycle,
                    neutron_integral_phase_start_ns=neutron_start_phase_ns,
                    neutron_integral_phase_end_ns=neutron_end_phase_ns,
                    neutron_integral_counts_per_accepted_cycle=neutron_integral_counts_per_cycle,
                )
            )

        reference_center = phase_centers_by_width[reference_half_width]
        for index in range(len(results) - len(args.half_widths_ns), len(results)):
            phase_center = results[index].phase_center_ns
            shift = signed_delta_ns(np.asarray([phase_center]), reference_center, period_ns)[0]
            results[index] = GammaWindowResult(
                run_key=results[index].run_key,
                position_m=results[index].position_m,
                duration_seconds=results[index].duration_seconds,
                accepted_exposure_seconds=results[index].accepted_exposure_seconds,
                accepted_window_count=results[index].accepted_window_count,
                total_window_count=results[index].total_window_count,
                validation_efficiency_window=results[index].validation_efficiency_window,
                rejected_window_fraction=results[index].rejected_window_fraction,
                assumed_trigger_window_width_ns=results[index].assumed_trigger_window_width_ns,
                accepted_rf_cycles=results[index].accepted_rf_cycles,
                peak_mu_ns=results[index].peak_mu_ns,
                half_width_ns=results[index].half_width_ns,
                gamma_tot_min_ns=results[index].gamma_tot_min_ns,
                gamma_tot_max_ns=results[index].gamma_tot_max_ns,
                phase_center_ns=results[index].phase_center_ns,
                phase_center_shift_from_1p25_ns=float(shift),
                gamma_total_counts=results[index].gamma_total_counts,
                gamma_total_rate_hz=results[index].gamma_total_rate_hz,
                neutron_tot_min_ns=results[index].neutron_tot_min_ns,
                neutron_tot_max_ns=results[index].neutron_tot_max_ns,
                gamma_integral_phase_start_ns=results[index].gamma_integral_phase_start_ns,
                gamma_integral_phase_end_ns=results[index].gamma_integral_phase_end_ns,
                gamma_integral_counts_per_accepted_cycle=results[index].gamma_integral_counts_per_accepted_cycle,
                neutron_integral_phase_start_ns=results[index].neutron_integral_phase_start_ns,
                neutron_integral_phase_end_ns=results[index].neutron_integral_phase_end_ns,
                neutron_integral_counts_per_accepted_cycle=results[index].neutron_integral_counts_per_accepted_cycle,
            )
        print(f"{run_key}: completed centered gamma-window sensitivity scan")

    results.sort(key=lambda item: (item.position_m, item.run_key, item.half_width_ns))
    write_results(results, output_dir / "gamma_window_sensitivity_summary.tsv")
    write_timewalk_fit_summary(timewalk_fits, output_dir / "timewalk_fit_summary.tsv")
    write_phase_profiles(profiles, output_dir / "phase_rate_slices.tsv")
    write_intervals(intervals, output_dir / "phase_rate_excess_intervals.tsv")
    draw_phase_center_summary(results, output_dir / "gamma_window_phase_center_summary.png")
    write_readme(output_dir, args)
    print(f"Wrote {output_dir}")


if __name__ == "__main__":
    main()