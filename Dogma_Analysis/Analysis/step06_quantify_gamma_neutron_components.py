from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, fields
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

from matplotlib import colors, patches
import matplotlib.pyplot as plt
import numpy as np

from step05_estimate_gamma_neutron_background import (
    SparseHistogram,
    circular_delta_ns,
    circular_smooth,
    count_in_mask,
    find_sparse_files,
    fold_projection_to_one_period,
    mask_between,
    parse_run_key,
    periodic_ranges,
    read_channel_histogram,
    safe_fraction,
)


C_M_PER_NS = 0.299792458
NEUTRON_REST_MEV = 939.56542052


@dataclass(frozen=True)
class PhaseStats:
    center_ns: float
    sigma_ns: float
    q16_ns: float
    q50_ns: float
    q84_ns: float
    q05_ns: float
    q95_ns: float
    q68_half_width_ns: float
    q90_half_width_ns: float
    selected_counts: float


@dataclass(frozen=True)
class LinearStats:
    mean_ns: float
    sigma_ns: float
    q16_ns: float
    q50_ns: float
    q84_ns: float
    q05_ns: float
    q95_ns: float
    selected_counts: float


@dataclass(frozen=True)
class RunMeasurement:
    run_key: str
    position_m: float
    input_file: Path
    duration_seconds: float
    period_ns: float
    gamma_center_ns: float
    gamma_phase_sigma_ns: float
    gamma_phase_q68_half_width_ns: float
    gamma_phase_q90_half_width_ns: float
    gamma_phase_roi_half_width_ns: float
    gamma_tot_mean_ns: float
    gamma_tot_sigma_ns: float
    gamma_tot_q16_ns: float
    gamma_tot_q84_ns: float
    neutron_30mev_delay_ns: float
    neutron_30mev_delay_mod_ns: float
    neutron_center_30mev_ns: float
    gamma_neutron_phase_separation_ns: float
    neutron_raw_phase_sigma_ns: float
    neutron_raw_phase_q68_half_width_ns: float
    neutron_raw_phase_q90_half_width_ns: float


@dataclass(frozen=True)
class ComponentResult:
    run_key: str
    position_m: float
    duration_seconds: float
    period_ns: float
    gamma_center_ns: float
    gamma_phase_roi_half_width_ns: float
    gamma_phase_sigma_ns: float
    gamma_tot_mean_ns: float
    gamma_tot_sigma_ns: float
    neutron_center_30mev_ns: float
    neutron_phase_roi_half_width_ns: float
    neutron_raw_phase_q68_half_width_ns: float
    neutron_width_source: str
    gamma_neutron_phase_separation_ns: float
    gamma_neutron_overlap_ns: float
    component_status: str
    local_background_counts: float
    local_background_bins: int
    local_background_density_per_bin: float
    local_background_rate_per_bin_hz: float
    gamma_raw_counts: float
    gamma_roi_bins: int
    gamma_background_counts: float
    gamma_background_fraction: float
    gamma_net_counts: float
    gamma_rate_hz: float
    neutron_raw_counts: float
    neutron_roi_bins: int
    neutron_background_counts: float
    neutron_background_fraction: float
    neutron_net_counts: float
    neutron_rate_hz: float
    neutron_tot_mean_ns: float
    neutron_tot_sigma_ns: float
    neutron_tot_q16_ns: float
    neutron_tot_q84_ns: float


@dataclass(frozen=True)
class BackgroundControlResult:
    run_key: str
    position_m: float
    duration_seconds: float
    period_ns: float
    gamma_center_ns: float
    neutron_center_30mev_ns: float
    gamma_exclusion_half_width_ns: float
    neutron_exclusion_half_width_ns: float
    selected_counts: float
    selected_bins: int
    density_per_bin: float
    rate_per_bin_hz: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Quantify prompt-gamma, background-control, and neutron-candidate regions "
            "from cleaned Folded_RF_3x phase-vs-ToT sparse histograms."
        )
    )
    parser.add_argument("--results-root", default="/data6/Dogma_analysis_by_Dachi/Results")
    parser.add_argument("--output-dir")
    parser.add_argument("--channel", type=int, default=2)
    parser.add_argument("--run-glob", default="NCAL_20us_Pos_*")
    parser.add_argument("--gamma-tot-min-ns", type=float, default=10.0)
    parser.add_argument("--gamma-tot-max-ns", type=float, default=18.0)
    parser.add_argument("--signal-tot-min-ns", type=float, default=8.0)
    parser.add_argument("--signal-tot-max-ns", type=float, default=30.0)
    parser.add_argument("--gamma-fit-half-window-ns", type=float, default=6.0)
    parser.add_argument("--gamma-phase-sigma-scale", type=float, default=2.0)
    parser.add_argument("--gamma-phase-min-half-width-ns", type=float, default=1.5)
    parser.add_argument("--gamma-phase-max-half-width-ns", type=float, default=4.0)
    parser.add_argument("--neutron-energy-mev", type=float, default=30.0)
    parser.add_argument("--neutron-fit-half-window-ns", type=float, default=16.0)
    parser.add_argument("--neutron-width-reference-positions", type=float, nargs="*", default=[5.0, 5.8, 6.6])
    parser.add_argument("--background-control-positions", type=float, nargs="*", default=[3.4])
    parser.add_argument("--background-exclusion-scale", type=float, default=1.25)
    parser.add_argument("--background-min-gamma-exclusion-ns", type=float, default=6.0)
    parser.add_argument("--background-min-neutron-exclusion-ns", type=float, default=7.5)
    parser.add_argument("--neutron-width-baseline-quantile", type=float, default=0.20)
    parser.add_argument("--neutron-roi-min-half-width-ns", type=float, default=4.0)
    parser.add_argument("--neutron-roi-max-half-width-ns", type=float, default=16.0)
    parser.add_argument("--neutron-roi-gamma-margin-ns", type=float, default=1.0)
    parser.add_argument("--template-min-separation-ns", type=float, default=8.0)
    parser.add_argument("--display-tot-max-ns", type=float, default=55.0)
    return parser.parse_args()


def signed_delta_ns(values_ns: np.ndarray, center_ns: float, period_ns: float) -> np.ndarray:
    return (values_ns - center_ns + 0.5 * period_ns) % period_ns - 0.5 * period_ns


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    if values.size == 0 or float(np.sum(weights)) <= 0.0:
        return float("nan")
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights) - 0.5 * sorted_weights
    total = float(np.sum(sorted_weights))
    if total <= 0.0:
        return float("nan")
    return float(np.interp(quantile * total, cumulative, sorted_values, left=sorted_values[0], right=sorted_values[-1]))


def weighted_linear_stats(values: np.ndarray, weights: np.ndarray) -> LinearStats:
    selected = weights > 0.0
    if not np.any(selected):
        return LinearStats(*(float("nan") for _ in range(7)), selected_counts=0.0)
    selected_values = values[selected]
    selected_weights = weights[selected]
    total = float(np.sum(selected_weights))
    mean = float(np.average(selected_values, weights=selected_weights))
    variance = float(np.average((selected_values - mean) ** 2, weights=selected_weights))
    q16 = weighted_quantile(selected_values, selected_weights, 0.16)
    q50 = weighted_quantile(selected_values, selected_weights, 0.50)
    q84 = weighted_quantile(selected_values, selected_weights, 0.84)
    q05 = weighted_quantile(selected_values, selected_weights, 0.05)
    q95 = weighted_quantile(selected_values, selected_weights, 0.95)
    return LinearStats(mean, float(np.sqrt(max(0.0, variance))), q16, q50, q84, q05, q95, total)


def weighted_phase_stats(
    phase_centers_ns: np.ndarray,
    weights: np.ndarray,
    center_guess_ns: float,
    period_ns: float,
    half_window_ns: float,
) -> PhaseStats:
    deltas = signed_delta_ns(phase_centers_ns, center_guess_ns, period_ns)
    selected = (np.abs(deltas) <= half_window_ns) & (weights > 0.0)
    if not np.any(selected):
        return PhaseStats(*(float("nan") for _ in range(9)), selected_counts=0.0)
    selected_deltas = deltas[selected]
    selected_weights = weights[selected]
    total = float(np.sum(selected_weights))
    mean_delta = float(np.average(selected_deltas, weights=selected_weights))
    variance = float(np.average((selected_deltas - mean_delta) ** 2, weights=selected_weights))
    q16 = weighted_quantile(selected_deltas, selected_weights, 0.16)
    q50 = weighted_quantile(selected_deltas, selected_weights, 0.50)
    q84 = weighted_quantile(selected_deltas, selected_weights, 0.84)
    q05 = weighted_quantile(selected_deltas, selected_weights, 0.05)
    q95 = weighted_quantile(selected_deltas, selected_weights, 0.95)
    q68_half_width = max(abs(q16), abs(q84))
    q90_half_width = max(abs(q05), abs(q95))
    center = float((center_guess_ns + mean_delta) % period_ns)
    return PhaseStats(
        center,
        float(np.sqrt(max(0.0, variance))),
        q16,
        q50,
        q84,
        q05,
        q95,
        q68_half_width,
        q90_half_width,
        total,
    )


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
    if "covered_duration_seconds" in summary:
        return float(summary["covered_duration_seconds"])
    return float("nan")


def neutron_gamma_delay_ns(distance_m: float, neutron_energy_mev: float) -> float:
    gamma_rel = 1.0 + neutron_energy_mev / NEUTRON_REST_MEV
    beta = float(np.sqrt(1.0 - 1.0 / (gamma_rel * gamma_rel)))
    gamma_tof_ns = distance_m / C_M_PER_NS
    neutron_tof_ns = distance_m / (beta * C_M_PER_NS)
    return neutron_tof_ns - gamma_tof_ns


def close_to_any(value: float, candidates: list[float], tolerance: float = 1.0e-6) -> bool:
    return any(abs(value - candidate) <= tolerance for candidate in candidates)


def phase_mask_3x(histogram: SparseHistogram, center_ns: float, half_width_ns: float, period_ns: float) -> np.ndarray:
    x_phase_ns = histogram.x_centers_ns % period_ns
    return circular_delta_ns(x_phase_ns, center_ns, period_ns) <= half_width_ns


def count_bins(y_mask: np.ndarray, x_mask: np.ndarray) -> int:
    return int(np.count_nonzero(y_mask) * np.count_nonzero(x_mask))


def measure_run(histogram: SparseHistogram, args: argparse.Namespace) -> RunMeasurement:
    run_key = histogram.metadata.get("run_key", parse_run_key(histogram.path))
    position_m = float("nan")
    if "Pos_" in run_key:
        position_text = run_key.split("Pos_", 1)[1].split("m_", 1)[0]
        position_m = float(position_text)
    duration_seconds = run_duration_seconds(histogram, run_key)
    period_ns = float(histogram.metadata["deduced_period_ns"])
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns

    gamma_tot_mask = mask_between(histogram.y_centers_ns, args.gamma_tot_min_ns, args.gamma_tot_max_ns)
    signal_tot_mask = mask_between(histogram.y_centers_ns, args.signal_tot_min_ns, args.signal_tot_max_ns)
    gamma_projection_3x = histogram.counts[gamma_tot_mask, :].sum(axis=0)
    gamma_projection_1x = fold_projection_to_one_period(gamma_projection_3x, period_bins)
    gamma_smoothed = circular_smooth(gamma_projection_1x, half_width_bins=2)
    gamma_guess_ns = float(phase_centers_ns[int(np.argmax(gamma_smoothed))])
    gamma_phase = weighted_phase_stats(phase_centers_ns, gamma_projection_1x, gamma_guess_ns, period_ns, args.gamma_fit_half_window_ns)
    gamma_half_width = float(
        np.clip(
            args.gamma_phase_sigma_scale * gamma_phase.sigma_ns,
            args.gamma_phase_min_half_width_ns,
            args.gamma_phase_max_half_width_ns,
        )
    )
    if not np.isfinite(gamma_half_width):
        gamma_half_width = args.gamma_phase_min_half_width_ns

    gamma_phase_mask = phase_mask_3x(histogram, gamma_phase.center_ns, gamma_half_width, period_ns)
    gamma_tot_projection = histogram.counts[:, gamma_phase_mask].sum(axis=1)
    gamma_tot_weights = np.where(gamma_tot_mask, gamma_tot_projection, 0.0)
    gamma_tot = weighted_linear_stats(histogram.y_centers_ns, gamma_tot_weights)

    delay_ns = neutron_gamma_delay_ns(position_m, args.neutron_energy_mev)
    delay_mod_ns = float(delay_ns % period_ns)
    neutron_center_ns = float((gamma_phase.center_ns + delay_mod_ns) % period_ns)
    separation_ns = float(circular_delta_ns(np.asarray([neutron_center_ns]), gamma_phase.center_ns, period_ns)[0])

    signal_projection_3x = histogram.counts[signal_tot_mask, :].sum(axis=0)
    signal_projection_1x = fold_projection_to_one_period(signal_projection_3x, period_bins)
    gamma_exclusion = max(args.background_min_gamma_exclusion_ns, gamma_half_width * args.background_exclusion_scale)
    neutron_deltas = signed_delta_ns(phase_centers_ns, neutron_center_ns, period_ns)
    gamma_deltas = circular_delta_ns(phase_centers_ns, gamma_phase.center_ns, period_ns)
    neutron_fit_weights = np.where(
        (np.abs(neutron_deltas) <= args.neutron_fit_half_window_ns) & (gamma_deltas > gamma_exclusion),
        signal_projection_1x,
        0.0,
    )
    positive_neutron_weights = neutron_fit_weights[neutron_fit_weights > 0.0]
    if positive_neutron_weights.size:
        baseline = float(np.quantile(positive_neutron_weights, args.neutron_width_baseline_quantile))
        neutron_fit_weights = np.maximum(0.0, neutron_fit_weights - baseline)
    neutron_phase = weighted_phase_stats(
        phase_centers_ns,
        neutron_fit_weights,
        neutron_center_ns,
        period_ns,
        args.neutron_fit_half_window_ns,
    )

    return RunMeasurement(
        run_key=run_key,
        position_m=position_m,
        input_file=histogram.path,
        duration_seconds=duration_seconds,
        period_ns=period_ns,
        gamma_center_ns=gamma_phase.center_ns,
        gamma_phase_sigma_ns=gamma_phase.sigma_ns,
        gamma_phase_q68_half_width_ns=gamma_phase.q68_half_width_ns,
        gamma_phase_q90_half_width_ns=gamma_phase.q90_half_width_ns,
        gamma_phase_roi_half_width_ns=gamma_half_width,
        gamma_tot_mean_ns=gamma_tot.mean_ns,
        gamma_tot_sigma_ns=gamma_tot.sigma_ns,
        gamma_tot_q16_ns=gamma_tot.q16_ns,
        gamma_tot_q84_ns=gamma_tot.q84_ns,
        neutron_30mev_delay_ns=delay_ns,
        neutron_30mev_delay_mod_ns=delay_mod_ns,
        neutron_center_30mev_ns=neutron_center_ns,
        gamma_neutron_phase_separation_ns=separation_ns,
        neutron_raw_phase_sigma_ns=neutron_phase.sigma_ns,
        neutron_raw_phase_q68_half_width_ns=neutron_phase.q68_half_width_ns,
        neutron_raw_phase_q90_half_width_ns=neutron_phase.q90_half_width_ns,
    )


def fit_neutron_width_model(measurements: list[RunMeasurement], args: argparse.Namespace) -> tuple[float, float, str]:
    reference = [
        measurement
        for measurement in measurements
        if close_to_any(measurement.position_m, args.neutron_width_reference_positions)
        and np.isfinite(measurement.neutron_raw_phase_q68_half_width_ns)
        and measurement.gamma_neutron_phase_separation_ns >= args.template_min_separation_ns
    ]
    unique_positions = sorted({round(measurement.position_m, 6) for measurement in reference})
    if len(reference) >= 2 and len(unique_positions) >= 2:
        positions = np.asarray([measurement.position_m for measurement in reference], dtype=np.float64)
        widths = np.asarray([measurement.neutron_raw_phase_q68_half_width_ns for measurement in reference], dtype=np.float64)
        slope, intercept = np.polyfit(positions, widths, deg=1)
        if slope < 0.0:
            width = float(np.nanmedian(widths))
            return width, 0.0, "constant_median_separated_runs_nonmonotonic_widths"
        return float(intercept), float(slope), "linear_fit_separated_runs"
    if reference:
        width = float(np.nanmedian([measurement.neutron_raw_phase_q68_half_width_ns for measurement in reference]))
        return width, 0.0, "constant_median_separated_runs"
    return args.neutron_roi_min_half_width_ns, 0.0, "default_minimum"


def neutron_model_width(position_m: float, intercept: float, slope: float, args: argparse.Namespace) -> float:
    width = intercept + slope * position_m
    return float(np.clip(width, args.neutron_roi_min_half_width_ns, args.neutron_roi_max_half_width_ns))


def analyze_components(
    histogram: SparseHistogram,
    measurement: RunMeasurement,
    neutron_width_ns: float,
    width_source: str,
    args: argparse.Namespace,
) -> tuple[ComponentResult, BackgroundControlResult | None, dict[str, np.ndarray]]:
    period_ns = measurement.period_ns
    signal_tot_mask = mask_between(histogram.y_centers_ns, args.signal_tot_min_ns, args.signal_tot_max_ns)
    gamma_tot_mask = mask_between(histogram.y_centers_ns, args.gamma_tot_min_ns, args.gamma_tot_max_ns)
    gamma_phase_mask = phase_mask_3x(histogram, measurement.gamma_center_ns, measurement.gamma_phase_roi_half_width_ns, period_ns)
    max_nonoverlap_width = measurement.gamma_neutron_phase_separation_ns - measurement.gamma_phase_roi_half_width_ns - args.neutron_roi_gamma_margin_ns
    if max_nonoverlap_width >= args.neutron_roi_min_half_width_ns:
        neutron_roi_width_ns = min(neutron_width_ns, max_nonoverlap_width)
    else:
        neutron_roi_width_ns = neutron_width_ns
    neutron_phase_mask = phase_mask_3x(histogram, measurement.neutron_center_30mev_ns, neutron_roi_width_ns, period_ns)

    gamma_exclusion = max(args.background_min_gamma_exclusion_ns, measurement.gamma_phase_roi_half_width_ns * args.background_exclusion_scale)
    neutron_exclusion = max(args.background_min_neutron_exclusion_ns, neutron_width_ns * args.background_exclusion_scale)
    gamma_exclusion_mask = phase_mask_3x(histogram, measurement.gamma_center_ns, gamma_exclusion, period_ns)
    neutron_exclusion_mask = phase_mask_3x(histogram, measurement.neutron_center_30mev_ns, neutron_exclusion, period_ns)
    background_phase_mask = ~(gamma_exclusion_mask | neutron_exclusion_mask)

    local_background_counts = count_in_mask(histogram.counts, signal_tot_mask, background_phase_mask)
    local_background_bins = count_bins(signal_tot_mask, background_phase_mask)
    local_density = safe_fraction(local_background_counts, local_background_bins)
    local_rate_density = safe_fraction(local_density, measurement.duration_seconds)

    gamma_counts = count_in_mask(histogram.counts, gamma_tot_mask, gamma_phase_mask)
    gamma_bins = count_bins(gamma_tot_mask, gamma_phase_mask)
    gamma_background = local_density * gamma_bins
    gamma_net = gamma_counts - gamma_background
    gamma_rate = safe_fraction(gamma_net, measurement.duration_seconds)

    neutron_counts = count_in_mask(histogram.counts, signal_tot_mask, neutron_phase_mask)
    neutron_bins = count_bins(signal_tot_mask, neutron_phase_mask)
    neutron_background = local_density * neutron_bins
    neutron_net = neutron_counts - neutron_background
    neutron_rate = safe_fraction(neutron_net, measurement.duration_seconds)

    neutron_tot_projection = histogram.counts[:, neutron_phase_mask].sum(axis=1)
    neutron_tot_weights = np.where(signal_tot_mask, neutron_tot_projection, 0.0)
    neutron_tot = weighted_linear_stats(histogram.y_centers_ns, neutron_tot_weights)

    overlap_ns = max(0.0, measurement.gamma_phase_roi_half_width_ns + neutron_roi_width_ns - measurement.gamma_neutron_phase_separation_ns)
    component_status = "separated"
    if overlap_ns > 0.0:
        component_status = "gamma_neutron_roi_overlap"
    elif measurement.gamma_neutron_phase_separation_ns < gamma_exclusion + neutron_width_ns:
        component_status = "close_to_gamma_tail"

    result = ComponentResult(
        run_key=measurement.run_key,
        position_m=measurement.position_m,
        duration_seconds=measurement.duration_seconds,
        period_ns=period_ns,
        gamma_center_ns=measurement.gamma_center_ns,
        gamma_phase_roi_half_width_ns=measurement.gamma_phase_roi_half_width_ns,
        gamma_phase_sigma_ns=measurement.gamma_phase_sigma_ns,
        gamma_tot_mean_ns=measurement.gamma_tot_mean_ns,
        gamma_tot_sigma_ns=measurement.gamma_tot_sigma_ns,
        neutron_center_30mev_ns=measurement.neutron_center_30mev_ns,
        neutron_phase_roi_half_width_ns=neutron_roi_width_ns,
        neutron_raw_phase_q68_half_width_ns=measurement.neutron_raw_phase_q68_half_width_ns,
        neutron_width_source=width_source,
        gamma_neutron_phase_separation_ns=measurement.gamma_neutron_phase_separation_ns,
        gamma_neutron_overlap_ns=overlap_ns,
        component_status=component_status,
        local_background_counts=local_background_counts,
        local_background_bins=local_background_bins,
        local_background_density_per_bin=local_density,
        local_background_rate_per_bin_hz=local_rate_density,
        gamma_raw_counts=gamma_counts,
        gamma_roi_bins=gamma_bins,
        gamma_background_counts=gamma_background,
        gamma_background_fraction=safe_fraction(gamma_background, gamma_counts),
        gamma_net_counts=gamma_net,
        gamma_rate_hz=gamma_rate,
        neutron_raw_counts=neutron_counts,
        neutron_roi_bins=neutron_bins,
        neutron_background_counts=neutron_background,
        neutron_background_fraction=safe_fraction(neutron_background, neutron_counts),
        neutron_net_counts=neutron_net,
        neutron_rate_hz=neutron_rate,
        neutron_tot_mean_ns=neutron_tot.mean_ns,
        neutron_tot_sigma_ns=neutron_tot.sigma_ns,
        neutron_tot_q16_ns=neutron_tot.q16_ns,
        neutron_tot_q84_ns=neutron_tot.q84_ns,
    )

    control: BackgroundControlResult | None = None
    if close_to_any(measurement.position_m, args.background_control_positions):
        control = BackgroundControlResult(
            run_key=measurement.run_key,
            position_m=measurement.position_m,
            duration_seconds=measurement.duration_seconds,
            period_ns=period_ns,
            gamma_center_ns=measurement.gamma_center_ns,
            neutron_center_30mev_ns=measurement.neutron_center_30mev_ns,
            gamma_exclusion_half_width_ns=gamma_exclusion,
            neutron_exclusion_half_width_ns=neutron_exclusion,
            selected_counts=local_background_counts,
            selected_bins=local_background_bins,
            density_per_bin=local_density,
            rate_per_bin_hz=local_rate_density,
        )

    masks = {
        "signal_tot_mask": signal_tot_mask,
        "gamma_tot_mask": gamma_tot_mask,
        "gamma_phase_mask": gamma_phase_mask,
        "neutron_phase_mask": neutron_phase_mask,
        "background_phase_mask": background_phase_mask,
    }
    return result, control, masks


def write_dataclass_table(rows: list[object], output_path: Path) -> None:
    if not rows:
        return
    fieldnames = [field.name for field in fields(rows[0])]  # type: ignore[arg-type]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row_object in rows:
            row = {field_name: getattr(row_object, field_name) for field_name in fieldnames}
            if "input_file" in row:
                row["input_file"] = str(row["input_file"])
            writer.writerow(row)


def write_width_model_table(measurements: list[RunMeasurement], intercept: float, slope: float, source: str, output_path: Path, args: argparse.Namespace) -> None:
    fieldnames = [
        "run_key",
        "position_m",
        "width_reference_used",
        "gamma_neutron_phase_separation_ns",
        "neutron_raw_phase_q68_half_width_ns",
        "model_neutron_phase_half_width_ns",
        "model_intercept_ns",
        "model_slope_ns_per_m",
        "model_source",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for measurement in measurements:
            writer.writerow(
                {
                    "run_key": measurement.run_key,
                    "position_m": measurement.position_m,
                    "width_reference_used": close_to_any(measurement.position_m, args.neutron_width_reference_positions)
                    and measurement.gamma_neutron_phase_separation_ns >= args.template_min_separation_ns,
                    "gamma_neutron_phase_separation_ns": measurement.gamma_neutron_phase_separation_ns,
                    "neutron_raw_phase_q68_half_width_ns": measurement.neutron_raw_phase_q68_half_width_ns,
                    "model_neutron_phase_half_width_ns": neutron_model_width(measurement.position_m, intercept, slope, args),
                    "model_intercept_ns": intercept,
                    "model_slope_ns_per_m": slope,
                    "model_source": source,
                }
            )


def draw_roi_overlay(
    histogram: SparseHistogram,
    measurement: RunMeasurement,
    component: ComponentResult,
    masks: dict[str, np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    display_tot_mask = histogram.y_centers_ns <= args.display_tot_max_ns
    display_counts = histogram.counts[display_tot_mask, :]
    positive_counts = display_counts[display_counts > 0.0]
    vmax = float(np.quantile(positive_counts, 0.995)) if positive_counts.size else 1.0
    masked_counts = np.ma.masked_less_equal(display_counts, 0.0)
    figure, axis = plt.subplots(figsize=(13.5, 5.4), constrained_layout=True)
    axis.imshow(
        masked_counts,
        aspect="auto",
        origin="lower",
        extent=(
            float(histogram.x_centers_ns[0] - 0.5 * histogram.x_bin_width_ns),
            float(histogram.x_centers_ns[-1] + 0.5 * histogram.x_bin_width_ns),
            0.0,
            args.display_tot_max_ns,
        ),
        interpolation="nearest",
        cmap="magma",
        norm=colors.LogNorm(vmin=1.0, vmax=max(1.0, vmax)),
    )

    fold_factor = int(histogram.metadata.get("fold_factor", "3"))
    for start_ns, end_ns in periodic_ranges(measurement.gamma_center_ns, component.gamma_phase_roi_half_width_ns, measurement.period_ns, fold_factor):
        axis.add_patch(
            patches.Rectangle(
                (start_ns, args.gamma_tot_min_ns),
                end_ns - start_ns,
                args.gamma_tot_max_ns - args.gamma_tot_min_ns,
                edgecolor="#00e5ff",
                facecolor="none",
                linewidth=1.4,
            )
        )
    for start_ns, end_ns in periodic_ranges(measurement.neutron_center_30mev_ns, component.neutron_phase_roi_half_width_ns, measurement.period_ns, fold_factor):
        axis.add_patch(
            patches.Rectangle(
                (start_ns, args.signal_tot_min_ns),
                end_ns - start_ns,
                args.signal_tot_max_ns - args.signal_tot_min_ns,
                edgecolor="#ffb000",
                facecolor="none",
                linewidth=1.4,
            )
        )

    background_phase_mask = masks["background_phase_mask"]
    period_bins = max(1, round(measurement.period_ns / histogram.x_bin_width_ns))
    background_one_period = background_phase_mask[:period_bins]
    for index, selected in enumerate(background_one_period):
        if not selected:
            continue
        phase_start = index * histogram.x_bin_width_ns
        phase_end = phase_start + histogram.x_bin_width_ns
        for fold_index in range(fold_factor):
            axis.add_patch(
                patches.Rectangle(
                    (fold_index * measurement.period_ns + phase_start, args.signal_tot_min_ns),
                    phase_end - phase_start,
                    args.signal_tot_max_ns - args.signal_tot_min_ns,
                    edgecolor="none",
                    facecolor="#5ec2ff",
                    alpha=0.035,
                )
            )

    axis.set_xlim(0.0, float(histogram.metadata["x_max_ns"]))
    axis.set_ylim(0.0, args.display_tot_max_ns)
    axis.set_xlabel("RF phase, three periods (ns)")
    axis.set_ylabel("ToT (ns)")
    axis.set_title(f"{measurement.run_key} | Ch{args.channel:02d} step06 ROIs | {component.component_status}")
    axis.text(
        0.01,
        0.98,
        "cyan: prompt gamma ROI\norange: 30 MeV neutron-candidate phase ROI\nblue shade: local background sideband",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none", "pad": 4},
    )
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_projection_table(
    histogram: SparseHistogram,
    measurement: RunMeasurement,
    masks: dict[str, np.ndarray],
    output_path: Path,
) -> None:
    duration = measurement.duration_seconds
    period_bins = max(1, round(measurement.period_ns / histogram.x_bin_width_ns))
    phase_centers = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    signal_tot_mask = masks["signal_tot_mask"]
    gamma_tot_mask = masks["gamma_tot_mask"]
    gamma_phase_mask = masks["gamma_phase_mask"]
    neutron_phase_mask = masks["neutron_phase_mask"]

    gamma_phase_projection = fold_projection_to_one_period(histogram.counts[gamma_tot_mask, :].sum(axis=0), period_bins)
    neutron_phase_projection = fold_projection_to_one_period(histogram.counts[signal_tot_mask, :].sum(axis=0), period_bins)
    gamma_tot_projection = histogram.counts[:, gamma_phase_mask].sum(axis=1)
    neutron_tot_projection = histogram.counts[:, neutron_phase_mask].sum(axis=1)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["run_key", "component", "axis", "center_ns", "counts", "rate_hz"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for center_ns, counts in zip(phase_centers, gamma_phase_projection):
            writer.writerow({"run_key": measurement.run_key, "component": "gamma", "axis": "rf_phase", "center_ns": center_ns, "counts": counts, "rate_hz": safe_fraction(float(counts), duration)})
        for center_ns, counts in zip(phase_centers, neutron_phase_projection):
            writer.writerow({"run_key": measurement.run_key, "component": "signal_band", "axis": "rf_phase", "center_ns": center_ns, "counts": counts, "rate_hz": safe_fraction(float(counts), duration)})
        for center_ns, counts in zip(histogram.y_centers_ns, gamma_tot_projection):
            writer.writerow({"run_key": measurement.run_key, "component": "gamma_phase_roi", "axis": "tot", "center_ns": center_ns, "counts": counts, "rate_hz": safe_fraction(float(counts), duration)})
        for center_ns, counts in zip(histogram.y_centers_ns, neutron_tot_projection):
            writer.writerow({"run_key": measurement.run_key, "component": "neutron_phase_roi", "axis": "tot", "center_ns": center_ns, "counts": counts, "rate_hz": safe_fraction(float(counts), duration)})


def draw_projection_plot(histogram: SparseHistogram, measurement: RunMeasurement, component: ComponentResult, masks: dict[str, np.ndarray], output_path: Path) -> None:
    duration = measurement.duration_seconds
    period_bins = max(1, round(measurement.period_ns / histogram.x_bin_width_ns))
    phase_centers = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    signal_tot_mask = masks["signal_tot_mask"]
    gamma_tot_mask = masks["gamma_tot_mask"]
    gamma_phase_mask = masks["gamma_phase_mask"]
    neutron_phase_mask = masks["neutron_phase_mask"]

    gamma_phase_rate = fold_projection_to_one_period(histogram.counts[gamma_tot_mask, :].sum(axis=0), period_bins) / duration
    signal_phase_rate = fold_projection_to_one_period(histogram.counts[signal_tot_mask, :].sum(axis=0), period_bins) / duration
    gamma_tot_rate = histogram.counts[:, gamma_phase_mask].sum(axis=1) / duration
    neutron_tot_rate = histogram.counts[:, neutron_phase_mask].sum(axis=1) / duration

    figure, axes = plt.subplots(2, 1, figsize=(9.4, 7.2), constrained_layout=True)
    axes[0].plot(phase_centers, signal_phase_rate, color="#6b7280", linewidth=1.0, label="signal ToT band")
    axes[0].plot(phase_centers, gamma_phase_rate, color="#0891b2", linewidth=1.0, label="gamma ToT band")
    axes[0].axvline(measurement.gamma_center_ns, color="#00bcd4", linestyle="--", linewidth=1.0, label="gamma center")
    axes[0].axvline(measurement.neutron_center_30mev_ns, color="#f59e0b", linestyle="--", linewidth=1.0, label="30 MeV neutron expected")
    axes[0].set_xlabel("RF phase folded to one period (ns)")
    axes[0].set_ylabel("rate per phase bin (Hz)")
    axes[0].grid(True, linewidth=0.3, alpha=0.5)
    axes[0].legend(fontsize=8)

    axes[1].plot(histogram.y_centers_ns, gamma_tot_rate, color="#0891b2", linewidth=1.0, label="inside gamma phase ROI")
    axes[1].plot(histogram.y_centers_ns, neutron_tot_rate, color="#f59e0b", linewidth=1.0, label="inside neutron-candidate phase ROI")
    axes[1].set_xlim(0.0, 55.0)
    axes[1].set_xlabel("ToT (ns)")
    axes[1].set_ylabel("rate per ToT bin (Hz)")
    axes[1].grid(True, linewidth=0.3, alpha=0.5)
    axes[1].legend(fontsize=8)
    figure.suptitle(f"{measurement.run_key} | Ch02 component projections | {component.component_status}")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_readme(output_dir: Path, args: argparse.Namespace, width_source: str) -> None:
    text = f"""# Step06 Gamma/Neutron Component Quantification

This directory contains the first quantitative gamma/background/neutron-candidate protocol built from cleaned `Folded_RF_3x` phase-vs-ToT sparse histograms.

Channel: Ch{args.channel:02d}

Purpose:
- detect the compact prompt-gamma centroid and width independently for each run;
- predict the neutron-candidate phase from the prompt gamma plus the gamma-neutron TOF delay for a {args.neutron_energy_mev:g} MeV reference neutron;
- derive a position-dependent neutron-candidate phase width from separated reference positions when possible;
- define background sidebands by excluding both the gamma ROI and the expected neutron-candidate ROI;
- report raw counts, background estimates, net counts, and rates for gamma and neutron-candidate ROIs.

Important interpretation notes:
- The 3.4 m runs are treated as background-control candidates because the neutron-candidate band is expected to overlap or hide near the prompt-gamma region, not because neutrons are absent.
- The neutron-candidate ROI is not a final neutron identification. It is a TOF-guided candidate region for later template or mixture fitting.
- Neutron-candidate width is allowed to vary with detector distance. Current width model source: `{width_source}`.
- Runs flagged as `gamma_neutron_roi_overlap` or `close_to_gamma_tail` should not be used for direct neutron-yield extraction without a template fit.

Default regions:
- Gamma ToT band: [{args.gamma_tot_min_ns}, {args.gamma_tot_max_ns}) ns.
- Signal/neutron-candidate ToT band: [{args.signal_tot_min_ns}, {args.signal_tot_max_ns}) ns.
- Background-control positions: {', '.join(f'{position:g} m' for position in args.background_control_positions)}.
- Neutron-width reference positions: {', '.join(f'{position:g} m' for position in args.neutron_width_reference_positions)}.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    results_root = Path(args.results_root).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else results_root / "Ch02_Position_Comparison" / "gamma_neutron_components"
    overlay_dir = output_dir / "overlays"
    projection_dir = output_dir / "projections"
    overlay_dir.mkdir(parents=True, exist_ok=True)
    projection_dir.mkdir(parents=True, exist_ok=True)

    sparse_files = find_sparse_files(results_root, args.run_glob)
    if not sparse_files:
        raise FileNotFoundError(f"No Folded_RF_3x phase-ToT sparse files found under {results_root}")

    histograms = [read_channel_histogram(path, args.channel) for path in sparse_files]
    measurements = [measure_run(histogram, args) for histogram in histograms]
    measurements.sort(key=lambda measurement: (measurement.position_m, measurement.run_key))
    histogram_by_run = {
        histogram.metadata.get("run_key", parse_run_key(histogram.path)): histogram
        for histogram in histograms
    }

    intercept, slope, width_source = fit_neutron_width_model(measurements, args)

    components: list[ComponentResult] = []
    controls: list[BackgroundControlResult] = []
    for measurement in measurements:
        histogram = histogram_by_run[measurement.run_key]
        width = neutron_model_width(measurement.position_m, intercept, slope, args)
        component, control, masks = analyze_components(histogram, measurement, width, width_source, args)
        components.append(component)
        if control is not None:
            controls.append(control)
        draw_roi_overlay(histogram, measurement, component, masks, overlay_dir / f"{measurement.run_key}_ch{args.channel:02d}_step06_rois.png", args)
        write_projection_table(histogram, measurement, masks, projection_dir / f"{measurement.run_key}_ch{args.channel:02d}_component_projections.tsv")
        draw_projection_plot(histogram, measurement, component, masks, projection_dir / f"{measurement.run_key}_ch{args.channel:02d}_component_projections.png")
        print(
            f"{measurement.run_key}: status={component.component_status}, "
            f"gamma_rate={component.gamma_rate_hz:.3f} Hz, neutron_candidate_rate={component.neutron_rate_hz:.3f} Hz"
        )

    write_dataclass_table(measurements, output_dir / "gamma_neutron_measurements.tsv")
    write_dataclass_table(components, output_dir / "gamma_neutron_component_summary.tsv")
    write_dataclass_table(controls, output_dir / "background_control_regions.tsv")
    write_width_model_table(measurements, intercept, slope, width_source, output_dir / "neutron_width_model.tsv", args)
    write_readme(output_dir, args, width_source)
    print(f"Wrote {output_dir}")


if __name__ == "__main__":
    main()
