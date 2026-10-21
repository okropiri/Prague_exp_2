from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import TypedDict
from analysis_provenance import record_arguments

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import uproot

from step02_run_cleaned_all_channel_rates import compile_cpp, parse_metadata, run_command
from step04_write_cleaned_all_channel_rf_period_outputs import (
    build_dense,
    read_matrix,
    read_sparse_histogram,
    write_profile_outputs,
    write_sparse_2d_outputs,
)


@dataclass
class ClockModel:
    periods: np.ndarray
    errors: np.ndarray
    residuals: np.ndarray
    blocks: list[dict[str, float | int]]


class Sensitivity(TypedDict):
    half_span: int
    max_phase_change_ns: float
    p999_phase_change_ns: float


class ClockAssessment(TypedDict, total=False):
    verdict: str
    windows: int
    prescaled_ticks: int
    missed_prescaled_ticks: int
    prescaler: int
    half_span: int
    mean_period_ns: float
    endpoint_period_ns: float
    min_local_period_ns: float
    max_local_period_ns: float
    prescaled_interval_residual_rms_ns: float
    prescaled_interval_residual_p999_abs_ns: float
    max_half_window_phase_disagreement_ns: float
    max_3sigma_conditional_phase_error_ns: float
    max_folded_delay_ns: float
    phase_tolerance_ns: float
    interval_span_sensitivity: list[Sensitivity]
    phase_precision_status: str
    assumed_ch0_edge_precision_ns: float
    limitations: list[str]
    board: str
    tdc_clock_mhz: float | None
    coarse_quantum_ns: float
    reference_gate_rise_ns: list[float]
    reference_gate_tot_ns: list[float]
    phase_origin: str
    input_state: dict[str, tuple[int, int]]
    source_sha256: dict[str, str]
    folded_output_verification: str
    png_files: int
    pdf_files: int
    tdc_ordinal: int
    allow_reference_gaps: bool
    decoded_selected_board_windows: int
    decoded_board0_headers: int
    represented_cleaned_windows: int
    reference_gap_windows: int
    empty_cleaned_windows: int
    mean_frequency_mhz: float
    detailed_channels: list[int]
    timing_axes_ns: dict[str, list[float]]
    tot_axis_ns: list[float]
    folded_channel_counts: list[int]
    unreferenced_channel_counts: list[int]
    unresolved_clock_channel_counts: list[int]
    clock_model: str
    calibrated_reference_windows: int
    unresolved_reference_windows: int
    cadence_regions: list[dict[str, float | int]]
    cadence_regimes: list[dict[str, float | int]]
    clock_regions: list[dict[str, float | int]]
    clock_gap_tolerance_ns: float
    clock_tick_scope: str
    unbridged_clock_gaps: int
    clock_support_reason_counts: dict[str, int]
    burst_clock: bool
    fit_reference_windows: int
    validation_regions: int
    maximum_regional_phase_bound_ns: float
    minimum_fit_references: int
    minimum_validation_references: int
    include_all_ch0_hits: bool
    clock_interval_bounds_ns: list[float]
    noncore_clock_intervals: int
    longest_noncore_interval_run: int


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def tick_numbers(times: np.ndarray, prescaler: int) -> np.ndarray:
    require(prescaler > 0, "Prescaler must be positive")
    require(times.ndim == 1 and times.size >= 32 and np.isfinite(times).all(), "Invalid ch0 timeline")
    intervals = np.diff(times)
    require(np.all(intervals > 0), "Ch0 timestamps must increase")
    nominal = float(np.median(intervals))
    require(35 <= nominal / prescaler <= 45, "Clock interval outside the declared 35-45 ns RF range")
    gaps = np.rint(intervals / nominal).astype(np.int64)
    require(np.all(gaps >= 1), "Ambiguous short/duplicate ch0 interval")
    require(np.all((intervals / gaps / prescaler >= 35) & (intervals / gaps / prescaler <= 45)),
            "Ch0 intervals are incompatible with integer prescaler ticks")
    return np.r_[0, np.cumsum(gaps)]


def interval_periods(times: np.ndarray, ticks: np.ndarray, prescaler: int,
                     half_span: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return interval_block(times, ticks, prescaler, half_span, 0, times.size)


def interval_block(times: np.ndarray, ticks: np.ndarray, prescaler: int,
                   half_span: int, start: int, stop: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    require(half_span >= 1, "Clock interval half-span must be positive")
    indexes = np.arange(start, stop)
    low = np.maximum(0, indexes - half_span)
    high = np.minimum(times.size - 1, indexes + half_span)
    elapsed_ticks = ticks[high] - ticks[low]
    require(np.all(elapsed_ticks > 0), "Empty clock interval")
    periods = (times[high] - times[low]) / elapsed_ticks / prescaler
    require(np.isfinite(periods).all() and np.all((periods >= 35) & (periods <= 45)),
            "Measured RF interval average outside 35-45 ns")
    return periods, low, high


def calibrate_clock(times: np.ndarray, prescaler: int, max_delay_ns: float,
                    phase_tolerance_ns: float, half_span: int = 4,
                    edge_precision_ns: float = .1,
                    window_ids: np.ndarray | None = None) -> tuple[ClockModel, ClockAssessment]:
    require(max_delay_ns > 0 and phase_tolerance_ns > 0, "Invalid phase-error requirement")
    require(half_span >= 2, "Use at least two neighboring clock intervals on each side")
    require(np.isfinite(edge_precision_ns) and edge_precision_ns > 0, "Invalid ch0 edge precision")
    if window_ids is not None:
        require(window_ids.shape == times.shape and np.all(np.diff(window_ids) > 0),
                "Invalid original reference-window IDs")
    ticks = tick_numbers(times, prescaler)
    periods, errors = np.empty(times.size), np.empty(times.size)
    blocks: list[dict[str, float | int]] = []
    maximum_disagreement = maximum_statistical_phase = 0.0
    chunk = 65536
    for start in range(0, times.size, chunk):
        stop = min(times.size, start + chunk)
        local, low, high = interval_block(times, ticks, prescaler, half_span, start, stop)
        edge_errors = np.sqrt(2) * edge_precision_ns / (ticks[high] - ticks[low]) / prescaler
        midpoint = (low + high) // 2
        first_half = (times[midpoint] - times[low]) / (ticks[midpoint] - ticks[low]) / prescaler
        second_half = (times[high] - times[midpoint]) / (ticks[high] - ticks[midpoint]) / prescaler
        require(np.isfinite(first_half).all() and np.isfinite(second_half).all(),
                "Insufficient half-interval support")
        disagreement = max_delay_ns * abs(first_half - second_half) / local
        maximum_disagreement = max(maximum_disagreement, float(disagreement.max()))
        maximum_statistical_phase = max(maximum_statistical_phase,
                                        float((max_delay_ns * 3 * edge_errors / local).max()))
        periods[start:stop], errors[start:stop] = local, edge_errors
        for offset in range(0, stop - start, 512):
            end = min(stop - start, offset + 512)
            first, last = start + offset, start + end
            blocks.append({
                "first_window": int(window_ids[first]) if window_ids is not None else first + 1,
                "last_window": int(window_ids[last - 1]) if window_ids is not None else last,
                "start_seconds": float(times[first] * 1e-9),
                "end_seconds": float(times[last - 1] * 1e-9),
                "mean_local_period_ns": float(local[offset:end].mean()),
                "min_local_period_ns": float(local[offset:end].min()),
                "max_local_period_ns": float(local[offset:end].max()),
                "max_conditional_edge_error_ns": float(edge_errors[offset:end].max()),
                "max_half_interval_phase_difference_ns": float(disagreement[offset:end].max()),
            })
    sensitivity: list[Sensitivity] = []
    for span in sorted({max(1, half_span // 2), half_span, 2 * half_span}):
        phase_change = np.empty(times.size)
        for start in range(0, times.size, chunk):
            stop = min(times.size, start + chunk)
            alternative, _, _ = interval_block(times, ticks, prescaler, span, start, stop)
            phase_change[start:stop] = max_delay_ns * abs(alternative - periods[start:stop]) / periods[start:stop]
        sensitivity.append({"half_span": span, "max_phase_change_ns": float(phase_change.max()),
                            "p999_phase_change_ns": float(np.quantile(phase_change, .999))})
    del phase_change
    require(maximum_statistical_phase <= phase_tolerance_ns,
            "Assumed edge-measurement error exceeds phase budget; use a longer interval")
    endpoint_period = float((times[-1] - times[0]) / (ticks[-1] - ticks[0]) / prescaler)
    require(abs(float(periods.mean()) - endpoint_period) < 1e-5,
            "Local interval averages disagree with full-span clock measurement")
    residuals = np.empty(times.size - 1)
    for start in range(0, residuals.size, chunk):
        stop = min(residuals.size, start + chunk)
        residuals[start:stop] = (np.diff(times[start:stop + 1]) - np.diff(ticks[start:stop + 1])
                                * prescaler * (periods[start + 1:stop + 1] + periods[start:stop]) / 2)
    result: ClockAssessment = {
        "verdict": "PASS_CH0_PERIOD_AND_RECONSTRUCTION_CONDITIONAL_ON_CLOCK_AND_PRESCALER",
        "phase_precision_status": "INSTANTANEOUS_RF_PHASE_NOT_CERTIFIED_BETWEEN_PRESCALED_EDGES",
        "windows": int(times.size), "prescaled_ticks": int(ticks[-1]),
        "missed_prescaled_ticks": int(ticks[-1] - (times.size - 1)),
        "prescaler": prescaler, "half_span": half_span,
        "assumed_ch0_edge_precision_ns": edge_precision_ns,
        "mean_period_ns": float(periods.mean()), "endpoint_period_ns": endpoint_period,
        "min_local_period_ns": float(periods.min()), "max_local_period_ns": float(periods.max()),
        "prescaled_interval_residual_rms_ns": float(np.sqrt(np.mean(residuals ** 2))),
        "prescaled_interval_residual_p999_abs_ns": float(np.quantile(abs(residuals), .999)),
        "max_half_window_phase_disagreement_ns": maximum_disagreement,
        "max_3sigma_conditional_phase_error_ns": maximum_statistical_phase,
        "max_folded_delay_ns": max_delay_ns, "phase_tolerance_ns": phase_tolerance_ns,
        "interval_span_sensitivity": sensitivity,
        "limitations": [
            "Absolute frequency scale is conditional on the supplied decoder clock and physical prescaler.",
            "Edge-error propagation excludes clock-scale/calibration systematics and RF phase noise.",
            "Boards sharing a DAQ timing chain do not provide independent absolute-frequency traceability.",
            "Local periods are averages over neighboring prescaled pulses, not observations of individual RF cycles.",
            "Span sensitivity and half-interval disagreement are reported, not hidden or treated as edge-measurement noise.",
        ],
    }
    return ClockModel(periods, errors, residuals, blocks), result


def write_table(path: Path, rows: list[dict[str, float | int]] | list[Sensitivity]) -> None:
    require(bool(rows), f"Empty output table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def calibrate_cadence_clock(times: np.ndarray, prescaler: int, max_delay_ns: float,
                            phase_tolerance_ns: float, factors: list[int],
                            half_span: int = 4, edge_precision_ns: float = .1,
                            window_ids: np.ndarray | None = None) -> tuple[ClockModel, ClockAssessment]:
    require(window_ids is not None and window_ids.shape == times.shape
            and np.all(np.diff(window_ids) > 0), "Cadence mode requires original reference IDs")
    require(bool(factors) and factors == sorted(set(factors)) and factors[0] == 1
            and all(1 <= value <= 16 for value in factors), "Invalid cadence factors")
    require(times.ndim == 1 and times.size >= 32 and np.isfinite(times).all()
            and prescaler > 0, "Invalid cadence timeline or prescaler")
    delta = np.diff(times)
    require(np.all(delta > 0), "Ch0 timestamps must increase")
    nominal = float(np.median(delta))
    require(35 <= nominal / prescaler <= 45, "Dominant cadence outside 35-45 ns baseline")
    labels = np.zeros(delta.size, dtype=np.int8)
    interval_ticks = np.zeros(delta.size, dtype=np.int64)
    for factor in factors:
        counts = np.rint(delta / (nominal / factor)).astype(np.int64)
        compatible = ((labels == 0) & (counts >= 1)
                      & (np.abs(delta / (nominal / factor) - counts) <= .01 * np.maximum(counts, 1)))
        labels[compatible] = factor
        interval_ticks[compatible] = counts[compatible]
    require(np.all(labels > 0), "Intervals incompatible with declared cadence factors")
    # A longer interval across missing references can be a missed tick within a fast region.
    for factor in factors:
        middle = ((labels[:-2] == factor) & (labels[2:] == factor)
                  & (labels[1:-1] < factor) & (np.diff(window_ids)[1:-1] > 1))
        positions = np.flatnonzero(middle) + 1
        counts = np.rint(delta[positions] / (nominal / factor)).astype(np.int64)
        require(np.all(np.abs(delta[positions] / (nominal / factor) - counts) <= .01 * counts),
                "Ambiguous missing tick at a cadence boundary")
        labels[positions] = factor
        interval_ticks[positions] = counts
    starts = np.r_[0, np.flatnonzero(np.diff(labels)) + 1]
    stops = np.r_[starts[1:], delta.size]
    model = ClockModel(np.full(times.size, np.nan), np.full(times.size, np.nan),
                       np.full(delta.size, np.nan), [])
    regions: list[dict[str, float | int]] = []
    supported: list[tuple[int, int, int]] = []
    maximum_disagreement = maximum_statistical_phase = 0.0
    for start, stop in zip(starts, stops, strict=True):
        first, last, factor = int(start), int(stop), int(labels[start])
        tick_count = int(interval_ticks[first:last].sum())
        region: dict[str, float | int] = {
            "first_reference_window": int(window_ids[first]),
            "last_reference_window": int(window_ids[last]),
            "cadence_factor": factor, "intervals": last - first,
            "inferred_prescaled_ticks": tick_count,
            "endpoint_period_ns": float((times[last] - times[first]) / tick_count / prescaler),
            "calibrated_reference_windows": 0,
        }
        if last - first >= 31:
            scaled = (times[first:last + 1] - times[first]) * factor
            regional, result = calibrate_clock(
                scaled, prescaler, max_delay_ns, phase_tolerance_ns,
                half_span, edge_precision_ns * factor, window_ids[first:last + 1])
            low = int(first != 0)
            high = regional.periods.size - int(last != times.size - 1)
            model.periods[first + low:first + high] = regional.periods[low:high] / factor
            model.errors[first + low:first + high] = regional.errors[low:high] / factor
            model.residuals[first:last] = regional.residuals / factor
            for block in regional.blocks:
                block["start_seconds"] = float(block["start_seconds"]) / factor + times[first] * 1e-9
                block["end_seconds"] = float(block["end_seconds"]) / factor + times[first] * 1e-9
                for field in ("mean_local_period_ns", "min_local_period_ns", "max_local_period_ns",
                              "max_conditional_edge_error_ns"):
                    block[field] = float(block[field]) / factor
                block["cadence_factor"] = factor
            model.blocks.extend(regional.blocks)
            maximum_disagreement = max(maximum_disagreement, result["max_half_window_phase_disagreement_ns"])
            maximum_statistical_phase = max(maximum_statistical_phase, result["max_3sigma_conditional_phase_error_ns"])
            region["calibrated_reference_windows"] = high - low
            supported.append((first, last, factor))
            del regional, scaled
        regions.append(region)
    resolved = np.isfinite(model.periods)
    require(np.count_nonzero(resolved) >= 32, "Insufficient supported cadence references")
    sensitivity: list[Sensitivity] = []
    for span in sorted({max(1, half_span // 2), half_span, 2 * half_span}):
        changes = np.full(times.size, np.nan)
        for first, last, factor in supported:
            scaled = (times[first:last + 1] - times[first]) * factor
            ticks = tick_numbers(scaled, prescaler)
            low, high = int(first != 0), len(scaled) - int(last != times.size - 1)
            for start in range(low, high, 65536):
                stop = min(high, start + 65536)
                alternative, _, _ = interval_block(scaled, ticks, prescaler, span, start, stop)
                base = model.periods[first + start:first + stop]
                changes[first + start:first + stop] = max_delay_ns * abs(alternative / factor - base) / base
        valid_changes = changes[resolved]
        sensitivity.append({"half_span": span, "max_phase_change_ns": float(valid_changes.max()),
                            "p999_phase_change_ns": float(np.quantile(valid_changes, .999))})
        del changes, valid_changes
    regimes: list[dict[str, float | int]] = []
    for factor in factors:
        selected = [row for row in regions if row["cadence_factor"] == factor]
        if not selected:
            continue
        ticks = sum(int(row["inferred_prescaled_ticks"]) for row in selected)
        elapsed = sum(float(row["endpoint_period_ns"]) * int(row["inferred_prescaled_ticks"]) * prescaler
                      for row in selected)
        measured = elapsed / ticks / prescaler
        regimes.append({"cadence_factor": factor, "intervals": sum(int(row["intervals"]) for row in selected),
                        "inferred_prescaled_ticks": ticks, "period_ns": measured,
                        "frequency_mhz": 1000 / measured,
                        "calibrated_reference_windows": sum(int(row["calibrated_reference_windows"]) for row in selected)})
    residuals = model.residuals[np.isfinite(model.residuals)]
    result = {
        "verdict": "PASS_CADENCE_AWARE_CH0_RECONSTRUCTION_WITH_EXPLICIT_UNRESOLVED_WINDOWS",
        "phase_precision_status": "INSTANTANEOUS_RF_PHASE_NOT_CERTIFIED_BETWEEN_PRESCALED_EDGES",
        "clock_model": "Fixed physical divider; declared frequency-cadence factors; no cross-boundary averaging",
        "windows": int(times.size), "calibrated_reference_windows": int(resolved.sum()),
        "unresolved_reference_windows": int((~resolved).sum()),
        "prescaled_ticks": int(interval_ticks.sum()),
        "missed_prescaled_ticks": int(interval_ticks.sum() - delta.size),
        "prescaler": prescaler, "half_span": half_span,
        "assumed_ch0_edge_precision_ns": edge_precision_ns,
        "mean_period_ns": float(model.periods[resolved].mean()),
        "endpoint_period_ns": float((times[-1] - times[0]) / interval_ticks.sum() / prescaler),
        "min_local_period_ns": float(model.periods[resolved].min()),
        "max_local_period_ns": float(model.periods[resolved].max()),
        "prescaled_interval_residual_rms_ns": float(np.sqrt(np.mean(residuals ** 2))),
        "prescaled_interval_residual_p999_abs_ns": float(np.quantile(abs(residuals), .999)),
        "max_half_window_phase_disagreement_ns": maximum_disagreement,
        "max_3sigma_conditional_phase_error_ns": maximum_statistical_phase,
        "max_folded_delay_ns": max_delay_ns, "phase_tolerance_ns": phase_tolerance_ns,
        "interval_span_sensitivity": sensitivity, "cadence_regions": regions, "cadence_regimes": regimes,
        "limitations": [
            "Physical divider remains explicitly assumed, not hardware-confirmed; divider changes remain an alternative interpretation.",
            "Cadence factors describe measured interval ratios, not forced oscillator frequencies.",
            "References at cadence transitions or in regions with fewer than 32 references have no assigned period.",
            "Frequency/phase estimates are conditional on decoder scale; unobserved cycle phase noise is not certified.",
            "The overall mean/display period is a mixed-regime normalization, not a single oscillator period.",
        ],
    }
    return model, result


def calibrate_gap_clock(times: np.ndarray, prescaler: int, max_delay_ns: float,
                        phase_tolerance_ns: float, gap_tolerance_ns: float,
                        half_span: int = 4, edge_precision_ns: float = .1,
                        window_ids: np.ndarray | None = None) -> tuple[ClockModel, ClockAssessment]:
    require(times.ndim == 1 and times.size >= 32 and np.isfinite(times).all()
            and prescaler > 0 and half_span >= 2, "Invalid gap-aware clock inputs")
    require(window_ids is not None and window_ids.shape == times.shape
            and np.all(np.diff(window_ids) > 0), "Gap-aware mode requires original reference IDs")
    require(max_delay_ns > 0 and phase_tolerance_ns > 0
            and np.isfinite(gap_tolerance_ns) and gap_tolerance_ns > 0
            and np.isfinite(edge_precision_ns) and edge_precision_ns > 0,
            "Invalid gap-aware error/support requirements")
    delta = np.diff(times)
    require(np.all(delta > 0), "Ch0 timestamps must increase")
    nominal = float(np.median(delta))
    require(35 <= nominal / prescaler <= 45, "Dominant interval outside declared RF range")
    connected = (np.diff(window_ids) == 1) & (abs(delta - nominal) <= gap_tolerance_ns)
    starts = np.r_[0, np.flatnonzero(~connected) + 1]
    stops = np.r_[starts[1:], len(times)]
    sizes = stops - starts
    supported = sizes >= max(32, 4 * half_span + 1)
    model = ClockModel(np.full(times.size, np.nan), np.full(times.size, np.nan),
                       np.full(times.size - 1, np.nan), [])
    spans = sorted({max(1, half_span // 2), half_span, 2 * half_span})
    changes = {span: [] for span in spans}
    maximum_disagreement = maximum_statistical_phase = 0.0
    excluded = {"unsupported_or_boundary": 0, "phase_consistency": 0, "resolved": 0}
    for start in range(0, times.size, 65536):
        positions = np.arange(start, min(times.size, start + 65536))
        region = np.searchsorted(stops, positions, side="right")
        eligible = (supported[region] & (positions - starts[region] >= 2 * half_span)
                    & (stops[region] - 1 - positions >= 2 * half_span))
        excluded["unsupported_or_boundary"] += int((~eligible).sum())
        selected = positions[eligible]
        local = (times[selected + half_span] - times[selected - half_span]) / (2 * half_span * prescaler)
        first = (times[selected] - times[selected - half_span]) / (half_span * prescaler)
        second = (times[selected + half_span] - times[selected]) / (half_span * prescaler)
        disagreement = max_delay_ns * abs(first - second) / local
        error = np.sqrt(2) * edge_precision_ns / (2 * half_span * prescaler)
        statistical = max_delay_ns * 3 * error / local
        valid = ((local >= 35) & (local <= 45) & (disagreement <= phase_tolerance_ns)
                 & (statistical <= phase_tolerance_ns))
        alternatives = {}
        for span in spans:
            alternative = (times[selected + span] - times[selected - span]) / (2 * span * prescaler)
            alternatives[span] = max_delay_ns * abs(alternative - local) / local
            valid &= alternatives[span] <= phase_tolerance_ns
        excluded["phase_consistency"] += int((~valid).sum())
        excluded["resolved"] += int(valid.sum())
        chosen = selected[valid]
        model.periods[chosen], model.errors[chosen] = local[valid], error
        if chosen.size:
            maximum_disagreement = max(maximum_disagreement, float(disagreement[valid].max()))
            maximum_statistical_phase = max(maximum_statistical_phase, float(statistical[valid].max()))
        for span in spans:
            changes[span].append(alternatives[span][valid])
        for offset in range(start, min(times.size, start + 65536), 512):
            end = min(times.size, offset + 512, start + 65536)
            values = model.periods[offset:end]
            finite = values[np.isfinite(values)]
            if finite.size:
                model.blocks.append({"first_window": int(window_ids[offset]),
                                     "last_window": int(window_ids[end - 1]),
                                     "mean_local_period_ns": float(finite.mean()),
                                     "min_local_period_ns": float(finite.min()),
                                     "max_local_period_ns": float(finite.max())})
    resolved = np.isfinite(model.periods)
    require(int(resolved.sum()) >= 32, "Insufficient supported gap-aware references")
    for start in range(0, len(delta), 65536):
        stop = min(len(delta), start + 65536)
        mask = resolved[start:stop] & resolved[start + 1:stop + 1] & connected[start:stop]
        residual = delta[start:stop] - prescaler * (model.periods[start:stop] + model.periods[start + 1:stop + 1]) / 2
        model.residuals[start:stop][mask] = residual[mask]
    residuals = model.residuals[np.isfinite(model.residuals)]
    require(residuals.size > 0, "No adjacent supported references for residual checks")
    sensitivity = []
    for span in spans:
        values = np.concatenate(changes[span])
        sensitivity.append({"half_span": span, "max_phase_change_ns": float(values.max()),
                            "p999_phase_change_ns": float(np.quantile(values, .999))})
    regions = [{"first_reference_window": int(window_ids[first]),
                "last_reference_window": int(window_ids[last - 1]),
                "references": int(last - first), "sufficient_support": int(allowed)}
               for first, last, allowed in zip(starts, stops, supported, strict=True)]
    endpoint = float(np.sum(times[stops[supported] - 1] - times[starts[supported]])
                     / np.sum(sizes[supported] - 1) / prescaler)
    result: ClockAssessment = {
        "verdict": "PASS_GAP_AWARE_CH0_FOLD_WITH_EXPLICIT_UNRESOLVED_WINDOWS",
        "clock_model": "Unit-cadence contiguous regions; symmetric local periods; no inter-region tick inference",
        "windows": len(times), "calibrated_reference_windows": int(resolved.sum()),
        "unresolved_reference_windows": int((~resolved).sum()),
        "prescaled_ticks": int(np.sum(sizes[supported] - 1)), "missed_prescaled_ticks": 0,
        "clock_tick_scope": "Supported contiguous regions only; missing ticks across gaps are unknown, not zero.",
        "unbridged_clock_gaps": len(starts) - 1, "clock_gap_tolerance_ns": gap_tolerance_ns,
        "clock_regions": regions, "clock_support_reason_counts": excluded,
        "prescaler": prescaler, "half_span": half_span, "assumed_ch0_edge_precision_ns": edge_precision_ns,
        "mean_period_ns": float(model.periods[resolved].mean()), "endpoint_period_ns": endpoint,
        "min_local_period_ns": float(model.periods[resolved].min()),
        "max_local_period_ns": float(model.periods[resolved].max()),
        "prescaled_interval_residual_rms_ns": float(np.sqrt(np.mean(residuals ** 2))),
        "prescaled_interval_residual_p999_abs_ns": float(np.quantile(abs(residuals), .999)),
        "max_half_window_phase_disagreement_ns": maximum_disagreement,
        "max_3sigma_conditional_phase_error_ns": maximum_statistical_phase,
        "max_folded_delay_ns": max_delay_ns, "phase_tolerance_ns": phase_tolerance_ns,
        "interval_span_sensitivity": sensitivity,
        "phase_precision_status": "INSTANTANEOUS_RF_PHASE_NOT_CERTIFIED_BETWEEN_PRESCALED_EDGES",
        "limitations": [
            "Opt-in regional variant; the original single-clock rejection and default settings are unchanged.",
            "A gap is any missing reference ID or interval differing from the dominant interval by more than one declared coarse quantum.",
            "At least 32 regional references and full symmetric 2x averaging-span support are required.",
            "Half-interval disagreement, averaging-span sensitivity and assumed edge error must each meet the phase budget.",
            "Unsupported references and all their pulses remain counted but are not folded.",
            "The endpoint estimate combines supported regional endpoints, not the gap-contaminated wall-clock span.",
            "Clock scale, physical divider and unobserved intra-interval RF phase remain unverified assumptions.",
        ],
    }
    return model, result


def calibrate_burst_clock(times: np.ndarray, prescaler: int, max_delay_ns: float,
                          phase_tolerance_ns: float, edge_precision_ns: float = .1,
                          window_ids: np.ndarray | None = None,
                          coarse_quantum_ns: float = 5) -> tuple[ClockModel, ClockAssessment]:
    require(times.ndim == 1 and times.size >= 32 and np.isfinite(times).all()
            and prescaler > 0, "Invalid burst-clock timeline or prescaler")
    require(window_ids is not None and window_ids.shape == times.shape
            and np.all(np.diff(window_ids) > 0), "Burst clock requires original reference IDs")
    require(all(np.isfinite(value) and value > 0 for value in
                (max_delay_ns, phase_tolerance_ns, edge_precision_ns, coarse_quantum_ns)),
            "Invalid burst-clock precision requirements")
    delta = np.diff(times)
    require(np.all(delta > 0), "Ch0 timestamps must increase")
    nominal = float(np.median(delta))
    require(35 <= nominal / prescaler <= 45, "Dominant cadence outside 35-45 ns")
    connected = (np.diff(window_ids) == 1) & (abs(delta - nominal) <= coarse_quantum_ns)
    noncore = (np.diff(window_ids) == 1) & ~connected & (delta < 5 * nominal)
    transitions = np.diff(np.r_[False, noncore, False].astype(np.int8))
    noncore_lengths = np.flatnonzero(transitions == -1) - np.flatnonzero(transitions == 1)
    longest_noncore = int(noncore_lengths.max()) if noncore_lengths.size else 0
    require(longest_noncore < 15,
            "Shared RF-period model fails regional stability budget: sustained non-core cadence")
    starts = np.r_[0, np.flatnonzero(~connected) + 1]
    stops = np.r_[starts[1:], len(times)]
    sizes = stops - starts
    region = np.repeat(np.arange(len(starts)), sizes)
    x = np.arange(len(times)) - starts[region] - (sizes[region] - 1) / 2
    y = times - times[starts[region]]
    sxx = np.add.reduceat(x * x, starts)
    sxy = np.add.reduceat(x * y, starts)
    slopes = np.divide(sxy, sxx, out=np.full(len(starts), np.nan), where=sxx > 0)
    means = np.add.reduceat(y, starts) / sizes
    residuals = y - means[region] - slopes[region] * x
    sse = np.add.reduceat(residuals ** 2, starts)
    fit = sizes >= 32
    validation = sizes >= 16
    require(np.any(fit), "Insufficient within-burst support for shared RF period")
    period = float(sxy[fit].sum() / sxx[fit].sum() / prescaler)
    require(35 <= period <= 45, "Shared RF period outside 35-45 ns")
    variance = np.divide(sse, sizes - 2, out=np.full(len(starts), np.nan), where=sizes > 2)
    observed_errors = np.sqrt(variance[validation] / sxx[validation]) / prescaler
    assumed_errors = edge_precision_ns / np.sqrt(sxx[validation]) / prescaler
    regional_errors = np.maximum(observed_errors, assumed_errors)
    differences = max_delay_ns * abs(slopes[validation] / prescaler - period) / period
    bounds = differences + max_delay_ns * 3 * regional_errors / period
    require(np.isfinite(bounds).all() and float(bounds.max()) <= phase_tolerance_ns,
            f"Shared RF-period model fails regional stability budget: {float(bounds.max()):.6g} ns")
    pooled_variance = float(sse[fit].sum() / (sizes[fit] - 2).sum())
    error = max(np.sqrt(pooled_variance), edge_precision_ns) / np.sqrt(sxx[fit].sum()) / prescaler
    model = ClockModel(np.full(len(times), period), np.full(len(times), error),
                       np.where(connected, delta - prescaler * period, np.nan), [])
    sensitivity: list[Sensitivity] = []
    for minimum in (16, 24, 32):
        selected = sizes >= minimum
        alternative = float(sxy[selected].sum() / sxx[selected].sum() / prescaler)
        change = max_delay_ns * abs(alternative - period) / period
        require(change <= phase_tolerance_ns, "Shared RF period depends on calibration support threshold")
        sensitivity.append({"half_span": minimum, "max_phase_change_ns": change,
                            "p999_phase_change_ns": change})
    regions = []
    for index, (first, last) in enumerate(zip(starts, stops, strict=True)):
        regions.append({
            "first_reference_window": int(window_ids[first]),
            "last_reference_window": int(window_ids[last - 1]),
            "references": int(sizes[index]),
            "used_for_period_fit": int(fit[index]),
            "used_for_stability_validation": int(validation[index]),
            "regional_period_ns": float(slopes[index] / prescaler),
        })
    for first in range(0, len(times), 65536):
        last = min(len(times), first + 65536)
        model.blocks.append({"first_window": int(window_ids[first]),
                             "last_window": int(window_ids[last - 1]),
                             "mean_local_period_ns": period, "min_local_period_ns": period,
                             "max_local_period_ns": period})
    valid_residuals = model.residuals[np.isfinite(model.residuals)]
    require(valid_residuals.size > 0, "No measured within-burst clock intervals")
    endpoint = float(np.sum(times[stops[fit] - 1] - times[starts[fit]])
                     / np.sum(sizes[fit] - 1) / prescaler)
    result: ClockAssessment = {
        "verdict": "PASS_SHARED_BURST_RF_PERIOD_ALL_REFERENCED_WINDOWS",
        "burst_clock": True,
        "clock_model": "Shared RF period from within-region least squares with independent intercepts",
        "windows": len(times), "calibrated_reference_windows": len(times),
        "unresolved_reference_windows": 0, "fit_reference_windows": int(sizes[fit].sum()),
        "validation_regions": int(validation.sum()),
        "maximum_regional_phase_bound_ns": float(bounds.max()),
        "minimum_fit_references": 32, "minimum_validation_references": 16,
        "prescaled_ticks": int((sizes[fit] - 1).sum()), "missed_prescaled_ticks": 0,
        "clock_tick_scope": "Measured within-region intervals only; inter-region cycle counts unknown.",
        "unbridged_clock_gaps": len(starts) - 1,
        "clock_interval_bounds_ns": [35 * prescaler, 45 * prescaler],
        "clock_gap_tolerance_ns": coarse_quantum_ns,
        "noncore_clock_intervals": int(noncore.sum()),
        "longest_noncore_interval_run": longest_noncore,
        "clock_regions": regions,
        "clock_support_reason_counts": {"resolved": len(times), "unsupported_or_boundary": 0,
                                       "phase_consistency": 0},
        "prescaler": prescaler, "half_span": 0,
        "assumed_ch0_edge_precision_ns": edge_precision_ns,
        "mean_period_ns": period, "endpoint_period_ns": endpoint,
        "min_local_period_ns": period, "max_local_period_ns": period,
        "prescaled_interval_residual_rms_ns": float(np.sqrt(np.mean(valid_residuals ** 2))),
        "prescaled_interval_residual_p999_abs_ns": float(np.quantile(abs(valid_residuals), .999)),
        "max_half_window_phase_disagreement_ns": float(differences.max()),
        "max_3sigma_conditional_phase_error_ns": float(max_delay_ns * 3 * error / period),
        "max_folded_delay_ns": max_delay_ns, "phase_tolerance_ns": phase_tolerance_ns,
        "interval_span_sensitivity": sensitivity,
        "phase_precision_status": "SHARED_PERIOD_ASSUMPTION_VALIDATED;INSTANTANEOUS_RF_PHASE_NOT_CERTIFIED",
        "limitations": [
            "No cycle counting or phase propagation through inter-burst gaps.",
            "Every passing ch0 reference anchors its own window; all such windows receive the shared period.",
            "Minimum regional lengths select calibration/validation evidence, never folding populations.",
            "Non-core timing intervals are not fitted; their windows/hits are retained under the shared-period assumption.",
            "A sustained non-core cadence of 15 intervals rejects the model rather than excluding its windows.",
            "Short/isolated regions use the explicitly assumed shared period, not an independently measured local period.",
            "Stability bounds are conditional diagnostics, not a measurement of instantaneous intra-cycle RF phase.",
            "Windows without passing ch0 remain in cleaned data and window/pulse trees with invalid phase.",
            "Clock scale and physical prescaler remain declared hardware assumptions.",
        ],
    }
    return model, result


def clock_plots(folder: Path, board: str, times: np.ndarray, model: ClockModel,
                result: ClockAssessment) -> None:
    specifications = [
        ("period_vs_time", "Local period averages from ch0", "Recorded time (s)", "Period (ns)"),
        ("clock_residuals", "Prescaled intervals minus local mean model", "Residual (ns)", "Intervals"),
        ("prescaled_intervals", "Consecutive selected ch0 intervals", "Interval (ns)", "Intervals"),
        ("interval_span_sensitivity", "Averaging-span sensitivity at maximum folded delay", "Neighboring intervals on each side", "Maximum phase change (ns)"),
    ]
    for name, title, xlabel, ylabel in specifications:
        if name == "interval_span_sensitivity" and result.get("burst_clock"):
            title = "Shared-period sensitivity to calibration-region length"
            xlabel = "Minimum references per calibration region"
        fig, ax = plt.subplots(figsize=(9, 4.5))
        if name == "period_vs_time":
            step = max(1, times.size // 10000)
            ax.plot(times[::step] * 1e-9, model.periods[::step], linewidth=.7)
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        elif name == "clock_residuals":
            ax.hist(model.residuals[np.isfinite(model.residuals)], bins=160)
        elif name == "prescaled_intervals":
            intervals = np.diff(times)
            ax.hist(intervals, bins=160)
            ax.ticklabel_format(axis="x", style="plain", useOffset=False)
        else:
            sensitivity = result["interval_span_sensitivity"]
            ax.plot([row["half_span"] for row in sensitivity],
                    [row["max_phase_change_ns"] for row in sensitivity], "o-")
            ax.axhline(float(result["phase_tolerance_ns"]), color="red", linestyle="--", label="100 ps edge-error reference; not a source-phase certificate")
            ax.legend()
        ax.set(title=f"{board}: {title}", xlabel=xlabel, ylabel=ylabel)
        fig.tight_layout()
        for suffix in ("png", "pdf"):
            fig.savefig(folder / f"{board}_{name}.{suffix}", dpi=150)
        plt.close(fig)


def write_folded_outputs(tables: Path, folder: Path, board: str = "Zebra",
                         tot_display_max_ns: float | None = None) -> None:
    for factor, profile_name, group in ((1, "folded_matrix.tsv", "FoldedRF"), (3, "folded3_matrix.tsv", "FoldedRF3x")):
        output = folder / ("Folded_RF" if factor == 1 else "Folded_RF_3x")
        output.mkdir()
        suffix = "folded_rf" if factor == 1 else "folded_rf_3x"
        root_suffix = "Folded_RF" if factor == 1 else "Folded_RF_3x"
        write_profile_outputs(
            tables / profile_name, output, combined_png_name=f"combined_{suffix}.png",
            combined_pdf_name=f"combined_{suffix}.pdf", per_channel_suffix=suffix,
            root_name_suffix=root_suffix, color="#1f77b4", cmap="viridis",
            title_prefix=f"{board} ch0-clock phase | {factor}x | zero = selected ch0",
            x_label="Local RF cycles scaled by board reference period (ns)", root_dir_name=group,
        )
        for quantity, label, hist_group in (
            ("tot", "ToT (ns)", "ToT"), ("time", "Time relative to selected ch0 (ns)", "Ch0RefTime"),
            ("trigger", "Time relative to trigger window (ns)", "TriggerTime"),
        ):
            prefix = "phase" if factor == 1 else "phase3"
            stem = f"rf_{prefix}_vs_{quantity}"
            write_sparse_2d_outputs(
                tables / f"{prefix}_{quantity}.tsv", output,
                combined_filename_suffix=stem, per_channel_filename_suffix=stem,
                root_name_suffix=stem, hist_group=f"RFPhase{'' if factor == 1 else '3x'}{hist_group}",
                x_label="Local RF cycles scaled by board reference period (ns)", y_label=label,
                title_prefix=f"{board} ch0-clock phase vs {quantity} | {factor}x | zero = selected ch0",
                display_y_max_ns=tot_display_max_ns if quantity == "tot" else None,
            )


def verify_folded(tables: Path, expected_input: np.ndarray, windows: int,
                  include_all_ch0: bool = False) -> list[dict[str, float | int]]:
    expected = expected_input.copy()
    expected[0] = windows
    with (tables / "channel_counts.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(len(rows) == 32, "Missing channel accounting")
    missing = np.asarray([int(row["without_reference_pulses"]) for row in rows])
    unresolved = np.asarray([int(row["without_clock_pulses"]) for row in rows])
    expected[1:] -= missing[1:] + unresolved[1:]
    if include_all_ch0:
        expected[0] = expected_input[0] - missing[0] - unresolved[0]
    for ch, row in enumerate(rows):
        require(int(row["channel"]) == ch and int(row["input_pulses"]) == expected_input[ch]
                and int(row["folded_pulses"]) == expected[ch], f"Ch{ch}: input/folded count mismatch")
        require(int(row["folded_pulses"]) + int(row["without_reference_pulses"])
                + int(row["nonreference_ch0_pulses"]) + int(row["without_clock_pulses"]) == int(row["input_pulses"]),
                f"Ch{ch}: excluded-pulse accounting mismatch")
        require(ch == 0 or int(row["nonreference_ch0_pulses"]) == 0,
                f"Ch{ch}: unexpected nonreference-ch0 exclusion")
    for factor, filename in ((1, "folded_matrix.tsv"), (3, "folded3_matrix.tsv")):
        _, _, counts, total = read_matrix(tables / filename)
        require(np.array_equal(counts.sum(axis=0), expected), f"{factor}x: folded profile count mismatch")
        require(np.array_equal(counts.sum(axis=1), total), f"{factor}x: folded row totals mismatch")
        for quantity in ("tot", "time", "trigger"):
            prefix = "phase" if factor == 1 else "phase3"
            metadata, entries, _, _ = read_sparse_histogram(tables / f"{prefix}_{quantity}.tsv")
            detailed = {int(value) for value in metadata["stored_channels"].split(",")}
            require(set(entries).issubset(detailed), "Unexpected detailed channel set")
            for ch in detailed:
                cells = entries.get(ch, [])
                target = expected[ch] - (int(rows[ch]["tot_outside_display"]) if quantity == "tot" else 0)
                require(sum(cell[2] for cell in cells) == target, f"{factor}x {quantity} ch{ch}: lost pulses")
                require(all(0 <= x < int(metadata["x_bins"]) and 0 <= y < int(metadata["y_bins"])
                            and count > 0 for x, y, count in cells), "Invalid sparse bin")
    return [{key: int(value) for key, value in row.items()} for row in rows]


def verify_products(tables: Path, folder: Path) -> None:
    for factor, profile_name, group in ((1, "folded_matrix.tsv", "FoldedRF"), (3, "folded3_matrix.tsv", "FoldedRF3x")):
        directory = folder / ("Folded_RF" if factor == 1 else "Folded_RF_3x")
        metadata, _, counts, total = read_matrix(tables / profile_name)
        key = metadata["run_key"]
        suffix = "Folded_RF" if factor == 1 else "Folded_RF_3x"
        with uproot.open(directory / f"{key}_{suffix}.root") as root:
            require(np.array_equal(root[f"{group}/hTotal"].values(), total), "ROOT total profile mismatch")
            for ch in np.flatnonzero(counts.sum(axis=0)):
                require(np.array_equal(root[f"{group}/hCh{ch:02d}"].values(), counts[:, ch]), "ROOT profile bins mismatch")
        for quantity, hist_group in (("tot", "ToT"), ("time", "Ch0RefTime"), ("trigger", "TriggerTime")):
            prefix = "phase" if factor == 1 else "phase3"
            meta, entries, aggregate, _ = read_sparse_histogram(tables / f"{prefix}_{quantity}.tsv")
            group_name = f"RFPhase{'' if factor == 1 else '3x'}{hist_group}"
            with uproot.open(directory / f"{key}_rf_{prefix}_vs_{quantity}.root") as root:
                for ch, cells in entries.items():
                    require(np.array_equal(root[f"{group_name}/hCh{ch:02d}"].values().T,
                                           build_dense(cells, int(meta["x_bins"]), int(meta["y_bins"]))),
                            "ROOT 2D bin mismatch")
                summed = build_dense([(x, y, count) for (x, y), count in aggregate.items()],
                                     int(meta["x_bins"]), int(meta["y_bins"]))
                require(np.array_equal(root[f"{group_name}/hAllChannels"].values().T, summed), "ROOT aggregate mismatch")
    for path in folder.rglob("*.png"):
        with Image.open(path) as image:
            image.verify()
    for path in folder.rglob("*.pdf"):
        with path.open("rb") as handle:
            require(handle.read(5) == b"%PDF-", f"Invalid PDF: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Parallel ch0-only clock calibration and folding for explicit decoder/board settings.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--decoded", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--board", required=True)
    parser.add_argument("--tdc-ordinal", type=int)
    scale = parser.add_mutually_exclusive_group(required=True)
    scale.add_argument("--tdc-clock-mhz", type=float)
    scale.add_argument("--coarse-quantum-ns", type=float)
    parser.add_argument("--prescaler", type=int, required=True)
    parser.add_argument("--clock-half-span", type=int, default=4)
    parser.add_argument("--ch0-edge-precision-ns", type=float, default=.1)
    parser.add_argument("--allow-reference-gaps", action="store_true")
    parser.add_argument("--expected-windows", type=int)
    parser.add_argument("--run-key")
    parser.add_argument("--ch0-valid-rise-min-ns", type=float, default=-392)
    parser.add_argument("--ch0-valid-rise-max-ns", type=float, default=-374)
    parser.add_argument("--ch0-valid-tot-min-ns", type=float, default=16.5)
    parser.add_argument("--ch0-valid-tot-max-ns", type=float, default=19.5)
    parser.add_argument("--trigger-ref-min-ns", type=float, default=-11000)
    parser.add_argument("--trigger-ref-max-ns", type=float, default=2000)
    parser.add_argument("--trigger-ref-bins", type=int, default=2600)
    parser.add_argument("--ch0ref-min-ns", type=float, default=-11000)
    parser.add_argument("--ch0ref-max-ns", type=float, default=2000)
    parser.add_argument("--ch0ref-bins", type=int, default=2600)
    parser.add_argument("--tot-min-ns", type=float, default=0)
    parser.add_argument("--tot-max-ns", type=float, default=80)
    parser.add_argument("--tot-bins", type=int, default=800)
    parser.add_argument("--phase-bin-width-ns", type=float, default=.25)
    parser.add_argument("--cadence-factors", help="Explicit frequency-rate factors relative to the dominant cadence; e.g. 1,2,4")
    parser.add_argument("--split-clock-gaps", action="store_true",
                        help="Opt-in unit-cadence regions; never infer ticks across acquisition gaps")
    parser.add_argument("--burst-clock", action="store_true",
                        help="Shared period fitted inside bursts; fold every referenced window and every ch0 pulse")
    parser.add_argument("--calibration-only", action="store_true",
                        help="Save reference inventory and calibration without folding")
    parser.add_argument("--tot-display-max-ns", type=float)
    parser.add_argument("--histogram-channels", default="0,1,2,3,6")
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--build-dir", type=Path, required=True)
    args = record_arguments(parser.parse_args())
    require(sum((args.split_clock_gaps, args.burst_clock, args.cadence_factors is not None)) <= 1,
            "Gap-aware, burst-clock and multi-cadence modes are mutually exclusive")
    require(not args.split_clock_gaps or args.allow_reference_gaps,
            "Gap-aware mode requires explicit unresolved-window accounting")
    require(not args.burst_clock or args.allow_reference_gaps,
            "Burst-clock mode requires explicit reference-gap accounting")
    if args.coarse_quantum_ns is None:
        require(np.isfinite(args.tdc_clock_mhz) and args.tdc_clock_mhz > 0, "Invalid decoder clock")
    quantum = args.coarse_quantum_ns if args.coarse_quantum_ns is not None else 2000 / args.tdc_clock_mhz
    require(np.isfinite(quantum) and 0 < quantum < 100, "Invalid decoder coarse quantum")
    ordinal = args.tdc_ordinal
    if ordinal is None:
        require(args.board in ("AE", "AF"), "Non-Zebra boards require an explicit --tdc-ordinal")
        ordinal = 1 if args.board == "AE" else 2
    require(ordinal >= 0, "Invalid TDC ordinal")
    for low, high, bins in ((args.trigger_ref_min_ns, args.trigger_ref_max_ns, args.trigger_ref_bins),
                            (args.ch0ref_min_ns, args.ch0ref_max_ns, args.ch0ref_bins),
                            (args.tot_min_ns, args.tot_max_ns, args.tot_bins)):
        require(np.isfinite(low) and np.isfinite(high) and low < high and bins > 0,
                "Invalid histogram bounds")
    detailed = [int(value) for value in args.histogram_channels.split(",")]
    require(len(detailed) == len(set(detailed)) and bool(detailed)
            and all(0 <= ch < 32 for ch in detailed), "Invalid detailed channels")
    require(not args.output.exists(), f"Refusing to replace output: {args.output}")
    args.output.mkdir(parents=True)
    args.build_dir.mkdir(parents=True, exist_ok=True)
    inputs = [args.input.resolve(), args.decoded.resolve()]
    before = {str(path): (path.stat().st_size, path.stat().st_mtime_ns) for path in inputs}
    source = Path(__file__).with_name("zebra_ch0_clock.cpp")
    executable = args.build_dir / "zebra_ch0_clock"
    compile_cpp(source, executable, args.compiler)
    clock_dir = args.output / "Clock_Calibration"
    clock_dir.mkdir()
    references = clock_dir / "ch0_reference_timestamps.tsv"
    gate_options = [
        "--ch0-valid-rise-min-ns", str(args.ch0_valid_rise_min_ns),
        "--ch0-valid-rise-max-ns", str(args.ch0_valid_rise_max_ns),
        "--ch0-valid-tot-min-ns", str(args.ch0_valid_tot_min_ns),
        "--ch0-valid-tot-max-ns", str(args.ch0_valid_tot_max_ns),
        "--allow-reference-gaps", "1" if args.allow_reference_gaps else "0",
    ]
    window_options = [] if args.expected_windows is None else ["--expected-windows", str(args.expected_windows)]
    run_command([str(executable), "extract", "--input", str(inputs[0]), "--decoded", str(inputs[1]),
                 "--ordinal", str(ordinal), "--hardware-id", args.board,
                 "--coarse-quantum-ns", str(quantum), "--output", str(references),
                 *gate_options, *window_options])
    data = np.loadtxt(references, comments="#", usecols=(0, 4, 5))
    window_ids = data[:, 0].astype(np.int64)
    require(np.all(data[:, 0] == window_ids) and np.all(np.diff(window_ids) > 0),
            "Reference sequence mismatch")
    if not args.allow_reference_gaps:
        require(np.array_equal(window_ids, np.arange(1, data.shape[0] + 1)), "Reference sequence mismatch")
    require(args.burst_clock or np.all(data[:, 2] == 1),
            "Multiple valid ch0 candidates require explicit assessment")
    times = np.ascontiguousarray(data[:, 1])
    del data
    max_delay = max(abs(args.ch0ref_min_ns), abs(args.ch0ref_max_ns))
    factors = None if args.cadence_factors is None else [int(value) for value in args.cadence_factors.split(",")]
    if args.burst_clock:
        model, result = calibrate_burst_clock(times, args.prescaler, max_delay, .1,
                                             args.ch0_edge_precision_ns, window_ids, quantum)
    elif args.split_clock_gaps:
        model, result = calibrate_gap_clock(times, args.prescaler, max_delay, .1, quantum,
                                           args.clock_half_span, args.ch0_edge_precision_ns, window_ids)
    elif factors is None:
        model, result = calibrate_clock(times, args.prescaler, max_delay, .1,
                                       args.clock_half_span, args.ch0_edge_precision_ns, window_ids)
    else:
        require(args.allow_reference_gaps, "Cadence mode requires explicit unresolved-window accounting")
        model, result = calibrate_cadence_clock(times, args.prescaler, max_delay, .1, factors,
                                               args.clock_half_span, args.ch0_edge_precision_ns, window_ids)
    reference_metadata = parse_metadata(references)
    result.update({"board": args.board, "tdc_clock_mhz": args.tdc_clock_mhz,
                   "tdc_ordinal": ordinal, "allow_reference_gaps": args.allow_reference_gaps,
                   "coarse_quantum_ns": quantum,
                   "reference_gate_rise_ns": [args.ch0_valid_rise_min_ns, args.ch0_valid_rise_max_ns],
                   "reference_gate_tot_ns": [args.ch0_valid_tot_min_ns, args.ch0_valid_tot_max_ns],
                   "phase_origin": "first passing ch0 in file order; no detector-derived alignment",
                   "mean_frequency_mhz": 1000 / (result["mean_period_ns"] if args.burst_clock
                                                else result["endpoint_period_ns"]),
                   "include_all_ch0_hits": args.burst_clock,
                   "detailed_channels": detailed,
                   "timing_axes_ns": {"ch0": [args.ch0ref_min_ns, args.ch0ref_max_ns],
                                      "trigger": [args.trigger_ref_min_ns, args.trigger_ref_max_ns]},
                   "tot_axis_ns": [args.tot_min_ns, args.tot_max_ns],
                   "decoded_selected_board_windows": int(reference_metadata["decoded_selected_board_windows"]),
                   "decoded_board0_headers": int(reference_metadata["decoded_board0_headers"]),
                   "represented_cleaned_windows": int(reference_metadata["represented_cleaned_windows"]),
                   "reference_gap_windows": int(reference_metadata["reference_gap_windows"]),
                   "empty_cleaned_windows": int(reference_metadata["empty_cleaned_windows"]),
                   "input_state": before,
                   "source_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                     for path in (source, Path(__file__))}})
    write_table(clock_dir / "clock_blocks.tsv", model.blocks)
    write_table(clock_dir / "interval_span_sensitivity.tsv", result["interval_span_sensitivity"])
    resolved = np.isfinite(model.periods)
    unresolved_file = clock_dir / "unresolved_clock_windows.tsv"
    with unresolved_file.open("w") as handle:
        handle.write("# window_index; cadence transition or fewer than 32 same-cadence references\n")
        for index in window_ids[~resolved]:
            handle.write(f"{index}\n")
    if factors is not None:
        write_table(clock_dir / "cadence_regions.tsv", result["cadence_regions"])
        write_table(clock_dir / "cadence_regimes.tsv", result["cadence_regimes"])
    if args.split_clock_gaps or args.burst_clock:
        write_table(clock_dir / "clock_regions.tsv", result["clock_regions"])
        del result["clock_regions"]
    period_file = clock_dir / "window_periods.tsv"
    with period_file.open("w", encoding="utf-8") as handle:
        handle.write("# columns=window_index local_period_ns\n")
        for index, value in zip(window_ids, model.periods, strict=True):
            if np.isfinite(value):
                handle.write(f"{index}\t{value:.12f}\n")
    if args.calibration_only:
        with (args.output / "assessment.json").open("w") as handle:
            json.dump(result, handle, indent=2)
        return
    clock_plots(clock_dir, args.board, times, model, result)
    tables = args.output / "Fold_Tables"
    key = args.run_key or parse_metadata(inputs[0])["run_key"].removesuffix("_rawRefined")
    unresolved_options = ([] if factors is None and not args.split_clock_gaps else
                          ["--unresolved-clock-windows", str(unresolved_file)])
    run_command([str(executable), "fold", "--input", str(inputs[0]), "--periods", str(period_file),
                 "--output", str(tables), "--run-key", key, "--ordinal", str(ordinal),
                 *gate_options, *window_options,
                 "--trigger-ref-min-ns", str(args.trigger_ref_min_ns),
                 "--trigger-ref-max-ns", str(args.trigger_ref_max_ns),
                 "--trigger-ref-bins", str(args.trigger_ref_bins),
                 "--ch0ref-min-ns", str(args.ch0ref_min_ns), "--ch0ref-max-ns", str(args.ch0ref_max_ns),
                 "--ch0ref-bins", str(args.ch0ref_bins), "--tot-min-ns", str(args.tot_min_ns),
                 "--tot-max-ns", str(args.tot_max_ns), "--tot-bins", str(args.tot_bins),
                 "--phase-bin-width-ns", str(args.phase_bin_width_ns),
                 "--period-min-ns", str(35 if factors is None else 35 / max(factors)),
                 "--period-max-ns", "45", *unresolved_options,
                 "--include-all-ch0-hits", "1" if args.burst_clock else "0",
                 "--period-model", "shared_burst_rf_period" if args.burst_clock
                 else "independent_ch0_neighbor_interval_average",
                 "--histogram-channels", args.histogram_channels])
    expected = np.asarray([int(value) for value in parse_metadata(references)["input_channel_counts"].split(",")])
    accounting = verify_folded(tables, expected, int(resolved.sum()), args.burst_clock)
    write_table(args.output / "channel_count_verification.tsv", accounting)
    write_folded_outputs(tables, args.output, args.board, args.tot_display_max_ns)
    verify_products(tables, args.output)
    require(before == {str(path): (path.stat().st_size, path.stat().st_mtime_ns) for path in inputs},
            "Input files changed during processing")
    result["folded_output_verification"] = "PASS: all 32 profile counts, detailed sparse/ROOT bins, PNG/PDF integrity"
    result["folded_channel_counts"] = [int(row["folded_pulses"]) for row in accounting]
    result["unreferenced_channel_counts"] = [int(row["without_reference_pulses"]) for row in accounting]
    result["unresolved_clock_channel_counts"] = [int(row["without_clock_pulses"]) for row in accounting]
    result["png_files"] = len(list(args.output.rglob("*.png")))
    result["pdf_files"] = len(list(args.output.rglob("*.pdf")))
    with (args.output / "assessment.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
