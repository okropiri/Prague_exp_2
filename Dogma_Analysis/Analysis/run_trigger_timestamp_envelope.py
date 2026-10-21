from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import socket

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(chunk)
    return sha.hexdigest()


def quantiles(values: np.ndarray) -> dict[str, float]:
    require(values.size > 0 and bool(np.all(np.isfinite(values))), "Empty/nonfinite measurement")
    return {name: float(value) for name, value in zip(
        ("min", "q001", "q01", "q16", "median", "q84", "q99", "q999", "max"),
        np.quantile(values, (0, .001, .01, .16, .5, .84, .99, .999, 1)), strict=True)}


def burst_groups(ticks: np.ndarray, gap_ns: float) -> tuple[np.ndarray, np.ndarray]:
    require(ticks.ndim == 1 and len(ticks) >= 2 and bool(np.all(ticks[1:] > ticks[:-1])),
            "Trigger timestamps must strictly increase")
    require(np.isfinite(gap_ns) and gap_ns > 0, "Invalid grouping threshold")
    starts = np.r_[0, np.flatnonzero(np.diff(ticks) * 5 > gap_ns) + 1]
    stops = np.r_[starts[1:], len(ticks)]
    require(int(np.sum(stops - starts)) == len(ticks), "Grouping lost trigger windows")
    return starts, stops


def fit_repetition(start_ns: np.ndarray) -> dict[str, float | int | dict[str, float]]:
    require(len(start_ns) >= 20 and bool(np.all(np.diff(start_ns) > 0)), "Insufficient burst starts")
    delta = np.diff(start_ns)
    centre = float(np.median(delta))
    connected = abs(delta / centre - 1) <= .1
    first = np.r_[0, np.flatnonzero(~connected) + 1]
    stop = np.r_[first[1:], len(start_ns)]
    first = np.concatenate([np.arange(low, high, 256) for low, high in zip(first, stop, strict=True)])
    stop = np.r_[first[1:], len(start_ns)]
    lengths = stop - first
    supported = lengths >= 20
    require(bool(np.any(supported)), "No uninterrupted burst-start sequence with at least 20 bursts")
    group = np.repeat(np.arange(len(first)), lengths)
    positions = np.arange(len(start_ns)) - first[group]
    centred_x = positions - (lengths[group] - 1) / 2
    elapsed = start_ns - start_ns[first[group]]
    means = np.add.reduceat(elapsed, first) / lengths
    centred_y = elapsed - means[group]
    mask = supported[group]
    period = float(np.sum(centred_x[mask] * centred_y[mask]) / np.sum(centred_x[mask] ** 2))
    endpoint = float(np.sum(start_ns[stop[supported] - 1] - start_ns[first[supported]])
                     / np.sum(lengths[supported] - 1))
    residual = centred_y[mask] - period * centred_x[mask]
    numerator = np.add.reduceat(centred_x * centred_y, first)[supported]
    denominator = np.add.reduceat(centred_x ** 2, first)[supported]
    return {"median_start_spacing_ns": centre, "fitted_period_ns": period,
            "frequency_hz": 1e9 / period, "regional_endpoint_period_ns": endpoint,
            "fitted_bursts": int(mask.sum()), "fit_regions": int(supported.sum()),
            "fit_region_maximum_bursts": 256,
            "block_frequency_hz": quantiles(1e9 / (numerator / denominator)),
            "excluded_start_intervals": int((~connected).sum()),
            "start_quantization_residual_rms_ns": float(np.sqrt(np.mean(residual ** 2))),
            "maximum_abs_fit_residual_ns": float(abs(residual).max())}


def measure(ticks: np.ndarray, gap_ns: float, minimum_dense_triggers: int = 8) -> tuple[dict, np.ndarray, np.ndarray]:
    starts, stops = burst_groups(ticks, gap_ns)
    all_counts = stops - starts
    dense = all_counts >= minimum_dense_triggers
    sparse_triggers = int(all_counts[~dense].sum())
    starts, stops = starts[dense], stops[dense]
    intervals = np.diff(ticks).astype(np.int64) * 5
    within = intervals[intervals <= gap_ns]
    spacing = float(np.median(within))
    start_ns = (ticks[starts] - ticks[0]).astype(np.float64) * 5
    spans = (ticks[stops - 1] - ticks[starts]).astype(np.float64) * 5
    counts = stops - starts
    repeat = np.diff(start_ns)
    fit = fit_repetition(start_ns)
    regular = abs(repeat / float(fit["median_start_spacing_ns"]) - 1) <= .1
    quiet = (ticks[starts[1:]] - ticks[stops[:-1] - 1]).astype(np.float64) * 5
    occupancy = (spans[:-1][regular] + spacing) / repeat[regular]
    values, numbers = np.unique(counts, return_counts=True)
    all_values, all_numbers = np.unique(all_counts, return_counts=True)
    result = {
        "grouping_gap_threshold_ns": gap_ns, "trigger_windows": len(ticks), "bursts": len(starts),
        "minimum_dense_group_triggers": minimum_dense_triggers,
        "all_group_count": len(all_counts), "sparse_group_count": int((~dense).sum()),
        "dense_group_trigger_windows": int(counts.sum()), "sparse_group_trigger_windows": sparse_triggers,
        "all_trigger_accounting_verified": int(counts.sum()) + sparse_triggers == len(ticks),
        "all_group_trigger_count_inventory": {str(value): int(number)
                                             for value, number in zip(all_values, all_numbers, strict=True)},
        "within_burst_interval_ns": quantiles(within.astype(np.float64)),
        "regular_burst_start_spacing_ns": quantiles(repeat[regular]),
        "all_burst_start_spacing_ns": quantiles(repeat),
        "regular_burst_observed_first_to_last_trigger_span_ns": quantiles(spans[:-1][regular]),
        "dense_burst_last_to_next_dense_first_interval_ns": quantiles(quiet[regular]),
        "regular_burst_trigger_count": quantiles(counts[:-1][regular].astype(np.float64)),
        "regular_burst_span_plus_one_trigger_interval_duty_proxy": quantiles(occupancy),
        "duty_proxy_definition": "(last_trigger - first_trigger + median_within_burst_spacing) / next_burst_start_spacing",
        "duty_proxy_limitation": "Trigger occupancy estimate, not an independently measured RF-on gate or beam-current duty.",
        "dense_burst_trigger_count_inventory": {str(value): int(number)
                                         for value, number in zip(values, numbers, strict=True)},
        "regular_start_intervals": int(regular.sum()), "nonregular_start_intervals": int((~regular).sum()),
        "gaps_above_100ms": int(np.count_nonzero(intervals > 1e8)),
        "largest_intertrigger_gap_ns": int(intervals.max()),
        "repetition_fit": fit,
    }
    return result, starts, stops


def save_figure(fig: plt.Figure, folder: Path, name: str) -> None:
    for extension in ("png", "pdf"):
        fig.savefig(folder / f"{name}.{extension}", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Unbiased native-trigger envelope measurements; no ch0 or pulse cuts.")
    parser.add_argument("--timestamps", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(socket.gethostname().split(".")[0] == "trbsrv-b", "DOGMA analysis runs only on padiwa")
    require(not args.output.exists(), "Refusing existing envelope output")
    args.output.mkdir(parents=True)
    ticks_path = args.timestamps / "trigger_ticks_5ns.bin"
    summary = dict(line.split("=", 1) for line in
                   (args.timestamps / "extraction_summary.txt").read_text().splitlines())
    expected = int(summary["trigger_timestamps"])
    require(ticks_path.stat().st_size == expected * 8 and summary["status"] == "PASS",
            "Unverified/incomplete timestamp extraction")
    ticks = np.memmap(ticks_path, dtype="<u8", mode="r", shape=(expected,))
    require(int(ticks[0]) == int(summary["first_tick"]) and int(ticks[-1]) == int(summary["last_tick"]),
            "Extraction endpoint mismatch")
    samples = 0
    with (args.timestamps / "timestamp_samples.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            actual = int((Decimal(row["original_seconds"]) * Decimal(200000000)).to_integral_value())
            require(actual == int(row["coarse_tick_5ns"])
                    == int(ticks[int(row["window_index"]) - 1]), "Independent Decimal timestamp mismatch")
            samples += 1
    measurements = []
    for threshold in (100000., 200000., 500000., 1000000.):
        measurement, starts, stops = measure(ticks, threshold)
        measurements.append(measurement)
    result, starts, stops = measure(ticks, 200000.)
    density_sensitivity = [measure(ticks, 200000., minimum)[0] for minimum in (4, 8, 16, 24)]
    require(all(item["all_trigger_accounting_verified"] for item in measurements + density_sensitivity),
            "Sparse/dense accounting lost windows")
    variants = measurements + density_sensitivity
    frequency_variation = max(abs(item["repetition_fit"]["frequency_hz"]
                                  / result["repetition_fit"]["frequency_hz"] - 1) for item in variants)
    duty_variation = max(abs(item["regular_burst_span_plus_one_trigger_interval_duty_proxy"]["median"]
                            - result["regular_burst_span_plus_one_trigger_interval_duty_proxy"]["median"])
                         for item in variants)
    require(frequency_variation <= .01 and duty_variation <= .01,
            "Dominant envelope is not stable under classification choices")
    before = (ticks_path.stat().st_size, ticks_path.stat().st_mtime_ns)
    time = (ticks - ticks[0]).astype(np.float64) * 5e-9
    first_seconds = int(ticks[0]) * 5e-9
    start_seconds = time[starts]
    spans_ns = (ticks[stops - 1] - ticks[starts]).astype(np.float64) * 5
    following_ns = np.diff(start_seconds) * 1e9
    segmentation = []
    for index in range(10):
        mask = (start_seconds >= time[-1] * index / 10) & (start_seconds < time[-1] * (index + 1) / 10)
        begin = int(np.searchsorted(time, time[-1] * index / 10))
        end = len(ticks) if index == 9 else int(np.searchsorted(time, time[-1] * (index + 1) / 10))
        local = measure(ticks[begin:end], 200000.)[0]
        fit = fit_repetition(start_seconds[mask] * 1e9)
        segmentation.append({"segment": index + 1, "start_elapsed_seconds": time[-1] * index / 10,
                             "stop_elapsed_seconds": time[-1] * (index + 1) / 10,
                             "trigger_windows": end - begin,
                             "frequency_hz": fit["frequency_hz"], "fitted_period_ns": fit["fitted_period_ns"],
                             "fitted_bursts": fit["fitted_bursts"],
                             "duty_proxy_median": local["regular_burst_span_plus_one_trigger_interval_duty_proxy"]["median"],
                             "trigger_count_median": local["regular_burst_trigger_count"]["median"]})
    require(sum(item["trigger_windows"] for item in segmentation) == expected,
            "Recording segments lost trigger windows")
    fig, panels = plt.subplots(4, 1, figsize=(12, 11), constrained_layout=True)
    edges = np.arange(0, time[-1] + .1, .1)
    count, _ = np.histogram(time, edges)
    require(int(count.sum()) == expected, "Absolute distribution lost timestamps")
    panels[0].stairs(count, edges + first_seconds)
    panels[0].set(title=f"All {expected:,} board1049 trigger windows; no hit/reference cuts",
                  xlabel="Absolute DAQ trigger timestamp (s; not UTC)", ylabel="Windows / 0.1 s")
    example_rows = []
    for ax, fraction, label in zip(panels[1:], (0., .5, .99), ("Early", "Middle", "Late"), strict=True):
        burst = min(np.searchsorted(start_seconds, time[-1] * fraction), len(starts) - 10)
        low = start_seconds[burst]
        begin, end = np.searchsorted(time, (low, low + .05))
        selected = time[begin:end]
        ax.vlines((selected - low) * 1000, 0, 1, linewidth=.45)
        ax.set(title=f"{label}: absolute timestamps {first_seconds + low:.6f} to {first_seconds + low + .05:.6f} s",
               xlabel="Milliseconds after labelled absolute start", ylabel="Recorded trigger", ylim=(0, 1.1))
        example_rows.extend({"example": label, "window_index": int(position + 1),
                             "absolute_seconds": f"{int(ticks[position]) * 5e-9:.9f}",
                             "offset_ms": f"{(time[position] - low) * 1000:.6f}"}
                            for position in range(begin, end))
    save_figure(fig, args.output, "absolute_trigger_timestamp_distribution")
    with (args.output / "absolute_timestamp_examples.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(example_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(example_rows)
    with (args.output / "absolute_timestamp_distribution.tsv").open("w") as handle:
        handle.write("absolute_left_seconds\tabsolute_right_seconds\ttrigger_windows\n")
        for low, high, value in zip(edges[:-1] + first_seconds, edges[1:] + first_seconds, count, strict=True):
            handle.write(f"{low:.9f}\t{high:.9f}\t{value}\n")
    interval_ns = np.diff(ticks).astype(np.int64) * 5
    regular = abs(following_ns / np.median(following_ns) - 1) <= .1
    fig, panels = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    panels[0, 0].hist(interval_ns * 1e-6, bins=np.geomspace(interval_ns.min() * 1e-6, interval_ns.max() * 1e-6, 200))
    panels[0, 0].set(xscale="log", yscale="log", xlabel="Consecutive trigger spacing (ms)", ylabel="Intervals",
                      title="Every consecutive recorded trigger interval")
    panels[0, 1].hist(interval_ns[interval_ns < 200000] / 1000, bins=np.linspace(30, 50, 201))
    panels[0, 1].set(xlabel="Within-burst trigger spacing (microseconds)", ylabel="Intervals",
                      title="Zoom; no nominal RF period imposed")
    panels[1, 0].hist(following_ns * 1e-6, bins=np.linspace(0, 20, 401))
    panels[1, 0].set(xlabel="Consecutive burst-start spacing (ms)", ylabel="Intervals",
                      title="Measured envelope repetition; no 150 Hz template imposed")
    panels[1, 1].hist((spans_ns[:-1][regular] + result["within_burst_interval_ns"]["median"]) * 1e-6,
                      bins=np.linspace(0, 2, 161))
    panels[1, 1].set(xlabel="Observed trigger span + one spacing (ms)", ylabel="Bursts",
                      title="Trigger-occupancy proxy, not a measured RF gate")
    save_figure(fig, args.output, "trigger_intervals_and_burst_structure")
    fig, panels = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True)
    panels[0].plot([item["segment"] for item in segmentation], [item["frequency_hz"] for item in segmentation], "o-")
    panels[0].set(xlabel="Equal-duration recording segment", ylabel="Envelope repetition (Hz)")
    panels[0].ticklabel_format(axis="y", useOffset=False)
    panels[1].plot([item["segment"] for item in segmentation], [100 * item["duty_proxy_median"] for item in segmentation], "o-")
    panels[1].set(xlabel="Equal-duration recording segment", ylabel="Median trigger-duty proxy (%)")
    save_figure(fig, args.output, "envelope_stability_across_recording")
    output = {
        "status": "PASS_UNBIASED_TRIGGER_MEASUREMENT",
        "board": "TDC_001049", "selection": "Every decoded selected-board header; no edge, ToT, ch0-reference or local-clock cut",
        "absolute_time_definition": "Original DAQ GlobalTriggerTime on native 5 ns lattice, not wall-clock UTC",
        "trigger_timestamps": expected, "first_absolute_seconds": first_seconds,
        "last_absolute_seconds": int(ticks[-1]) * 5e-9, "elapsed_seconds": float(time[-1]),
        "independent_decimal_samples_verified": samples, "absolute_histogram_integral": int(count.sum()),
        "measurements": result, "threshold_sensitivity": measurements,
        "density_classification_sensitivity": density_sensitivity,
        "classification_sensitivity_maximum_relative_frequency_change": frequency_variation,
        "classification_sensitivity_maximum_absolute_duty_proxy_change": duty_variation,
        "equal_duration_segments": segmentation,
        "timestamp_binary_sha256": digest(ticks_path),
        "source_sha256": digest(Path(__file__)),
        "hardware_causality_certified": False,
        "limitations": [
            "A periodic recorded-trigger envelope can be established; RF suppression, discriminator threshold, trigger gating and beam-current duty are not distinguished by timestamps alone.",
            "Trigger occupancy uses one-spacing endpoint correction; this is not an exact measured hardware duty cycle.",
            "Long pauses and nonregular burst intervals are explicitly reported, not silently assigned missing modulation cycles.",
            "Sparse groups remain in absolute distributions; the dense-group classification is used only to measure the dominant envelope.",
            "No trigger windows or detector hits are discarded by this diagnostic, and existing folded products are unchanged.",
        ],
    }
    require(before == (ticks_path.stat().st_size, ticks_path.stat().st_mtime_ns), "Timestamp source changed")
    for path in args.output.glob("*.png"):
        with Image.open(path) as image:
            image.verify()
        require(path.with_suffix(".pdf").read_bytes().startswith(b"%PDF-"), "Invalid PDF plot")
    (args.output / "measurement.json").write_text(json.dumps(output, indent=2) + "\n")
    with (args.output / "segment_measurements.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(segmentation[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(segmentation)
    print(json.dumps({"status": output["status"], "trigger_windows": expected,
                      "bursts": result["bursts"], "repetition_fit": result["repetition_fit"],
                      "duty_proxy": result["regular_burst_span_plus_one_trigger_interval_duty_proxy"],
                      "sparse_trigger_windows": result["sparse_group_trigger_windows"],
                      "all_trigger_accounting_verified": True}, indent=2), flush=True)


if __name__ == "__main__":
    main()
