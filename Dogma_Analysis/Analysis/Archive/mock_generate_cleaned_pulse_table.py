from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from statistics import NormalDist


class RunningStats:
    def __init__(self) -> None:
        self.count = 0
        self.sum = 0.0
        self.sum_sq = 0.0
        self.minimum = math.inf
        self.maximum = -math.inf

    def add(self, value: float) -> None:
        self.count += 1
        self.sum += value
        self.sum_sq += value * value
        self.minimum = min(self.minimum, value)
        self.maximum = max(self.maximum, value)

    def as_dict(self) -> dict[str, float | int]:
        if self.count == 0:
            return {"count": 0}
        mean = self.sum / self.count
        variance = max(0.0, self.sum_sq / self.count - mean * mean)
        return {
            "count": self.count,
            "mean": mean,
            "sigma": math.sqrt(variance),
            "min": self.minimum,
            "max": self.maximum,
        }


def format_float(value: float) -> str:
    return f"{value:.6f}".rstrip("0").rstrip(".")


def build_gaussian_profile(size: int, mean: float, sigma: float, seed: int) -> list[float]:
    if size <= 0:
        raise ValueError("Profile size must be positive")
    distribution = NormalDist(mu=mean, sigma=sigma)
    values = [distribution.inv_cdf((index + 0.5) / size) for index in range(size)]
    random.Random(seed).shuffle(values)
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate deterministic mock cleaned DOGMA pulse data for testing the cleaned-data analysis chain."
    )
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--run-key", default="Mock_data")
    parser.add_argument("--trigger-count", type=int, default=10_000)
    parser.add_argument("--trigger-rate-hz", type=float, default=25_000_000.0 / (2**10))
    parser.add_argument("--window-width-ns", type=float, default=20_000.0)
    parser.add_argument("--ch0-channel", type=int, default=0)
    parser.add_argument("--ch2-channel", type=int, default=2)
    parser.add_argument("--ch0-rise-ns", type=float, default=-400.0)
    parser.add_argument("--ch0-tot-mean-ns", type=float, default=17.0)
    parser.add_argument("--ch0-tot-sigma-ns", type=float, default=0.12)
    parser.add_argument("--rf-period-ns", type=float, default=40.0)
    parser.add_argument("--rf-phase-center-ns", type=float, default=20.0)
    parser.add_argument("--ch2-time-jitter-sigma-ns", type=float, default=2.0)
    parser.add_argument("--ch2-tot-mean-ns", type=float, default=13.5)
    parser.add_argument("--ch2-tot-sigma-ns", type=float, default=1.2)
    parser.add_argument("--trigger-populated-min-ns", type=float, default=-5880.0)
    parser.add_argument("--trigger-populated-max-ns", type=float, default=5880.0)
    parser.add_argument("--ch2-pulses-per-bucket", type=int, default=1)
    parser.add_argument("--profile-size", type=int, default=8192)
    parser.add_argument("--seed", type=int, default=20260515)
    parser.add_argument("--period-tolerance-ns", type=float, default=0.01)
    parser.add_argument("--phase-origin-tolerance-ns", type=float, default=0.25)
    parser.add_argument("--sigma-tolerance-ns", type=float, default=0.1)
    parser.add_argument("--selected-fraction-min", type=float, default=0.999)
    return parser.parse_args()


def write_mock_spec_markdown(path: Path, spec: dict[str, object]) -> None:
    parameters = spec["parameters"]
    expected = spec["expected"]
    actual = spec["actual_samples"]
    lines = [
        "# Mock Data Specification",
        "",
        "## Overview",
        f"- Run key: {spec['run_key']}",
        f"- Pulse table: {spec['pulse_table_path']}",
        f"- Generation algorithm: {spec['generation_algorithm']}",
        "",
        "## Parameters",
        f"- Trigger count: {parameters['trigger_count']}",
        f"- Trigger rate: {parameters['trigger_rate_hz']} Hz",
        f"- Trigger period: {parameters['trigger_period_ns']} ns",
        f"- Trigger window width: {parameters['window_width_ns']} ns",
        f"- Ch0 rise time: {parameters['ch0_rise_ns']} ns",
        f"- Ch0 ToT target: mean {parameters['ch0_tot_mean_ns']} ns, sigma {parameters['ch0_tot_sigma_ns']} ns",
        f"- RF period target: {parameters['rf_period_ns']} ns",
        f"- RF phase center target: {parameters['rf_phase_center_ns']} ns",
        f"- Ch2 time jitter target: sigma {parameters['ch2_time_jitter_sigma_ns']} ns",
        f"- Ch2 ToT target: mean {parameters['ch2_tot_mean_ns']} ns, sigma {parameters['ch2_tot_sigma_ns']} ns",
        f"- Populated trigger-time range: [{parameters['trigger_populated_min_ns']}, {parameters['trigger_populated_max_ns']}] ns",
        f"- RF bucket index range: [{parameters['rf_bucket_index_min']}, {parameters['rf_bucket_index_max']}], bucket count {parameters['rf_bucket_count']}",
        f"- Ch2 pulses per bucket per trigger: {parameters['ch2_pulses_per_bucket']}",
        f"- Fixed profile size: {parameters['profile_size']}",
        f"- Seed: {parameters['seed']}",
        "",
        "## Expected Counts",
        f"- Ch0 pulses: {expected['total_ch0_pulses']}",
        f"- Ch2 pulses: {expected['total_ch2_pulses']}",
        f"- Total pulse rows: {expected['total_rows']}",
        f"- Expected reconstructed RF period: {expected['expected_rf_period_ns']} ns",
        f"- Expected reconstructed phase origin: {expected['expected_phase_origin_ns']} ns",
        f"- Expected reconstructed sigma: {expected['expected_sigma_ns']} ns",
        f"- Expected minimum selected fraction: {expected['selected_fraction_min']}",
        "",
        "## Realized Sample Statistics",
        f"- Ch0 ToT: mean {format_float(actual['ch0_tot']['mean'])} ns, sigma {format_float(actual['ch0_tot']['sigma'])} ns, range [{format_float(actual['ch0_tot']['min'])}, {format_float(actual['ch0_tot']['max'])}] ns",
        f"- Ch2 ToT: mean {format_float(actual['ch2_tot']['mean'])} ns, sigma {format_float(actual['ch2_tot']['sigma'])} ns, range [{format_float(actual['ch2_tot']['min'])}, {format_float(actual['ch2_tot']['max'])}] ns",
        f"- Ch2 phase jitter: mean {format_float(actual['ch2_phase_jitter']['mean'])} ns, sigma {format_float(actual['ch2_phase_jitter']['sigma'])} ns, range [{format_float(actual['ch2_phase_jitter']['min'])}, {format_float(actual['ch2_phase_jitter']['max'])}] ns",
        f"- Ch2 absolute trigger time: mean {format_float(actual['ch2_trigger_time']['mean'])} ns, sigma {format_float(actual['ch2_trigger_time']['sigma'])} ns, range [{format_float(actual['ch2_trigger_time']['min'])}, {format_float(actual['ch2_trigger_time']['max'])}] ns",
        "",
        "## QC Features",
        "- Deterministic shuffled-quantile Gaussian profiles are used instead of ad-hoc random draws, so the same mock data can be regenerated exactly.",
        "- Ch0 has exactly one valid reference pulse in every trigger window, with no missing or ambiguous candidates expected.",
        "- Ch2 timing peaks are placed at RF phase centers separated by exact 40 ns bucket spacing, so period reconstruction can be checked against a known truth.",
        "- The mock specification stores explicit tolerances for RF period, phase origin, sigma, and selected-fraction checks.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.trigger_count <= 0:
        raise ValueError("trigger-count must be positive")
    if args.trigger_rate_hz <= 0.0:
        raise ValueError("trigger-rate-hz must be positive")
    if args.window_width_ns <= 0.0:
        raise ValueError("window-width-ns must be positive")
    if args.rf_period_ns <= 0.0:
        raise ValueError("rf-period-ns must be positive")
    if args.ch2_pulses_per_bucket <= 0:
        raise ValueError("ch2-pulses-per-bucket must be positive")
    if args.profile_size <= 0:
        raise ValueError("profile-size must be positive")
    if not args.trigger_populated_min_ns < args.trigger_populated_max_ns:
        raise ValueError("trigger-populated-min-ns must be smaller than trigger-populated-max-ns")

    output_root = Path(args.output_root).resolve()
    output_dir = output_root / args.run_key
    output_dir.mkdir(parents=True, exist_ok=True)

    pulse_table_path = output_dir / f"{args.run_key}_pulses.tsv"
    variables_json_path = output_dir / f"{args.run_key}_variables.json"
    spec_markdown_path = output_dir / "00_mock_data_spec.md"

    trigger_period_seconds = 1.0 / args.trigger_rate_hz
    trigger_period_ns = trigger_period_seconds * 1.0e9

    base_absolute_phase_center_ns = args.ch0_rise_ns + args.rf_phase_center_ns
    bucket_index_min = math.ceil((args.trigger_populated_min_ns - base_absolute_phase_center_ns) / args.rf_period_ns)
    bucket_index_max = math.floor((args.trigger_populated_max_ns - base_absolute_phase_center_ns) / args.rf_period_ns)
    if bucket_index_min > bucket_index_max:
        raise ValueError("No RF buckets fit inside the populated trigger-time range")
    bucket_indices = list(range(bucket_index_min, bucket_index_max + 1))
    bucket_count = len(bucket_indices)

    ch0_tot_profile = build_gaussian_profile(args.trigger_count, args.ch0_tot_mean_ns, args.ch0_tot_sigma_ns, args.seed + 11)
    ch2_time_profile = build_gaussian_profile(args.profile_size, 0.0, args.ch2_time_jitter_sigma_ns, args.seed + 23)
    ch2_tot_profile = build_gaussian_profile(args.profile_size, args.ch2_tot_mean_ns, args.ch2_tot_sigma_ns, args.seed + 37)

    ch0_tot_stats = RunningStats()
    ch2_tot_stats = RunningStats()
    ch2_phase_jitter_stats = RunningStats()
    ch2_trigger_time_stats = RunningStats()

    total_ch2_pulses = args.trigger_count * bucket_count * args.ch2_pulses_per_bucket
    total_rows = args.trigger_count + total_ch2_pulses

    with pulse_table_path.open("w", encoding="utf-8", buffering=1024 * 1024) as handle:
        handle.write(f"# run_key={args.run_key}\n")
        handle.write(f"# input_file={args.run_key}.synthetic.dld.dat\n")
        handle.write("# source_kind=mock_cleaned_pulse_table\n")
        handle.write("# mock_profile_algorithm=deterministic_shuffled_quantiles\n")
        handle.write(f"# trigger_count={args.trigger_count}\n")
        handle.write(f"# trigger_rate_hz={format_float(args.trigger_rate_hz)}\n")
        handle.write(f"# trigger_period_ns={format_float(trigger_period_ns)}\n")
        handle.write(f"# trigger_window_width_ns={format_float(args.window_width_ns)}\n")
        handle.write(f"# ch0_rise_ns={format_float(args.ch0_rise_ns)}\n")
        handle.write(f"# ch0_tot_mean_ns={format_float(args.ch0_tot_mean_ns)}\n")
        handle.write(f"# ch0_tot_sigma_ns={format_float(args.ch0_tot_sigma_ns)}\n")
        handle.write(f"# rf_period_ns={format_float(args.rf_period_ns)}\n")
        handle.write(f"# rf_phase_center_ns={format_float(args.rf_phase_center_ns)}\n")
        handle.write(f"# ch2_time_jitter_sigma_ns={format_float(args.ch2_time_jitter_sigma_ns)}\n")
        handle.write(f"# ch2_tot_mean_ns={format_float(args.ch2_tot_mean_ns)}\n")
        handle.write(f"# ch2_tot_sigma_ns={format_float(args.ch2_tot_sigma_ns)}\n")
        handle.write(f"# trigger_populated_min_ns={format_float(args.trigger_populated_min_ns)}\n")
        handle.write(f"# trigger_populated_max_ns={format_float(args.trigger_populated_max_ns)}\n")
        handle.write(f"# rf_bucket_index_min={bucket_index_min}\n")
        handle.write(f"# rf_bucket_index_max={bucket_index_max}\n")
        handle.write(f"# rf_bucket_count={bucket_count}\n")
        handle.write(f"# ch2_pulses_per_bucket={args.ch2_pulses_per_bucket}\n")
        handle.write(f"# seed={args.seed}\n")
        handle.write("window_index tdc_ordinal global_trigger_seconds channel rise_tdc_index_one_based fall_tdc_index_one_based rise_ns fall_ns tot_ns\n")

        pulse_counter = 0
        for window_index in range(args.trigger_count):
            global_trigger_seconds = window_index * trigger_period_seconds
            ch0_tot_ns = ch0_tot_profile[window_index]
            ch0_tot_stats.add(ch0_tot_ns)

            entries: list[tuple[int, float, float]] = [(args.ch0_channel, args.ch0_rise_ns, ch0_tot_ns)]
            for bucket_offset, bucket_index in enumerate(bucket_indices):
                center_rise_ns = base_absolute_phase_center_ns + bucket_index * args.rf_period_ns
                for repeat_index in range(args.ch2_pulses_per_bucket):
                    profile_index = pulse_counter % args.profile_size
                    time_jitter_ns = ch2_time_profile[(profile_index * 2053 + 17) % args.profile_size]
                    ch2_tot_ns = ch2_tot_profile[(profile_index * 4093 + 97) % args.profile_size]
                    rise_ns = center_rise_ns + time_jitter_ns
                    entries.append((args.ch2_channel, rise_ns, ch2_tot_ns))
                    ch2_phase_jitter_stats.add(time_jitter_ns)
                    ch2_tot_stats.add(ch2_tot_ns)
                    ch2_trigger_time_stats.add(rise_ns)
                    pulse_counter += 1

            entries.sort(key=lambda item: item[1])
            for channel, rise_ns, tot_ns in entries:
                fall_ns = rise_ns + tot_ns
                handle.write(
                    f"{window_index}\t1\t{global_trigger_seconds:.9f}\t{channel}\t1\t2\t"
                    f"{rise_ns:.6f}\t{fall_ns:.6f}\t{tot_ns:.6f}\n"
                )

    spec = {
        "run_key": args.run_key,
        "pulse_table_path": str(pulse_table_path),
        "generation_algorithm": "deterministic shuffled Gaussian quantiles with fixed seeds",
        "parameters": {
            "trigger_count": args.trigger_count,
            "trigger_rate_hz": args.trigger_rate_hz,
            "trigger_period_ns": trigger_period_ns,
            "window_width_ns": args.window_width_ns,
            "ch0_channel": args.ch0_channel,
            "ch2_channel": args.ch2_channel,
            "ch0_rise_ns": args.ch0_rise_ns,
            "ch0_tot_mean_ns": args.ch0_tot_mean_ns,
            "ch0_tot_sigma_ns": args.ch0_tot_sigma_ns,
            "rf_period_ns": args.rf_period_ns,
            "rf_phase_center_ns": args.rf_phase_center_ns,
            "ch2_time_jitter_sigma_ns": args.ch2_time_jitter_sigma_ns,
            "ch2_tot_mean_ns": args.ch2_tot_mean_ns,
            "ch2_tot_sigma_ns": args.ch2_tot_sigma_ns,
            "trigger_populated_min_ns": args.trigger_populated_min_ns,
            "trigger_populated_max_ns": args.trigger_populated_max_ns,
            "rf_bucket_index_min": bucket_index_min,
            "rf_bucket_index_max": bucket_index_max,
            "rf_bucket_count": bucket_count,
            "ch2_pulses_per_bucket": args.ch2_pulses_per_bucket,
            "profile_size": args.profile_size,
            "seed": args.seed,
        },
        "expected": {
            "total_ch0_pulses": args.trigger_count,
            "total_ch2_pulses": total_ch2_pulses,
            "total_rows": total_rows,
            "expected_rf_period_ns": args.rf_period_ns,
            "expected_phase_origin_ns": args.rf_phase_center_ns,
            "expected_sigma_ns": args.ch2_time_jitter_sigma_ns,
            "expected_valid_ch0_windows": args.trigger_count,
            "expected_missing_valid_ch0_windows": 0,
            "selected_fraction_min": args.selected_fraction_min,
            "tolerance_period_ns": args.period_tolerance_ns,
            "tolerance_phase_origin_ns": args.phase_origin_tolerance_ns,
            "tolerance_sigma_ns": args.sigma_tolerance_ns,
        },
        "actual_samples": {
            "ch0_tot": ch0_tot_stats.as_dict(),
            "ch2_tot": ch2_tot_stats.as_dict(),
            "ch2_phase_jitter": ch2_phase_jitter_stats.as_dict(),
            "ch2_trigger_time": ch2_trigger_time_stats.as_dict(),
        },
    }

    variables_json_path.write_text(json.dumps(spec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_mock_spec_markdown(spec_markdown_path, spec)

    print(f"Pulse table: {pulse_table_path}")
    print(f"Variables JSON: {variables_json_path}")
    print(f"Specification Markdown: {spec_markdown_path}")
    print(f"Trigger count: {args.trigger_count}")
    print(f"Ch2 bucket count: {bucket_count}")
    print(f"Ch2 total pulses: {total_ch2_pulses}")
    print(f"Expected RF period (ns): {format_float(args.rf_period_ns)}")
    print(f"Expected RF phase center (ns): {format_float(args.rf_phase_center_ns)}")


if __name__ == "__main__":
    main()