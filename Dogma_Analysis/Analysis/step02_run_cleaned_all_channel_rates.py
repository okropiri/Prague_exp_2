from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import subprocess
import sys


SAMPLE_LIMIT_PER_CATEGORY = 12


def strip_known_suffixes(text: str) -> str:
    result = text
    changed = True
    while changed:
        changed = False
        for suffix in ("_rawRefined_pulses.tsv", ".tsv", ".txt"):
            if result.endswith(suffix):
                result = result[: -len(suffix)]
                changed = True
                break
    return result


def parse_metadata(input_path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    with input_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if "=" not in line:
                continue
            key, value = line[1:].split("=", 1)
            metadata[key.strip()] = value.strip()
    return metadata


def derive_run_key(input_path: Path, metadata: dict[str, str], explicit_run_key: str | None) -> str:
    if explicit_run_key:
        return explicit_run_key
    if metadata.get("run_key"):
        return metadata["run_key"]
    return strip_known_suffixes(input_path.name)


def run_command(command: list[str]) -> None:
    subprocess.run(command, check=True)


def compile_cpp(source_path: Path, output_path: Path, compiler: str) -> None:
    if output_path.exists() and output_path.stat().st_mtime_ns >= source_path.stat().st_mtime_ns:
        print(f"Reusing compiled binary: {output_path}")
        return
    run_command([
        compiler,
        "-std=c++17",
        "-O2",
        "-Wall",
        "-Wextra",
        "-pedantic",
        str(source_path),
        "-o",
        str(output_path),
    ])


def quantile(values: list[float], fraction: float) -> float:
    if not values:
        return float("nan")
    index = int(round((len(values) - 1) * fraction))
    index = min(len(values) - 1, max(0, index))
    return values[index]


def summarize(values: list[float]) -> dict[str, float | int]:
    sorted_values = sorted(values)
    if not sorted_values:
        return {"count": 0}
    return {
        "count": len(sorted_values),
        "min": sorted_values[0],
        "q01": quantile(sorted_values, 0.01),
        "q05": quantile(sorted_values, 0.05),
        "median": quantile(sorted_values, 0.50),
        "q95": quantile(sorted_values, 0.95),
        "q99": quantile(sorted_values, 0.99),
        "max": sorted_values[-1],
    }


def format_value(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    if value != value:
        return "nan"
    return f"{value:.6f}".rstrip("0").rstrip(".")


@dataclass
class Ch0Hit:
    rise_ns: float
    tot_ns: float


@dataclass
class SampleWindow:
    category: str
    window_index: int
    selected_candidate_ordinal: int
    valid_candidate_count: int
    hits: list[Ch0Hit]
    reasons: list[str]


def classify_ch0_hit(hit: Ch0Hit, args: argparse.Namespace) -> str:
    outside_time = hit.rise_ns < args.ch0_valid_rise_min_ns or hit.rise_ns > args.ch0_valid_rise_max_ns
    outside_tot = hit.tot_ns < args.ch0_valid_tot_min_ns or hit.tot_ns > args.ch0_valid_tot_max_ns
    if not outside_time and not outside_tot:
        return "valid"
    if outside_time and outside_tot:
        return "outside_time_and_tot"
    if outside_time:
        return "outside_time"
    return "outside_tot"


def write_odd_ch0_stats(input_path: Path, output_dir: Path, run_key: str, metadata: dict[str, str], args: argparse.Namespace) -> tuple[Path, Path]:
    stats: Counter[str] = Counter()
    selected_candidate_ordinal_counts: Counter[int] = Counter()
    rejected_rise_ns: list[float] = []
    rejected_tot_ns: list[float] = []
    leading_rejected_rise_ns: list[float] = []
    leading_rejected_tot_ns: list[float] = []
    sample_windows: list[SampleWindow] = []
    sample_counts: Counter[str] = Counter()

    def add_sample(category: str, window_index: int, selected_candidate_ordinal: int, valid_candidate_count: int, hits: list[Ch0Hit], reasons: list[str]) -> None:
        if sample_counts[category] >= SAMPLE_LIMIT_PER_CATEGORY:
            return
        sample_windows.append(SampleWindow(category, window_index, selected_candidate_ordinal, valid_candidate_count, list(hits), list(reasons)))
        sample_counts[category] += 1

    def flush_window(window_index: int | None, hits: list[Ch0Hit]) -> None:
        if window_index is None:
            return
        if not hits:
            stats["windows_without_ch0"] += 1
            return

        stats["windows_with_ch0"] += 1
        stats["total_ch0_candidates"] += len(hits)
        if len(hits) > 1:
            stats["windows_with_multiple_ch0"] += 1

        reasons = [classify_ch0_hit(hit, args) for hit in hits]
        valid_ordinals = [index for index, reason in enumerate(reasons, start=1) if reason == "valid"]
        for hit, reason in zip(hits, reasons):
            if reason == "valid":
                stats["total_valid_ch0_candidates"] += 1
                continue
            stats["total_rejected_ch0_candidates"] += 1
            rejected_rise_ns.append(hit.rise_ns)
            rejected_tot_ns.append(hit.tot_ns)
            stats[f"rejected_{reason}"] += 1

        if not valid_ordinals:
            stats["windows_without_valid_ch0_candidate"] += 1
            add_sample("no_valid_ch0_candidate", window_index, 0, 0, hits, reasons)
            return

        selected_candidate_ordinal = valid_ordinals[0]
        selected_candidate_ordinal_counts[selected_candidate_ordinal] += 1
        if selected_candidate_ordinal == 1:
            stats["windows_using_first_valid_ch0_candidate"] += 1
        else:
            stats["windows_rescued_by_later_valid_ch0_candidate"] += 1
            stats["rejected_leading_ch0_candidates"] += selected_candidate_ordinal - 1
            for hit in hits[: selected_candidate_ordinal - 1]:
                leading_rejected_rise_ns.append(hit.rise_ns)
                leading_rejected_tot_ns.append(hit.tot_ns)
            add_sample("rescued_by_later_valid_ch0_candidate", window_index, selected_candidate_ordinal, len(valid_ordinals), hits, reasons)

        if len(valid_ordinals) == 1:
            stats["windows_with_exactly_one_valid_ch0_candidate"] += 1
        else:
            stats["windows_with_multiple_valid_ch0_candidates"] += 1
            add_sample("multiple_valid_ch0_candidates", window_index, selected_candidate_ordinal, len(valid_ordinals), hits, reasons)

    current_window: int | None = None
    current_hits: list[Ch0Hit] = []
    with input_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#") or line.startswith("window_index"):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue
            window_index = int(parts[0])
            if current_window is None:
                current_window = window_index
            if window_index != current_window:
                flush_window(current_window, current_hits)
                current_window = window_index
                current_hits = []
            if parts[3] == "0":
                current_hits.append(Ch0Hit(rise_ns=float(parts[6]), tot_ns=float(parts[8])))
    flush_window(current_window, current_hits)

    stats_path = output_dir / f"{run_key}_odd_ch0_reference_stats.txt"
    with stats_path.open("w", encoding="utf-8") as handle:
        handle.write(f"input_pulse_table={input_path}\n")
        if metadata.get("input_file"):
            handle.write(f"input_raw_file={metadata['input_file']}\n")
        handle.write(f"run_key={run_key}\n")
        handle.write(f"ch0_valid_rise_min_ns={format_value(args.ch0_valid_rise_min_ns)}\n")
        handle.write(f"ch0_valid_rise_max_ns={format_value(args.ch0_valid_rise_max_ns)}\n")
        handle.write(f"ch0_valid_tot_min_ns={format_value(args.ch0_valid_tot_min_ns)}\n")
        handle.write(f"ch0_valid_tot_max_ns={format_value(args.ch0_valid_tot_max_ns)}\n")
        for key in (
            "windows_without_ch0",
            "windows_with_ch0",
            "windows_with_multiple_ch0",
            "windows_without_valid_ch0_candidate",
            "windows_using_first_valid_ch0_candidate",
            "windows_rescued_by_later_valid_ch0_candidate",
            "windows_with_exactly_one_valid_ch0_candidate",
            "windows_with_multiple_valid_ch0_candidates",
            "total_ch0_candidates",
            "total_valid_ch0_candidates",
            "total_rejected_ch0_candidates",
            "rejected_outside_time",
            "rejected_outside_tot",
            "rejected_outside_time_and_tot",
            "rejected_leading_ch0_candidates",
        ):
            handle.write(f"{key}={stats[key]}\n")
        for ordinal in sorted(selected_candidate_ordinal_counts):
            handle.write(f"selected_candidate_ordinal_{ordinal}={selected_candidate_ordinal_counts[ordinal]}\n")

        for label, values in (
            ("rejected_ch0_rise_ns", rejected_rise_ns),
            ("rejected_ch0_tot_ns", rejected_tot_ns),
            ("leading_rejected_ch0_rise_ns", leading_rejected_rise_ns),
            ("leading_rejected_ch0_tot_ns", leading_rejected_tot_ns),
        ):
            summary = summarize(values)
            for summary_key, summary_value in summary.items():
                handle.write(f"{label}_{summary_key}={format_value(summary_value)}\n")

    samples_path = output_dir / f"{run_key}_odd_ch0_reference_samples.tsv"
    with samples_path.open("w", encoding="utf-8") as handle:
        handle.write(
            "category\twindow_index\tselected_candidate_ordinal\tvalid_candidate_count\ttotal_ch0_hits\t"
            "hit_ordinal\trise_ns\ttot_ns\tclassification\n"
        )
        for sample in sample_windows:
            for hit_ordinal, (hit, reason) in enumerate(zip(sample.hits, sample.reasons), start=1):
                handle.write(
                    f"{sample.category}\t{sample.window_index}\t{sample.selected_candidate_ordinal}\t"
                    f"{sample.valid_candidate_count}\t{len(sample.hits)}\t{hit_ordinal}\t"
                    f"{format_value(hit.rise_ns)}\t{format_value(hit.tot_ns)}\t{reason}\n"
                )

    return stats_path, samples_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compile and run the cleaned all-channel absolute and ch0-referenced rate pipeline, including odd ch0 stats."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--run-key")
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--build-dir", default="/tmp")
    parser.add_argument("--abs-display-bin-width-seconds", type=float, default=0.1)
    parser.add_argument("--ch0ref-min-ns", type=float, default=-6000.0)
    parser.add_argument("--ch0ref-max-ns", type=float, default=6000.0)
    parser.add_argument("--ch0ref-bins", type=int, default=12000)
    parser.add_argument("--ch0-valid-rise-min-ns", type=float, default=-410.0)
    parser.add_argument("--ch0-valid-rise-max-ns", type=float, default=-395.0)
    parser.add_argument("--ch0-valid-tot-min-ns", type=float, default=16.5)
    parser.add_argument("--ch0-valid-tot-max-ns", type=float, default=19.5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    analysis_dir = Path(__file__).resolve().parent
    input_path = Path(args.input).resolve()
    output_root = Path(args.output_root).resolve()
    build_dir = Path(args.build_dir).resolve()
    build_dir.mkdir(parents=True, exist_ok=True)
    output_root.mkdir(parents=True, exist_ok=True)

    metadata = parse_metadata(input_path)
    run_key = derive_run_key(input_path, metadata, args.run_key)

    analyzer_source = analysis_dir / "step02_dogma_cleaned_all_channel_rates.cpp"
    analyzer_executable = build_dir / "dogma_cleaned_all_channel_rates"
    writer_path = analysis_dir / "step02_write_cleaned_all_channel_rates_outputs.py"

    compile_cpp(analyzer_source, analyzer_executable, args.compiler)
    run_command([sys.executable, "-m", "py_compile", str(writer_path)])

    run_command([
        str(analyzer_executable),
        "--input",
        str(input_path),
        "--output-root",
        str(output_root),
        "--run-key",
        run_key,
        "--abs-display-bin-width-seconds",
        str(args.abs_display_bin_width_seconds),
        "--ch0ref-min-ns",
        str(args.ch0ref_min_ns),
        "--ch0ref-max-ns",
        str(args.ch0ref_max_ns),
        "--ch0ref-bins",
        str(args.ch0ref_bins),
        "--ch0-valid-rise-min-ns",
        str(args.ch0_valid_rise_min_ns),
        "--ch0-valid-rise-max-ns",
        str(args.ch0_valid_rise_max_ns),
        "--ch0-valid-tot-min-ns",
        str(args.ch0_valid_tot_min_ns),
        "--ch0-valid-tot-max-ns",
        str(args.ch0_valid_tot_max_ns),
    ])

    run_dir = output_root / run_key
    abs_matrix = run_dir / "Abs_rates" / f"{run_key}_Abs_rates_matrix.tsv"
    ch0ref_matrix = run_dir / "Ch0_ref_Rates" / f"{run_key}_Ch0_ref_Rates_matrix.tsv"
    run_command([
        sys.executable,
        str(writer_path),
        "--abs-matrix-file",
        str(abs_matrix),
        "--ch0ref-matrix-file",
        str(ch0ref_matrix),
    ])

    odd_stats_path, odd_samples_path = write_odd_ch0_stats(input_path, run_dir, run_key, metadata, args)
    print(f"Odd ch0 stats: {odd_stats_path}")
    print(f"Odd ch0 samples: {odd_samples_path}")


if __name__ == "__main__":
    main()