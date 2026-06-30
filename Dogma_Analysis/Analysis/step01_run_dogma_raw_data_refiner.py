from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
import tempfile


def strip_known_suffixes(text: str) -> str:
    result = text
    changed = True
    while changed:
        changed = False
        for suffix in (".dld.dat", ".dld", ".dat"):
            if result.endswith(suffix):
                result = result[: -len(suffix)]
                changed = True
                break
    return result


def split_tokens(text: str) -> list[str]:
    return [token for token in text.split("_") if token]


def is_all_digits(text: str) -> bool:
    return bool(text) and text.isdigit()


def derive_run_key(input_path: Path) -> str:
    base_name = strip_known_suffixes(input_path.name)
    tokens = split_tokens(base_name)
    if not tokens:
        return "dogma_run"

    detector = tokens[0]

    window_token = ""
    for token in tokens:
        if token == "Pos":
            break
        if token.endswith(("ns", "us", "ms")):
            window_token = token

    position_token = ""
    for index in range(len(tokens) - 1):
        if tokens[index] == "Pos":
            position_token = f"Pos_{tokens[index + 1]}"
            break

    run_id = ""
    for token in reversed(tokens):
        if len(token) == 4 and is_all_digits(token):
            run_id = token
            break

    if detector and window_token and position_token and run_id:
        return f"{detector}_{window_token}_{position_token}_{run_id}"

    fallback_tokens = list(tokens)
    if len(fallback_tokens) >= 2:
        if len(fallback_tokens[-2]) == 8 and is_all_digits(fallback_tokens[-2]) and len(fallback_tokens[-1]) == 4 and is_all_digits(fallback_tokens[-1]):
            fallback_tokens.pop(-2)
    return "_".join(fallback_tokens)


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compile and run the step01 DOGMA raw-data refiner, then write a ROOT TTree from the cleaned pulse table."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--run-key")
    parser.add_argument("--analysis-suffix", default="rawRefined")
    parser.add_argument("--signal-tdc-ordinal", type=int)
    parser.add_argument("--trigger-tdc-ordinal", type=int)
    parser.add_argument("--trigger-reset-threshold-seconds", type=float, default=1.0)
    parser.add_argument("--max-accepted-hit-time-ns", type=float, default=50000.0)
    parser.add_argument("--tot-min-ns", type=float, default=1.0)
    parser.add_argument("--tot-max-ns", type=float, default=5000.0)
    parser.add_argument("--max-anomaly-samples", type=int, default=10)
    parser.add_argument("--progress-report-interval-seconds", type=float, default=15.0)
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--binary-path")
    parser.add_argument("--python-executable", default=sys.executable)
    parser.add_argument("--skip-root-tree", action="store_true")
    parser.add_argument("--root-tree-output")
    parser.add_argument("--tree-name", default="CleanedPulses")
    parser.add_argument("--root-tree-chunk-size", type=int, default=1_000_000)
    parser.add_argument("--root-tree-max-entries-per-file", type=int, default=20_000_000)
    parser.add_argument("--root-tree-compression", choices=("zlib", "none"), default="zlib")
    parser.add_argument("--root-tree-compression-level", type=int, default=1)
    parser.add_argument("--root-tree-progress-report-interval-seconds", type=float, default=15.0)
    parser.add_argument("--root-tree-ch0-valid-rise-min-ns", type=float, default=-410.0)
    parser.add_argument("--root-tree-ch0-valid-rise-max-ns", type=float, default=-395.0)
    parser.add_argument("--root-tree-ch0-valid-tot-min-ns", type=float, default=16.5)
    parser.add_argument("--root-tree-ch0-valid-tot-max-ns", type=float, default=19.5)
    parser.add_argument("--root-tree-rf-summary-file")
    parser.add_argument("--root-tree-rf-period-ns", type=float)
    parser.add_argument("--root-tree-rf-phase-origin-ns", type=float)
    return parser.parse_args()


def build_cleaner_command(args: argparse.Namespace, binary_path: Path, run_key: str) -> list[str]:
    command = [
        str(binary_path),
        "--input",
        args.input,
        "--output-root",
        args.output_root,
        "--run-key",
        run_key,
        "--analysis-suffix",
        args.analysis_suffix,
        "--trigger-reset-threshold-seconds",
        str(args.trigger_reset_threshold_seconds),
        "--max-accepted-hit-time-ns",
        str(args.max_accepted_hit_time_ns),
        "--tot-min-ns",
        str(args.tot_min_ns),
        "--tot-max-ns",
        str(args.tot_max_ns),
        "--max-anomaly-samples",
        str(args.max_anomaly_samples),
        "--progress-report-interval-seconds",
        str(args.progress_report_interval_seconds),
    ]
    if args.signal_tdc_ordinal is not None:
        command.extend(["--signal-tdc-ordinal", str(args.signal_tdc_ordinal)])
    if args.trigger_tdc_ordinal is not None:
        command.extend(["--trigger-tdc-ordinal", str(args.trigger_tdc_ordinal)])
    return command


def build_tree_command(
    args: argparse.Namespace,
    converter_path: Path,
    pulse_table_path: Path,
    summary_path: Path,
    root_tree_output: Path,
) -> list[str]:
    return [
        args.python_executable,
        str(converter_path),
        "--pulse-table",
        str(pulse_table_path),
        "--summary-file",
        str(summary_path),
        "--output-root-file",
        str(root_tree_output),
        "--tree-name",
        args.tree_name,
        "--chunk-size",
        str(args.root_tree_chunk_size),
        "--max-entries-per-file",
        str(args.root_tree_max_entries_per_file),
        "--compression",
        args.root_tree_compression,
        "--compression-level",
        str(args.root_tree_compression_level),
        "--progress-report-interval-seconds",
        str(args.root_tree_progress_report_interval_seconds),
        "--ch0-valid-rise-min-ns",
        str(args.root_tree_ch0_valid_rise_min_ns),
        "--ch0-valid-rise-max-ns",
        str(args.root_tree_ch0_valid_rise_max_ns),
        "--ch0-valid-tot-min-ns",
        str(args.root_tree_ch0_valid_tot_min_ns),
        "--ch0-valid-tot-max-ns",
        str(args.root_tree_ch0_valid_tot_max_ns),
    ]

    if args.root_tree_rf_summary_file:
        command.extend(["--rf-summary-file", args.root_tree_rf_summary_file])
    if args.root_tree_rf_period_ns is not None:
        command.extend(["--rf-period-ns", str(args.root_tree_rf_period_ns)])
    if args.root_tree_rf_phase_origin_ns is not None:
        command.extend(["--rf-phase-origin-ns", str(args.root_tree_rf_phase_origin_ns)])
    return command


def read_summary_value(summary_path: Path, key: str) -> str | None:
    if not summary_path.exists():
        return None
    with summary_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            current_key, value = line.split("=", 1)
            if current_key.strip() == key:
                return value.strip().strip('"')
    return None


def main() -> int:
    try:
        args = parse_args()
        script_dir = Path(__file__).resolve().parent
        source_path = script_dir / "step01_dogma_raw_data_refiner.cpp"
        converter_path = script_dir / "step01_write_cleaned_root_tree.py"
        input_path = Path(args.input).expanduser().resolve()
        output_root = Path(args.output_root).expanduser().resolve()
        run_key = args.run_key or derive_run_key(input_path)
        prefix = f"{run_key}_{args.analysis_suffix}"
        output_dir = output_root / prefix
        pulse_table_path = output_dir / f"{prefix}_pulses.tsv"
        summary_path = output_dir / f"{prefix}_summary.txt"
        root_tree_output = (
            Path(args.root_tree_output).expanduser().resolve()
            if args.root_tree_output
            else output_dir / f"{prefix}_pulses.root"
        )

        if args.binary_path:
            binary_path = Path(args.binary_path).expanduser().resolve()
            binary_path.parent.mkdir(parents=True, exist_ok=True)
        else:
            binary_path = Path(tempfile.gettempdir()) / f"step01_dogma_raw_data_refiner_{run_key}"

        print(f"Compiling: {source_path}")
        compile_cpp(source_path, binary_path, args.compiler)

        cleaner_command = build_cleaner_command(args, binary_path, run_key)
        print("Running cleaner...")
        run_command(cleaner_command)

        if not args.skip_root_tree:
            print("Writing ROOT tree...")
            run_command(build_tree_command(args, converter_path, pulse_table_path, summary_path, root_tree_output))

        print(f"Run key: {run_key}")
        print(f"Pulse table: {pulse_table_path}")
        print(f"Summary: {summary_path}")
        if not args.skip_root_tree:
            actual_root_tree_output = read_summary_value(summary_path, "root_tree_output")
            print(f"ROOT tree file: {actual_root_tree_output or root_tree_output}")
    except subprocess.CalledProcessError as error:
        print(f"Command failed with exit code {error.returncode}: {error.cmd}", file=sys.stderr)
        return error.returncode
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())