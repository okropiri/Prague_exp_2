from __future__ import annotations

import argparse
from pathlib import Path
import sys

from step02_run_cleaned_all_channel_rates import (
    compile_cpp,
    derive_run_key,
    parse_metadata,
    run_command,
    write_odd_ch0_stats,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compile and run the cleaned all-channel time-vs-ToT pipeline, including odd ch0 stats."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--run-key")
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--build-dir", default="/tmp")
    parser.add_argument("--trigger-ref-min-ns", type=float, default=-6000.0)
    parser.add_argument("--trigger-ref-max-ns", type=float, default=6000.0)
    parser.add_argument("--trigger-ref-bins", type=int, default=2400)
    parser.add_argument("--ch0ref-min-ns", type=float, default=-6000.0)
    parser.add_argument("--ch0ref-max-ns", type=float, default=6000.0)
    parser.add_argument("--ch0ref-bins", type=int, default=2400)
    parser.add_argument("--tot-min-ns", type=float, default=0.0)
    parser.add_argument("--tot-max-ns", type=float, default=128.0)
    parser.add_argument("--tot-bins", type=int, default=1280)
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

    analyzer_source = analysis_dir / "step03_dogma_cleaned_all_channel_time_tot.cpp"
    analyzer_executable = build_dir / "dogma_cleaned_all_channel_time_tot"
    writer_path = analysis_dir / "step03_write_cleaned_all_channel_time_tot_outputs.py"

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
        "--trigger-ref-min-ns",
        str(args.trigger_ref_min_ns),
        "--trigger-ref-max-ns",
        str(args.trigger_ref_max_ns),
        "--trigger-ref-bins",
        str(args.trigger_ref_bins),
        "--ch0ref-min-ns",
        str(args.ch0ref_min_ns),
        "--ch0ref-max-ns",
        str(args.ch0ref_max_ns),
        "--ch0ref-bins",
        str(args.ch0ref_bins),
        "--tot-min-ns",
        str(args.tot_min_ns),
        "--tot-max-ns",
        str(args.tot_max_ns),
        "--tot-bins",
        str(args.tot_bins),
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
    trigger_ref_file = run_dir / "Trigger_ref_ToT" / f"{run_key}_Trigger_ref_ToT_sparse.tsv"
    ch0ref_file = run_dir / "Ch0_ref_TOT" / f"{run_key}_Ch0_ref_TOT_sparse.tsv"
    tot_distrib_dir = run_dir / "TOT_distrib"
    tot_distrib_file = tot_distrib_dir / f"{run_key}_TOT_distrib.tsv"
    run_command([
        sys.executable,
        str(writer_path),
        "--trigger-ref-file",
        str(trigger_ref_file),
        "--ch0ref-file",
        str(ch0ref_file),
        "--tot-distrib-file",
        str(tot_distrib_file),
        "--tot-distrib-output-dir",
        str(tot_distrib_dir),
    ])

    odd_stats_path, odd_samples_path = write_odd_ch0_stats(input_path, run_dir, run_key, metadata, args)
    print(f"Odd ch0 stats: {odd_stats_path}")
    print(f"Odd ch0 samples: {odd_samples_path}")


if __name__ == "__main__":
    main()