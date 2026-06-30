from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


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


def derive_run_key(input_path: Path) -> str:
    return strip_known_suffixes(input_path.name)


def run_command(command: list[str]) -> None:
    subprocess.run(command, check=True)


def compile_cpp(source_path: Path, output_path: Path, compiler: str) -> None:
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
        description="Compile and run the absolute-rate and ch0-referenced rate/plot analyses."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--run-key")
    parser.add_argument("--mode", choices=("abs", "ch0ref", "both"), default="both")
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--build-dir", default="/tmp")

    parser.add_argument("--abs-bin-width-ns", type=float, default=1000.0)
    parser.add_argument("--abs-trigger-reset-threshold-seconds", type=float, default=1.0)

    parser.add_argument("--ch0-folded-min-ns", type=float, default=-6000.0)
    parser.add_argument("--ch0-folded-max-ns", type=float, default=6000.0)
    parser.add_argument("--ch0-folded-bins", type=int, default=12000)
    parser.add_argument("--ch0-tot-min-ns", type=float, default=0.0)
    parser.add_argument("--ch0-tot-max-ns", type=float, default=128.0)
    parser.add_argument("--ch0-tot-bins", type=int, default=512)
    parser.add_argument("--ch0-rate-bin-width-ns", type=float, default=1000.0)
    parser.add_argument("--ch0-asymmetry-block-size-windows", type=int, default=1000)
    parser.add_argument("--ch0-trigger-reset-threshold-seconds", type=float, default=1.0)
    return parser.parse_args()


def run_abs_pipeline(args: argparse.Namespace, analysis_dir: Path, build_dir: Path, input_path: Path, output_root: Path, run_key: str) -> None:
    executable_path = build_dir / "dogma_absolute_rate_profile"
    source_path = analysis_dir / "dogma_absolute_rate_profile.cpp"
    writer_path = analysis_dir / "write_absolute_rate_outputs.py"
    compile_cpp(source_path, executable_path, args.compiler)

    output_dir = output_root / f"{run_key}_absRate"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_dir / f"{run_key}_absRate"

    run_command([
        str(executable_path),
        "--input",
        str(input_path),
        "--output-prefix",
        str(output_prefix),
        "--bin-width-ns",
        str(args.abs_bin_width_ns),
        "--trigger-reset-threshold-seconds",
        str(args.abs_trigger_reset_threshold_seconds),
    ])

    run_command([
        sys.executable,
        str(writer_path),
        "--rates-file",
        str(output_prefix) + "_rates.txt",
        "--summary-file",
        str(output_prefix) + "_summary.txt",
        "--root-output",
        str(output_prefix) + ".root",
        "--pdf-output",
        str(output_prefix) + ".pdf",
        "--png-output",
        str(output_prefix) + ".png",
    ])



def run_ch0ref_pipeline(args: argparse.Namespace, analysis_dir: Path, build_dir: Path, input_path: Path, output_root: Path, run_key: str) -> None:
    executable_path = build_dir / "dogma_ch0ref_full_analysis"
    source_path = analysis_dir / "dogma_ch0ref_full_analysis.cpp"
    writer_path = analysis_dir / "write_ch0ref_full_analysis_outputs.py"
    compile_cpp(source_path, executable_path, args.compiler)

    output_dir = output_root / f"{run_key}_ch0Ref_rawdata"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = output_dir / f"{run_key}_ch0Ref_rawdata"

    run_command([
        str(executable_path),
        "--input",
        str(input_path),
        "--output-prefix",
        str(output_prefix),
        "--folded-min-ns",
        str(args.ch0_folded_min_ns),
        "--folded-max-ns",
        str(args.ch0_folded_max_ns),
        "--folded-bins",
        str(args.ch0_folded_bins),
        "--tot-min-ns",
        str(args.ch0_tot_min_ns),
        "--tot-max-ns",
        str(args.ch0_tot_max_ns),
        "--tot-bins",
        str(args.ch0_tot_bins),
        "--rate-bin-width-ns",
        str(args.ch0_rate_bin_width_ns),
        "--asymmetry-block-size-windows",
        str(args.ch0_asymmetry_block_size_windows),
        "--trigger-reset-threshold-seconds",
        str(args.ch0_trigger_reset_threshold_seconds),
    ])

    run_command([
        sys.executable,
        str(writer_path),
        "--input-prefix",
        str(output_prefix),
        "--root-output",
        str(output_prefix) + ".root",
        "--pdf-prefix",
        str(output_prefix),
        "--png-prefix",
        str(output_prefix),
    ])



def main() -> None:
    args = parse_args()
    analysis_dir = Path(__file__).resolve().parent
    input_path = Path(args.input).resolve()
    output_root = Path(args.output_root).resolve()
    build_dir = Path(args.build_dir).resolve()
    build_dir.mkdir(parents=True, exist_ok=True)
    output_root.mkdir(parents=True, exist_ok=True)

    run_key = args.run_key or derive_run_key(input_path)

    if args.mode in ("abs", "both"):
        run_abs_pipeline(args, analysis_dir, build_dir, input_path, output_root, run_key)
    if args.mode in ("ch0ref", "both"):
        run_ch0ref_pipeline(args, analysis_dir, build_dir, input_path, output_root, run_key)


if __name__ == "__main__":
    main()
