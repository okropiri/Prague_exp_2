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


def read_key_value_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip().strip('"')
    return values


def can_reuse_odd_ch0_stats(stats_path: Path, samples_path: Path, input_path: Path, args: argparse.Namespace) -> bool:
    if not stats_path.exists() or not samples_path.exists():
        return False
    values = read_key_value_file(stats_path)
    expected_values = {
        "input_pulse_table": str(input_path),
        "ch0_valid_rise_min_ns": str(args.ch0_valid_rise_min_ns),
        "ch0_valid_rise_max_ns": str(args.ch0_valid_rise_max_ns),
        "ch0_valid_tot_min_ns": str(args.ch0_valid_tot_min_ns),
        "ch0_valid_tot_max_ns": str(args.ch0_valid_tot_max_ns),
    }
    return all(values.get(key) == expected for key, expected in expected_values.items())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compile and run the cleaned all-channel RF-period pipeline, including odd ch0 stats."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--run-key")
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--build-dir", default="/tmp")
    parser.add_argument("--ch0ref-min-ns", type=float, default=-6000.0)
    parser.add_argument("--ch0ref-max-ns", type=float, default=6000.0)
    parser.add_argument("--ch0ref-bins", type=int, default=12000)
    parser.add_argument("--trigger-ref-min-ns", type=float, default=-6000.0)
    parser.add_argument("--trigger-ref-max-ns", type=float, default=6000.0)
    parser.add_argument("--trigger-ref-bins", type=int, default=2400)
    parser.add_argument("--tot-min-ns", type=float, default=0.0)
    parser.add_argument("--tot-max-ns", type=float, default=128.0)
    parser.add_argument("--tot-bins", type=int, default=1280)
    parser.add_argument("--score-channel", type=int, default=2)
    parser.add_argument("--score-tot-min-ns", type=float, default=0.0)
    parser.add_argument("--score-tot-max-ns", type=float, default=128.0)
    parser.add_argument("--initial-period-min-ns", type=float, default=37.0)
    parser.add_argument("--initial-period-max-ns", type=float, default=41.0)
    parser.add_argument("--initial-step-ns", type=float, default=0.05)
    parser.add_argument("--refine-rounds", type=int, default=3)
    parser.add_argument("--refine-half-span-steps", type=int, default=25)
    parser.add_argument("--refine-factor", type=float, default=10.0)
    parser.add_argument("--phase-bin-width-ns", type=float, default=0.25)
    parser.add_argument("--peak-window-ns", type=float, default=8.0)
    parser.add_argument("--min-selected-pulses", type=int, default=1000)
    parser.add_argument("--score-stride", type=int, default=1)
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

    analyzer_source = analysis_dir / "step04_dogma_cleaned_all_channel_rf_period.cpp"
    analyzer_executable = build_dir / "dogma_cleaned_all_channel_rf_period"
    writer_path = analysis_dir / "step04_write_cleaned_all_channel_rf_period_outputs.py"
    scan_writer_path = analysis_dir / "step04_write_rf_period_scan_outputs.py"

    compile_cpp(analyzer_source, analyzer_executable, args.compiler)
    run_command([sys.executable, "-m", "py_compile", str(writer_path)])
    run_command([sys.executable, "-m", "py_compile", str(scan_writer_path)])

    run_command([
        str(analyzer_executable),
        "--input",
        str(input_path),
        "--output-root",
        str(output_root),
        "--run-key",
        run_key,
        "--ch0ref-min-ns",
        str(args.ch0ref_min_ns),
        "--ch0ref-max-ns",
        str(args.ch0ref_max_ns),
        "--ch0ref-bins",
        str(args.ch0ref_bins),
        "--trigger-ref-min-ns",
        str(args.trigger_ref_min_ns),
        "--trigger-ref-max-ns",
        str(args.trigger_ref_max_ns),
        "--trigger-ref-bins",
        str(args.trigger_ref_bins),
        "--tot-min-ns",
        str(args.tot_min_ns),
        "--tot-max-ns",
        str(args.tot_max_ns),
        "--tot-bins",
        str(args.tot_bins),
        "--score-channel",
        str(args.score_channel),
        "--score-tot-min-ns",
        str(args.score_tot_min_ns),
        "--score-tot-max-ns",
        str(args.score_tot_max_ns),
        "--initial-period-min-ns",
        str(args.initial_period_min_ns),
        "--initial-period-max-ns",
        str(args.initial_period_max_ns),
        "--initial-step-ns",
        str(args.initial_step_ns),
        "--refine-rounds",
        str(args.refine_rounds),
        "--refine-half-span-steps",
        str(args.refine_half_span_steps),
        "--refine-factor",
        str(args.refine_factor),
        "--phase-bin-width-ns",
        str(args.phase_bin_width_ns),
        "--peak-window-ns",
        str(args.peak_window_ns),
        "--min-selected-pulses",
        str(args.min_selected_pulses),
        "--score-stride",
        str(args.score_stride),
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
    ch0ref_file = run_dir / "Ch0_ref_Rates" / f"{run_key}_Ch0_ref_Rates_scan_matrix.tsv"
    folded_file = run_dir / "Folded_RF" / f"{run_key}_Folded_RF_matrix.tsv"
    folded3x_file = run_dir / "Folded_RF_3x" / f"{run_key}_Folded_RF_3x_matrix.tsv"
    folded_phase_tot_file = run_dir / "Folded_RF" / f"{run_key}_Folded_RF_phase_tot_sparse.tsv"
    folded3x_phase_tot_file = run_dir / "Folded_RF_3x" / f"{run_key}_Folded_RF_3x_phase_tot_sparse.tsv"
    folded_phase_ch0_time_file = run_dir / "Folded_RF" / f"{run_key}_Folded_RF_phase_vs_ch0_ref_time_sparse.tsv"
    folded_phase_trigger_time_file = run_dir / "Folded_RF" / f"{run_key}_Folded_RF_phase_vs_trigger_time_sparse.tsv"
    folded3x_phase_ch0_time_file = run_dir / "Folded_RF_3x" / f"{run_key}_Folded_RF_3x_phase_vs_ch0_ref_time_sparse.tsv"
    folded3x_phase_trigger_time_file = run_dir / "Folded_RF_3x" / f"{run_key}_Folded_RF_3x_phase_vs_trigger_time_sparse.tsv"
    scan_dir = run_dir / "RF_period_scan"
    scan_prefix = scan_dir / f"{run_key}_rf_period_scan"
    run_command([
        sys.executable,
        str(writer_path),
        "--ch0ref-matrix-file",
        str(ch0ref_file),
        "--folded-matrix-file",
        str(folded_file),
        "--folded3x-matrix-file",
        str(folded3x_file),
        "--folded-phase-tot-file",
        str(folded_phase_tot_file),
        "--folded3x-phase-tot-file",
        str(folded3x_phase_tot_file),
        "--folded-phase-ch0-time-file",
        str(folded_phase_ch0_time_file),
        "--folded-phase-trigger-time-file",
        str(folded_phase_trigger_time_file),
        "--folded3x-phase-ch0-time-file",
        str(folded3x_phase_ch0_time_file),
        "--folded3x-phase-trigger-time-file",
        str(folded3x_phase_trigger_time_file),
    ])
    run_command([
        sys.executable,
        str(scan_writer_path),
        "--input-prefix",
        str(scan_prefix),
        "--pdf-prefix",
        str(scan_prefix),
        "--png-prefix",
        str(scan_prefix),
        "--root-output",
        str(scan_dir / f"{run_key}_rf_period_scan.root"),
    ])

    odd_stats_path = run_dir / f"{run_key}_odd_ch0_reference_stats.txt"
    odd_samples_path = run_dir / f"{run_key}_odd_ch0_reference_samples.tsv"
    if not can_reuse_odd_ch0_stats(odd_stats_path, odd_samples_path, input_path, args):
        odd_stats_path, odd_samples_path = write_odd_ch0_stats(input_path, run_dir, run_key, metadata, args)
    print(f"Odd ch0 stats: {odd_stats_path}")
    print(f"Odd ch0 samples: {odd_samples_path}")


if __name__ == "__main__":
    main()