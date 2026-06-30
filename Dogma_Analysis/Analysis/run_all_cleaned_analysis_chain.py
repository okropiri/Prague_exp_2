from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, fields
from pathlib import Path
import shutil
import subprocess
import sys
import time

from step01_run_dogma_raw_data_refiner import derive_run_key


RAW_REFINED_SUFFIX = "rawRefined"
COMPARISON_SOURCE_DIRS = (
    "Ch0_ref_Rates",
    "Ch0_ref_TOT",
    "Folded_RF",
    "Folded_RF_3x",
)
PLOT_EXTENSIONS = (".png", ".pdf")


@dataclass
class RunReportRow:
    run_key: str
    raw_file: str
    status: str
    failed_step: str
    error: str
    cleaned_pulse_table: str
    results_dir: str
    cleaned_hits_summary: str
    comparison_files_copied: int
    elapsed_seconds: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Clean every DOGMA raw run, run the full cleaned-data analysis chain, "
            "write the final RF-annotated CleanedHits ROOT tree, and collect Ch02 comparison plots."
        )
    )
    parser.add_argument("--input-dir", default="/data6/Data")
    parser.add_argument("--pattern", default="*.dld.dat")
    parser.add_argument("--raw-files", nargs="*", default=[])
    parser.add_argument("--refined-root", default="/data6/Data/Refined_Data")
    parser.add_argument("--results-root", default="/data6/Dogma_analysis_by_Dachi/Results")
    parser.add_argument("--comparison-dir")
    parser.add_argument("--report-file")
    parser.add_argument("--analysis-dir")
    parser.add_argument("--python-executable", default=sys.executable)
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--build-dir", default="/tmp/dogma_cleaned_batch")
    parser.add_argument("--step04-score-channel", type=int, default=2)
    parser.add_argument("--max-runs", type=int)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--stop-on-error", action="store_true")
    return parser.parse_args()


def discover_raw_files(args: argparse.Namespace) -> list[Path]:
    if args.raw_files:
        raw_files = [Path(raw_file).expanduser().resolve() for raw_file in args.raw_files]
    else:
        input_dir = Path(args.input_dir).expanduser().resolve()
        raw_files = sorted(path.resolve() for path in input_dir.glob(args.pattern) if path.is_file())
    if args.max_runs is not None:
        raw_files = raw_files[: args.max_runs]
    return raw_files


def ensure_exists(path: Path, description: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing {description}: {path}")


def run_command(command: list[str], step_name: str) -> None:
    print(f"[{step_name}] {' '.join(command)}", flush=True)
    subprocess.run(command, check=True)


def write_report(report_path: Path, rows: list[RunReportRow]) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[field.name for field in fields(RunReportRow)],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field.name: getattr(row, field.name) for field in fields(RunReportRow)})


def remove_existing_outputs(run_key: str, refined_dir: Path, results_dir: Path, comparison_dir: Path) -> None:
    if refined_dir.exists():
        shutil.rmtree(refined_dir)
    if results_dir.exists():
        shutil.rmtree(results_dir)
    if comparison_dir.exists():
        for path in comparison_dir.glob(f"{run_key}__*"):
            path.unlink()


def collect_comparison_plots(run_key: str, results_run_dir: Path, comparison_dir: Path) -> int:
    comparison_dir.mkdir(parents=True, exist_ok=True)
    copied = 0
    for subdir_name in COMPARISON_SOURCE_DIRS:
        source_dir = results_run_dir / subdir_name
        if not source_dir.is_dir():
            continue
        for extension in PLOT_EXTENSIONS:
            for source_path in sorted(source_dir.glob(f"ch02_*{extension}")):
                destination = comparison_dir / f"{run_key}__{subdir_name}__{source_path.name}"
                shutil.copy2(source_path, destination)
                copied += 1
    return copied


def refined_dir_for_run(refined_root: Path, run_key: str) -> Path:
    return refined_root / f"{run_key}_{RAW_REFINED_SUFFIX}"


def pulse_table_for_run(refined_root: Path, run_key: str) -> Path:
    refined_dir = refined_dir_for_run(refined_root, run_key)
    prefix = f"{run_key}_{RAW_REFINED_SUFFIX}"
    return refined_dir / f"{prefix}_pulses.tsv"


def refined_summary_for_run(refined_root: Path, run_key: str) -> Path:
    refined_dir = refined_dir_for_run(refined_root, run_key)
    prefix = f"{run_key}_{RAW_REFINED_SUFFIX}"
    return refined_dir / f"{prefix}_summary.txt"


def run_chain_for_file(
    raw_file: Path,
    analysis_dir: Path,
    refined_root: Path,
    results_root: Path,
    comparison_dir: Path,
    python_executable: str,
    compiler: str,
    build_dir: Path,
    step04_score_channel: int,
    fresh: bool,
) -> RunReportRow:
    run_key = derive_run_key(raw_file)
    refined_dir = refined_dir_for_run(refined_root, run_key)
    pulse_table = pulse_table_for_run(refined_root, run_key)
    refined_summary = refined_summary_for_run(refined_root, run_key)
    results_run_dir = results_root / run_key
    rf_summary = results_run_dir / f"{run_key}_cleaned_rf_period_summary.txt"
    cleaned_hits_summary = results_run_dir / f"{run_key}_cleaned_hits_summary.txt"
    cleaned_hits_root = results_run_dir / f"{run_key}_cleaned_hits.root"
    started_at = time.monotonic()
    current_step = "setup"

    try:
        if fresh:
            remove_existing_outputs(run_key, refined_dir, results_run_dir, comparison_dir)

        refined_root.mkdir(parents=True, exist_ok=True)
        results_root.mkdir(parents=True, exist_ok=True)
        build_dir.mkdir(parents=True, exist_ok=True)

        current_step = "step01_clean"
        run_command(
            [
                python_executable,
                str(analysis_dir / "step01_run_dogma_raw_data_refiner.py"),
                "--input",
                str(raw_file),
                "--output-root",
                str(refined_root),
                "--run-key",
                run_key,
                "--analysis-suffix",
                RAW_REFINED_SUFFIX,
                "--compiler",
                compiler,
                "--binary-path",
                str(build_dir / "step01_dogma_raw_data_refiner"),
                "--python-executable",
                python_executable,
                "--skip-root-tree",
            ],
            current_step,
        )
        ensure_exists(pulse_table, "cleaned pulse table")
        ensure_exists(refined_summary, "cleaner summary")

        current_step = "step02_rates"
        run_command(
            [
                python_executable,
                str(analysis_dir / "step02_run_cleaned_all_channel_rates.py"),
                "--input",
                str(pulse_table),
                "--output-root",
                str(results_root),
                "--run-key",
                run_key,
                "--compiler",
                compiler,
                "--build-dir",
                str(build_dir),
            ],
            current_step,
        )
        ensure_exists(results_run_dir / f"{run_key}_cleaned_rates_summary.txt", "step02 summary")

        current_step = "step03_time_tot"
        run_command(
            [
                python_executable,
                str(analysis_dir / "step03_run_cleaned_all_channel_time_tot.py"),
                "--input",
                str(pulse_table),
                "--output-root",
                str(results_root),
                "--run-key",
                run_key,
                "--compiler",
                compiler,
                "--build-dir",
                str(build_dir),
            ],
            current_step,
        )
        ensure_exists(results_run_dir / f"{run_key}_cleaned_time_tot_summary.txt", "step03 summary")
        ensure_exists(results_run_dir / "Ch0_ref_TOT" / "ch02_ch0_ref_tot.png", "step03 Ch02 ch0-ref PNG")

        current_step = "step04_rf_period"
        run_command(
            [
                python_executable,
                str(analysis_dir / "step04_run_cleaned_all_channel_rf_period.py"),
                "--input",
                str(pulse_table),
                "--output-root",
                str(results_root),
                "--run-key",
                run_key,
                "--compiler",
                compiler,
                "--build-dir",
                str(build_dir),
                "--score-channel",
                str(step04_score_channel),
            ],
            current_step,
        )
        ensure_exists(rf_summary, "step04 RF summary")
        ensure_exists(results_run_dir / "Folded_RF" / "ch02_rf_phase_vs_tot.png", "step04 Ch02 RF phase-vs-ToT PNG")

        current_step = "cleaned_hits_root"
        run_command(
            [
                python_executable,
                str(analysis_dir / "step01_write_cleaned_root_tree.py"),
                "--pulse-table",
                str(pulse_table),
                "--summary-file",
                str(cleaned_hits_summary),
                "--output-root-file",
                str(cleaned_hits_root),
                "--rf-summary-file",
                str(rf_summary),
            ],
            current_step,
        )
        ensure_exists(cleaned_hits_summary, "cleaned hits summary")
        if not any(results_run_dir.glob(f"{run_key}_cleaned_hits*.root")):
            raise FileNotFoundError(f"Missing cleaned hits ROOT files in {results_run_dir}")

        current_step = "collect_comparison_plots"
        copied = collect_comparison_plots(run_key, results_run_dir, comparison_dir)
        if copied == 0:
            raise FileNotFoundError(f"No Ch02 comparison plots found for {run_key}")

        elapsed_seconds = time.monotonic() - started_at
        return RunReportRow(
            run_key=run_key,
            raw_file=str(raw_file),
            status="ok",
            failed_step="",
            error="",
            cleaned_pulse_table=str(pulse_table),
            results_dir=str(results_run_dir),
            cleaned_hits_summary=str(cleaned_hits_summary),
            comparison_files_copied=copied,
            elapsed_seconds=elapsed_seconds,
        )
    except Exception as error:
        elapsed_seconds = time.monotonic() - started_at
        return RunReportRow(
            run_key=run_key,
            raw_file=str(raw_file),
            status="failed",
            failed_step=current_step,
            error=str(error),
            cleaned_pulse_table=str(pulse_table),
            results_dir=str(results_run_dir),
            cleaned_hits_summary=str(cleaned_hits_summary),
            comparison_files_copied=0,
            elapsed_seconds=elapsed_seconds,
        )


def main() -> int:
    args = parse_args()
    analysis_dir = (
        Path(args.analysis_dir).expanduser().resolve()
        if args.analysis_dir
        else Path(__file__).resolve().parent
    )
    refined_root = Path(args.refined_root).expanduser().resolve()
    results_root = Path(args.results_root).expanduser().resolve()
    comparison_dir = (
        Path(args.comparison_dir).expanduser().resolve()
        if args.comparison_dir
        else results_root / "Ch02_Position_Comparison"
    )
    report_file = (
        Path(args.report_file).expanduser().resolve()
        if args.report_file
        else comparison_dir / "overnight_batch_report.tsv"
    )
    build_dir = Path(args.build_dir).expanduser().resolve()

    raw_files = discover_raw_files(args)
    if not raw_files:
        print("No raw files matched the requested selection.", file=sys.stderr)
        return 1

    report_rows: list[RunReportRow] = []
    for raw_file in raw_files:
        row = run_chain_for_file(
            raw_file=raw_file,
            analysis_dir=analysis_dir,
            refined_root=refined_root,
            results_root=results_root,
            comparison_dir=comparison_dir,
            python_executable=args.python_executable,
            compiler=args.compiler,
            build_dir=build_dir,
            step04_score_channel=args.step04_score_channel,
            fresh=args.fresh,
        )
        report_rows.append(row)
        write_report(report_file, report_rows)
        print(
            f"[{row.status}] {row.run_key} | step={row.failed_step or 'complete'} | "
            f"comparison_files={row.comparison_files_copied} | elapsed_s={row.elapsed_seconds:.1f}",
            flush=True,
        )
        if row.status != "ok" and args.stop_on_error:
            break

    failed = [row for row in report_rows if row.status != "ok"]
    print(f"Runs processed: {len(report_rows)}", flush=True)
    print(f"Runs failed: {len(failed)}", flush=True)
    print(f"Report: {report_file}", flush=True)
    print(f"Comparison dir: {comparison_dir}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())