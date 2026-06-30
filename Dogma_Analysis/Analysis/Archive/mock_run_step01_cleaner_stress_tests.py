from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
import shutil
import subprocess
import sys


DEFAULT_ANALYSIS_SUFFIX = "rawRefined"


def format_float(value: float) -> str:
    return f"{value:.6f}".rstrip("0").rstrip(".")


def pass_fail(condition: bool) -> str:
    return "PASS" if condition else "FAIL"


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


def to_int(values: dict[str, str], key: str) -> int:
    return int(float(values[key]))


def to_float(values: dict[str, str], key: str) -> float:
    return float(values[key])


@dataclass(frozen=True)
class FloatExpectation:
    expected: float
    tolerance: float = 1e-6


@dataclass(frozen=True)
class PulseExpectation:
    window_index: int
    tdc_ordinal: int
    global_trigger_seconds: float | None
    channel: int
    rise_tdc_index_one_based: int
    fall_tdc_index_one_based: int
    rise_ns: float
    fall_ns: float
    tot_ns: float


@dataclass(frozen=True)
class StressCase:
    name: str
    description: str
    input_lines: list[str]
    cleaner_args: list[str] = field(default_factory=list)
    expected_summary_ints: dict[str, int] = field(default_factory=dict)
    expected_summary_floats: dict[str, FloatExpectation] = field(default_factory=dict)
    expected_summary_strings: dict[str, str] = field(default_factory=dict)
    expected_pulses: list[PulseExpectation] = field(default_factory=list)


@dataclass
class CheckResult:
    label: str
    status: str
    expected: str
    actual: str
    notes: str = ""


@dataclass
class CaseResult:
    name: str
    description: str
    status: str
    raw_input_path: Path
    output_dir: Path
    summary_path: Path
    pulse_table_path: Path
    validation_path: Path
    anomaly_path: Path
    checks: list[CheckResult]
    accepted_pulses: list[PulseExpectation]
    stdout: str
    stderr: str


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_output_root = script_dir.parent / "Results" / "Mock_step01_cleaner_stress"

    parser = argparse.ArgumentParser(
        description="Run step01 cleaner stress tests on synthetic raw DOGMA text streams and write a plain-text QC report."
    )
    parser.add_argument("--output-root", default=str(default_output_root))
    parser.add_argument("--report-name", default="00_step01_cleaner_stress_report.txt")
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--analysis-suffix", default=DEFAULT_ANALYSIS_SUFFIX)
    return parser.parse_args()


def build_cases() -> list[StressCase]:
    return [
        StressCase(
            name="nominal_pairs",
            description="Two clean windows with only valid rise/fall pairs across two channels.",
            input_lines=[
                "TDC 1 of total 1",
                "GlobalTriggerTime 0.001",
                "1 2 1 100",
                "2 2 0 112",
                "3 5 1 300",
                "4 5 0 316",
                "TDC 1 of total 1",
                "GlobalTriggerTime 0.002",
                "1 2 1 -50",
                "2 2 0 -38",
            ],
            expected_summary_ints={
                "parsed_headers": 2,
                "parsed_rows": 6,
                "processed_rows": 6,
                "processed_windows_seen": 2,
                "finalized_windows": 2,
                "refined_windows": 2,
                "accepted_pairs": 3,
                "total_rise_entries": 3,
                "total_fall_entries": 3,
                "orphaned_rises": 0,
                "orphaned_falls": 0,
                "tot_below_min": 0,
                "tot_above_max": 0,
                "invalid_hit_time_rows": 0,
                "duplicate_restart_count": 0,
            },
            expected_summary_strings={
                "duplicate_restart_detected": "false",
            },
            expected_pulses=[
                PulseExpectation(1, 1, 0.001, 2, 1, 2, 100.0, 112.0, 12.0),
                PulseExpectation(1, 1, 0.001, 5, 3, 4, 300.0, 316.0, 16.0),
                PulseExpectation(2, 1, 0.002, 2, 1, 2, -50.0, -38.0, 12.0),
            ],
        ),
        StressCase(
            name="anomaly_mix",
            description="Exercises orphan fall, replaced rise, ToT min/max rejection, invalid hit time, unexpected edge, and window-end orphaning.",
            cleaner_args=["--tot-max-ns", "20"],
            input_lines=[
                "TDC 1 of total 1",
                "GlobalTriggerTime 1.0",
                "1 2 0 40.0",
                "2 2 1 50.0",
                "3 2 1 55.0",
                "4 2 0 68.0",
                "5 3 1 100.0",
                "6 3 0 100.4",
                "7 4 1 200.0",
                "8 4 0 250.0",
                "9 5 1 1000000.0",
                "10 6 3 10.0",
                "11 7 1 400.0",
            ],
            expected_summary_ints={
                "parsed_headers": 1,
                "parsed_rows": 11,
                "processed_rows": 11,
                "processed_windows_seen": 1,
                "finalized_windows": 1,
                "refined_windows": 1,
                "accepted_pairs": 1,
                "total_rise_entries": 6,
                "total_fall_entries": 4,
                "orphaned_rises": 2,
                "orphaned_falls": 1,
                "tot_below_min": 1,
                "tot_above_max": 1,
                "invalid_hit_time_rows": 1,
                "duplicate_restart_count": 0,
            },
            expected_summary_strings={
                "duplicate_restart_detected": "false",
            },
            expected_pulses=[
                PulseExpectation(1, 1, 1.0, 2, 3, 4, 55.0, 68.0, 13.0),
            ],
        ),
        StressCase(
            name="all_rejected_window",
            description="A window with no accepted pairs; checks that refined output stays empty while summary counters remain correct.",
            input_lines=[
                "TDC 1 of total 1",
                "GlobalTriggerTime 3.0",
                "1 2 1 0.0",
                "2 2 0 0.5",
                "3 2 1 10.0",
            ],
            expected_summary_ints={
                "parsed_headers": 1,
                "parsed_rows": 3,
                "processed_rows": 3,
                "processed_windows_seen": 1,
                "finalized_windows": 1,
                "refined_windows": 0,
                "accepted_pairs": 0,
                "total_rise_entries": 2,
                "total_fall_entries": 1,
                "orphaned_rises": 1,
                "orphaned_falls": 0,
                "tot_below_min": 1,
                "tot_above_max": 0,
                "invalid_hit_time_rows": 0,
                "duplicate_restart_count": 0,
            },
            expected_summary_strings={
                "duplicate_restart_detected": "false",
            },
            expected_pulses=[],
        ),
        StressCase(
            name="restart_stitching",
            description="Backward trigger jump must be detected and stitched into a monotonic adjusted trigger axis.",
            input_lines=[
                "TDC 1 of total 1",
                "GlobalTriggerTime 10.0",
                "1 2 1 10.0",
                "2 2 0 22.0",
                "TDC 1 of total 1",
                "GlobalTriggerTime 8.0",
                "1 2 1 30.0",
                "2 2 0 44.0",
                "TDC 1 of total 1",
                "GlobalTriggerTime 8.5",
                "1 3 1 50.0",
                "2 3 0 60.0",
            ],
            expected_summary_ints={
                "parsed_headers": 3,
                "parsed_rows": 6,
                "processed_rows": 6,
                "processed_windows_seen": 3,
                "finalized_windows": 3,
                "refined_windows": 3,
                "accepted_pairs": 3,
                "total_rise_entries": 3,
                "total_fall_entries": 3,
                "orphaned_rises": 0,
                "orphaned_falls": 0,
                "tot_below_min": 0,
                "tot_above_max": 0,
                "invalid_hit_time_rows": 0,
                "duplicate_restart_count": 1,
            },
            expected_summary_floats={
                "duplicate_restart_time_offset_seconds": FloatExpectation(2.0),
                "restart_event_1_applied_offset_seconds": FloatExpectation(2.0),
                "restart_event_1_new_adjusted_trigger_seconds": FloatExpectation(10.0),
            },
            expected_summary_strings={
                "duplicate_restart_detected": "true",
            },
            expected_pulses=[
                PulseExpectation(1, 1, 10.0, 2, 1, 2, 10.0, 22.0, 12.0),
                PulseExpectation(2, 1, 10.0, 2, 1, 2, 30.0, 44.0, 14.0),
                PulseExpectation(3, 1, 10.5, 3, 1, 2, 50.0, 60.0, 10.0),
            ],
        ),
        StressCase(
            name="tdc_filtering",
            description="Only the selected signal TDC windows should be processed while trigger times are inherited from the trigger TDC blocks.",
            cleaner_args=["--signal-tdc-ordinal", "2", "--trigger-tdc-ordinal", "1"],
            input_lines=[
                "TDC 1 of total 2",
                "GlobalTriggerTime 1.0",
                "TDC 2 of total 2",
                "1 2 1 10.0",
                "2 2 0 20.0",
                "TDC 1 of total 2",
                "GlobalTriggerTime 2.0",
                "TDC 2 of total 2",
                "1 2 1 30.0",
                "2 2 0 41.0",
            ],
            expected_summary_ints={
                "parsed_headers": 4,
                "parsed_rows": 4,
                "processed_rows": 4,
                "processed_windows_seen": 2,
                "finalized_windows": 2,
                "refined_windows": 2,
                "accepted_pairs": 2,
                "total_rise_entries": 2,
                "total_fall_entries": 2,
                "orphaned_rises": 0,
                "orphaned_falls": 0,
                "tot_below_min": 0,
                "tot_above_max": 0,
                "invalid_hit_time_rows": 0,
                "duplicate_restart_count": 0,
            },
            expected_summary_strings={
                "duplicate_restart_detected": "false",
                "processed_tdc_ordinal": "2",
                "trigger_tdc_ordinal": "1",
            },
            expected_pulses=[
                PulseExpectation(1, 2, 1.0, 2, 1, 2, 10.0, 20.0, 10.0),
                PulseExpectation(2, 2, 2.0, 2, 1, 2, 30.0, 41.0, 11.0),
            ],
        ),
    ]


def compile_cleaner(source_path: Path, binary_path: Path, compiler: str) -> None:
    binary_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            compiler,
            "-std=c++17",
            "-O2",
            "-Wall",
            "-Wextra",
            "-pedantic",
            str(source_path),
            "-o",
            str(binary_path),
        ],
        check=True,
    )


def write_raw_input(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_pulse_table(path: Path) -> list[PulseExpectation]:
    pulses: list[PulseExpectation] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("window_index"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                raise ValueError(f"Expected 9 columns in pulse table row: {line}")
            trigger_text = parts[2].strip()
            pulses.append(
                PulseExpectation(
                    window_index=int(parts[0]),
                    tdc_ordinal=int(parts[1]),
                    global_trigger_seconds=float(trigger_text) if trigger_text else None,
                    channel=int(parts[3]),
                    rise_tdc_index_one_based=int(parts[4]),
                    fall_tdc_index_one_based=int(parts[5]),
                    rise_ns=float(parts[6]),
                    fall_ns=float(parts[7]),
                    tot_ns=float(parts[8]),
                )
            )
    return pulses


def pulse_to_text(pulse: PulseExpectation) -> str:
    trigger = "" if pulse.global_trigger_seconds is None else format_float(pulse.global_trigger_seconds)
    return (
        f"window={pulse.window_index}, tdc={pulse.tdc_ordinal}, trig={trigger}, ch={pulse.channel}, "
        f"rise_idx={pulse.rise_tdc_index_one_based}, fall_idx={pulse.fall_tdc_index_one_based}, "
        f"rise_ns={format_float(pulse.rise_ns)}, fall_ns={format_float(pulse.fall_ns)}, tot_ns={format_float(pulse.tot_ns)}"
    )


def compare_pulses(expected: list[PulseExpectation], actual: list[PulseExpectation]) -> tuple[bool, str, str, str]:
    if len(expected) != len(actual):
        return False, str(len(expected)), str(len(actual)), "accepted pulse count differs"
    tolerance = 1e-6
    for index, (exp, act) in enumerate(zip(expected, actual), start=1):
        numeric_match = (
            exp.window_index == act.window_index
            and exp.tdc_ordinal == act.tdc_ordinal
            and exp.channel == act.channel
            and exp.rise_tdc_index_one_based == act.rise_tdc_index_one_based
            and exp.fall_tdc_index_one_based == act.fall_tdc_index_one_based
            and ((exp.global_trigger_seconds is None and act.global_trigger_seconds is None) or (exp.global_trigger_seconds is not None and act.global_trigger_seconds is not None and abs(exp.global_trigger_seconds - act.global_trigger_seconds) <= tolerance))
            and abs(exp.rise_ns - act.rise_ns) <= tolerance
            and abs(exp.fall_ns - act.fall_ns) <= tolerance
            and abs(exp.tot_ns - act.tot_ns) <= tolerance
        )
        if not numeric_match:
            return False, pulse_to_text(exp), pulse_to_text(act), f"first differing pulse at index {index}"
    return True, str(len(expected)), str(len(actual)), "exact pulse-table match"


def build_checks(case: StressCase, summary: dict[str, str], pulses: list[PulseExpectation], required_paths: dict[str, Path]) -> list[CheckResult]:
    checks: list[CheckResult] = []

    for label, path in required_paths.items():
        checks.append(CheckResult(label, pass_fail(path.exists()), "present", "present" if path.exists() else "missing", "output artifact"))

    for key, expected in case.expected_summary_ints.items():
        if key not in summary:
            checks.append(CheckResult(key, "FAIL", str(expected), "missing", "summary key missing"))
            continue
        actual = to_int(summary, key)
        checks.append(CheckResult(key, pass_fail(actual == expected), str(expected), str(actual), "summary integer"))

    for key, expected in case.expected_summary_strings.items():
        actual = summary.get(key, "missing")
        checks.append(CheckResult(key, pass_fail(actual == expected), expected, actual, "summary string"))

    for key, expectation in case.expected_summary_floats.items():
        if key not in summary:
            checks.append(CheckResult(key, "FAIL", format_float(expectation.expected), "missing", "summary float missing"))
            continue
        actual = to_float(summary, key)
        condition = abs(actual - expectation.expected) <= expectation.tolerance
        checks.append(
            CheckResult(
                key,
                pass_fail(condition),
                f"{format_float(expectation.expected)} +/- {format_float(expectation.tolerance)}",
                format_float(actual),
                "summary float",
            )
        )

    pulse_match, expected_text, actual_text, note = compare_pulses(case.expected_pulses, pulses)
    checks.append(CheckResult("accepted_pulse_rows", pass_fail(pulse_match), expected_text, actual_text, note))
    return checks


def run_case(binary_path: Path, case: StressCase, case_root: Path, analysis_suffix: str) -> CaseResult:
    raw_input_path = case_root / f"{case.name}.dld.dat"
    write_raw_input(raw_input_path, case.input_lines)

    results_root = case_root / "results"
    output_dir = results_root / f"{case.name}_{analysis_suffix}"
    if results_root.exists():
        shutil.rmtree(results_root)
    results_root.mkdir(parents=True, exist_ok=True)

    command = [
        str(binary_path),
        "--input",
        str(raw_input_path),
        "--output-root",
        str(results_root),
        "--run-key",
        case.name,
        "--analysis-suffix",
        analysis_suffix,
        "--progress-report-interval-seconds",
        "3600",
    ] + case.cleaner_args

    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"Cleaner failed for case {case.name} with exit code {completed.returncode}\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )

    prefix = f"{case.name}_{analysis_suffix}"
    summary_path = output_dir / f"{prefix}_summary.txt"
    pulse_table_path = output_dir / f"{prefix}_pulses.tsv"
    validation_path = output_dir / f"{prefix}_validation.tsv"
    anomaly_path = output_dir / f"{prefix}_anomaly_samples.tsv"
    required_paths = {
        "summary_file": summary_path,
        "pulse_table": pulse_table_path,
        "validation_table": validation_path,
        "anomaly_table": anomaly_path,
    }

    summary = read_key_value_file(summary_path)
    pulses = read_pulse_table(pulse_table_path)
    checks = build_checks(case, summary, pulses, required_paths)
    status = "PASS" if all(check.status == "PASS" for check in checks) else "FAIL"
    return CaseResult(
        name=case.name,
        description=case.description,
        status=status,
        raw_input_path=raw_input_path,
        output_dir=output_dir,
        summary_path=summary_path,
        pulse_table_path=pulse_table_path,
        validation_path=validation_path,
        anomaly_path=anomaly_path,
        checks=checks,
        accepted_pulses=pulses,
        stdout=completed.stdout.strip(),
        stderr=completed.stderr.strip(),
    )


def render_summary_table(results: list[CaseResult]) -> list[str]:
    rows = [("Case", "Status", "Passed", "Total", "Accepted pulses")]
    for result in results:
        passed = sum(check.status == "PASS" for check in result.checks)
        total = len(result.checks)
        rows.append((result.name, result.status, str(passed), str(total), str(len(result.accepted_pulses))))

    widths = [max(len(row[index]) for row in rows) for index in range(len(rows[0]))]
    lines = []
    for row_index, row in enumerate(rows):
        lines.append("  ".join(value.ljust(widths[index]) for index, value in enumerate(row)))
        if row_index == 0:
            lines.append("  ".join("-" * width for width in widths))
    return lines


def render_case_checks(result: CaseResult) -> list[str]:
    rows = [("Check", "Status", "Expected", "Actual", "Notes")]
    for check in result.checks:
        rows.append((check.label, check.status, check.expected, check.actual, check.notes))

    widths = [max(len(row[index]) for row in rows) for index in range(len(rows[0]))]
    lines = []
    for row_index, row in enumerate(rows):
        lines.append("  ".join(value.ljust(widths[index]) for index, value in enumerate(row)))
        if row_index == 0:
            lines.append("  ".join("-" * width for width in widths))
    return lines


def write_report(output_path: Path, binary_path: Path, results: list[CaseResult]) -> None:
    overall_status = "PASS" if all(result.status == "PASS" for result in results) else "FAIL"
    lines = [
        "Step01 Cleaner Stress-Test Report",
        "=" * 33,
        "",
        f"Overall status: {overall_status}",
        f"Cleaner binary: {binary_path}",
        f"Cases run: {len(results)}",
        "",
        "Top-level summary:",
        *render_summary_table(results),
        "",
    ]

    for index, result in enumerate(results, start=1):
        lines.extend(
            [
                f"Case {index}: {result.name}",
                "-" * (len(result.name) + 8),
                f"Description: {result.description}",
                f"Status: {result.status}",
                f"Raw input: {result.raw_input_path}",
                f"Output dir: {result.output_dir}",
                f"Summary file: {result.summary_path}",
                f"Pulse table: {result.pulse_table_path}",
                f"Validation table: {result.validation_path}",
                f"Anomaly table: {result.anomaly_path}",
                "",
                *render_case_checks(result),
                "",
                "Accepted pulses:",
            ]
        )
        if result.accepted_pulses:
            lines.extend(f"  - {pulse_to_text(pulse)}" for pulse in result.accepted_pulses)
        else:
            lines.append("  - none")
        if result.stdout:
            lines.extend(["", "Cleaner stdout:", *[f"  {line}" for line in result.stdout.splitlines()]])
        if result.stderr:
            lines.extend(["", "Cleaner stderr:", *[f"  {line}" for line in result.stderr.splitlines()]])
        lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    source_path = script_dir / "step01_dogma_raw_data_refiner.cpp"
    output_root = Path(args.output_root).resolve()
    binary_path = output_root / "_bin" / "step01_dogma_raw_data_refiner_stress"

    cases = build_cases()
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    try:
        compile_cleaner(source_path, binary_path, args.compiler)
        results = []
        for index, case in enumerate(cases, start=1):
            case_root = output_root / f"{index:02d}_{case.name}"
            results.append(run_case(binary_path, case, case_root, args.analysis_suffix))
    except subprocess.CalledProcessError as error:
        print(f"Compile failed with exit code {error.returncode}: {error.cmd}", file=sys.stderr)
        return error.returncode
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    report_path = output_root / args.report_name
    write_report(report_path, binary_path, results)
    overall_status = "PASS" if all(result.status == "PASS" for result in results) else "FAIL"
    print(f"Report: {report_path}")
    print(f"Overall status: {overall_status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())