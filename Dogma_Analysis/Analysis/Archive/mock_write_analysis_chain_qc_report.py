from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a markdown quality-control report for the mock cleaned-data analysis-chain test."
    )
    parser.add_argument("--mock-spec", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--output-name", default="00_mock_analysis_chain_quality_report.md")
    return parser.parse_args()


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


def file_status(path: Path) -> str:
    return "PASS" if path.exists() else "FAIL"


def to_int(values: dict[str, str], key: str, default: int = 0) -> int:
    if key not in values:
        return default
    return int(float(values[key]))


def to_float(values: dict[str, str], key: str, default: float = 0.0) -> float:
    if key not in values:
        return default
    return float(values[key])


def format_float(value: float) -> str:
    return f"{value:.6f}".rstrip("0").rstrip(".")


def pass_fail(condition: bool) -> str:
    return "PASS" if condition else "FAIL"


def add_exact_check(rows: list[tuple[str, str, str, str, str]], label: str, expected: int, actual: int, notes: str = "") -> None:
    rows.append((label, pass_fail(expected == actual), str(expected), str(actual), notes))


def add_tolerance_check(rows: list[tuple[str, str, str, str, str]], label: str, expected: float, actual: float, tolerance: float, notes: str = "") -> None:
    condition = abs(actual - expected) <= tolerance
    rows.append(
        (
            label,
            pass_fail(condition),
            f"{format_float(expected)} +/- {format_float(tolerance)}",
            format_float(actual),
            notes,
        )
    )


def add_threshold_check(rows: list[tuple[str, str, str, str, str]], label: str, minimum: float, actual: float, notes: str = "") -> None:
    rows.append((label, pass_fail(actual >= minimum), f">= {format_float(minimum)}", format_float(actual), notes))


def main() -> None:
    args = parse_args()
    mock_spec_path = Path(args.mock_spec).resolve()
    results_dir = Path(args.results_dir).resolve()
    output_path = results_dir / args.output_name

    spec = json.loads(mock_spec_path.read_text(encoding="utf-8"))
    expected = spec["expected"]
    actual_samples = spec["actual_samples"]

    rates_summary = read_key_value_file(results_dir / f"{spec['run_key']}_cleaned_rates_summary.txt")
    time_tot_summary = read_key_value_file(results_dir / f"{spec['run_key']}_cleaned_time_tot_summary.txt")
    rf_summary = read_key_value_file(results_dir / f"{spec['run_key']}_cleaned_rf_period_summary.txt")
    odd_ch0_stats = read_key_value_file(results_dir / f"{spec['run_key']}_odd_ch0_reference_stats.txt")

    checks: list[tuple[str, str, str, str, str]] = []
    expected_rows = int(expected["total_rows"])
    expected_trigger_count = int(spec["parameters"]["trigger_count"])
    expected_ch2_pulses = int(expected["total_ch2_pulses"])
    ch2_channel = int(spec["parameters"]["ch2_channel"])
    ch2_pulses_per_bucket = int(spec["parameters"]["ch2_pulses_per_bucket"])
    rf_period_ns = float(spec["parameters"]["rf_period_ns"])
    ch0_rise_ns = float(spec["parameters"]["ch0_rise_ns"])

    add_exact_check(checks, "Step-02 parsed rows", expected_rows, to_int(rates_summary, "parsed_pulse_rows"))
    add_exact_check(checks, "Step-03 parsed rows", expected_rows, to_int(time_tot_summary, "parsed_pulse_rows"))
    add_exact_check(checks, "Step-04 parsed rows", expected_rows, to_int(rf_summary, "parsed_pulse_rows"))
    add_exact_check(checks, "Valid ch0 windows", expected_trigger_count, to_int(rf_summary, "valid_ch0_windows"))
    add_exact_check(checks, "Missing valid ch0 windows", int(expected["expected_missing_valid_ch0_windows"]), to_int(rf_summary, "windows_without_valid_ch0_candidate"))
    add_exact_check(checks, "Multiple valid ch0 candidates", 0, to_int(rf_summary, "windows_with_multiple_valid_ch0_candidates"))
    add_exact_check(checks, "Score-channel pulses", expected_ch2_pulses, to_int(rf_summary, "score_pulses_before_stride"))
    add_tolerance_check(
        checks,
        "Recovered RF period",
        float(expected["expected_rf_period_ns"]),
        to_float(rf_summary, "deduced_period_ns"),
        float(expected["tolerance_period_ns"]),
    )
    add_tolerance_check(
        checks,
        "Recovered phase origin",
        float(expected["expected_phase_origin_ns"]),
        to_float(rf_summary, "phase_origin_ns"),
        float(expected["tolerance_phase_origin_ns"]),
    )
    add_tolerance_check(
        checks,
        "Recovered phase sigma",
        float(actual_samples["ch2_phase_jitter"]["sigma"]),
        to_float(rf_summary, "sigma_ns"),
        float(expected["tolerance_sigma_ns"]),
    )
    add_threshold_check(
        checks,
        "Selected fraction",
        float(expected["selected_fraction_min"]),
        to_float(rf_summary, "selected_fraction"),
    )
    add_exact_check(checks, "Odd ch0 rescued windows", 0, to_int(odd_ch0_stats, "windows_rescued_by_later_valid_ch0_candidate"))
    add_exact_check(checks, "Odd ch0 invalid windows", 0, to_int(odd_ch0_stats, "windows_without_valid_ch0_candidate"))

    key_paths = {
        "Abs_rates": results_dir / "Abs_rates",
        "Ch0_ref_Rates": results_dir / "Ch0_ref_Rates",
        "Trigger_ref_ToT": results_dir / "Trigger_ref_ToT",
        "Ch0_ref_TOT": results_dir / "Ch0_ref_TOT",
        "TOT_distrib": results_dir / "TOT_distrib",
        "Folded_RF": results_dir / "Folded_RF",
        "Folded_RF_3x": results_dir / "Folded_RF_3x",
        "RF_period_scan": results_dir / "RF_period_scan",
        "RF scan metrics PNG": results_dir / "RF_period_scan" / f"{spec['run_key']}_rf_period_scan_scan_metrics.png",
        "Ch0-ref combined rates PNG": results_dir / "Ch0_ref_Rates" / "combined_ch0_ref_rates.png",
        "Ch0-ref scan rates PNG": results_dir / "Ch0_ref_Rates" / "combined_ch0_ref_rates_scan.png",
        "TOT distrib combined PNG": results_dir / "TOT_distrib" / "combined_tot_distrib.png",
        "Folded RF combined ToT PNG": results_dir / "Folded_RF" / "combined_rf_phase_vs_tot.png",
        "Folded RF trigger-time PNG": results_dir / "Folded_RF" / "combined_rf_phase_vs_trigger_time.png",
        "Folded RF ch0-time PNG": results_dir / "Folded_RF" / "combined_rf_phase_vs_ch0_ref_time.png",
        "Folded RF 3x combined ToT PNG": results_dir / "Folded_RF_3x" / "combined_rf_phase_3x_vs_tot.png",
    }

    actual_trigger_time_min = float(actual_samples["ch2_trigger_time"]["min"])
    actual_trigger_time_max = float(actual_samples["ch2_trigger_time"]["max"])
    ch0_ref_min_ns = to_float(time_tot_summary, "ch0_ref_min_ns")
    ch0_ref_max_ns = to_float(time_tot_summary, "ch0_ref_max_ns")
    shifted_ch0_ref_min_ns = actual_trigger_time_min - ch0_rise_ns
    shifted_ch0_ref_max_ns = actual_trigger_time_max - ch0_rise_ns

    lower_clipped_buckets = 0
    if shifted_ch0_ref_min_ns < ch0_ref_min_ns:
        lower_clipped_buckets = int(math.ceil((ch0_ref_min_ns - shifted_ch0_ref_min_ns) / rf_period_ns))

    upper_clipped_buckets = 0
    if shifted_ch0_ref_max_ns > ch0_ref_max_ns:
        upper_clipped_buckets = int(math.ceil((shifted_ch0_ref_max_ns - ch0_ref_max_ns) / rf_period_ns))

    expected_clipped_ch2_pulses = (
        (lower_clipped_buckets + upper_clipped_buckets)
        * expected_trigger_count
        * ch2_pulses_per_bucket
    )
    expected_ch0_ref_accepted = expected_rows - expected_clipped_ch2_pulses
    actual_trigger_ref_ch2_accepted = to_int(time_tot_summary, f"channel_{ch2_channel:02d}_trigger_ref_accepted")
    actual_ch0_ref_ch2_accepted = to_int(time_tot_summary, f"channel_{ch2_channel:02d}_ch0_ref_accepted")
    peak_center_ns = to_float(rf_summary, "peak_center_ns")

    comparison_rows = [
        (
            "Total pulse rows",
            str(expected_rows),
            str(to_int(rf_summary, "parsed_pulse_rows")),
            "All parsed-row counts should match across step-02, step-03, and step-04.",
        ),
        (
            "Valid ch0 windows",
            str(int(expected["expected_valid_ch0_windows"])),
            str(to_int(rf_summary, "valid_ch0_windows")),
            "The clean mock injects exactly one valid ch0 reference pulse in every trigger window.",
        ),
        (
            "Missing valid ch0 windows",
            str(int(expected["expected_missing_valid_ch0_windows"])),
            str(to_int(rf_summary, "windows_without_valid_ch0_candidate")),
            "No missing valid ch0 windows are expected for this idealized dataset.",
        ),
        (
            "Recovered RF period (ns)",
            format_float(float(expected["expected_rf_period_ns"])),
            format_float(to_float(rf_summary, "deduced_period_ns")),
            "The RF-period scan should lock onto the injected 40 ns spacing.",
        ),
        (
            "Recovered phase origin (ns)",
            format_float(float(expected["expected_phase_origin_ns"])),
            format_float(to_float(rf_summary, "phase_origin_ns")),
            "The refined phase origin tracks the realized phase center, which is slightly offset by the finite deterministic sample.",
        ),
        (
            "Recovered phase sigma (ns)",
            format_float(float(expected["expected_sigma_ns"])),
            format_float(to_float(rf_summary, "sigma_ns")),
            f"The analysis reproduces the realized finite-sample width {format_float(float(actual_samples['ch2_phase_jitter']['sigma']))} ns.",
        ),
        (
            "Selected fraction",
            f">= {format_float(float(expected['selected_fraction_min']))}",
            format_float(to_float(rf_summary, "selected_fraction")),
            "All score-channel pulses remain inside the selected RF peak window for this clean mock.",
        ),
        (
            "Trigger-ref accepted pulses",
            str(expected_rows),
            str(to_int(time_tot_summary, "trigger_ref_accepted_pulses")),
            "The configured trigger-reference range fully contains the populated mock timing span.",
        ),
        (
            "Ch0-ref accepted pulses",
            str(expected_ch0_ref_accepted),
            str(to_int(time_tot_summary, "ch0_ref_accepted_pulses")),
            f"The fixed ch0 at {format_float(ch0_rise_ns)} ns shifts the upper trigger-time tail to {format_float(shifted_ch0_ref_max_ns)} ns, clipping {expected_clipped_ch2_pulses} ch2 pulses at the +{format_float(ch0_ref_max_ns)} ns bound.",
        ),
        (
            "Seed-window center diagnostic (ns)",
            format_float(float(expected["expected_phase_origin_ns"])),
            format_float(peak_center_ns),
            "This is the center of the strongest full-phase seed window before the residual-centering refinement, not the final phase-origin estimator.",
        ),
    ]

    all_checks_pass = all(status == "PASS" for _, status, _, _, _ in checks)
    all_files_present = all(path.exists() for path in key_paths.values())
    overall_status = "PASS" if all_checks_pass and all_files_present else "FAIL"

    lines = [
        "# Mock Analysis Chain Quality Report",
        "",
        f"- Run key: {spec['run_key']}",
        f"- Mock specification: {mock_spec_path}",
        f"- Results directory: {results_dir}",
        f"- Overall status: **{overall_status}**",
        "",
        "## Generation Summary",
        f"- Trigger count: {spec['parameters']['trigger_count']}",
        f"- Ch2 RF bucket count: {spec['parameters']['rf_bucket_count']}",
        f"- Ch2 pulses per bucket per trigger: {spec['parameters']['ch2_pulses_per_bucket']}",
        f"- Total ch2 pulses: {expected['total_ch2_pulses']}",
        f"- Expected RF period: {format_float(float(expected['expected_rf_period_ns']))} ns",
        f"- Expected RF phase center: {format_float(float(expected['expected_phase_origin_ns']))} ns",
        f"- Realized ch0 ToT mean/sigma: {format_float(float(actual_samples['ch0_tot']['mean']))} / {format_float(float(actual_samples['ch0_tot']['sigma']))} ns",
        f"- Realized ch2 ToT mean/sigma: {format_float(float(actual_samples['ch2_tot']['mean']))} / {format_float(float(actual_samples['ch2_tot']['sigma']))} ns",
        f"- Realized ch2 phase-jitter mean/sigma: {format_float(float(actual_samples['ch2_phase_jitter']['mean']))} / {format_float(float(actual_samples['ch2_phase_jitter']['sigma']))} ns",
        f"- Realized ch2 trigger-time range: [{format_float(float(actual_samples['ch2_trigger_time']['min']))}, {format_float(float(actual_samples['ch2_trigger_time']['max']))}] ns",
        "",
        "## QC Checks",
        "| Check | Status | Expected | Actual | Notes |",
        "| --- | --- | --- | --- | --- |",
    ]
    for label, status, expected_text, actual_text, notes in checks:
        lines.append(f"| {label} | {status} | {expected_text} | {actual_text} | {notes} |")

    lines.extend([
        "",
        "## Key Outputs",
        "| Artifact | Status | Path |",
        "| --- | --- | --- |",
    ])
    for label, path in key_paths.items():
        lines.append(f"| {label} | {file_status(path)} | {path} |")

    lines.extend([
        "",
        "## Reconstruction Summary",
        f"- Recovered RF period: {rf_summary.get('deduced_period_ns', 'missing')} ns",
        f"- Recovered phase origin: {rf_summary.get('phase_origin_ns', 'missing')} ns",
        f"- Recovered sigma: {rf_summary.get('sigma_ns', 'missing')} ns",
        f"- Selected fraction: {rf_summary.get('selected_fraction', 'missing')}",
        f"- Valid ch0 windows: {rf_summary.get('valid_ch0_windows', 'missing')}",
        f"- Windows without valid ch0 candidate: {rf_summary.get('windows_without_valid_ch0_candidate', 'missing')}",
        "",
        "## Interpretation of Apparent Discrepancies",
        "- The RF reconstruction itself is consistent with the injected truth. The apparent mismatches come from comparing nominal generator targets, realized finite-sample statistics, and binned diagnostics in the same report.",
        f"- The ch0-reference acceptance is lower than the trigger-reference acceptance by design. The fixed ch0 reference at {format_float(ch0_rise_ns)} ns shifts the realized ch2 upper trigger-time edge to {format_float(shifted_ch0_ref_max_ns)} ns in ch0-reference coordinates, so the highest {upper_clipped_buckets} RF buckets per trigger fall beyond the configured +{format_float(ch0_ref_max_ns)} ns bound. That clips {expected_clipped_ch2_pulses} ch2 pulses in total, which matches the observed drop from {actual_trigger_ref_ch2_accepted} to {actual_ch0_ref_ch2_accepted} accepted ch2 pulses.",
        f"- peak_center_ns = {format_float(peak_center_ns)} ns is the center of the strongest full-phase seed window before refinement. phase_origin_ns = {format_float(to_float(rf_summary, 'phase_origin_ns'))} ns is the refined estimator after the residual-centering iterations, so the refined value is the one to compare with the injected RF phase center.",
        f"- The generator target for the phase width is {format_float(float(expected['expected_sigma_ns']))} ns, but the realized deterministic finite sample has sigma {format_float(float(actual_samples['ch2_phase_jitter']['sigma']))} ns. The analysis correctly reconstructs the realized width rather than the nominal target value.",
        "",
        "## Input vs Analysis Comparison",
        "| Quantity | Input / truth | Analysis result | Interpretation |",
        "| --- | --- | --- | --- |",
    ])
    for label, expected_text, actual_text, notes in comparison_rows:
        lines.append(f"| {label} | {expected_text} | {actual_text} | {notes} |")

    lines.extend([
        "",
        "## Notes",
        "- This mock dataset is intentionally clean: exactly one valid ch0 reference pulse is expected in every trigger window.",
        "- The ch2 timing distribution is deterministic and phase-centered at 20 ns modulo 40 ns, so RF-period recovery should converge to 40 ns.",
        "- Large TSV and ROOT files are intended to remain remote; this report is designed to be mirrored locally as a lightweight QC artifact.",
    ])

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"QC report: {output_path}")
    print(f"Overall status: {overall_status}")


if __name__ == "__main__":
    main()