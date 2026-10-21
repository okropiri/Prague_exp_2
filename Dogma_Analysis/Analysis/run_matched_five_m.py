from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import socket
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analysis_provenance import file_identity, record_arguments, record_command, write_json
from run_zebra_ch0_clock import require, verify_folded, verify_products, write_folded_outputs
from step02_run_cleaned_all_channel_rates import compile_cpp, parse_metadata
from step04_run_cleaned_all_channel_rf_period import read_key_value_file
from step04_write_cleaned_all_channel_rf_period_outputs import read_matrix


KEY = "NCAL_20us_Pos_5m_0000"
ACQUISITION = "NCAL_Test_Prescaler1024_TrigWindow_0x83e803e8_20us_Pos_5m_20260425_0000"
RAW_SHA256 = "f6d99216d9d452575286eead2859c08f7c4a7d8eda7177d5d57ed3ca02f20513"
DECODE_SHA256 = "5ce19b459accb789024be5192a21998e6769aaba2d87d62d9052d0b0620af579"
STEPS = {
    1: "01. Cleaning and validation", 2: "02. Channel rates", 3: "03. Timing and ToT",
    5: "05. Gamma-neutron background", 6: "06. Gamma-neutron components",
    7: "07. Phase-rate diagnostics", 8: "08. Gamma ToT peak fit",
    9: "09. Gamma-window sensitivity", 10: "10. Gamma-reference CFD",
}
COMMON = {
    "run_key": KEY, "board": "TDC_001049", "ordinal": 1,
    "ch0_valid_rise_min_ns": -410, "ch0_valid_rise_max_ns": -395,
    "ch0_valid_tot_min_ns": 16.5, "ch0_valid_tot_max_ns": 19.5,
    "reference_order": "first passing hit in file order",
    "cleaning_tot_ns": [1, 5000], "max_accepted_hit_time_ns": 50000,
    "ch0ref_axis_ns": [-6000, 6000], "trigger_axis_ns": [-6000, 6000],
    "time_map_bins": 2400, "step02_ch0ref_bins": 12000,
    "absolute_display_bin_seconds": 0.1,
    "tot_storage_ns": [0, 128], "tot_bins": 1280, "tot_display_ns": [0, 60],
    "phase_bin_width_target_ns": 0.25, "phase_binning": "ceil(T/0.25), 3P uses three times that bin count",
    "detailed_channels": list(range(9)), "include_all_ch0_hits": True,
    "phase_origin_ns": 0, "detector_phase_alignment": False,
    "assumed_prescaler": 1024, "coarse_quantum_ns": 5,
    "assumed_trigger_window_ns": 10000, "duty_factor_correction": False,
    "zebra_timewalk_correction": False,
}


def options(values: dict[str, object]) -> list[str]:
    return [item for key, value in values.items()
            for item in (f"--{key.replace('_', '-')}", str(value))]


def gate_options(prefix: str = "ch0-valid") -> list[str]:
    return options({key.replace("ch0_valid", prefix.replace("-", "_")): value
                    for key, value in COMMON.items() if key.startswith("ch0_valid_")})


def product(stage: Path) -> Path:
    return stage / "Results" / KEY


def inventory(folder: Path) -> list[dict[str, str | int]]:
    return [file_identity(path) for path in sorted(folder.rglob("*"))
            if path.is_file() and not path.is_symlink() and "Reports" not in path.relative_to(folder).parts]


def tracked(stage: Path, label: str, command: list[str]) -> None:
    reports = stage / "Reports"
    reports.mkdir(parents=True, exist_ok=True)
    report = reports / f"{label}.json"
    if report.exists():
        saved = json.loads(report.read_text())
        require(saved["status"] == "PASS" and saved["command"] == command,
                f"Existing failed/different command requires a new run: {report}")
        for entry in saved["outputs"]:
            path = Path(entry["path"])
            require(path.is_file() and file_identity(path)["sha256"] == entry["sha256"],
                    f"Previously verified output changed: {path}")
        print(f"REUSE verified stage {stage.name}: {label}", flush=True)
        return
    state = {"status": "RUNNING", "command": command, "started_utc": datetime.now(timezone.utc).isoformat()}
    write_json(report, state)
    env = os.environ.copy()
    env["DOGMA_STEP_REPORT_DIR"] = str(reports)
    print(f"START {stage.name}: {label}", flush=True)
    with (reports / f"{label}.log").open("x") as log:
        result = subprocess.run(["nice", "-n", "15", *command], env=env,
                                stdout=log, stderr=subprocess.STDOUT)
    state.update({"status": "PASS" if result.returncode == 0 else "FAIL",
                  "exit_code": result.returncode, "completed_utc": datetime.now(timezone.utc).isoformat(),
                  "outputs": inventory(stage) if result.returncode == 0 else []})
    write_json(report, state)
    require(result.returncode == 0, f"Failed {stage}: inspect {reports / (label + '.log')}")
    print(f"PASS {stage.name}: {label}", flush=True)


def fold_shared(output: Path, calibration: Path, pulse: Path, build: Path,
                period: float, method: str, expected_windows: int) -> None:
    require(method in ("ch0", "ch2") and np.isfinite(period) and 35 < period < 45,
            "Invalid shared RF model")
    tables = output / "Fold_Tables"
    clock = output / "Clock_Calibration"
    if clock != calibration:
        clock.mkdir()
        (clock / "ch0_reference_timestamps.tsv").symlink_to(calibration / "ch0_reference_timestamps.tsv")
    references = clock / "ch0_reference_timestamps.tsv"
    metadata = parse_metadata(references)
    ids = np.loadtxt(references, usecols=0, dtype=np.int64)
    periods = clock / "window_periods.tsv"
    if periods.exists():
        periods = clock / "folding_window_periods.tsv"
    with periods.open("x") as handle:
        handle.write("# columns=window_index local_period_ns\n")
        for index in ids:
            handle.write(f"{index}\t{period:.17g}\n")
    binary = build / "zebra_ch0_clock"
    compile_cpp(Path(__file__).with_name("zebra_ch0_clock.cpp"), binary, "g++")
    model = "shared_burst_rf_period" if method == "ch0" else "shared_ch2_detector_period"
    command = [str(binary), "fold", "--input", str(pulse), "--periods", str(periods),
               "--output", str(tables), "--run-key", KEY, "--ordinal", "1",
               *gate_options(), *options({
                   "allow_reference_gaps": 1, "expected_windows": expected_windows,
                   "period_min_ns": 35, "period_max_ns": 45, "include_all_ch0_hits": 1,
                   "period_model": model, "phase_bin_width_ns": .25,
                   "ch0ref_min_ns": -6000, "ch0ref_max_ns": 6000, "ch0ref_bins": 2400,
                   "trigger_ref_min_ns": -6000, "trigger_ref_max_ns": 6000, "trigger_ref_bins": 2400,
                   "tot_min_ns": 0, "tot_max_ns": 128, "tot_bins": 1280,
                   "histogram_channels": ",".join(map(str, range(9))),
               })]
    record_command(command)
    subprocess.run(command, check=True)
    counts = np.asarray([int(value) for value in metadata["input_channel_counts"].split(",")])
    accounting = verify_folded(tables, counts, len(ids), include_all_ch0=True)
    if method == "ch0":
        assessment = json.loads((calibration.parent / "assessment.json").read_text())
    else:
        assessment = {"board": "TDC_001049", "coarse_quantum_ns": 5,
                      "phase_origin": "first passing ch0 in file order; no detector-derived alignment",
                      "burst_clock": False, "prescaler": 1024}
    assessment.update({
        "mean_period_ns": period, "period_method": method, "include_all_ch0_hits": True,
        "decoded_selected_board_windows": expected_windows, "unresolved_reference_windows": 0,
        "folded_channel_counts": [row["folded_pulses"] for row in accounting],
        "unreferenced_channel_counts": [row["without_reference_pulses"] for row in accounting],
        "display_tot_max_ns": 60, "phase_origin_ns": 0,
    })
    write_json(output / "assessment.json", assessment)
    write_json(output / "channel_count_verification.json", accounting)
    write_folded_outputs(tables, output, f"{method}-derived T; unshifted ch0 origin", 60)
    verify_products(tables, output)
    for factor, prefix in ((1, "phase"), (3, "phase3")):
        target = output / ("Folded_RF" if factor == 1 else "Folded_RF_3x")
        suffix = "Folded_RF" if factor == 1 else "Folded_RF_3x"
        for name, quantity in (("phase_tot", "tot"), ("phase_vs_ch0_ref_time", "time"),
                               ("phase_vs_trigger_time", "trigger")):
            (target / f"{KEY}_{suffix}_{name}_sparse.tsv").symlink_to(tables / f"{prefix}_{quantity}.tsv")
    summary = {
        "deduced_period_ns": period, "phase_origin_ns": 0, "period_method": method,
        "valid_ch0_windows": len(ids), "windows_seen": expected_windows,
        "detector_phase_alignment": "false", "phase_origin_definition": assessment["phase_origin"],
        **{key: value for key, value in COMMON.items() if key.startswith("ch0_valid_")},
    }
    (output / f"{KEY}_cleaned_rf_period_summary.txt").write_text(
        "".join(f"{key}={value}\n" for key, value in summary.items()))


def verify_pair(ch0: Path, ch2: Path, output: Path) -> None:
    a, b = (json.loads((folder / "assessment.json").read_text()) for folder in (ch0, ch2))
    for key in ("folded_channel_counts", "unreferenced_channel_counts", "decoded_selected_board_windows"):
        require(a[key] == b[key], f"Matched population differs: {key}")
    output.mkdir(exist_ok=True)
    figure, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
    for folder, assessment, label in ((ch0, a, "ch0"), (ch2, b, "ch2")):
        meta, centers, counts, _ = read_matrix(folder / "Fold_Tables/folded_matrix.tsv")
        axes[0].step(centers, counts[:, 2], where="mid",
                     label=f"{label}: T={assessment['mean_period_ns']:.9f} ns")
        axes[1].step(centers / float(meta["x_max_ns"]), counts[:, 2] / counts[:, 2].sum(),
                     where="mid", label=label)
    axes[0].set(xlabel="Ch0-referenced RF phase (ns)", ylabel="NCAL counts / native bin")
    axes[1].set(xlabel="Ch0-referenced RF phase (cycles)", ylabel="NCAL fraction / native bin")
    for axis in axes:
        axis.legend()
    figure.savefig(output / "NCAL_period_method_profiles.png", dpi=180)
    figure.savefig(output / "NCAL_period_method_profiles.pdf")
    plt.close(figure)
    write_json(output / "population_comparison.json", {
        "status": "PASS_MATCHED_POPULATION", "ch0_period_ns": a["mean_period_ns"],
        "ch2_period_ns": b["mean_period_ns"],
        "period_difference_ns": b["mean_period_ns"] - a["mean_period_ns"],
        "folded_channel_counts": a["folded_channel_counts"],
        "unreferenced_channel_counts": a["unreferenced_channel_counts"],
        "alignment_applied": False,
        "profile_comparison": "Native bins labelled explicitly; cycle-axis overlay is not an equal-ns-bin statistical test.",
    })


def main() -> None:
    parser = argparse.ArgumentParser(description="Controlled 5 m Steps01–10 with frozen shared inputs and unshifted folding.")
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--verified-input-root", type=Path, required=True)
    parser.add_argument("--phase", choices=("pilot", "full"), required=True)
    parser.add_argument("--stop-after", type=int, choices=(3, 4, 10), default=10)
    parser.add_argument("--fold-only", choices=("ch0", "ch2"))
    parser.add_argument("--period-ns", type=float)
    args = record_arguments(parser.parse_args())
    require(socket.gethostname().split(".")[0] == "trbsrv-b", "DOGMA execution is restricted to padiwa")
    source = Path(__file__).resolve().parent
    stage = args.root / args.phase
    if args.fold_only:
        output = stage / f"04. {args.fold_only.capitalize()}-based RF folding"
        pulse = stage / STEPS[1] / f"Refined_Data/{KEY}_rawRefined/{KEY}_rawRefined_pulses.tsv"
        calibration = stage / "04. Ch0-based RF folding/Calibration/Clock_Calibration"
        require(args.period_ns is not None, "Explicit measured period required")
        fold_shared(output, calibration, pulse, args.root / "build",
                    args.period_ns, args.fold_only, expected_windows(args))
        return
    args.root.mkdir(parents=True, exist_ok=True)
    build = args.root / "build"
    build.mkdir(exist_ok=True)
    stage.mkdir(exist_ok=True)
    decode = args.verified_input_root / "step00_decode" / args.phase
    decoded = decode / f"{ACQUISITION}.dld.dat"
    os.environ.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MPLBACKEND="Agg")
    windows = expected_windows(args)
    manifest = stage / "Run_settings.json"
    settings = {
        "common": COMMON, "phase": args.phase, "source": [file_identity(p) for p in sorted(source.glob("*.py"))]
        + [file_identity(p) for p in sorted(source.glob("*.cpp"))],
        "decoded": file_identity(decoded), "python": sys.executable, "platform": platform.platform(),
        "packages": {name: importlib.metadata.version(name) for name in ("numpy", "matplotlib", "uproot", "pillow")},
        "compiler": subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
        "raw_tail_warning": "Final raw record lacks 292 bytes; all 10,521,678 complete events decoded.",
        "step00_validation": (decode / "ascii_validation.txt").read_text(),
        "source_control": json.loads((args.root / "source_control.json").read_text()),
        "verification_sources": [file_identity(path) for path in sorted((args.root / "verifiers").glob("*.py"))],
    }
    if args.phase == "full":
        raw = Path("/data6/Data") / f"{ACQUISITION}.dld"
        settings["raw"] = file_identity(raw)
        require(settings["raw"]["sha256"] == RAW_SHA256 and settings["decoded"]["sha256"] == DECODE_SHA256,
                "Original or decoded input identity differs; refuse legacy duplicated input")
        require((args.root / "pilot/Comparison/completion.json").exists(), "Verified pilot required")
    if manifest.exists():
        require(json.loads(manifest.read_text()) == settings, "Frozen run settings changed; use a new run")
    else:
        write_json(manifest, settings)
    shared = [stage / STEPS[index] for index in (1, 2, 3)]
    pulse_dir = shared[0] / f"Refined_Data/{KEY}_rawRefined"
    pulse = pulse_dir / f"{KEY}_rawRefined_pulses.tsv"
    gate = gate_options()
    settings_txt = {
        "run_key": KEY, "assumed_prescaler": 1024,
        "raw_input_warning": settings["raw_tail_warning"],
        "step01_dir": str(shared[0]), "step02_dir": str(shared[1]),
        "input_pulse_table": str(pulse), "ch0ref_min_ns": -6000, "ch0ref_max_ns": 6000,
        "trigger_ref_min_ns": -6000, "trigger_ref_max_ns": 6000,
        "trigger_ref_bins": 2400,
        "tot_min_ns": 0, "tot_max_ns": 128, "tot_bins": 1280,
        "tot_distrib_display_max_ns": 60, "abs_display_bin_width_seconds": .1,
        **{key: value for key, value in COMMON.items() if key.startswith("ch0_valid_")},
        **{key.replace("ch0_valid", "ch0_reference"): value for key, value in COMMON.items()
           if key.startswith("ch0_valid_")},
    }
    for index, directory in enumerate(shared, 1):
        directory.mkdir(exist_ok=True)
        if not (directory / "analysis_source").exists():
            (directory / "analysis_source").symlink_to(source, target_is_directory=True)
        (directory / "settings.txt").write_text("".join(f"{key}={value}\n" for key, value in
             {**settings_txt, "ch0ref_bins": 12000 if index == 2 else 2400}.items()))
    tracked(shared[0], "analysis", [
        sys.executable, str(source / "step01_run_dogma_raw_data_refiner.py"),
        "--input", str(decoded), "--output-root", str(shared[0] / "Refined_Data"),
        "--run-key", KEY, "--signal-tdc-ordinal", "1", "--trigger-tdc-ordinal", "1",
        "--reject-trigger-restarts", "--tot-min-ns", "1", "--tot-max-ns", "5000",
        "--max-accepted-hit-time-ns", "50000", "--root-tree-chunk-size", "250000",
        "--binary-path", str(build / "cleaner"), *gate_options("root-tree-ch0-valid")])
    tracked(shared[0], "verification", [sys.executable, str(args.root / "verifiers/verify_step01.py"),
                                      str(shared[0]), "--expected-windows", str(windows)])
    for index in (2, 3):
        directory = shared[index - 1]
        command = [sys.executable, str(source / ("step02_run_cleaned_all_channel_rates.py" if index == 2
                                                else "step03_run_cleaned_all_channel_time_tot.py")),
                   "--input", str(pulse), "--output-root", str(directory / "Results"), "--run-key", KEY,
                   "--build-dir", str(build), *gate, "--ch0ref-min-ns", "-6000", "--ch0ref-max-ns", "6000",
                   "--ch0ref-bins", str(12000 if index == 2 else 2400)]
        if index == 2:
            command += ["--abs-display-bin-width-seconds", ".1"]
        else:
            command += options({"trigger_ref_min_ns": -6000, "trigger_ref_max_ns": 6000,
                                "trigger_ref_bins": 2400, "tot_min_ns": 0, "tot_max_ns": 128,
                                "tot_bins": 1280, "tot_distrib_display_max_ns": 60,
                                "tot_map_display_max_ns": 60})
        tracked(directory, "analysis", command)
        tracked(directory, "verification", [sys.executable, str(args.root / f"verifiers/verify_step0{index}.py"),
                                            str(directory)])
    if args.stop_after == 3:
        return
    ch0 = stage / "04. Ch0-based RF folding"
    ch0.mkdir(exist_ok=True)
    tracked(ch0, "calibration", [sys.executable, str(source / "run_zebra_ch0_clock.py"),
        "--input", str(pulse), "--decoded", str(decoded), "--output", str(ch0 / "Calibration"),
        "--board", "TDC_001049", "--tdc-ordinal", "1", "--coarse-quantum-ns", "5",
        "--prescaler", "1024", "--expected-windows", str(windows),
        "--allow-reference-gaps", "--burst-clock", "--calibration-only",
        "--run-key", KEY, "--build-dir", str(build), *gate,
        "--ch0ref-min-ns", "-6000", "--ch0ref-max-ns", "6000", "--ch0ref-bins", "2400",
        "--trigger-ref-min-ns", "-6000", "--trigger-ref-max-ns", "6000", "--trigger-ref-bins", "2400",
        "--tot-min-ns", "0", "--tot-max-ns", "128", "--tot-bins", "1280",
        "--histogram-channels", "0,1,2,3,4,5,6,7,8"])
    ch2 = stage / "04. Ch2-based RF folding"
    ch2.mkdir(exist_ok=True)
    tracked(ch2, "period_search", [sys.executable, str(source / "step04_run_cleaned_all_channel_rf_period.py"),
        "--input", str(pulse), "--output-root", str(ch2 / "Period_Search"), "--run-key", KEY,
        "--build-dir", str(build), *gate, "--estimate-only", "--score-channel", "2",
        "--score-tot-min-ns", "0", "--score-tot-max-ns", "128", "--score-stride", "1",
        "--initial-period-min-ns", "37", "--initial-period-max-ns", "41", "--initial-step-ns", ".05",
        "--refine-rounds", "3", "--refine-half-span-steps", "25", "--refine-factor", "10",
        "--phase-bin-width-ns", ".25", "--peak-window-ns", "8", "--min-selected-pulses", "1000",
        "--single-rf-segment", "--skip-cycle-residual-diagnostics"])
    estimate = read_key_value_file(ch2 / f"Period_Search/{KEY}/period_estimate.txt")
    write_json(ch2 / "period_estimate.json", estimate)
    periods = {"ch0": json.loads((ch0 / "Calibration/assessment.json").read_text())["mean_period_ns"],
               "ch2": float(estimate["deduced_period_ns"])}
    for method, output in (("ch0", ch0), ("ch2", ch2)):
        tracked(output, "folding", [sys.executable, str(Path(__file__).resolve()),
            "--root", str(args.root), "--verified-input-root", str(args.verified_input_root),
            "--phase", args.phase, "--fold-only", method, "--period-ns", str(periods[method])])
        verify = [sys.executable, str(source / "verify_ch0_burst_fold.py"),
                  "--output", str(output), "--cleaned-parts", str(pulse_dir),
                  "--trigger-ticks", str(args.verified_input_root /
                    "Trigger_Timestamp_Envelope_20261015/Timestamps/trigger_ticks_5ns.bin"),
                  "--expected-windows", str(windows), "--period-method", method]
        if method == "ch2":
            verify += ["--pulse-source-folder", str(ch0 / "ROOT_Trees")]
        tracked(output, "verification", verify)
    comparison = stage / "Comparison"
    verify_pair(ch0, ch2, comparison)
    if args.stop_after == 4:
        return
    for method, output in (("ch0", ch0), ("ch2", ch2)):
        view = stage / ".downstream_inputs" / method
        run_view = view / KEY
        run_view.mkdir(parents=True, exist_ok=True)
        links = {path.name: path for path in output.iterdir() if path.name in
                 ("Folded_RF", "Folded_RF_3x", f"{KEY}_cleaned_rf_period_summary.txt")}
        for upstream in shared[1:]:
            for path in product(upstream).iterdir():
                if path.name in ("TOT_distrib", f"{KEY}_cleaned_rates_summary.txt"):
                    links[path.name] = path
        for name, path in links.items():
            if not (run_view / name).exists():
                (run_view / name).symlink_to(path)
        for index, script in (
            (5, "step05_estimate_gamma_neutron_background.py"),
            (6, "step06_quantify_gamma_neutron_components.py"),
            (7, "step07_phase_rate_slice_diagnostics.py"),
            (8, "step08_fit_gamma_tot_peak.py"),
            (9, "step09_gamma_window_sensitivity.py"),
            (10, "step10_gamma_reference_cfd.py"),
        ):
            if index == 8 and method == "ch2":
                continue
            target = stage / STEPS[index]
            if index != 8:
                target /= f"{method.capitalize()}-based"
            command = [sys.executable, str(source / script), "--results-root", str(view),
                       "--output-dir", str(target), "--channel", "2"]
            if index != 10:
                command += ["--run-glob", KEY]
            if index in (5, 6):
                command += ["--display-tot-max-ns", "60"]
            if index == 6:
                command += ["--phase-display", "raw", "--z-scale", "log"]
            if index in (7, 9, 10):
                command += ["--gamma-display-tot-max-ns", "60"]
            if index == 9:
                command += ["--gamma-fit-summary", str(stage / STEPS[8] / "gamma_tot_peak_fit_summary.tsv"),
                            "--assumed-trigger-window-width-ns", "10000"]
            if index == 10:
                command += ["--input-dir", str(stage / STEPS[9] / f"{method.capitalize()}-based")]
            tracked(target, "analysis", command)
    write_json(comparison / "completion.json", {"status": "PASS_STEPS01_TO10", "phase": args.phase,
        "shared_upstream": True, "step00_reused": True, "tot_display_ns": [0, 60],
        "detector_alignment": False, "periods_ns": periods})


def expected_windows(args: argparse.Namespace) -> int:
    text = (args.verified_input_root / "step00_decode" / args.phase / "ascii_validation.txt").read_text()
    import re
    match = re.search(r"Board 1 TDC_001049: blocks=(\d+)", text)
    require("ASCII_VALIDATION PASSED" in text and match is not None, "Verified selected-board decode required")
    return int(match[1])


if __name__ == "__main__":
    main()
