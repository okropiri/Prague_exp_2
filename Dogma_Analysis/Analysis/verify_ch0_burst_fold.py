from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import socket
from analysis_provenance import record_arguments

import numpy as np
import uproot

from run_zebra_ch0_clock import require
from step04_write_cleaned_all_channel_rf_period_outputs import (
    build_dense, read_matrix, read_sparse_histogram,
)


class PulseTrees:
    def __init__(self, folder: Path):
        self.folder = folder
        self.file = None
        self.tree = None
        self.rows = 0
        self.parts: list[Path] = []
        self.total = 0

    def extend(self, data: dict[str, np.ndarray]) -> None:
        position = 0
        length = len(data["trigger_window_id"])
        while position < length:
            if self.file is None:
                path = self.folder / f"NCAL_ch02_pulses_part{len(self.parts) + 1:03d}.root"
                self.file = uproot.recreate(path, compression=uproot.ZLIB(4))
                self.tree = self.file.mktree("Pulses", {key: value.dtype for key, value in data.items()})
                self.parts.append(path)
                self.rows = 0
            count = min(length - position, 5_000_000 - self.rows)
            require(self.tree is not None, "Missing pulse TTree writer")
            self.tree.extend({key: value[position:position + count] for key, value in data.items()})
            position += count
            self.rows += count
            self.total += count
            if self.rows == 5_000_000:
                self.close()

    def close(self) -> None:
        if self.file is not None:
            self.file.close()
            self.file = None
            self.tree = None
            require(self.parts[-1].stat().st_size <= 250 * 1024 ** 2, "Pulse ROOT part exceeds 250 MiB")


def write_replay(tables: Path, output: Path) -> None:
    with uproot.recreate(output, compression=uproot.ZLIB(4)) as root:
        for factor, filename in ((1, "folded_matrix.tsv"), (3, "folded3_matrix.tsv")):
            metadata, centers, counts, _ = read_matrix(tables / filename)
            channel, index = np.nonzero(counts.T)
            name = f"phase{factor}_profiles"
            payload = {"channel": channel.astype(np.uint8), "phase_ns": centers[index],
                       "x_bin": index.astype(np.uint16),
                       "count": counts[index, channel].astype(np.uint64)}
            tree = root.mktree(name, {key: value.dtype for key, value in payload.items()})
            tree.extend(payload)
            root[f"{name}_axes"] = json.dumps(metadata)
            for quantity in ("tot", "time", "trigger"):
                prefix = "phase" if factor == 1 else "phase3"
                meta, entries, _, _ = read_sparse_histogram(tables / f"{prefix}_{quantity}.tsv")
                rows = np.asarray([(ch, x, y, n) for ch, cells in entries.items()
                                   for x, y, n in cells], dtype=np.uint64).reshape(-1, 4)
                x = (rows[:, 1] + .5) * float(meta["x_max_ns"]) / int(meta["x_bins"])
                y = float(meta["y_min_ns"]) + (rows[:, 2] + .5) * (
                    float(meta["y_max_ns"]) - float(meta["y_min_ns"])) / int(meta["y_bins"])
                payload = {"channel": rows[:, 0].astype(np.uint8),
                           "phase_ns": x, "value_ns": y,
                           "x_bin": rows[:, 1].astype(np.uint16),
                           "y_bin": rows[:, 2].astype(np.uint16), "count": rows[:, 3]}
                name = f"phase{factor}_{quantity}"
                tree = root.mktree(name, {key: value.dtype for key, value in payload.items()})
                tree.extend(payload)
                root[f"{name}_axes"] = json.dumps(meta)
    with uproot.open(output) as root:
        for name, classname in root.classnames().items():
            if classname != "TTree":
                continue
            tree = root[name]
            data = tree.arrays(library="np")
            meta = json.loads(str(root[f"{name.split(';')[0]}_axes"]))
            require(tree.classname == "TTree", "Plot replay is not a ROOT TTree")
            if "y_bin" not in data:
                filename = "folded_matrix.tsv" if name.startswith("phase1") else "folded3_matrix.tsv"
                _, _, expected, _ = read_matrix(tables / filename)
                replay = np.zeros(expected.shape, dtype=np.uint64)
                replay[data["x_bin"], data["channel"]] = data["count"]
                require(np.array_equal(replay, expected), "Profile replay tree differs from plotted bins")
            else:
                prefix = "phase" if name.startswith("phase1") else "phase3"
                quantity = name.split(";")[0].rsplit("_", 1)[1]
                _, entries, _, _ = read_sparse_histogram(tables / f"{prefix}_{quantity}.tsv")
                for ch, cells in entries.items():
                    selected = data["channel"] == ch
                    replay = np.zeros((int(meta["y_bins"]), int(meta["x_bins"])), dtype=np.uint64)
                    replay[data["y_bin"][selected], data["x_bin"][selected]] = data["count"][selected]
                    require(np.array_equal(replay, build_dense(cells, replay.shape[1], replay.shape[0])),
                            "2D replay tree differs from plotted bins")


def main() -> None:
    parser = argparse.ArgumentParser(description="Independent full-population ch0 QC and reproducible ROOT TTrees.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cleaned-parts", type=Path, required=True)
    parser.add_argument("--trigger-ticks", type=Path, required=True)
    parser.add_argument("--expected-windows", type=int, required=True)
    parser.add_argument("--period-method", choices=("ch0", "ch2"), default="ch0")
    parser.add_argument("--pulse-source-folder", type=Path,
                        help="Reuse complete NCAL pulse trees; recompute phases from their unshifted times")
    args = record_arguments(parser.parse_args())
    require(socket.gethostname().split(".")[0] == "trbsrv-b", "Run DOGMA verification only on padiwa")
    folder = args.output / "ROOT_Trees"
    require(not folder.exists(), "Refusing previous verification/export")
    folder.mkdir()
    assessment = json.loads((args.output / "assessment.json").read_text())
    require(assessment["include_all_ch0_hits"]
            and (args.period_method == "ch2" or assessment.get("burst_clock")), "Wrong folding mode")
    require(args.expected_windows == assessment["decoded_selected_board_windows"],
            "Requested window inventory differs from the original decode")
    require(assessment["board"] == "TDC_001049" and assessment["coarse_quantum_ns"] == 5,
            "Native board1049 five-ns timestamps required")
    require(assessment["unresolved_reference_windows"] == 0, "Clock-support window exclusion remains")
    clock, tables = args.output / "Clock_Calibration", args.output / "Fold_Tables"
    raw = np.loadtxt(clock / "ch0_reference_timestamps.tsv", usecols=(0, 1, 2, 4), ndmin=2)
    ids, coarse, reference, times = raw[:, 0].astype(np.int64), raw[:, 1].astype(np.uint64), raw[:, 2], raw[:, 3]
    period_rows = np.loadtxt(clock / "window_periods.tsv", ndmin=2)
    require(np.array_equal(period_rows[:, 0], ids), "Some reference windows lack periods")
    require(np.all(period_rows[:, 1] == period_rows[0, 1]), "Expected explicitly shared period")
    w = args.expected_windows
    periods, reference_rise = np.zeros(w + 1), np.full(w + 1, np.nan)
    periods[ids], reference_rise[ids] = period_rows[:, 1], reference
    ticks = np.memmap(args.trigger_ticks, dtype="<u8", mode="r")
    require(len(ticks) >= w and np.array_equal(ticks[ids - 1], coarse), "Original trigger identity differs")
    delta = np.diff(times)
    joined = (np.diff(ids) == 1) & (abs(delta - np.median(delta)) <= assessment["coarse_quantum_ns"])
    starts = np.r_[0, np.flatnonzero(~joined) + 1]
    stops = np.r_[starts[1:], len(ids)]
    n = stops - starts
    region = np.repeat(np.arange(len(starts)), n)
    j = np.arange(1, len(ids)) - starts[region[:-1]]
    weights = j * (n[region[:-1]] - j) / 2
    select = joined & (n[region[:-1]] >= 32)
    independent_period = float(np.sum(delta[select] * weights[select])
                               / np.sum(n[n >= 32] * (n[n >= 32] ** 2 - 1) / 12) / assessment["prescaler"])
    if args.period_method == "ch0":
        require(abs(independent_period - assessment["mean_period_ns"]) < 1e-10,
                "Independent interval-weighted RF period differs")
    else:
        estimate = json.loads((args.output / "period_estimate.json").read_text())
        require(float(estimate["deduced_period_ns"]) == assessment["mean_period_ns"],
                "Ch2 fold period differs from independently recorded scan result")
    fit_support = np.zeros(w + 1, dtype=bool)
    if args.period_method == "ch0":
        fit_support[ids] = n[region] >= 32
    del raw, times, coarse, period_rows, delta, region, weights, select, joined, j
    profiles, expected_profiles, maps, expected_maps, axes = {}, {}, {}, {}, {}
    for factor, filename in ((1, "folded_matrix.tsv"), (3, "folded3_matrix.tsv")):
        _, _, expected, _ = read_matrix(tables / filename)
        expected_profiles[factor] = expected
        profiles[factor] = np.zeros_like(expected, dtype=np.uint64)
        for quantity in ("tot", "time", "trigger"):
            prefix = "phase" if factor == 1 else "phase3"
            meta, entries, _, _ = read_sparse_histogram(tables / f"{prefix}_{quantity}.tsv")
            require(max(entries, default=0) < 9, "Unexpected detailed channel")
            shape = (9, int(meta["y_bins"]), int(meta["x_bins"]))
            expected_maps[factor, quantity] = np.zeros(shape, dtype=np.uint64)
            for ch, cells in entries.items():
                expected_maps[factor, quantity][ch] = build_dense(cells, shape[2], shape[1])
            maps[factor, quantity] = np.zeros(shape, dtype=np.uint64)
            axes[factor, quantity] = meta
    input_counts, folded_counts, missing_counts = (np.zeros(32, dtype=np.int64) for _ in range(3))
    window_counts = np.zeros(w + 1, dtype=np.uint32)
    pulse_trees = PulseTrees(folder)
    reference_count = 0
    parts = sorted(args.cleaned_parts.glob("*_pulses*.root"))
    require(bool(parts), "No cleaned ROOT TTrees")
    part_state = {str(p): (p.stat().st_size, p.stat().st_mtime_ns) for p in parts}
    fields = ["trigger_window_id", "channel", "rise_time_trigger_ref_ns", "tot_ns",
              "rise_time_ch0_ref_ns", "is_ch0_reference_hit"]
    with ThreadPoolExecutor(max_workers=1) as executor:
        for part in parts:
            for data in uproot.iterate({str(part): "CleanedPulses"}, expressions=fields,
                                      step_size=250000, library="np", decompression_executor=executor):
                window, channel = data["trigger_window_id"], data["channel"]
                require(np.all((window >= 1) & (window <= w)), "Pulse outside original window range")
                rise, tot = data["rise_time_trigger_ref_ns"], data["tot_ns"]
                valid = np.isfinite(reference_rise[window])
                require(np.array_equal(valid, np.isfinite(data["rise_time_ch0_ref_ns"])),
                        "Native reference availability differs from cleaned ROOT")
                native_relative = rise - reference_rise[window]
                require(np.all(abs(native_relative[valid] - data["rise_time_ch0_ref_ns"][valid]) <= .0011),
                        "Ch0-relative time mismatch")
                chosen = data["is_ch0_reference_hit"]
                require(np.all(channel[chosen] == 0) and np.all(valid[chosen])
                        and np.all(rise[chosen] == reference_rise[window[chosen]]),
                        "First-passing reference identity differs")
                reference_count += int(chosen.sum())
                input_counts += np.bincount(channel, minlength=32)
                folded_counts += np.bincount(channel[valid], minlength=32)
                missing_counts += np.bincount(channel[~valid], minlength=32)
                unique, counts = np.unique(window, return_counts=True)
                window_counts[unique] += counts.astype(np.uint32)
                phase_arrays = {}
                for factor in (1, 3):
                    meta = axes[factor, "tot"]
                    phase = np.full(len(window), np.nan)
                    phase[valid] = np.mod(native_relative[valid] / periods[window[valid]], factor) * float(meta["deduced_period_ns"])
                    phase_arrays[factor] = phase
                    x_bins = int(meta["x_bins"])
                    x = np.floor(phase[valid] / float(meta["x_max_ns"]) * x_bins).astype(np.int64)
                    require(np.all((x >= 0) & (x < x_bins)), "Independent phase overflow")
                    flat = x * 32 + channel[valid]
                    profiles[factor] += np.bincount(flat, minlength=x_bins * 32).reshape(x_bins, 32).astype(np.uint64)
                    for quantity, values in (("tot", tot), ("time", native_relative), ("trigger", rise)):
                        axis = axes[factor, quantity]
                        keep = valid & (channel < 9)
                        y = np.floor((values[keep] - float(axis["y_min_ns"])) /
                                     (float(axis["y_max_ns"]) - float(axis["y_min_ns"])) * int(axis["y_bins"])).astype(np.int64)
                        xp = np.floor(phase[keep] / float(axis["x_max_ns"]) * int(axis["x_bins"])).astype(np.int64)
                        inside = (y >= 0) & (y < int(axis["y_bins"]))
                        require(quantity == "tot" or np.all(inside), "Timing display loses valid pulses")
                        flat = ((channel[keep][inside] * int(axis["y_bins"]) + y[inside]) *
                                int(axis["x_bins"]) + xp[inside])
                        values_flat, multiplicities = np.unique(flat, return_counts=True)
                        maps[factor, quantity].reshape(-1)[values_flat] += multiplicities.astype(np.uint64)
                selected = channel == 2
                if args.pulse_source_folder is None:
                    pulse_trees.extend({
                        "trigger_window_id": window[selected].astype(np.uint32),
                        "rise_time_trigger_ref_ns": rise[selected], "rise_time_ch0_ref_ns": native_relative[selected],
                        "tot_ns": tot[selected], "rf_phase_ns": phase_arrays[1][selected],
                        "rf_phase3_ns": phase_arrays[3][selected], "phase_valid": valid[selected],
                    })
    pulse_trees.close()
    require(reference_count == len(ids), "Reference TTree count differs")
    require(np.array_equal(input_counts, folded_counts + missing_counts), "Lost cleaned pulses")
    require(np.array_equal(folded_counts, assessment["folded_channel_counts"]), "Folded population differs")
    require(np.array_equal(missing_counts, assessment["unreferenced_channel_counts"]), "Missing-phase accounting differs")
    for factor in (1, 3):
        require(np.array_equal(profiles[factor], expected_profiles[factor]), f"{factor}P independent profile differs")
        for quantity in ("tot", "time", "trigger"):
            require(np.array_equal(maps[factor, quantity], expected_maps[factor, quantity]),
                    f"{factor}P {quantity}: independent plot bins differ")
    if args.pulse_source_folder is None:
        require(pulse_trees.total == int(input_counts[2]), "NCAL pulse TTree loses data")
    else:
        pulse_trees.parts = sorted(args.pulse_source_folder.glob("NCAL_ch02_pulses_part*.root"))
        require(bool(pulse_trees.parts), "Missing shared NCAL pulse TTrees")
    with uproot.recreate(folder / "Trigger_windows.root", compression=uproot.ZLIB(4)) as root:
        tree = root.mktree("TriggerWindows", {
            "trigger_window_id": "uint32", "coarse_tick_5ns": "uint64", "cleaned_pulses": "uint32",
            "ch0_reference_rise_ns": "float64", "phase_valid": "bool", "direct_period_fit_support": "bool",
        })
        for start in range(0, w, 250000):
            stop = min(w, start + 250000)
            tree.extend({"trigger_window_id": np.arange(start + 1, stop + 1, dtype=np.uint32),
                         "coarse_tick_5ns": np.asarray(ticks[start:stop]),
                         "cleaned_pulses": window_counts[start + 1:stop + 1],
                         "ch0_reference_rise_ns": reference_rise[start + 1:stop + 1],
                         "phase_valid": np.isfinite(reference_rise[start + 1:stop + 1]),
                         "direct_period_fit_support": fit_support[start + 1:stop + 1]})
        root["RFModel"] = json.dumps({"period_ns": assessment["mean_period_ns"], "phase_origin": assessment["phase_origin"],
                                      "absolute_seconds": "coarse_tick_5ns * 5e-9; DAQ counter, not UTC"})
    write_replay(tables, folder / "All_channel_plot_replay.root")
    with uproot.open(folder / "Trigger_windows.root") as root:
        require(root["TriggerWindows"].classname == "TTree" and root["TriggerWindows"].num_entries == w,
                "Trigger-window TTree lost windows")
        data = root["TriggerWindows"].arrays(["coarse_tick_5ns", "cleaned_pulses", "phase_valid"], library="np")
        require(np.array_equal(data["coarse_tick_5ns"], ticks[:w])
                and np.array_equal(data["cleaned_pulses"], window_counts[1:])
                and np.array_equal(data["phase_valid"], np.isfinite(reference_rise[1:])),
                "Persisted trigger-window TTree differs")
        del data
    persisted_profiles = {factor: np.zeros_like(profiles[factor][:, 2]) for factor in (1, 3)}
    persisted_maps = {key: np.zeros_like(value[2]) for key, value in maps.items()}
    persisted_pulses = persisted_invalid = 0
    for part in pulse_trees.parts:
        with uproot.open(part) as root:
            require(root["Pulses"].classname == "TTree", "NCAL pulse export is not a TTree")
            for data in root["Pulses"].iterate(step_size=250000, library="np"):
                valid = data["phase_valid"]
                require(np.array_equal(valid, np.isfinite(reference_rise[data["trigger_window_id"]]))
                        and np.allclose(data["rise_time_ch0_ref_ns"][valid],
                                        (data["rise_time_trigger_ref_ns"]
                                         - reference_rise[data["trigger_window_id"]])[valid],
                                        rtol=0, atol=1e-10),
                        "Persisted NCAL reference identity differs")
                persisted_pulses += len(valid)
                persisted_invalid += int((~valid).sum())
                require(np.all(np.isnan(data["rf_phase_ns"][~valid]))
                        and np.all(np.isnan(data["rf_phase3_ns"][~valid])), "Invalid phases were fabricated")
                for factor, field in ((1, "rf_phase_ns"), (3, "rf_phase3_ns")):
                    meta = axes[factor, "tot"]
                    phase = (data[field][valid] if args.pulse_source_folder is None else
                             np.mod(data["rise_time_ch0_ref_ns"][valid] /
                                    periods[data["trigger_window_id"][valid]], factor)
                             * float(meta["deduced_period_ns"]))
                    x = np.floor(phase / float(meta["x_max_ns"]) * int(meta["x_bins"])).astype(np.int64)
                    persisted_profiles[factor] += np.bincount(x, minlength=int(meta["x_bins"])).astype(np.uint64)
                    for quantity, value in (("tot", "tot_ns"), ("time", "rise_time_ch0_ref_ns"),
                                            ("trigger", "rise_time_trigger_ref_ns")):
                        axis = axes[factor, quantity]
                        values = data[value][valid]
                        y = np.floor((values - float(axis["y_min_ns"])) /
                                     (float(axis["y_max_ns"]) - float(axis["y_min_ns"])) *
                                     int(axis["y_bins"])).astype(np.int64)
                        inside = (y >= 0) & (y < int(axis["y_bins"]))
                        flat = y[inside] * int(axis["x_bins"]) + x[inside]
                        cells, counts = np.unique(flat, return_counts=True)
                        persisted_maps[factor, quantity].reshape(-1)[cells] += counts.astype(np.uint64)
    require(persisted_pulses == int(input_counts[2]) and persisted_invalid == int(missing_counts[2]),
            "Persisted NCAL TTree population differs")
    for factor in (1, 3):
        require(np.array_equal(persisted_profiles[factor], profiles[factor][:, 2]),
                "Persisted NCAL TTree cannot reproduce profile")
        for quantity in ("tot", "time", "trigger"):
            require(np.array_equal(persisted_maps[factor, quantity], maps[factor, quantity][2]),
                    "Persisted NCAL TTree cannot reproduce detailed plot")
    require(part_state == {str(p): (p.stat().st_size, p.stat().st_mtime_ns) for p in parts},
            "Cleaned ROOT inputs changed during verification")
    root_files = sorted(folder.glob("*.root"))
    require(all(p.stat().st_size <= 250 * 1024 ** 2 for p in root_files), "ROOT file exceeds local transfer limit")
    result = {
        "status": "PASS_INDEPENDENT_ALL_WINDOW_ALL_PULSE_ALL_PLOT_BIN_VERIFICATION",
        "trigger_windows": w, "reference_windows": len(ids), "clock_support_excluded_windows": 0,
        "empty_cleaned_windows": int(np.count_nonzero(window_counts[1:] == 0)),
        "input_pulses": int(input_counts.sum()), "input_channel_counts": input_counts.tolist(),
        "folded_channel_counts": folded_counts.tolist(), "unreferenced_channel_counts": missing_counts.tolist(),
        "period_method": args.period_method,
        "independent_ch0_period_ns": independent_period,
        "difference_from_ch0_period_ns": independent_period - assessment["mean_period_ns"],
        "all_1p_3p_profile_and_detailed_bins_equal": True,
        "ncal_pulse_tree_entries": persisted_pulses,
        "shared_pulse_source": str(args.pulse_source_folder) if args.pulse_source_folder else None,
        "persisted_ncal_pulse_tree_entries_verified": persisted_pulses,
        "root_files": [{"file": p.name, "bytes": p.stat().st_size} for p in root_files],
        "root_total_bytes": sum(p.stat().st_size for p in root_files),
        "entire_root_export_fits_local_1000mib_limit": sum(p.stat().st_size for p in root_files) <= 1000 * 1024 ** 2,
        "root_replay": "Weighted bin TTrees exactly reproduce all channel profiles and detailed maps. "
        "Shared NCAL trees retain every ch2 pulse, including invalid-phase rows; for the ch2-period method, "
        "recompute phase from rise_time_ch0_ref_ns and its recorded period, not the stored ch0-period phase branches.",
    }
    (folder / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
