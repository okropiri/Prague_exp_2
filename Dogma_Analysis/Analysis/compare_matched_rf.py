from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import uproot

from analysis_provenance import file_identity, record_arguments, write_json
from run_zebra_ch0_clock import require


def main() -> None:
    parser = argparse.ArgumentParser(description="Equal-bin, equal-population comparison of unshifted RF methods.")
    parser.add_argument("--ch0", type=Path, required=True)
    parser.add_argument("--ch2", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = record_arguments(parser.parse_args())
    args.output.mkdir(parents=True, exist_ok=True)
    models = {name: json.loads((folder / "assessment.json").read_text())
              for name, folder in (("ch0", args.ch0), ("ch2", args.ch2))}
    for field in ("folded_channel_counts", "unreferenced_channel_counts", "decoded_selected_board_windows"):
        require(models["ch0"][field] == models["ch2"][field], f"Unmatched {field}")
    periods = {name: float(model["mean_period_ns"]) for name, model in models.items()}
    edges = {factor: np.arange(int(np.ceil(factor * max(periods.values()) / .25)) + 1) * .25
             for factor in (1, 3)}
    tot_edges = np.linspace(0, 128, 1281)
    profiles = {(name, factor): np.zeros(len(edges[factor]) - 1, dtype=np.uint64)
                for name in models for factor in (1, 3)}
    maps = {key: np.zeros((1280, len(value)), dtype=np.uint64) for key, value in profiles.items()}
    parts = sorted((args.ch0 / "ROOT_Trees").glob("NCAL_ch02_pulses_part*.root"))
    require(bool(parts), "Complete shared NCAL pulse TTrees required")
    input_count = invalid_count = above_display = 0
    for part in parts:
        with uproot.open(part) as root:
            for data in root["Pulses"].iterate(step_size=250000, library="np"):
                valid = data["phase_valid"]
                input_count += len(valid)
                invalid_count += int((~valid).sum())
                relative, tot = data["rise_time_ch0_ref_ns"][valid], data["tot_ns"][valid]
                above_display += int((tot >= 60).sum())
                y = np.floor(tot * 10).astype(np.int64)
                inside = (y >= 0) & (y < 1280)
                for name, period in periods.items():
                    for factor in (1, 3):
                        phase = np.mod(relative / period, factor) * period
                        x = np.floor(phase / .25).astype(np.int64)
                        bins = len(edges[factor]) - 1
                        profiles[name, factor] += np.bincount(x, minlength=bins).astype(np.uint64)
                        flat = y[inside] * bins + x[inside]
                        maps[name, factor] += np.bincount(flat, minlength=1280 * bins).reshape(
                            1280, bins).astype(np.uint64)
    require(input_count == models["ch0"]["folded_channel_counts"][2]
            + models["ch0"]["unreferenced_channel_counts"][2], "Shared NCAL set is incomplete")
    require(invalid_count == models["ch0"]["unreferenced_channel_counts"][2], "Invalid-phase population differs")
    with uproot.recreate(args.output / "Equal_bin_NCAL_replay.root", compression=uproot.ZLIB(4)) as root:
        for factor in (1, 3):
            xedges = edges[factor]
            width = .25
            figure, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True, sharex=True, sharey=True)
            maximum = max(int(maps[name, factor][:600].max()) for name in models)
            norm = LogNorm(vmin=1, vmax=max(2, maximum))
            for axis, name in zip(axes, models, strict=True):
                counts = maps[name, factor]
                require(int(profiles[name, factor].sum()) == input_count - invalid_count,
                        "Equal-bin comparison lost phase-valid pulses")
                image = axis.pcolormesh(xedges, tot_edges, np.ma.masked_equal(counts, 0),
                                        cmap="magma", norm=norm, rasterized=True)
                axis.set(xlabel=f"Unshifted ch0-referenced phase, {factor}P (ns)",
                         title=f"{name}-derived T = {periods[name]:.9f} ns", ylim=(0, 60))
                x, y = np.nonzero(counts.T)
                payload = {"phase_ns": (x + .5) * width, "tot_ns": (y + .5) * .1,
                           "x_bin": x.astype(np.uint16), "y_bin": y.astype(np.uint16),
                           "count": counts[y, x]}
                key = f"{name}_phase{factor}_tot"
                root.mktree(key, {field: value.dtype for field, value in payload.items()}).extend(payload)
                root[f"{key}_axes"] = json.dumps({"x_edges_ns": xedges.tolist(), "y_min_ns": 0,
                    "y_max_ns": 128, "y_bins": 1280, "display_y_max_ns": 60})
            axes[0].set_ylabel("ToT (ns)")
            figure.colorbar(image, ax=axes, label="Counts / common 0.25 ns x 0.1 ns bin")
            figure.suptitle(f"Same {input_count - invalid_count:,} NCAL hits; no detector recentering")
            for extension in ("png", "pdf"):
                figure.savefig(args.output / f"NCAL_{factor}P_ToT_comparison.{extension}", dpi=180)
            plt.close(figure)
            figure, axis = plt.subplots(figsize=(12, 4.5), constrained_layout=True)
            for name in models:
                axis.stairs(profiles[name, factor], xedges, label=f"{name}-derived period")
            axis.set(xlabel=f"Unshifted ch0-referenced phase, {factor}P (ns)",
                     ylabel="NCAL counts / common 0.25 ns bin")
            axis.legend()
            for extension in ("png", "pdf"):
                figure.savefig(args.output / f"NCAL_{factor}P_equal_bin_profiles.{extension}", dpi=180)
            plt.close(figure)
    with uproot.open(args.output / "Equal_bin_NCAL_replay.root") as root:
        for name in models:
            for factor in (1, 3):
                data = root[f"{name}_phase{factor}_tot"].arrays(library="np")
                replay = np.zeros_like(maps[name, factor])
                replay[data["y_bin"], data["x_bin"]] = data["count"]
                require(np.array_equal(replay, maps[name, factor]), "Persisted comparison replay differs")
    scan = next((args.ch2 / "Period_Search").glob("*/period_scan.tsv"))
    raw = np.loadtxt(scan, dtype=str, ndmin=2)
    require(np.all(np.isin(raw[:, 13], ("true", "false"))), "Unknown period-scan validity encoding")
    rows, valid = raw[:, :13].astype(float), raw[:, 13] == "true"
    figure, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
    for round_index in np.unique(rows[:, 0]):
        selected = (rows[:, 0] == round_index) & valid
        axes[0].plot(rows[selected, 1], rows[selected, 12], ".-", label=f"Round {int(round_index)}")
    axes[0].set(xlabel="Trial RF period (ns)", ylabel="Detector-feature merit")
    refined = rows[(rows[:, 0] == rows[:, 0].max()) & valid]
    axes[1].plot(refined[:, 1], refined[:, 7], ".-")
    axes[1].set(xlabel="Refined trial RF period (ns)", ylabel="Selected-peak sigma (ns)")
    for axis in axes:
        axis.axvline(periods["ch2"], color="black", linestyle="--", label="Selected ch2 period")
        axis.legend()
    for extension in ("png", "pdf"):
        figure.savefig(args.output / f"Ch2_period_search_procedure.{extension}", dpi=180)
    plt.close(figure)
    write_json(args.output / "equal_bin_verification.json", {
        "status": "PASS_EQUAL_BIN_COMPLETE_NCAL_COMPARISON", "source": file_identity(Path(__file__)),
        "input_ncal_pulses": input_count, "phase_invalid_pulses": invalid_count,
        "phase_valid_pulses": input_count - invalid_count, "phase_valid_at_or_above_display_tot": above_display,
        "periods_ns": periods, "phase_bin_width_ns": .25, "tot_bin_width_ns": .1,
        "tot_stored_ns": [0, 128], "tot_display_ns": [0, 60],
        "phase_zero": "same selected guarded ch0; no fitted detector offset",
        "chi_square_claim": "None: the two histograms use the same hits and are statistically correlated.",
        "range_note": "Common axes extend to the larger fitted period (or 3P); the shorter-period tail is structurally empty.",
        "persisted_replay_bins_equal": True,
    })


if __name__ == "__main__":
    main()
