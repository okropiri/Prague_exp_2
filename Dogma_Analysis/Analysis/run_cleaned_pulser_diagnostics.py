from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import uproot

from run_zebra_ch0_clock import require
from step03_write_cleaned_all_channel_time_tot_outputs import build_heatmap_cmap, build_norm


FACTORS = (0, 1, 2, 4)
CHANNELS = range(9)
DELAY_EDGES = np.linspace(-300, 1100, 28001)
TOT_EDGES = np.linspace(0, 128, 2561)
BRANCHES = ["trigger_window_id", "channel", "rise_time_ch0_ref_ns", "tot_ns",
            "is_ch0_reference_hit"]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cadence_labels(regions: list[dict[str, str]], unresolved: np.ndarray,
                   windows: int) -> np.ndarray:
    require(windows > 0 and unresolved.ndim == 1
            and len(np.unique(unresolved)) == len(unresolved), "Invalid unresolved IDs")
    require(bool(np.all((unresolved >= 1) & (unresolved <= windows))),
            "Unresolved ID outside selected-board windows")
    labels = np.zeros(windows + 1, dtype=np.uint8)
    previous = 0
    for row in regions:
        first, last = int(row["first_reference_window"]), int(row["last_reference_window"])
        factor = int(row["cadence_factor"])
        require(1 <= first < last <= windows and first >= previous and factor in FACTORS[1:],
                "Invalid or overlapping cadence regions")
        require(first != previous or first in unresolved, "Shared cadence boundary must be unresolved")
        previous = last
        if int(row["calibrated_reference_windows"]) > 0:
            labels[first:last + 1] = factor
    labels[unresolved] = 0
    return labels


def histogram_quantile(counts: np.ndarray, edges: np.ndarray, q: float) -> float:
    require(counts.ndim == 1 and len(edges) == len(counts) + 1
            and bool(np.all(counts >= 0)) and counts.sum() > 0 and 0 < q < 1,
            "Cannot estimate quantile from invalid/empty histogram")
    cumulative = np.cumsum(counts)
    target = q * cumulative[-1]
    index = int(np.searchsorted(cumulative, target))
    below = 0 if index == 0 else cumulative[index - 1]
    return float(edges[index] + (target - below) / counts[index]
                 * (edges[index + 1] - edges[index]))


def zoom_edges(counts: np.ndarray, edges: np.ndarray) -> np.ndarray:
    low, high = (histogram_quantile(counts, edges, q) for q in (.001, .999))
    padding = max(2.0, .15 * (high - low))
    low, high = max(float(edges[0]), low - padding), min(float(edges[-1]), high + padding)
    bins = min(384, max(1, int(np.ceil((high - low) / .05))))
    return np.linspace(low, high, bins + 1)


def eligibility(channel: np.ndarray, relative: np.ndarray,
                reference: np.ndarray) -> np.ndarray:
    return np.isfinite(relative) & ((channel != 0) | reference)


def crop_histogram(counts: np.ndarray, x_edges: np.ndarray, y_edges: np.ndarray,
                   x_limits: np.ndarray, y_limits: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    indices = []
    for edges, limits in ((x_edges, x_limits), (y_edges, y_limits)):
        first = max(0, int(np.searchsorted(edges, limits[0], side="right")) - 1)
        last = min(len(edges) - 1, int(np.searchsorted(edges, limits[-1], side="left")))
        require(first < last, "Empty histogram crop")
        indices.append((first, last))
    (x_first, x_last), (y_first, y_last) = indices
    return (counts[x_first:x_last, y_first:y_last],
            x_edges[x_first:x_last + 1], y_edges[y_first:y_last + 1])


def fingerprints(paths: list[Path]) -> dict[str, dict[str, int | str]]:
    result = {}
    for path in paths:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        stat = path.stat()
        result[str(path)] = {"bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                             "sha256": digest.hexdigest()}
    return result


def save_figure(figure: plt.Figure, path: Path) -> None:
    for suffix in ("png", "pdf"):
        figure.savefig(path.with_suffix(f".{suffix}"), dpi=180)
    plt.close(figure)
    with Image.open(path.with_suffix(".png")) as image:
        image.verify()
    require(path.with_suffix(".pdf").read_bytes().startswith(b"%PDF-"), "Invalid PDF")


def heatmap(axis: plt.Axes, counts: np.ndarray, x_edges: np.ndarray,
            y_edges: np.ndarray, title: str) -> None:
    axis.set(title=title, ylabel="ToT (ns)")
    if counts.sum() == 0:
        axis.text(.5, .5, "No pulses", ha="center", va="center", transform=axis.transAxes)
        return
    image = axis.pcolormesh(x_edges, y_edges, np.ma.masked_less_equal(counts.T, 0),
                           cmap=build_heatmap_cmap(),
                           norm=build_norm(float(counts.max()), logarithmic=True),
                           shading="flat", rasterized=True)
    axis.figure.colorbar(image, ax=axis, label="Counts")


def write_core_views(samples: dict[int, tuple[np.ndarray, np.ndarray]],
                     plots: Path, reports: Path) -> list[dict[str, int | float]]:
    rows: list[dict[str, int | float]] = []
    for ch, (delay, tot) in samples.items():
        selected = (delay >= 40) & (delay <= 80)
        require(bool(np.any(selected)), f"Ch{ch}: no pulses in the explicit 40-80 ns core-display gate")
        core_delay, core_tot = delay[selected], tot[selected]
        x_edges = zoom_edges(np.histogram(core_delay, DELAY_EDGES)[0], DELAY_EDGES)
        y_edges = zoom_edges(np.histogram(core_tot, TOT_EDGES)[0], TOT_EDGES)
        counts = np.histogram2d(core_delay, core_tot, (x_edges, y_edges))[0]
        delay_q = np.quantile(core_delay, [.16, .5, .84])
        tot_q = np.quantile(core_tot, [.16, .5, .84])
        row: dict[str, int | float] = {
            "channel": ch, "all_referenced_pulses": len(delay), "core_display_pulses": int(selected.sum()),
            "outside_core_display_gate": int((~selected).sum()),
            "core_zoom_outside": int(selected.sum() - counts.sum()),
            "median_delay_ns": float(delay_q[1]),
            "delay_central68_halfwidth_ns": float((delay_q[2] - delay_q[0]) / 2),
            "median_tot_ns": float(tot_q[1]),
            "tot_central68_halfwidth_ns": float((tot_q[2] - tot_q[0]) / 2),
        }
        rows.append(row)
        fig, panels = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
        heatmap(panels[0], counts, x_edges, y_edges,
                f"Core display: {int(selected.sum()):,} pulses\n"
                f"Outside this zoom: {row['core_zoom_outside']:,}")
        panels[0].set_xlabel("Measured delay from selected ch0 (ns)")
        histogram = np.histogram(core_delay, DELAY_EDGES)[0]
        panels[1].stairs(histogram, DELAY_EDGES)
        panels[1].set(xlim=(x_edges[0], x_edges[-1]), xlabel="Measured delay from selected ch0 (ns)",
                      ylabel="Pulses / 0.05 ns bin",
                      title=f"Median = {delay_q[1]:.3f} ns\n"
                      f"Central-68% halfwidth = {(delay_q[2] - delay_q[0]) / 2:.3f} ns")
        fig.suptitle(f"Prague pulser | Ch{ch:02d} | explicit display gate: 40 <= delay <= 80 ns\n"
                     f"No ToT gate; {int((~selected).sum()):,} other pulses remain in full-population views")
        save_figure(fig, plots / f"ch{ch:02d}_main_pulse_zoom")
    with (reports / "main_pulse_display_summary.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Measured pulser delays, not an RF reconstruction")
    parser.add_argument("--cleaned-root-dir", type=Path, required=True)
    parser.add_argument("--clock-result-dir", type=Path, required=True)
    parser.add_argument("--verification", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    folder, output = args.clock_result_dir, args.output
    require(not output.exists(), f"Refusing to overwrite {output}")
    roots = sorted(args.cleaned_root_dir.glob("*_pulses_part*.root"))
    require(bool(roots), "No cleaned ROOT parts")
    clock = folder / "Clock_Calibration"
    assessment = json.loads((folder / "assessment.json").read_text())
    verified = json.loads(args.verification.read_text())
    require(verified["status"] == "PASS", "Step04 independent verification must pass first")
    regions = read_rows(clock / "cadence_regions.tsv")
    regimes = read_rows(clock / "cadence_regimes.tsv")
    unresolved = np.loadtxt(clock / "unresolved_clock_windows.tsv", dtype=np.int64, ndmin=1)
    labels = cadence_labels(regions, unresolved, verified["selected_board_windows"])
    folded_roots = list((folder / "Folded_RF_3x").glob("*_rf_phase3_vs_tot.root"))
    require(len(folded_roots) == 1, "Expected one existing 3x phase-ToT ROOT")
    inputs = roots + folded_roots + [folder / "assessment.json", args.verification,
                                     clock / "cadence_regions.tsv", clock / "cadence_regimes.tsv",
                                     clock / "unresolved_clock_windows.tsv"]
    before = fingerprints(inputs)
    output.mkdir(parents=True)
    plots, reports = output / "Plots", output / "Reports"
    plots.mkdir()
    reports.mkdir()
    input_counts = np.zeros(32, dtype=np.int64)
    no_reference = np.zeros(32, dtype=np.int64)
    finite_counts = np.zeros((5, 32), dtype=np.int64)
    eligible_counts = np.zeros_like(finite_counts)
    references = np.zeros(5, dtype=np.int64)
    delay_hist = np.zeros((5, 9, len(DELAY_EDGES) - 1), dtype=np.int64)
    tot_hist = np.zeros((5, 9, len(TOT_EDGES) - 1), dtype=np.int64)
    core_batches: dict[int, list[tuple[np.ndarray, np.ndarray]]] = {ch: [] for ch in range(2, 9)}
    last_reference = 0

    def batches():
        for root in roots:
            yield from uproot.iterate({str(root): "CleanedPulses"}, expressions=BRANCHES,
                                     step_size=250000, library="np")

    for data in batches():
        window, channel = data["trigger_window_id"], data["channel"]
        delay, tot = data["rise_time_ch0_ref_ns"], data["tot_ns"]
        reference = data["is_ch0_reference_hit"]
        require(bool(np.all((window >= 1) & (window < len(labels))))
                and bool(np.all((channel >= 0) & (channel < 32)))
                and bool(np.all(np.isfinite(tot) & (tot >= 0))), "Invalid cleaned pulse")
        require(bool(np.all((channel[reference] == 0) & (delay[reference] == 0))),
                "Invalid selected ch0 reference")
        ref_ids = window[reference]
        if len(ref_ids):
            require(int(ref_ids[0]) > last_reference and bool(np.all(ref_ids[1:] > ref_ids[:-1])),
                    "Duplicate/out-of-order selected references")
            last_reference = int(ref_ids[-1])
        group = labels[window]
        finite = np.isfinite(delay)
        eligible = eligibility(channel, delay, reference)
        input_counts += np.bincount(channel, minlength=32)
        no_reference += np.bincount(channel[~finite], minlength=32)
        references += np.bincount(group[reference], minlength=5)
        for factor in FACTORS:
            selected = (group == factor) & finite
            finite_counts[factor] += np.bincount(channel[selected], minlength=32)
            eligible_counts[factor] += np.bincount(channel[eligible & (group == factor)], minlength=32)
            for ch in CHANNELS:
                mask = selected & (channel == ch)
                delay_hist[factor, ch] += np.histogram(delay[mask], DELAY_EDGES)[0]
                tot_hist[factor, ch] += np.histogram(tot[mask], TOT_EDGES)[0]
                if ch in core_batches and np.any(mask):
                    core_batches[ch].append((delay[mask], tot[mask]))
    supported = np.asarray(FACTORS[1:])
    require(int(input_counts.sum()) == verified["cleaned_pulses"], "Input pulse count mismatch")
    require(np.array_equal(no_reference, assessment["unreferenced_channel_counts"]),
            "No-reference counts differ from Step04")
    require(np.array_equal(eligible_counts[supported].sum(axis=0), verified["folded_channel_counts"]),
            "Cadence eligibility differs from verified folding")
    require(np.array_equal(finite_counts[0], assessment["unresolved_clock_channel_counts"])
            and np.array_equal(eligible_counts[0], verified["eligible_clock_uncertainty_channel_counts"]),
            "Unresolved-clock accounting differs from Step04")
    require(np.array_equal(no_reference + finite_counts.sum(axis=0), input_counts),
            "Diagnostic categories do not conserve every cleaned pulse")
    for row in regimes:
        require(int(references[int(row["cadence_factor"])]) == int(row["calibrated_reference_windows"]),
                "Cadence reference membership mismatch")
    require(references[0] == verified["unresolved_clock_reference_windows"],
            "Unresolved reference membership mismatch")
    require(np.array_equal(delay_hist.sum(axis=2), finite_counts[:, :9])
            and np.array_equal(tot_hist.sum(axis=2), finite_counts[:, :9]),
            "Detailed channels have pulses outside diagnostic axes")
    print("PASS: all cleaned pulses and cadence/exclusion memberships reconcile", flush=True)

    axes = [(zoom_edges(delay_hist[:, ch].sum(axis=0), DELAY_EDGES),
             zoom_edges(tot_hist[:, ch].sum(axis=0), TOT_EDGES)) for ch in CHANNELS]
    maps = [np.zeros((5, len(x) - 1, len(y) - 1), dtype=np.int64) for x, y in axes]
    outside = np.zeros((5, 9), dtype=np.int64)
    for data in batches():
        channel, delay, tot = data["channel"], data["rise_time_ch0_ref_ns"], data["tot_ns"]
        group = labels[data["trigger_window_id"]]
        for ch, (x_edges, y_edges) in enumerate(axes):
            for factor in FACTORS:
                mask = (channel == ch) & np.isfinite(delay) & (group == factor)
                counts = np.histogram2d(delay[mask], tot[mask], (x_edges, y_edges))[0].astype(np.int64)
                maps[ch][factor] += counts
                outside[factor, ch] += int(mask.sum()) - int(counts.sum())
    for ch in CHANNELS:
        require(np.array_equal(maps[ch].sum(axis=(1, 2)) + outside[:, ch], finite_counts[:, ch]),
                f"Ch{ch}: zoomed maps lost pulses")
    period_us = {int(row["cadence_factor"]): float(row["period_ns"]) * assessment["prescaler"] / 1000
                 for row in regimes}
    names = {factor: f"Ch0 arrival cadence ~{period_us[factor]:.3f} us" for factor in FACTORS[1:]}
    summary: list[dict[str, int | float]] = []
    for ch in CHANNELS:
        for factor in FACTORS:
            count = int(finite_counts[factor, ch])
            if not count:
                continue
            row: dict[str, int | float] = {"channel": ch, "cadence_factor": factor,
                                         "referenced_pulses": count,
                                         "fold_eligible_pulses": int(eligible_counts[factor, ch]),
                                         "zoom_inside": int(maps[ch][factor].sum()),
                                         "zoom_outside": int(outside[factor, ch])}
            for name, hist, edges in (("delay", delay_hist[factor, ch], DELAY_EDGES),
                                      ("tot", tot_hist[factor, ch], TOT_EDGES)):
                row[f"{name}_mode_ns"] = float((edges[np.argmax(hist)] + edges[np.argmax(hist) + 1]) / 2)
                for q in (16, 50, 84):
                    row[f"{name}_q{q}_ns"] = histogram_quantile(hist, edges, q / 100)
                row[f"{name}_central68_halfwidth_ns"] = (float(row[f"{name}_q84_ns"])
                                                        - float(row[f"{name}_q16_ns"])) / 2
            summary.append(row)
        x_edges, y_edges = axes[ch]
        fig, panels = plt.subplots(2, 3, figsize=(15, 7), constrained_layout=True)
        for column, factor in enumerate(FACTORS[1:]):
            n = int(finite_counts[factor, ch])
            heatmap(panels[0, column], maps[ch][factor], x_edges, y_edges,
                    f"{names[factor]}\nN={n:,}; outside zoom={outside[factor, ch]:,}")
            panels[0, column].set_xlabel("Measured delay from selected ch0 (ns)")
            ax = panels[1, column]
            if n:
                ax.stairs(delay_hist[factor, ch] / n, DELAY_EDGES)
            else:
                ax.text(.5, .5, "No pulses in this cadence", transform=ax.transAxes, ha="center")
            ax.set(xlim=(x_edges[0], x_edges[-1]), xlabel="Measured delay from selected ch0 (ns)",
                   ylabel="Fraction of pulses / 0.05 ns bin")
        fig.suptitle(f"Prague pulser | Ch{ch:02d} | measured delays, no modulo or period rescaling\n"
                     f"Clock-unresolved pulses: {finite_counts[0, ch]:,}; "
                     + ("both cleaned ch0 families retained" if ch == 0 else "all referenced cleaned pulses"))
        save_figure(fig, plots / f"ch{ch:02d}_delay_tot_by_cadence")

    with uproot.open(folded_roots[0]) as root:
        for ch, (x_edges, y_edges) in enumerate(axes):
            phase, phase_edges, phase_tot_edges = root[f"RFPhase3xToT/hCh{ch:02d}"].to_numpy()
            require(int(phase.sum()) == verified["folded_channel_counts"][ch],
                    f"Ch{ch}: legacy folded ROOT count mismatch")
            phase_zoom = zoom_edges(phase.sum(axis=1), phase_edges)
            cropped, cropped_x, cropped_y = crop_histogram(
                phase, phase_edges, phase_tot_edges, phase_zoom, y_edges)
            phase_outside = int(phase.sum() - cropped.sum())
            fig, panels = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
            heatmap(panels[0], maps[ch].sum(axis=0), x_edges, y_edges,
                    f"Measured delay | N={finite_counts[:, ch].sum():,}\n"
                    f"Outside zoom={outside[:, ch].sum():,}")
            panels[0].set_xlabel("Measured delay from selected ch0 (ns)")
            heatmap(panels[1], cropped, cropped_x, cropped_y,
                    f"Existing conditional 3P fold | outside zoom={phase_outside:,}\n"
                    "Mixed normalization; NOT certified physical RF phase")
            panels[1].set_xlabel("Conditional cycle coordinate (display ns)")
            fig.suptitle(f"Prague pulser | Ch{ch:02d} | zoomed comparison; full-range originals preserved")
            save_figure(fig, plots / f"ch{ch:02d}_measured_vs_conditional_fold")

    fig, panels = plt.subplots(3, 3, figsize=(14, 10), constrained_layout=True)
    for ch, ax in enumerate(panels.flat):
        heatmap(ax, maps[ch].sum(axis=0), *axes[ch], f"Ch{ch:02d} | N={finite_counts[:, ch].sum():,}")
        ax.set_xlabel("Measured delay from ch0 (ns)")
    fig.suptitle("Prague pulser | measured delay-ToT zooms | no RF assumption\n"
                 "Each channel has its own display limits; excluded display tails are listed in the summary")
    save_figure(fig, plots / "measured_delay_tot_overview")
    fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
    for factor in FACTORS[1:]:
        ax.plot(list(CHANNELS), finite_counts[factor, :9], "o-", label=names[factor])
    ax.set_yscale("symlog", linthresh=1)
    ax.set(xticks=list(CHANNELS), xlabel="Channel",
           ylabel="Referenced cleaned pulses (both ch0 families)",
           title="Pulse membership by observed ch0 arrival cadence\n"
                 "Arrival cadence does not establish oscillator frequency or hardware divider")
    ax.legend()
    save_figure(fig, plots / "channel_counts_by_cadence")
    core_rows = write_core_views(
        {ch: (np.concatenate([part[0] for part in parts]), np.concatenate([part[1] for part in parts]))
         for ch, parts in core_batches.items()}, plots, reports)

    with (reports / "pulser_delay_tot_summary.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    result = {"status": "PASS_MEASURED_DIAGNOSTICS_NOT_RF_CERTIFICATION",
              "cleaned_pulses": int(input_counts.sum()),
              "cadence_factors": list(FACTORS), "referenced_channel_counts_by_factor": finite_counts.tolist(),
              "fold_eligible_channel_counts_by_factor": eligible_counts.tolist(),
              "references_by_factor": references.tolist(), "without_reference_channel_counts": no_reference.tolist(),
              "zoom_outside_by_factor": outside.tolist(), "statistics_bin_width_ns": .05,
              "legacy_fold_zoom_policy": "Crop outwards to whole original ROOT bins; list outside counts",
              "main_pulse_display_gate_ns": [40, 80], "main_pulse_display_summary": core_rows,
              "png_files": len(list(plots.glob("*.png"))), "pdf_files": len(list(plots.glob("*.pdf")))}
    require(result["png_files"] == result["pdf_files"] == 27, "Incomplete diagnostic image inventory")
    for path, state in before.items():
        stat = Path(path).stat()
        require(stat.st_size == state["bytes"] and stat.st_mtime_ns == state["mtime_ns"],
                f"Input changed during diagnostics: {path}")
    (reports / "input_fingerprints.json").write_text(json.dumps(before, indent=2) + "\n")
    (reports / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    text = (
        "PULSER DIAGNOSTICS - MEASURED DELAYS, NOT PHYSICAL RF CERTIFICATION\n"
        "STATUS: PASS for pulse accounting, cadence membership, zoom accounting and image integrity.\n\n"
        "These plots correct the broad-axis presentation; they do not repair or certify a physical RF model.\n"
        "Delay = cleaned rising time relative to selected first-valid ch0; no modulo/rescaling.\n"
        "All cleaned pulses are counted. Missing-reference and clock-unresolved categories remain explicit.\n"
        "Both ch0 families are included in measured diagnostics; only the selected reference enters ch0 folding.\n"
        "Zoom limits use histogram q0.001..q0.999 plus >=2 ns padding. Outside-zoom counts are listed,\n"
        "not deleted. Quantiles use within-bin interpolation in 0.05 ns histograms; half central-68%\n"
        "width is descriptive, not a Gaussian fit, electronics jitter or detector resolution.\n"
        "Ch0 contains two populations; its combined width must not be interpreted as reference jitter.\n"
        "Selected-ch0 relative delay is zero by construction, not a measurement of reference stability.\n\n"
        "Additional main-pulse views explicitly select 40 <= measured delay <= 80 ns for display only,\n"
        "with no ToT gate. Their exact quantiles describe that cohort, not all pulses or detector resolution.\n"
        "Late/low-ToT populations remain in the full-population views and input; all gate exclusions\n"
        "and zoom exclusions are counted. The core-display gate is not applied to cleaning or folding.\n\n"
        "The comparison's right panel only zooms the existing mixed-normalization 3P histogram.\n"
        "Its coordinate is conditional on divider 4096 and decoder timing scale; not certified RF phase.\n"
        "3P means modulo three cycles, not three pasted peaks. Timing cadence alone cannot distinguish\n"
        "an oscillator change from a divider/trigger/acquisition change. No later physics step is authorized.\n\n"
    )
    for row in regimes:
        factor = int(row["cadence_factor"])
        text += (f"Cadence factor {factor}: ~{period_us[factor]:.9f} us between nominal ch0 arrivals; "
                 f"{references[factor]:,} supported references; "
                 f"ch2-8 pulses={int(finite_counts[factor, 2:9].sum()):,}.\n")
    text += f"\nAll {int(input_counts.sum()):,} cleaned pulses reconcile. 27 PNG/PDF pairs verified.\n"
    for row in summary:
        if row["channel"] in range(2, 9) and row["cadence_factor"] == 1:
            text += (f"Ch{row['channel']:02d}: N={row['referenced_pulses']:,}; "
                     f"median delay={row['delay_q50_ns']:.4f} ns; "
                     f"central68 delay halfwidth={row['delay_central68_halfwidth_ns']:.4f} ns; "
                     f"median ToT={row['tot_q50_ns']:.4f} ns.\n")
    text += "\nMAIN-PULSE DISPLAY COHORT (40-80 ns, no ToT gate; not a physics cut):\n"
    for row in core_rows:
        text += (f"Ch{row['channel']:02d}: N={row['core_display_pulses']:,}; "
                 f"median delay={row['median_delay_ns']:.4f} ns; "
                 f"central68 halfwidth={row['delay_central68_halfwidth_ns']:.4f} ns; "
                 f"other referenced pulses={row['outside_core_display_gate']:,}.\n")
    (reports / "pulser_diagnostics_report.txt").write_text(text)
    print(text, flush=True)


if __name__ == "__main__":
    main()
