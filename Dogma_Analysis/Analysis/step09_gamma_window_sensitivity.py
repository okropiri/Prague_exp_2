from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from step05_estimate_gamma_neutron_background import find_sparse_files, mask_between, parse_run_key, read_channel_histogram
from step07_phase_rate_slice_diagnostics import (
    SliceProfile,
    draw_profiles,
    find_intervals,
    parse_position_m,
    run_duration_seconds,
    signed_delta_ns,
    write_intervals,
    write_phase_profiles,
)


@dataclass(frozen=True)
class GammaWindowResult:
    run_key: str
    position_m: float
    duration_seconds: float
    peak_mu_ns: float
    half_width_ns: float
    gamma_tot_min_ns: float
    gamma_tot_max_ns: float
    phase_center_ns: float
    phase_center_shift_from_1p25_ns: float
    gamma_total_counts: float
    gamma_total_rate_hz: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Test the sensitivity of the prompt-gamma phase reference to fixed-width gamma ToT windows "
            "centered on the fitted gamma ToT peak."
        )
    )
    parser.add_argument("--results-root", default="/data6/Dogma_analysis_by_Dachi/Results")
    parser.add_argument("--gamma-fit-summary", default=None)
    parser.add_argument("--output-dir")
    parser.add_argument("--channel", type=int, default=2)
    parser.add_argument("--run-glob", default="NCAL_20us_Pos_5m_*")
    parser.add_argument("--broad-gamma-tot-min-ns", type=float, default=10.0)
    parser.add_argument("--broad-gamma-tot-max-ns", type=float, default=18.0)
    parser.add_argument("--half-widths-ns", type=float, nargs="*", default=[1.0, 1.25, 1.5])
    parser.add_argument("--tot-slice-edges-ns", type=float, nargs="*", default=[8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 22.0, 26.0, 30.0])
    parser.add_argument("--gamma-exclusion-half-width-ns", type=float, default=4.0)
    parser.add_argument("--gamma-display-tot-max-ns", type=float, default=30.0)
    parser.add_argument("--baseline-quantile", type=float, default=0.20)
    parser.add_argument("--z-threshold", type=float, default=5.0)
    parser.add_argument("--min-contiguous-bins", type=int, default=3)
    return parser.parse_args()


def read_peak_summary(path: Path) -> dict[str, dict[str, float]]:
    summary: dict[str, dict[str, float]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            summary[row["run_key"]] = {
                "peak_mu_ns": float(row["peak_mu_ns"]),
                "peak_mu_err_ns": float(row["peak_mu_err_ns"]),
                "peak_sigma_ns": float(row["peak_sigma_ns"]),
            }
    return summary


def detect_gamma_center_for_window(histogram, gamma_tot_min_ns: float, gamma_tot_max_ns: float):
    period_ns = float(histogram.metadata["deduced_period_ns"])
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    gamma_tot_mask = mask_between(histogram.y_centers_ns, gamma_tot_min_ns, gamma_tot_max_ns)
    projection_3x = histogram.counts[gamma_tot_mask, :].sum(axis=0)
    projection_1x = np.zeros(period_bins, dtype=np.float64)
    for index, count in enumerate(projection_3x):
        projection_1x[index % period_bins] += count
    smoothed = np.convolve(np.r_[projection_1x[-2:], projection_1x, projection_1x[:2]], np.ones(5) / 5.0, mode="same")[2:-2]
    gamma_center_ns = float(phase_centers_ns[int(np.argmax(smoothed))])
    gamma_relative_phase_ns = signed_delta_ns(phase_centers_ns, gamma_center_ns, period_ns)
    order = np.argsort(gamma_relative_phase_ns)
    return gamma_center_ns, gamma_relative_phase_ns[order], projection_1x[order], period_ns


def build_slice_profile(histogram, run_key, position_m, duration_seconds, period_ns, gamma_center_ns, sorted_relative_phase_ns, tot_min_ns, tot_max_ns, args) -> SliceProfile:
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    tot_mask = mask_between(histogram.y_centers_ns, tot_min_ns, tot_max_ns)
    projection_3x = histogram.counts[tot_mask, :].sum(axis=0)
    projection_1x = np.zeros(period_bins, dtype=np.float64)
    for index, count in enumerate(projection_3x):
        projection_1x[index % period_bins] += count
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    relative_phase_ns = signed_delta_ns(phase_centers_ns, gamma_center_ns, period_ns)
    order = np.argsort(relative_phase_ns)
    sorted_counts = projection_1x[order]
    sorted_rate = sorted_counts / duration_seconds
    baseline_candidates = sorted_counts[np.abs(sorted_relative_phase_ns) > args.gamma_exclusion_half_width_ns]
    baseline_counts = float(np.quantile(baseline_candidates if baseline_candidates.size else sorted_counts, args.baseline_quantile))
    z_score = (sorted_counts - baseline_counts) / np.sqrt(np.maximum(baseline_counts, 1.0))
    return SliceProfile(
        run_key=run_key,
        position_m=position_m,
        tot_min_ns=tot_min_ns,
        tot_max_ns=tot_max_ns,
        phase_center_ns=gamma_center_ns,
        gamma_relative_phase_ns=sorted_relative_phase_ns,
        rate_hz=sorted_rate,
        counts=sorted_counts,
        baseline_rate_hz=baseline_counts / duration_seconds,
        baseline_counts_per_bin=baseline_counts,
        z_score=z_score,
    )


def build_phase_profile(histogram, duration_seconds, period_ns, gamma_center_ns, tot_min_ns, tot_max_ns):
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    tot_mask = mask_between(histogram.y_centers_ns, tot_min_ns, tot_max_ns)
    projection_3x = histogram.counts[tot_mask, :].sum(axis=0)
    projection_1x = np.zeros(period_bins, dtype=np.float64)
    for index, count in enumerate(projection_3x):
        projection_1x[index % period_bins] += count
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * histogram.x_bin_width_ns
    relative_phase_ns = signed_delta_ns(phase_centers_ns, gamma_center_ns, period_ns)
    order = np.argsort(relative_phase_ns)
    return relative_phase_ns[order], projection_1x[order], projection_1x[order] / duration_seconds


def draw_gamma_cut_comparison(
    histogram,
    run_key: str,
    position_m: float,
    duration_seconds: float,
    broad_gamma_center_ns: float,
    cut_gamma_center_ns: float,
    period_ns: float,
    broad_tot_min_ns: float,
    broad_tot_max_ns: float,
    cut_tot_min_ns: float,
    cut_tot_max_ns: float,
    output_path: Path,
    profile_path: Path,
    args: argparse.Namespace,
) -> None:
    broad_ref_relative_ns = float(signed_delta_ns(np.asarray([broad_gamma_center_ns]), cut_gamma_center_ns, period_ns)[0])
    broad_phase_ns, broad_counts, broad_rate = build_phase_profile(
        histogram, duration_seconds, period_ns, cut_gamma_center_ns, broad_tot_min_ns, broad_tot_max_ns
    )
    cut_phase_ns, cut_counts, cut_rate = build_phase_profile(
        histogram, duration_seconds, period_ns, cut_gamma_center_ns, cut_tot_min_ns, cut_tot_max_ns
    )

    with profile_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "run_key",
                "phase_relative_ns",
                "broad_gamma_counts",
                "broad_gamma_rate_hz",
                "cut_gamma_counts",
                "cut_gamma_rate_hz",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for phase_ns, broad_c, broad_r, cut_c, cut_r in zip(broad_phase_ns, broad_counts, broad_rate, cut_counts, cut_rate):
            writer.writerow(
                {
                    "run_key": run_key,
                    "phase_relative_ns": phase_ns,
                    "broad_gamma_counts": broad_c,
                    "broad_gamma_rate_hz": broad_r,
                    "cut_gamma_counts": cut_c,
                    "cut_gamma_rate_hz": cut_r,
                }
            )

    display_tot_mask = histogram.y_centers_ns <= args.gamma_display_tot_max_ns
    display_counts = histogram.counts[display_tot_mask, :]
    vmax = float(np.quantile(display_counts[display_counts > 0.0], 0.995)) if np.any(display_counts > 0.0) else 1.0

    figure, axes = plt.subplots(1, 2, figsize=(13.0, 4.8), constrained_layout=True)
    image = axes[0].imshow(
        np.ma.masked_less_equal(display_counts, 0.0),
        aspect="auto",
        origin="lower",
        extent=(
            float(histogram.x_centers_ns[0] - 0.5 * histogram.x_bin_width_ns),
            float(histogram.x_centers_ns[-1] + 0.5 * histogram.x_bin_width_ns),
            0.0,
            args.gamma_display_tot_max_ns,
        ),
        interpolation="nearest",
        cmap="magma",
        vmin=0.0,
        vmax=max(1.0, vmax),
    )
    figure.colorbar(image, ax=axes[0], label="counts/bin, linear scale, clipped at q0.995")
    axes[0].axhspan(broad_tot_min_ns, broad_tot_max_ns, color="#2563eb", alpha=0.10)
    axes[0].axhspan(cut_tot_min_ns, cut_tot_max_ns, color="#f59e0b", alpha=0.20)
    for fold_index in range(int(histogram.metadata.get("fold_factor", "3"))):
        axes[0].axvline(fold_index * period_ns + broad_gamma_center_ns, color="#2563eb", linestyle="--", linewidth=1.1)
        axes[0].axvline(fold_index * period_ns + cut_gamma_center_ns, color="#00e5ff", linestyle="--", linewidth=1.2)
    axes[0].set_xlabel("RF phase, three periods (ns)")
    axes[0].set_ylabel("ToT (ns)")
    axes[0].set_title(f"{run_key} | broad vs centered cut on 2D")

    axes[1].plot(broad_phase_ns, broad_rate, color="#0891b2", linewidth=1.2, label="broad gamma band rate")
    axes[1].plot(cut_phase_ns, cut_rate, color="#dc2626", linewidth=1.2, label="centered cut gamma rate")
    axes[1].axvline(broad_ref_relative_ns, color="#2563eb", linestyle="--", linewidth=1.1, label="broad-band reference")
    axes[1].axvline(0.0, color="#00e5ff", linestyle="--", linewidth=1.2, label="cut-band reference")
    axes[1].set_xlabel("RF phase relative to selected gamma reference (ns)")
    axes[1].set_ylabel("rate per phase bin (Hz)")
    axes[1].set_title(f"{run_key} | cut vs no-cut comparison")
    axes[1].grid(True, linewidth=0.3, alpha=0.5)
    axes[1].legend(fontsize=8)

    figure.suptitle(f"Prompt-gamma cut comparison | detector face {position_m:g} m")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_results(results: list[GammaWindowResult], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[field.name for field in GammaWindowResult.__dataclass_fields__.values()], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for result in results:
            writer.writerow({field: getattr(result, field) for field in GammaWindowResult.__dataclass_fields__})


def draw_phase_center_summary(results: list[GammaWindowResult], output_path: Path) -> None:
    if not results:
        return
    figure, axis = plt.subplots(figsize=(8.8, 5.2), constrained_layout=True)
    by_run: dict[str, list[GammaWindowResult]] = {}
    for result in results:
        by_run.setdefault(result.run_key, []).append(result)
    colors = ["#2563eb", "#dc2626", "#059669", "#d97706"]
    for color, (run_key, group) in zip(colors, sorted(by_run.items())):
        group = sorted(group, key=lambda item: item.half_width_ns)
        axis.plot([item.half_width_ns for item in group], [item.phase_center_shift_from_1p25_ns for item in group], "o-", color=color, label=run_key)
    axis.axhline(0.0, color="black", linewidth=0.8)
    axis.set_xlabel("Gamma ToT half-width around fitted mu (ns)")
    axis.set_ylabel("phase-center shift relative to ±1.25 ns (ns)")
    axis.set_title("5 m gamma-reference sensitivity to centered fixed-width ToT gate")
    axis.grid(True, linewidth=0.3, alpha=0.5)
    axis.legend(fontsize=8)
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_readme(output_dir: Path, args: argparse.Namespace) -> None:
    text = f"""# Step09 Gamma Window Sensitivity

This directory tests the 5 m prompt-gamma reference against centered fixed-width gamma ToT windows.

Method:
- use the fitted gamma ToT peak position `mu_gamma` from step08;
- compare the original broad gamma band [{args.broad_gamma_tot_min_ns}, {args.broad_gamma_tot_max_ns}] ns against centered cut gates `[mu_gamma - half_width, mu_gamma + half_width]`;
- scan half widths: {', '.join(f'{width:g} ns' for width in args.half_widths_ns)};
- rebuild the prompt-gamma phase reference and the 1D phase-rate slice diagnostics for each width.

Interpretation:
- the spread in `phase_center_ns` across the width scan is a direct measure of gamma-reference sensitivity to the ToT gate choice;
- the `*_gamma_reference.png` plots show both references explicitly: broad-band reference and centered-cut reference;
- if those two lines overlap, the broad and cut selections pick the same prompt-gamma phase bin;
- this is a systematic-selection study, not a statistical fit error.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    results_root = Path(args.results_root).expanduser().resolve()
    gamma_fit_summary = (
        Path(args.gamma_fit_summary).expanduser().resolve()
        if args.gamma_fit_summary
        else results_root / "Ch02_Position_Comparison" / "gamma_tot_peak_fits" / "gamma_tot_peak_fit_summary.tsv"
    )
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else results_root / "Ch02_Position_Comparison" / "gamma_window_sensitivity_5m"
    output_dir.mkdir(parents=True, exist_ok=True)

    peak_summary = read_peak_summary(gamma_fit_summary)
    sparse_files = find_sparse_files(results_root, args.run_glob)
    if not sparse_files:
        raise FileNotFoundError(f"No Folded_RF_3x phase-ToT sparse files found under {results_root}")

    if len(args.half_widths_ns) == 0:
        raise ValueError("Need at least one half-width")
    reference_half_width = min(args.half_widths_ns, key=lambda value: abs(value - 1.25))

    results: list[GammaWindowResult] = []
    intervals: list[object] = []
    profiles: list[SliceProfile] = []

    for sparse_path in sparse_files:
        histogram = read_channel_histogram(sparse_path, args.channel)
        run_key = histogram.metadata.get("run_key", parse_run_key(histogram.path))
        if run_key not in peak_summary:
            continue
        position_m = parse_position_m(run_key)
        duration_seconds = run_duration_seconds(histogram, run_key)
        mu_gamma = peak_summary[run_key]["peak_mu_ns"]
        broad_gamma_center_ns, _, _, _ = detect_gamma_center_for_window(histogram, args.broad_gamma_tot_min_ns, args.broad_gamma_tot_max_ns)
        phase_centers_by_width: dict[float, float] = {}

        for half_width_ns in args.half_widths_ns:
            gamma_tot_min_ns = mu_gamma - half_width_ns
            gamma_tot_max_ns = mu_gamma + half_width_ns
            gamma_center_ns, sorted_relative_phase_ns, sorted_counts, period_ns = detect_gamma_center_for_window(histogram, gamma_tot_min_ns, gamma_tot_max_ns)
            phase_centers_by_width[half_width_ns] = gamma_center_ns
            args.gamma_tot_min_ns = gamma_tot_min_ns
            args.gamma_tot_max_ns = gamma_tot_max_ns

            width_tag = f"hw{half_width_ns:.2f}".replace(".", "p")
            draw_gamma_cut_comparison(
                histogram,
                run_key,
                position_m,
                duration_seconds,
                broad_gamma_center_ns,
                gamma_center_ns,
                period_ns,
                args.broad_gamma_tot_min_ns,
                args.broad_gamma_tot_max_ns,
                gamma_tot_min_ns,
                gamma_tot_max_ns,
                output_dir / f"{run_key}_ch{args.channel:02d}_{width_tag}_gamma_reference.png",
                output_dir / f"{run_key}_ch{args.channel:02d}_{width_tag}_gamma_reference_profile.tsv",
                args,
            )

            run_profiles: list[SliceProfile] = []
            for tot_min_ns, tot_max_ns in zip(args.tot_slice_edges_ns[:-1], args.tot_slice_edges_ns[1:]):
                profile = build_slice_profile(
                    histogram,
                    run_key,
                    position_m,
                    duration_seconds,
                    period_ns,
                    gamma_center_ns,
                    sorted_relative_phase_ns,
                    tot_min_ns,
                    tot_max_ns,
                    args,
                )
                run_profiles.append(profile)
                profiles.append(profile)
                intervals.extend(find_intervals(profile, args))
            draw_profiles(run_key, position_m, run_profiles, output_dir / f"{run_key}_ch{args.channel:02d}_{width_tag}_phase_rate_slices.png", args.gamma_exclusion_half_width_ns)

            results.append(
                GammaWindowResult(
                    run_key=run_key,
                    position_m=position_m,
                    duration_seconds=duration_seconds,
                    peak_mu_ns=mu_gamma,
                    half_width_ns=half_width_ns,
                    gamma_tot_min_ns=gamma_tot_min_ns,
                    gamma_tot_max_ns=gamma_tot_max_ns,
                    phase_center_ns=gamma_center_ns,
                    phase_center_shift_from_1p25_ns=0.0,
                    gamma_total_counts=float(np.sum(sorted_counts)),
                    gamma_total_rate_hz=float(np.sum(sorted_counts) / duration_seconds),
                )
            )

        reference_center = phase_centers_by_width[reference_half_width]
        for index in range(len(results) - len(args.half_widths_ns), len(results)):
            phase_center = results[index].phase_center_ns
            shift = signed_delta_ns(np.asarray([phase_center]), reference_center, period_ns)[0]
            results[index] = GammaWindowResult(
                run_key=results[index].run_key,
                position_m=results[index].position_m,
                duration_seconds=results[index].duration_seconds,
                peak_mu_ns=results[index].peak_mu_ns,
                half_width_ns=results[index].half_width_ns,
                gamma_tot_min_ns=results[index].gamma_tot_min_ns,
                gamma_tot_max_ns=results[index].gamma_tot_max_ns,
                phase_center_ns=results[index].phase_center_ns,
                phase_center_shift_from_1p25_ns=float(shift),
                gamma_total_counts=results[index].gamma_total_counts,
                gamma_total_rate_hz=results[index].gamma_total_rate_hz,
            )
        print(f"{run_key}: completed centered gamma-window sensitivity scan")

    results.sort(key=lambda item: (item.position_m, item.run_key, item.half_width_ns))
    write_results(results, output_dir / "gamma_window_sensitivity_summary.tsv")
    write_phase_profiles(profiles, output_dir / "phase_rate_slices.tsv")
    write_intervals(intervals, output_dir / "phase_rate_excess_intervals.tsv")
    draw_phase_center_summary(results, output_dir / "gamma_window_phase_center_summary.png")
    write_readme(output_dir, args)
    print(f"Wrote {output_dir}")


if __name__ == "__main__":
    main()