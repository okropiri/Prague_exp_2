from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from step05_estimate_gamma_neutron_background import read_channel_histogram


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit a simple linear RF-phase walk correction phase_corrected = phase - const * ToT "
            "from a folded phase-ToT histogram and plot the corrected RF-phase projection."
        )
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--channel", type=int, default=2)
    parser.add_argument("--fit-tot-min-ns", type=float, default=10.0)
    parser.add_argument("--fit-tot-max-ns", type=float, default=18.0)
    parser.add_argument("--projection-tot-min-ns", type=float, default=8.0)
    parser.add_argument("--projection-tot-max-ns", type=float, default=30.0)
    parser.add_argument("--phase-window-half-width-ns", type=float, default=6.0)
    parser.add_argument("--min-row-counts", type=float, default=200.0)
    parser.add_argument("--display-tot-max-ns", type=float, default=30.0)
    return parser.parse_args()


def positive_mod(value: np.ndarray | float, period_ns: float) -> np.ndarray | float:
    wrapped = np.mod(value, period_ns)
    if isinstance(wrapped, np.ndarray):
        return np.where(wrapped < 0.0, wrapped + period_ns, wrapped)
    return wrapped + period_ns if wrapped < 0.0 else wrapped


def signed_delta_ns(values_ns: np.ndarray, center_ns: float, period_ns: float) -> np.ndarray:
    return (values_ns - center_ns + 0.5 * period_ns) % period_ns - 0.5 * period_ns


def fold_to_one_period(histogram, period_ns: float) -> np.ndarray:
    period_bins = max(1, round(period_ns / histogram.x_bin_width_ns))
    folded = np.zeros((histogram.counts.shape[0], period_bins), dtype=np.float64)
    for x_bin in range(histogram.counts.shape[1]):
        folded[:, x_bin % period_bins] += histogram.counts[:, x_bin]
    return folded


def weighted_window_center(
    phase_centers_ns: np.ndarray,
    counts: np.ndarray,
    seed_center_ns: float,
    period_ns: float,
    half_width_ns: float,
) -> tuple[float, float]:
    relative_phase_ns = signed_delta_ns(phase_centers_ns, seed_center_ns, period_ns)
    mask = np.abs(relative_phase_ns) <= half_width_ns
    selected_counts = counts[mask]
    total_counts = float(np.sum(selected_counts))
    if total_counts <= 0.0:
        return float("nan"), 0.0
    mean_relative_ns = float(np.sum(relative_phase_ns[mask] * selected_counts) / total_counts)
    return mean_relative_ns, total_counts


def fit_walk_constant(
    folded_counts: np.ndarray,
    phase_centers_ns: np.ndarray,
    tot_centers_ns: np.ndarray,
    period_ns: float,
    fit_tot_min_ns: float,
    fit_tot_max_ns: float,
    phase_window_half_width_ns: float,
    min_row_counts: float,
) -> tuple[float, float, np.ndarray, np.ndarray, np.ndarray]:
    fit_tot_mask = (tot_centers_ns >= fit_tot_min_ns) & (tot_centers_ns < fit_tot_max_ns)
    prompt_projection = np.sum(folded_counts[fit_tot_mask, :], axis=0)
    seed_center_ns = float(phase_centers_ns[int(np.argmax(prompt_projection))])

    used_tot_ns: list[float] = []
    used_phase_offsets_ns: list[float] = []
    used_weights: list[float] = []
    for tot_index, tot_ns in enumerate(tot_centers_ns):
        if not fit_tot_mask[tot_index]:
            continue
        row_counts = folded_counts[tot_index, :]
        phase_offset_ns, total_counts = weighted_window_center(
            phase_centers_ns,
            row_counts,
            seed_center_ns,
            period_ns,
            phase_window_half_width_ns,
        )
        if not np.isfinite(phase_offset_ns) or total_counts < min_row_counts:
            continue
        used_tot_ns.append(float(tot_ns))
        used_phase_offsets_ns.append(phase_offset_ns)
        used_weights.append(total_counts)

    if len(used_tot_ns) < 2:
        return 0.0, seed_center_ns, np.asarray([]), np.asarray([]), np.asarray([])

    used_tot = np.asarray(used_tot_ns, dtype=np.float64)
    used_offsets = np.asarray(used_phase_offsets_ns, dtype=np.float64)
    weights = np.asarray(used_weights, dtype=np.float64)
    tot_reference_ns = float(np.average(used_tot, weights=weights))
    coeffs = np.polyfit(used_tot - tot_reference_ns, used_offsets, deg=1, w=np.sqrt(weights))
    walk_const = float(coeffs[0])
    return walk_const, seed_center_ns, used_tot, used_offsets, weights


def build_corrected_histogram(
    folded_counts: np.ndarray,
    phase_centers_ns: np.ndarray,
    tot_centers_ns: np.ndarray,
    period_ns: float,
    walk_const_ns_per_ns: float,
) -> np.ndarray:
    corrected = np.zeros_like(folded_counts, dtype=np.float64)
    phase_bin_width_ns = period_ns / folded_counts.shape[1]
    for tot_index, tot_ns in enumerate(tot_centers_ns):
        corrected_phase_ns = positive_mod(phase_centers_ns - walk_const_ns_per_ns * tot_ns, period_ns)
        corrected_bins = np.floor(corrected_phase_ns / phase_bin_width_ns).astype(int)
        corrected_bins = np.clip(corrected_bins, 0, folded_counts.shape[1] - 1)
        for phase_index, corrected_index in enumerate(corrected_bins):
            corrected[tot_index, corrected_index] += folded_counts[tot_index, phase_index]
    return corrected


def write_row_fit_summary(
    output_path: Path,
    tot_ns: np.ndarray,
    phase_offsets_ns: np.ndarray,
    weights: np.ndarray,
) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["tot_ns", "prompt_phase_offset_ns", "row_counts"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for tot_value, offset_value, weight_value in zip(tot_ns, phase_offsets_ns, weights):
            writer.writerow(
                {
                    "tot_ns": float(tot_value),
                    "prompt_phase_offset_ns": float(offset_value),
                    "row_counts": float(weight_value),
                }
            )


def main() -> None:
    args = parse_args()
    histogram = read_channel_histogram(Path(args.input), args.channel)
    run_key = histogram.metadata.get("run_key", Path(args.input).stem)
    period_ns = float(histogram.metadata["deduced_period_ns"])
    one_period_counts = fold_to_one_period(histogram, period_ns)
    period_bins = one_period_counts.shape[1]
    phase_bin_width_ns = period_ns / period_bins
    phase_centers_ns = (np.arange(period_bins, dtype=np.float64) + 0.5) * phase_bin_width_ns
    tot_centers_ns = histogram.y_centers_ns.copy()

    walk_const, seed_center_ns, used_tot, used_offsets, weights = fit_walk_constant(
        one_period_counts,
        phase_centers_ns,
        tot_centers_ns,
        period_ns,
        args.fit_tot_min_ns,
        args.fit_tot_max_ns,
        args.phase_window_half_width_ns,
        args.min_row_counts,
    )
    corrected_counts = build_corrected_histogram(one_period_counts, phase_centers_ns, tot_centers_ns, period_ns, walk_const)

    projection_tot_mask = (tot_centers_ns >= args.projection_tot_min_ns) & (tot_centers_ns < args.projection_tot_max_ns)
    raw_projection = np.sum(one_period_counts[projection_tot_mask, :], axis=0)
    corrected_projection = np.sum(corrected_counts[projection_tot_mask, :], axis=0)
    corrected_seed_center_ns = float(phase_centers_ns[int(np.argmax(corrected_projection))])
    raw_relative_phase_ns = signed_delta_ns(phase_centers_ns, seed_center_ns, period_ns)
    corrected_relative_phase_ns = signed_delta_ns(phase_centers_ns, corrected_seed_center_ns, period_ns)
    raw_order = np.argsort(raw_relative_phase_ns)
    corrected_order = np.argsort(corrected_relative_phase_ns)

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    write_row_fit_summary(output_prefix.with_name(output_prefix.name + "_row_fit.tsv"), used_tot, used_offsets, weights)
    with output_prefix.with_name(output_prefix.name + "_summary.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"run_key={run_key}\n")
        handle.write(f"period_ns={period_ns:.12f}\n")
        handle.write(f"fit_tot_min_ns={args.fit_tot_min_ns:.6f}\n")
        handle.write(f"fit_tot_max_ns={args.fit_tot_max_ns:.6f}\n")
        handle.write(f"projection_tot_min_ns={args.projection_tot_min_ns:.6f}\n")
        handle.write(f"projection_tot_max_ns={args.projection_tot_max_ns:.6f}\n")
        handle.write(f"phase_window_half_width_ns={args.phase_window_half_width_ns:.6f}\n")
        handle.write(f"walk_const_ns_per_ns={walk_const:.12f}\n")
        handle.write(f"raw_seed_center_ns={seed_center_ns:.12f}\n")
        handle.write(f"corrected_seed_center_ns={corrected_seed_center_ns:.12f}\n")
        handle.write(f"rows_used_for_fit={len(used_tot)}\n")

    display_tot_mask = tot_centers_ns <= args.display_tot_max_ns
    vmax_raw = float(np.quantile(one_period_counts[display_tot_mask, :][one_period_counts[display_tot_mask, :] > 0.0], 0.995)) if np.any(one_period_counts[display_tot_mask, :] > 0.0) else 1.0
    vmax_corrected = float(np.quantile(corrected_counts[display_tot_mask, :][corrected_counts[display_tot_mask, :] > 0.0], 0.995)) if np.any(corrected_counts[display_tot_mask, :] > 0.0) else 1.0
    vmax = max(1.0, min(vmax_raw, vmax_corrected))

    figure, axes = plt.subplots(1, 3, figsize=(17.0, 4.8), constrained_layout=True)
    extent = (0.0, period_ns, 0.0, args.display_tot_max_ns)
    image0 = axes[0].imshow(
        np.ma.masked_less_equal(one_period_counts[display_tot_mask, :], 0.0),
        aspect="auto",
        origin="lower",
        extent=extent,
        interpolation="nearest",
        cmap="magma",
        vmin=0.0,
        vmax=vmax,
    )
    figure.colorbar(image0, ax=axes[0], label="counts / 2D bin")
    axes[0].set_xlabel("RF phase (ns)")
    axes[0].set_ylabel("ToT (ns)")
    axes[0].set_title(f"{run_key} | raw one-period phase-ToT")

    image1 = axes[1].imshow(
        np.ma.masked_less_equal(corrected_counts[display_tot_mask, :], 0.0),
        aspect="auto",
        origin="lower",
        extent=extent,
        interpolation="nearest",
        cmap="magma",
        vmin=0.0,
        vmax=vmax,
    )
    figure.colorbar(image1, ax=axes[1], label="counts / 2D bin")
    axes[1].set_xlabel("RF phase corrected (ns)")
    axes[1].set_ylabel("ToT (ns)")
    axes[1].set_title(f"{run_key} | corrected with const={walk_const:.4f} ns/ns")

    axes[2].plot(raw_relative_phase_ns[raw_order], raw_projection[raw_order], color="#2563eb", linewidth=1.4, label="raw projection")
    axes[2].plot(corrected_relative_phase_ns[corrected_order], corrected_projection[corrected_order], color="#dc2626", linewidth=1.4, label="walk-corrected projection")
    axes[2].set_xlabel("RF phase relative to prompt ridge (ns)")
    axes[2].set_ylabel("counts / phase bin")
    axes[2].set_title(
        f"projection for ToT {args.projection_tot_min_ns:g}-{args.projection_tot_max_ns:g} ns\n"
        f"fit band {args.fit_tot_min_ns:g}-{args.fit_tot_max_ns:g} ns"
    )
    axes[2].grid(True, linewidth=0.3, alpha=0.5)
    axes[2].legend(fontsize=8)

    figure.suptitle("RF-phase walk-correction diagnostic")
    figure.savefig(output_prefix.with_suffix(".png"), dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()