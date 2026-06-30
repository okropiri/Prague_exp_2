from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["agg.path.chunksize"] = 10000

import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.lines import Line2D
from uproot.writing.identify import to_TAxis, to_TH1x, to_TH2x


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write PDF and PNG diagnostics for the NCAL-only ch0-referenced RF scan."
    )
    parser.add_argument("--input-prefix", required=True)
    parser.add_argument("--pdf-prefix", required=True)
    parser.add_argument("--png-prefix", required=True)
    parser.add_argument("--root-output")
    return parser.parse_args()


def parse_metadata_line(line: str) -> tuple[str, str] | None:
    if not line.startswith("#") or "=" not in line:
        return None
    key, value = line[1:].split("=", 1)
    return key.strip(), value.strip()


def read_hist1(path: Path) -> dict[str, object]:
    metadata: dict[str, object] = {}
    rows: list[tuple[int, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            parsed = parse_metadata_line(line)
            if parsed is not None:
                key, value = parsed
                metadata[key] = value
                continue
            if not line or line.startswith("#"):
                continue
            index_text, count_text = line.split()
            rows.append((int(index_text), float(count_text)))

    bins = int(metadata["bins"])
    storage = np.zeros(bins + 2, dtype=np.float64)
    for index, count in rows:
        storage[index] = count

    metadata["bins"] = bins
    metadata["min"] = float(metadata["min"])
    metadata["max"] = float(metadata["max"])
    metadata["storage"] = storage
    return metadata


def read_hist2(path: Path) -> dict[str, object]:
    metadata: dict[str, object] = {}
    rows: list[tuple[int, int, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            parsed = parse_metadata_line(line)
            if parsed is not None:
                key, value = parsed
                metadata[key] = value
                continue
            if not line or line.startswith("#"):
                continue
            x_text, y_text, count_text = line.split()
            rows.append((int(x_text), int(y_text), float(count_text)))

    x_bins = int(metadata["x_bins"])
    y_bins = int(metadata["y_bins"])
    storage = np.zeros((y_bins + 2, x_bins + 2), dtype=np.float64)
    for x_index, y_index, count in rows:
        storage[y_index, x_index] = count

    metadata["x_bins"] = x_bins
    metadata["x_min"] = float(metadata["x_min"])
    metadata["x_max"] = float(metadata["x_max"])
    metadata["y_bins"] = y_bins
    metadata["y_min"] = float(metadata["y_min"])
    metadata["y_max"] = float(metadata["y_max"])
    metadata["storage"] = storage
    return metadata


def read_scan(path: Path) -> tuple[dict[str, str], list[dict[str, object]]]:
    metadata: dict[str, str] = {}
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            parsed = parse_metadata_line(line)
            if parsed is not None:
                key, value = parsed
                metadata[key] = value
                continue
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            rows.append(
                {
                    "round_index": int(parts[0]),
                    "period_ns": float(parts[1]),
                    "phase_origin_ns": float(parts[2]),
                    "peak_center_ns": float(parts[3]),
                    "peak_height": int(parts[4]),
                    "selected_pulses": int(parts[5]),
                    "selected_fraction": float(parts[6]),
                    "sigma_ns": float(parts[7]),
                    "mean_residual_ns": float(parts[8]),
                    "coherence": float(parts[9]),
                    "drift_slope_ns_per_cycle": float(parts[10]),
                    "drift_intercept_ns": float(parts[11]),
                    "merit": float(parts[12]),
                    "valid": parts[13] == "true",
                }
            )
    return metadata, rows


def parse_range_text(value: str | None) -> tuple[float, float] | None:
    if not value:
        return None
    text = value.strip()
    if not (text.startswith("[") and text.endswith("]")):
        return None
    left_text, separator, right_text = text[1:-1].partition(",")
    if not separator:
        return None
    left = float(left_text.strip())
    right = float(right_text.strip())
    if left >= right:
        return None
    return left, right


def best_valid_row(rows: list[dict[str, object]]) -> dict[str, object]:
    valid_rows = [row for row in rows if bool(row["valid"])]
    if not valid_rows:
        raise RuntimeError("Scan table does not contain a valid candidate")
    return max(valid_rows, key=lambda row: float(row["merit"]))


def best_selected_fraction_row(rows: list[dict[str, object]]) -> dict[str, object]:
    valid_rows = [row for row in rows if bool(row["valid"])]
    if not valid_rows:
        raise RuntimeError("Scan table does not contain a valid candidate")
    return max(valid_rows, key=lambda row: float(row["selected_fraction"]))


def best_sigma_row(rows: list[dict[str, object]]) -> dict[str, object]:
    valid_rows = [row for row in rows if bool(row["valid"])]
    if not valid_rows:
        raise RuntimeError("Scan table does not contain a valid candidate")
    return min(valid_rows, key=lambda row: float(row["sigma_ns"]))


def round_step_sizes(periods: np.ndarray, rounds: np.ndarray) -> dict[int, float | None]:
    step_sizes: dict[int, float | None] = {}
    for round_index in sorted(set(int(value) for value in rounds.tolist())):
        round_periods = np.unique(periods[rounds == round_index])
        if round_periods.size < 2:
            step_sizes[round_index] = None
            continue
        deltas = np.diff(np.sort(round_periods))
        step_sizes[round_index] = float(np.median(deltas))
    return step_sizes


def round_label(round_index: int, step_size: float | None) -> str:
    if round_index == 0:
        name = "coarse"
    elif round_index == 1:
        name = "fine"
    elif round_index == 2:
        name = "finest"
    else:
        name = f"refine {round_index}"
    if step_size is None:
        return name
    return f"{name}: {step_size:.4f} ns"


def read_cycle_residuals(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    cycle_index: list[int] = []
    counts: list[int] = []
    means: list[float] = []
    rms_values: list[float] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            cycle_text, count_text, mean_text, rms_text = line.split()
            cycle_index.append(int(cycle_text))
            counts.append(int(count_text))
            means.append(float(mean_text))
            rms_values.append(float(rms_text))
    return (
        np.asarray(cycle_index, dtype=np.int64),
        np.asarray(counts, dtype=np.int64),
        np.asarray(means, dtype=np.float64),
        np.asarray(rms_values, dtype=np.float64),
    )


def read_cycle_residual_points(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cycle_index: list[int] = []
    residuals: list[float] = []
    global_time_seconds: list[float] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            cycle_text, residual_text = parts[0], parts[1]
            cycle_index.append(int(cycle_text))
            residuals.append(float(residual_text))
            if len(parts) >= 3:
                global_time_seconds.append(float(parts[2]))
            else:
                global_time_seconds.append(np.nan)
    return (
        np.asarray(cycle_index, dtype=np.int64),
        np.asarray(residuals, dtype=np.float64),
        np.asarray(global_time_seconds, dtype=np.float64),
    )


def hist1_centers(hist: dict[str, object]) -> np.ndarray:
    bins = int(hist["bins"])
    x_min = float(hist["min"])
    x_max = float(hist["max"])
    return x_min + (np.arange(bins, dtype=np.float64) + 0.5) * ((x_max - x_min) / bins)


def hist1_values(hist: dict[str, object]) -> np.ndarray:
    return np.asarray(hist["storage"], dtype=np.float64)[1:-1]


def hist2_values(hist: dict[str, object]) -> np.ndarray:
    return np.asarray(hist["storage"], dtype=np.float64)[1:-1, 1:-1]


def hist1_edges(hist: dict[str, object]) -> np.ndarray:
    bins = int(hist["bins"])
    x_min = float(hist["min"])
    x_max = float(hist["max"])
    return np.linspace(x_min, x_max, bins + 1, dtype=np.float64)


def hist2_x_edges(hist: dict[str, object]) -> np.ndarray:
    x_bins = int(hist["x_bins"])
    x_min = float(hist["x_min"])
    x_max = float(hist["x_max"])
    return np.linspace(x_min, x_max, x_bins + 1, dtype=np.float64)


def hist2_y_edges(hist: dict[str, object]) -> np.ndarray:
    y_bins = int(hist["y_bins"])
    y_min = float(hist["y_min"])
    y_max = float(hist["y_max"])
    return np.linspace(y_min, y_max, y_bins + 1, dtype=np.float64)


def downsample_for_scatter(x_values: np.ndarray, y_values: np.ndarray, max_points: int = 400_000) -> tuple[np.ndarray, np.ndarray]:
    if x_values.size <= max_points:
        return x_values, y_values
    indices = np.linspace(0, x_values.size - 1, max_points, dtype=np.int64)
    return x_values[indices], y_values[indices]


def write_root_outputs(root_output: Path, phase_profile: dict[str, object], phase_tot: dict[str, object]) -> None:
    profile_storage = np.asarray(phase_profile["storage"], dtype=np.float64)
    profile_values = profile_storage[1:-1]
    profile_centers = hist1_centers(phase_profile)
    profile_edges = hist1_edges(phase_profile)
    profile_entries = float(profile_values.sum())

    profile_hist = to_TH1x(
        fName=None,
        fTitle=str(phase_profile["title"]),
        data=profile_storage,
        fEntries=profile_entries,
        fTsumw=profile_entries,
        fTsumw2=profile_entries,
        fTsumwx=float((profile_values * profile_centers).sum()),
        fTsumwx2=float((profile_values * profile_centers * profile_centers).sum()),
        fSumw2=None,
        fXaxis=to_TAxis(
            "xaxis",
            str(phase_profile["x_title"]),
            int(phase_profile["bins"]),
            float(profile_edges[0]),
            float(profile_edges[-1]),
        ),
    )

    phase_storage = np.asarray(phase_tot["storage"], dtype=np.float64)
    phase_values = phase_storage[1:-1, 1:-1]
    x_centers = hist1_centers(
        {
            "bins": phase_tot["x_bins"],
            "min": phase_tot["x_min"],
            "max": phase_tot["x_max"],
        }
    )
    y_centers = hist1_centers(
        {
            "bins": phase_tot["y_bins"],
            "min": phase_tot["y_min"],
            "max": phase_tot["y_max"],
        }
    )
    x_edges = hist2_x_edges(phase_tot)
    y_edges = hist2_y_edges(phase_tot)
    x_centers_2d = x_centers[np.newaxis, :]
    y_centers_2d = y_centers[:, np.newaxis]
    phase_entries = float(phase_values.sum())

    phase_hist = to_TH2x(
        fName=None,
        fTitle=str(phase_tot["title"]),
        data=phase_storage.T.ravel(),
        fEntries=phase_entries,
        fTsumw=phase_entries,
        fTsumw2=phase_entries,
        fTsumwx=float((phase_values * x_centers_2d).sum()),
        fTsumwx2=float((phase_values * x_centers_2d * x_centers_2d).sum()),
        fTsumwy=float((phase_values * y_centers_2d).sum()),
        fTsumwy2=float((phase_values * y_centers_2d * y_centers_2d).sum()),
        fTsumwxy=float((phase_values * x_centers_2d * y_centers_2d).sum()),
        fSumw2=None,
        fXaxis=to_TAxis(
            "xaxis",
            str(phase_tot["x_title"]),
            int(phase_tot["x_bins"]),
            float(x_edges[0]),
            float(x_edges[-1]),
        ),
        fYaxis=to_TAxis(
            "yaxis",
            str(phase_tot["y_title"]),
            int(phase_tot["y_bins"]),
            float(y_edges[0]),
            float(y_edges[-1]),
        ),
    )

    with uproot.recreate(root_output) as root_file:
        root_file["BestPhaseProfile"] = profile_hist
        root_file["BestPhaseTOT"] = phase_hist


def write_scan_plot(
    rows: list[dict[str, object]],
    pdf_path: Path,
    png_path: Path,
    x_limits: tuple[float, float] | None = None,
) -> None:
    periods = np.asarray([row["period_ns"] for row in rows], dtype=np.float64)
    selected_fractions = np.asarray([row["selected_fraction"] for row in rows], dtype=np.float64)
    sigmas = np.asarray([row["sigma_ns"] for row in rows], dtype=np.float64)
    merits = np.asarray([row["merit"] for row in rows], dtype=np.float64)
    drifts = np.asarray([row["drift_slope_ns_per_cycle"] for row in rows], dtype=np.float64)
    rounds = np.asarray([row["round_index"] for row in rows], dtype=np.int64)
    round_steps = round_step_sizes(periods, rounds)

    best_period_row = best_valid_row(rows)
    best_selected_row = best_selected_fraction_row(rows)
    best_sigma_period_row = best_sigma_row(rows)
    best_period = float(best_period_row["period_ns"])
    norm = mpl.colors.Normalize(vmin=float(np.min(rounds)), vmax=float(np.max(rounds)))
    cmap = mpl.colormaps["viridis"]

    figure, axes = plt.subplots(4, 1, figsize=(10, 14), sharex=True)
    annotation_kwargs = {
        "color": "#d62728",
        "fontsize": 9,
        "bbox": {"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "#d62728", "alpha": 0.8},
    }

    def mark_optimum(axis: plt.Axes, x_value: float, y_value: float, label: str, dx: float, dy: float, va: str = "bottom") -> None:
        axis.scatter(
            [x_value],
            [y_value],
            s=70,
            marker="o",
            facecolors="white",
            edgecolors="#d62728",
            linewidths=1.5,
            zorder=5,
        )
        axis.annotate(
            label,
            xy=(x_value, y_value),
            xytext=(dx, dy),
            textcoords="offset points",
            ha="left",
            va=va,
            **annotation_kwargs,
        )

    axes[0].scatter(periods, selected_fractions, c=rounds, cmap="viridis", s=20)
    axes[0].axvline(best_period, color="#d62728", linewidth=1.0, linestyle="--")
    axes[0].set_ylabel("selected fraction")
    axes[0].set_title("NCAL RF-period scan metrics")
    axes[0].grid(True, linewidth=0.3, alpha=0.4)
    mark_optimum(
        axes[0],
        float(best_selected_row["period_ns"]),
        float(best_selected_row["selected_fraction"]),
        f"sel. max {float(best_selected_row['period_ns']):.4f} ns",
        10,
        -10,
        va="top",
    )
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=cmap(norm(round_index)),
            markeredgecolor=cmap(norm(round_index)),
            markersize=6,
            label=round_label(round_index, round_steps.get(round_index)),
        )
        for round_index in sorted(round_steps)
    ]
    axes[0].legend(
        handles=legend_handles,
        loc="upper right",
        framealpha=0.9,
        facecolor="white",
        edgecolor="#cccccc",
        fontsize=9,
        title="scan rounds",
        title_fontsize=9,
    )

    axes[1].scatter(periods, sigmas, c=rounds, cmap="viridis", s=20)
    axes[1].axvline(best_period, color="#d62728", linewidth=1.0, linestyle="--")
    axes[1].set_ylabel("sigma (ns)")
    axes[1].grid(True, linewidth=0.3, alpha=0.4)
    mark_optimum(
        axes[1],
        float(best_sigma_period_row["period_ns"]),
        float(best_sigma_period_row["sigma_ns"]),
        f"sigma min {float(best_sigma_period_row['period_ns']):.4f} ns",
        10,
        12,
    )

    axes[2].scatter(periods, merits, c=rounds, cmap="viridis", s=20)
    axes[2].axvline(best_period, color="#d62728", linewidth=1.0, linestyle="--")
    axes[2].set_ylabel("merit")
    axes[2].grid(True, linewidth=0.3, alpha=0.4)
    mark_optimum(
        axes[2],
        float(best_period_row["period_ns"]),
        float(best_period_row["merit"]),
        f"merit max {best_period:.4f} ns",
        10,
        -10,
        va="top",
    )

    axes[3].scatter(periods, drifts, c=rounds, cmap="viridis", s=20)
    axes[3].axhline(0.0, color="#444444", linewidth=0.8, linestyle="--")
    axes[3].axvline(best_period, color="#d62728", linewidth=1.0, linestyle="--")
    axes[3].set_xlabel("trial period (ns)")
    axes[3].set_ylabel("drift slope (ns/cycle)")
    axes[3].grid(True, linewidth=0.3, alpha=0.4)

    if x_limits is not None:
        for axis in axes:
            axis.set_xlim(*x_limits)

    figure.tight_layout()
    figure.savefig(pdf_path)
    figure.savefig(png_path, dpi=180)
    plt.close(figure)


def write_phase_profile_plot(hist: dict[str, object], pdf_path: Path, png_path: Path) -> None:
    centers = hist1_centers(hist)
    values = hist1_values(hist)

    figure, axis = plt.subplots(figsize=(10, 5))
    axis.step(centers, values, where="mid", linewidth=1.0, color="#1f77b4")
    axis.set_xlabel(str(hist["x_title"]))
    axis.set_ylabel("NCAL counts / bin")
    axis.set_title(str(hist["title"]))
    axis.grid(True, linewidth=0.3, alpha=0.4)
    figure.tight_layout()
    figure.savefig(pdf_path)
    figure.savefig(png_path, dpi=180)
    plt.close(figure)


def write_phase_tot_plot(hist: dict[str, object], pdf_path: Path, png_path: Path) -> None:
    values = hist2_values(hist)
    x_min = float(hist["x_min"])
    x_max = float(hist["x_max"])
    y_min = float(hist["y_min"])
    y_max = float(hist["y_max"])

    figure, axis = plt.subplots(figsize=(10, 6))
    image = axis.imshow(
        np.log10(values + 1.0),
        origin="lower",
        aspect="auto",
        extent=(x_min, x_max, y_min, y_max),
        interpolation="nearest",
        cmap="magma",
    )
    axis.set_xlabel(str(hist["x_title"]))
    axis.set_ylabel(str(hist["y_title"]))
    axis.set_title(str(hist["title"]))
    figure.colorbar(image, ax=axis, label="log10(count + 1)")
    figure.tight_layout()
    figure.savefig(pdf_path)
    figure.savefig(png_path, dpi=180)
    plt.close(figure)


def write_phase_tot_linear_plot(hist: dict[str, object], pdf_path: Path, png_path: Path) -> None:
    values = hist2_values(hist)
    x_min = float(hist["x_min"])
    x_max = float(hist["x_max"])
    y_min = float(hist["y_min"])
    y_max = float(hist["y_max"])

    figure, axis = plt.subplots(figsize=(10, 6))
    image = axis.imshow(
        values,
        origin="lower",
        aspect="auto",
        extent=(x_min, x_max, y_min, y_max),
        interpolation="nearest",
        cmap="magma",
    )
    axis.set_xlabel(str(hist["x_title"]))
    axis.set_ylabel(str(hist["y_title"]))
    axis.set_title(f"{hist['title']} (linear Z)")
    figure.colorbar(image, ax=axis, label="count")
    figure.tight_layout()
    figure.savefig(pdf_path)
    figure.savefig(png_path, dpi=180)
    plt.close(figure)


def write_cycle_residual_plot(cycle_index: np.ndarray, counts: np.ndarray, means: np.ndarray, rms_values: np.ndarray, pdf_path: Path, png_path: Path) -> None:
    figure, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    axes[0].plot(cycle_index, means, linewidth=0.9, color="#d62728", label="cycle mean residual")
    axes[0].axhline(0.0, color="#444444", linewidth=0.8, linestyle="--", label="zero")
    axes[0].set_ylabel("mean residual (ns)")
    axes[0].set_title("Best-period residual stability by cycle index")
    axes[0].grid(True, linewidth=0.3, alpha=0.4)
    axes[0].legend(loc="upper left", framealpha=0.9)

    axes[1].plot(cycle_index, rms_values, linewidth=0.9, color="#1f77b4", label="RMS residual")
    axes[1].set_xlabel("cycle index")
    axes[1].set_ylabel("RMS residual (ns)")
    axes[1].grid(True, linewidth=0.3, alpha=0.4)

    twin = axes[1].twinx()
    twin.plot(cycle_index, counts, linewidth=0.8, color="#2ca02c", alpha=0.7, label="count per cycle")
    twin.set_ylabel("count per cycle")

    lines_left, labels_left = axes[1].get_legend_handles_labels()
    lines_right, labels_right = twin.get_legend_handles_labels()
    axes[1].legend(lines_left + lines_right, labels_left + labels_right, loc="upper left", framealpha=0.9)

    figure.tight_layout()
    figure.savefig(pdf_path)
    figure.savefig(png_path, dpi=180)
    plt.close(figure)


def compute_local_cycle_slopes(cycle_index: np.ndarray, counts: np.ndarray, means: np.ndarray, half_window: int = 6) -> tuple[np.ndarray, np.ndarray]:
    slopes: list[float] = []
    centers: list[int] = []
    if cycle_index.size < 2 * half_window + 1:
        return np.asarray([], dtype=np.int64), np.asarray([], dtype=np.float64)

    for index in range(half_window, cycle_index.size - half_window):
        left = index - half_window
        right = index + half_window + 1
        x = cycle_index[left:right].astype(np.float64)
        y = means[left:right]
        w = counts[left:right].astype(np.float64)
        if x.size < 3 or np.unique(x).size < 3:
            continue
        w_sum = np.sum(w)
        x_mean = np.sum(w * x) / w_sum
        y_mean = np.sum(w * y) / w_sum
        denominator = np.sum(w * (x - x_mean) ** 2)
        if denominator <= 0.0:
            continue
        slope = np.sum(w * (x - x_mean) * (y - y_mean)) / denominator
        centers.append(int(cycle_index[index]))
        slopes.append(float(slope))
    return np.asarray(centers, dtype=np.int64), np.asarray(slopes, dtype=np.float64)


def write_cycle_residual_detail_plot(point_cycle_index: np.ndarray,
                                     point_residuals: np.ndarray,
                                     cycle_index: np.ndarray,
                                     counts: np.ndarray,
                                     means: np.ndarray,
                                     pdf_path: Path,
                                     png_path: Path) -> None:
    local_cycle_index, local_slopes = compute_local_cycle_slopes(cycle_index, counts, means)
    scatter_cycle_index, scatter_residuals = downsample_for_scatter(point_cycle_index, point_residuals)

    figure, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    axes[0].scatter(
        scatter_cycle_index,
        scatter_residuals,
        s=4,
        alpha=0.12,
        color="#1f77b4",
        edgecolors="none",
        rasterized=True,
    )
    axes[0].plot(cycle_index, means, linewidth=1.2, color="#d62728", label="cycle mean")
    axes[0].axhline(0.0, color="#444444", linewidth=0.8, linestyle="--")
    axes[0].set_ylabel("residual (ns)")
    axes[0].set_title("Best-period residuals vs cycle index")
    axes[0].grid(True, linewidth=0.3, alpha=0.4)
    axes[0].legend(loc="upper left", framealpha=0.9)

    axes[1].plot(local_cycle_index, local_slopes, linewidth=1.1, color="#9467bd")
    axes[1].axhline(0.0, color="#444444", linewidth=0.8, linestyle="--")
    axes[1].set_xlabel("cycle index")
    axes[1].set_ylabel("local slope (ns/cycle)")
    axes[1].set_title("Local residual drift by cycle index")
    axes[1].grid(True, linewidth=0.3, alpha=0.4)

    figure.tight_layout()
    figure.savefig(pdf_path)
    figure.savefig(png_path, dpi=180)
    plt.close(figure)


def compute_local_time_slopes(minutes: np.ndarray, residuals: np.ndarray, half_window: int = 4000) -> tuple[np.ndarray, np.ndarray]:
    valid = np.isfinite(minutes) & np.isfinite(residuals)
    if np.count_nonzero(valid) < 2 * half_window + 1:
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)

    x_all = minutes[valid]
    y_all = residuals[valid]
    order = np.argsort(x_all)
    x_sorted = x_all[order]
    y_sorted = y_all[order]

    window = 2 * half_window + 1
    left = np.arange(0, x_sorted.size - window + 1, dtype=np.int64)
    right = left + window

    prefix_x = np.concatenate(([0.0], np.cumsum(x_sorted, dtype=np.float64)))
    prefix_y = np.concatenate(([0.0], np.cumsum(y_sorted, dtype=np.float64)))
    prefix_xx = np.concatenate(([0.0], np.cumsum(x_sorted * x_sorted, dtype=np.float64)))
    prefix_xy = np.concatenate(([0.0], np.cumsum(x_sorted * y_sorted, dtype=np.float64)))

    count = float(window)
    sum_x = prefix_x[right] - prefix_x[left]
    sum_y = prefix_y[right] - prefix_y[left]
    sum_xx = prefix_xx[right] - prefix_xx[left]
    sum_xy = prefix_xy[right] - prefix_xy[left]

    denominator = sum_xx - (sum_x * sum_x) / count
    numerator = sum_xy - (sum_x * sum_y) / count
    valid_denominator = denominator > 0.0
    if not np.any(valid_denominator):
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64)

    centers = x_sorted[left + half_window][valid_denominator]
    slopes = numerator[valid_denominator] / denominator[valid_denominator]
    return centers.astype(np.float64), slopes.astype(np.float64)


def write_global_time_validation_plot(global_time_seconds: np.ndarray,
                                      point_residuals: np.ndarray,
                                      pdf_path: Path,
                                      png_path: Path) -> None:
    valid = np.isfinite(global_time_seconds) & np.isfinite(point_residuals)
    if not np.any(valid):
        return

    time_minutes = (global_time_seconds[valid] - np.nanmin(global_time_seconds[valid])) / 60.0
    residuals = point_residuals[valid]
    local_minutes, local_slopes = compute_local_time_slopes(time_minutes, residuals)
    scatter_minutes, scatter_residuals = downsample_for_scatter(time_minutes, residuals)

    figure, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    axes[0].scatter(
        scatter_minutes,
        scatter_residuals,
        s=4,
        alpha=0.1,
        color="#1f77b4",
        edgecolors="none",
        rasterized=True,
    )
    axes[0].axhline(0.0, color="#444444", linewidth=0.8, linestyle="--")
    axes[0].set_ylabel("residual (ns)")
    axes[0].set_title("Best-period residuals vs global time")
    axes[0].grid(True, linewidth=0.3, alpha=0.4)

    axes[1].plot(local_minutes, local_slopes, linewidth=1.1, color="#9467bd")
    axes[1].axhline(0.0, color="#444444", linewidth=0.8, linestyle="--")
    axes[1].set_xlabel("global time since run start (min)")
    axes[1].set_ylabel("local slope (ns/min)")
    axes[1].set_title("Local residual drift vs global time")
    axes[1].grid(True, linewidth=0.3, alpha=0.4)

    figure.tight_layout()
    figure.savefig(pdf_path)
    figure.savefig(png_path, dpi=180)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    prefix = Path(args.input_prefix)

    scan_metadata, scan_rows = read_scan(prefix.with_name(prefix.name + "_scan.txt"))
    phase_profile = read_hist1(prefix.with_name(prefix.name + "_best_phase_profile_hist.txt"))
    phase_tot = read_hist2(prefix.with_name(prefix.name + "_best_phase_tot_hist.txt"))
    cycle_index, counts, means, rms_values = read_cycle_residuals(
        prefix.with_name(prefix.name + "_best_cycle_residuals.txt")
    )
    point_cycle_index, point_residuals, point_global_time_seconds = read_cycle_residual_points(
        prefix.with_name(prefix.name + "_best_cycle_residual_points.txt")
    )

    pdf_prefix = Path(args.pdf_prefix)
    png_prefix = Path(args.png_prefix)
    best_row = best_valid_row(scan_rows)

    write_scan_plot(
        scan_rows,
        pdf_prefix.with_name(pdf_prefix.name + "_scan_metrics.pdf"),
        png_prefix.with_name(png_prefix.name + "_scan_metrics.png"),
        x_limits=parse_range_text(scan_metadata.get("initial_period_range_ns")),
    )
    write_phase_profile_plot(
        phase_profile,
        pdf_prefix.with_name(pdf_prefix.name + "_best_phase_profile.pdf"),
        png_prefix.with_name(png_prefix.name + "_best_phase_profile.png"),
    )
    write_phase_tot_plot(
        phase_tot,
        pdf_prefix.with_name(pdf_prefix.name + "_best_phase_tot.pdf"),
        png_prefix.with_name(png_prefix.name + "_best_phase_tot.png"),
    )
    write_phase_tot_linear_plot(
        phase_tot,
        pdf_prefix.with_name(pdf_prefix.name + "_best_phase_tot_linearZ.pdf"),
        png_prefix.with_name(png_prefix.name + "_best_phase_tot_linearZ.png"),
    )
    write_cycle_residual_plot(
        cycle_index,
        counts,
        means,
        rms_values,
        pdf_prefix.with_name(pdf_prefix.name + "_best_cycle_residuals.pdf"),
        png_prefix.with_name(png_prefix.name + "_best_cycle_residuals.png"),
    )
    write_cycle_residual_detail_plot(
        point_cycle_index,
        point_residuals,
        cycle_index,
        counts,
        means,
        pdf_prefix.with_name(pdf_prefix.name + "_best_cycle_residual_detail.pdf"),
        png_prefix.with_name(png_prefix.name + "_best_cycle_residual_detail.png"),
    )
    write_global_time_validation_plot(
        point_global_time_seconds,
        point_residuals,
        pdf_prefix.with_name(pdf_prefix.name + "_best_global_time_residual_validation.pdf"),
        png_prefix.with_name(png_prefix.name + "_best_global_time_residual_validation.png"),
    )

    if args.root_output:
        write_root_outputs(Path(args.root_output), phase_profile, phase_tot)

    print(f"Best-by-merit period: {best_row['period_ns']:.6f} ns")
    print(f"Best-by-merit sigma: {best_row['sigma_ns']:.6f} ns")
    print(f"Best-by-merit drift slope: {best_row['drift_slope_ns_per_cycle']:.6f} ns/cycle")
    if args.root_output:
        print(f"ROOT output: {args.root_output}")


if __name__ == "__main__":
    main()