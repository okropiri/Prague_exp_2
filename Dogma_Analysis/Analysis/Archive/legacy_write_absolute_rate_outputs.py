from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import uproot
from uproot.writing.identify import to_TAxis, to_TH1x


mpl.rcParams["agg.path.chunksize"] = 10000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write ROOT, PDF, and PNG outputs for sparse absolute-time 1 us detector rates."
    )
    parser.add_argument("--rates-file", required=True)
    parser.add_argument("--summary-file", required=True)
    parser.add_argument("--root-output", required=True)
    parser.add_argument("--pdf-output", required=True)
    parser.add_argument("--png-output", required=True)
    return parser.parse_args()


def parse_metadata_line(line: str) -> tuple[str, str] | None:
    if not line.startswith("#") or "=" not in line:
        return None
    key, value = line[1:].split("=", 1)
    return key.strip(), value.strip()


def read_summary(path: Path) -> dict[str, str]:
    summary: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            summary[key.strip()] = value.strip()
    return summary


def read_sparse_rates(path: Path) -> tuple[dict[str, str], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    metadata: dict[str, str] = {}
    rows: list[tuple[int, int, int, int]] = []
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
            bin_text, ncal_text, l_text, s_text = line.split()
            rows.append((int(bin_text), int(ncal_text), int(l_text), int(s_text)))

    if rows:
        values = np.asarray(rows, dtype=np.int64)
        bin_indices = values[:, 0]
        ncal_counts = values[:, 1].astype(np.int32)
        lstilbene_counts = values[:, 2].astype(np.int32)
        sstilbene_counts = values[:, 3].astype(np.int32)
    else:
        bin_indices = np.asarray([], dtype=np.int64)
        ncal_counts = np.asarray([], dtype=np.int32)
        lstilbene_counts = np.asarray([], dtype=np.int32)
        sstilbene_counts = np.asarray([], dtype=np.int32)

    return metadata, bin_indices, ncal_counts, lstilbene_counts, sstilbene_counts


def choose_time_unit(max_seconds: float) -> tuple[float, str]:
    if max_seconds >= 3600.0:
        return 3600.0, "h"
    if max_seconds >= 120.0:
        return 60.0, "min"
    return 1.0, "s"


def determine_total_run_seconds(summary: dict[str, str], time_seconds: np.ndarray, bin_width_ns: float) -> float:
    candidates: list[float] = []
    run_start = float(summary.get("run_start_trigger_seconds", "0"))
    duplicate_previous = summary.get("duplicate_restart_previous_trigger_seconds")
    if duplicate_previous is not None:
        candidates.append(max(0.0, float(duplicate_previous) - run_start))
    last_hit = summary.get("last_hit_time_since_run_start_seconds")
    if last_hit is not None:
        candidates.append(max(0.0, float(last_hit)))
    if time_seconds.size > 0:
        candidates.append(float(time_seconds.max()) + bin_width_ns * 1e-9)
    return max(candidates, default=0.0)


def build_display_series(
    bin_indices: np.ndarray,
    counts: np.ndarray,
    total_run_seconds: float,
    bin_width_ns: float,
    max_display_points: int = 50000,
) -> tuple[np.ndarray, np.ndarray, float]:
    if total_run_seconds <= 0.0:
        return np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64), 1.0

    total_microsecond_bins = max(1, int(np.ceil(total_run_seconds * 1e6)))
    display_bin_width_us = max(1, int(np.ceil(total_microsecond_bins / max_display_points)))
    display_bin_count = int(np.ceil(total_microsecond_bins / display_bin_width_us))

    aggregated = np.zeros(display_bin_count, dtype=np.float64)
    if bin_indices.size > 0:
        display_indices = np.floor_divide(bin_indices, display_bin_width_us)
        valid = (display_indices >= 0) & (display_indices < display_bin_count)
        if np.any(valid):
            np.add.at(aggregated, display_indices[valid], counts[valid].astype(np.float64))

    aggregated /= float(display_bin_width_us)

    display_bin_width_seconds = display_bin_width_us * 1e-6
    time_axis = (np.arange(display_bin_count, dtype=np.float64) + 0.5) * display_bin_width_seconds
    return time_axis, aggregated, display_bin_width_seconds


def write_plot(
    bin_indices: np.ndarray,
    bin_width_ns: float,
    total_run_seconds: float,
    ncal_counts: np.ndarray,
    lstilbene_counts: np.ndarray,
    sstilbene_counts: np.ndarray,
    pdf_output: Path,
    png_output: Path,
) -> None:
    if bin_indices.size == 0:
        figure, axis = plt.subplots(figsize=(12, 6))
        axis.set_title("Absolute detector counts vs time")
        axis.set_xlabel("Time since run start (s)")
        axis.set_ylabel("counts / 1 us bin")
        axis.text(0.5, 0.5, "No valid detector hits", transform=axis.transAxes, ha="center", va="center")
        figure.tight_layout()
        figure.savefig(pdf_output)
        figure.savefig(png_output, dpi=180)
        plt.close(figure)
        return

    ncal_x, ncal_y, display_bin_width_seconds = build_display_series(
        bin_indices, ncal_counts, total_run_seconds, bin_width_ns
    )
    lstilbene_x, lstilbene_y, _ = build_display_series(
        bin_indices, lstilbene_counts, total_run_seconds, bin_width_ns
    )
    sstilbene_x, sstilbene_y, _ = build_display_series(
        bin_indices, sstilbene_counts, total_run_seconds, bin_width_ns
    )
    unit_scale, unit_label = choose_time_unit(total_run_seconds)

    figure, axis = plt.subplots(figsize=(13, 6))
    axis.plot(
        ncal_x / unit_scale,
        ncal_y,
        color="#1f77b4",
        linewidth=0.35,
        alpha=0.85,
        label="Ncal1",
        rasterized=True,
    )
    axis.plot(
        lstilbene_x / unit_scale,
        lstilbene_y,
        color="#d62728",
        linewidth=0.35,
        alpha=0.8,
        label="Lstilbene",
        rasterized=True,
    )
    axis.plot(
        sstilbene_x / unit_scale,
        sstilbene_y,
        color="#2ca02c",
        linewidth=0.45,
        alpha=0.9,
        label="Sstilbene",
        rasterized=True,
    )
    axis.set_xlabel(f"Time since run start ({unit_label})")
    axis.set_ylabel("average counts / 1 us bin")
    axis.set_title(
        f"Absolute detector rates over full run (display averaged over {display_bin_width_seconds:.6f} s bins)"
    )
    axis.set_xlim(0.0, max(total_run_seconds / unit_scale, 0.0))
    axis.set_ylim(bottom=0.0)
    axis.grid(True, linewidth=0.3, alpha=0.4)
    axis.legend(loc="upper right")
    figure.tight_layout()
    figure.savefig(pdf_output)
    figure.savefig(png_output, dpi=180)
    plt.close(figure)


def build_root_display_series(
    bin_indices: np.ndarray,
    bin_width_ns: float,
    total_run_seconds: float,
    ncal_counts: np.ndarray,
    lstilbene_counts: np.ndarray,
    sstilbene_counts: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    display_time_s, ncal_display, display_bin_width_seconds = build_display_series(
        bin_indices, ncal_counts, total_run_seconds, bin_width_ns
    )
    _, lstilbene_display, _ = build_display_series(
        bin_indices, lstilbene_counts, total_run_seconds, bin_width_ns
    )
    _, sstilbene_display, _ = build_display_series(
        bin_indices, sstilbene_counts, total_run_seconds, bin_width_ns
    )
    return display_time_s, ncal_display, lstilbene_display, sstilbene_display, display_bin_width_seconds


def make_rate_histogram(title: str, values: np.ndarray, total_run_seconds: float, entries: float):
    bins = int(values.size)
    storage = np.zeros(bins + 2, dtype=np.float64)
    if bins > 0:
        storage[1:-1] = values.astype(np.float64)

    if bins > 0 and total_run_seconds > 0.0:
        bin_width = total_run_seconds / bins
        centers = (np.arange(bins, dtype=np.float64) + 0.5) * bin_width
    else:
        centers = np.asarray([], dtype=np.float64)

    sumw = float(values.sum(dtype=np.float64))
    sumw2 = float(np.square(values, dtype=np.float64).sum(dtype=np.float64))
    sumwx = float((values.astype(np.float64) * centers).sum(dtype=np.float64)) if centers.size else 0.0
    sumwx2 = float((values.astype(np.float64) * np.square(centers)).sum(dtype=np.float64)) if centers.size else 0.0

    return to_TH1x(
        fName=None,
        fTitle=title,
        data=storage,
        fEntries=float(entries),
        fTsumw=sumw,
        fTsumw2=sumw2,
        fTsumwx=sumwx,
        fTsumwx2=sumwx2,
        fSumw2=None,
        fXaxis=to_TAxis("xaxis", "Time since run start (s)", bins, 0.0, float(total_run_seconds) if total_run_seconds > 0.0 else 1.0),
    )


def main() -> None:
    args = parse_args()
    rates_file = Path(args.rates_file)
    summary = read_summary(Path(args.summary_file))
    metadata, bin_indices, ncal_counts, lstilbene_counts, sstilbene_counts = read_sparse_rates(rates_file)

    bin_width_ns = float(metadata["bin_width_ns"])
    time_seconds = bin_indices.astype(np.float64) * (bin_width_ns * 1e-9)
    total_run_seconds = determine_total_run_seconds(summary, time_seconds, bin_width_ns)
    display_time_s, ncal_display, lstilbene_display, sstilbene_display, display_bin_width_seconds = build_root_display_series(
        bin_indices=bin_indices,
        bin_width_ns=bin_width_ns,
        total_run_seconds=total_run_seconds,
        ncal_counts=ncal_counts,
        lstilbene_counts=lstilbene_counts,
        sstilbene_counts=sstilbene_counts,
    )

    root_output = Path(args.root_output)
    with uproot.recreate(root_output) as root_file:
        root_file["Rates/hNcal1_display"] = make_rate_histogram(
            f"Ncal1 average counts per 1 us bin (display averaged over {display_bin_width_seconds:.6f} s bins)",
            ncal_display,
            total_run_seconds,
            float(ncal_counts.sum(dtype=np.int64)),
        )
        root_file["Rates/hLstilbene_display"] = make_rate_histogram(
            f"Lstilbene average counts per 1 us bin (display averaged over {display_bin_width_seconds:.6f} s bins)",
            lstilbene_display,
            total_run_seconds,
            float(lstilbene_counts.sum(dtype=np.int64)),
        )
        root_file["Rates/hSstilbene_display"] = make_rate_histogram(
            f"Sstilbene average counts per 1 us bin (display averaged over {display_bin_width_seconds:.6f} s bins)",
            sstilbene_display,
            total_run_seconds,
            float(sstilbene_counts.sum(dtype=np.int64)),
        )
        root_file["Rates/display_series"] = {
            "time_since_run_start_s": display_time_s.astype(np.float64),
            "ncal1_avg_count_per_1us": ncal_display.astype(np.float32),
            "lstilbene_avg_count_per_1us": lstilbene_display.astype(np.float32),
            "sstilbene_avg_count_per_1us": sstilbene_display.astype(np.float32),
        }
        root_file["Rates/summary"] = {
            "bin_width_ns": np.asarray([bin_width_ns], dtype=np.float64),
            "total_run_duration_s": np.asarray([total_run_seconds], dtype=np.float64),
            "display_bin_width_s": np.asarray([display_bin_width_seconds], dtype=np.float64),
            "min_bin_index": np.asarray([bin_indices.min(initial=0)], dtype=np.int64),
            "max_bin_index": np.asarray([bin_indices.max(initial=0)], dtype=np.int64),
            "nonzero_bins": np.asarray([bin_indices.size], dtype=np.int64),
            "total_ncal1_hits": np.asarray([int(ncal_counts.sum(dtype=np.int64))], dtype=np.int64),
            "total_lstilbene_hits": np.asarray([int(lstilbene_counts.sum(dtype=np.int64))], dtype=np.int64),
            "total_sstilbene_hits": np.asarray([int(sstilbene_counts.sum(dtype=np.int64))], dtype=np.int64),
        }

    write_plot(
        bin_indices=bin_indices,
        bin_width_ns=bin_width_ns,
        total_run_seconds=total_run_seconds,
        ncal_counts=ncal_counts,
        lstilbene_counts=lstilbene_counts,
        sstilbene_counts=sstilbene_counts,
        pdf_output=Path(args.pdf_output),
        png_output=Path(args.png_output),
    )

    print(f"ROOT output: {root_output}")
    print(f"Display points: {display_time_s.size}")
    print(f"Sparse bins: {bin_indices.size}")
    print(f"Ncal1 hits: {int(ncal_counts.sum(dtype=np.int64))}")
    print(f"Lstilbene hits: {int(lstilbene_counts.sum(dtype=np.int64))}")
    print(f"Sstilbene hits: {int(sstilbene_counts.sum(dtype=np.int64))}")


if __name__ == "__main__":
    main()