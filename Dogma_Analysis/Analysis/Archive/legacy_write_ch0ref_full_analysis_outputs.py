from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["agg.path.chunksize"] = 10000

import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.colors import LogNorm
from uproot.writing.identify import to_TAxis, to_TH1x, to_TH2x


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write ROOT, PDF, and PNG outputs for ch0-referenced DOGMA analysis products."
    )
    parser.add_argument("--input-prefix", required=True)
    parser.add_argument("--root-output", required=True)
    parser.add_argument("--pdf-prefix", required=True)
    parser.add_argument("--png-prefix", required=True)
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
    metadata["entries"] = float(metadata["entries"])
    metadata["sumw"] = float(metadata["sumw"])
    metadata["sumw2"] = float(metadata["sumw2"])
    metadata["sumwx"] = float(metadata["sumwx"])
    metadata["sumwx2"] = float(metadata["sumwx2"])
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
    metadata["entries"] = float(metadata["entries"])
    metadata["sumw"] = float(metadata["sumw"])
    metadata["sumw2"] = float(metadata["sumw2"])
    metadata["sumwx"] = float(metadata["sumwx"])
    metadata["sumwx2"] = float(metadata["sumwx2"])
    metadata["sumwy"] = float(metadata["sumwy"])
    metadata["sumwy2"] = float(metadata["sumwy2"])
    metadata["sumwxy"] = float(metadata["sumwxy"])
    metadata["storage"] = storage
    return metadata


def make_hist1_root_object(hist: dict[str, object]):
    bins = int(hist["bins"])
    min_value = float(hist["min"])
    max_value = float(hist["max"])
    storage = np.asarray(hist["storage"], dtype=np.float64)
    return to_TH1x(
        fName=None,
        fTitle=str(hist["title"]),
        data=storage,
        fEntries=float(hist["entries"]),
        fTsumw=float(hist["sumw"]),
        fTsumw2=float(hist["sumw2"]),
        fTsumwx=float(hist["sumwx"]),
        fTsumwx2=float(hist["sumwx2"]),
        fSumw2=None,
        fXaxis=to_TAxis("xaxis", str(hist["x_title"]), bins, min_value, max_value),
    )


def make_hist2_root_object(hist: dict[str, object]):
    x_bins = int(hist["x_bins"])
    x_min = float(hist["x_min"])
    x_max = float(hist["x_max"])
    y_bins = int(hist["y_bins"])
    y_min = float(hist["y_min"])
    y_max = float(hist["y_max"])
    storage = np.asarray(hist["storage"], dtype=np.float64).reshape(-1)
    return to_TH2x(
        fName=None,
        fTitle=str(hist["title"]),
        data=storage,
        fEntries=float(hist["entries"]),
        fTsumw=float(hist["sumw"]),
        fTsumw2=float(hist["sumw2"]),
        fTsumwx=float(hist["sumwx"]),
        fTsumwx2=float(hist["sumwx2"]),
        fTsumwy=float(hist["sumwy"]),
        fTsumwy2=float(hist["sumwy2"]),
        fTsumwxy=float(hist["sumwxy"]),
        fSumw2=None,
        fXaxis=to_TAxis("xaxis", str(hist["x_title"]), x_bins, x_min, x_max),
        fYaxis=to_TAxis("yaxis", str(hist["y_title"]), y_bins, y_min, y_max),
    )


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


def read_sparse_rates(path: Path) -> tuple[dict[str, str], np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    metadata: dict[str, str] = {}
    rows: list[tuple[int, int, int, int, int]] = []
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
            bin_text, ch0_text, n_text, l_text, s_text = line.split()
            rows.append((int(bin_text), int(ch0_text), int(n_text), int(l_text), int(s_text)))

    if rows:
        values = np.asarray(rows, dtype=np.int64)
        bin_indices = values[:, 0]
        ch0_counts = values[:, 1].astype(np.int32)
        ncal_counts = values[:, 2].astype(np.int32)
        lstilbene_counts = values[:, 3].astype(np.int32)
        sstilbene_counts = values[:, 4].astype(np.int32)
    else:
        bin_indices = np.asarray([], dtype=np.int64)
        ch0_counts = np.asarray([], dtype=np.int32)
        ncal_counts = np.asarray([], dtype=np.int32)
        lstilbene_counts = np.asarray([], dtype=np.int32)
        sstilbene_counts = np.asarray([], dtype=np.int32)
    return metadata, bin_indices, ch0_counts, ncal_counts, lstilbene_counts, sstilbene_counts


def read_asymmetry_blocks(path: Path) -> dict[str, object]:
    metadata: dict[str, object] = {}
    rows: list[list[int]] = []
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
            rows.append([int(value) for value in line.split()])

    metadata["block_size_windows"] = int(metadata["block_size_windows"])
    metadata["folded_min_ns"] = float(metadata["folded_min_ns"])
    metadata["folded_max_ns"] = float(metadata["folded_max_ns"])
    metadata["data"] = np.asarray(rows, dtype=np.int64) if rows else np.zeros((0, 11), dtype=np.int64)
    return metadata


def ratio(after: np.ndarray, total: np.ndarray) -> np.ndarray:
    values = np.full(after.shape, np.nan, dtype=np.float64)
    mask = total > 0
    values[mask] = after[mask] / total[mask]
    return values


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


def make_series_hist(values: np.ndarray, x_min: float, x_max: float, title: str, x_title: str):
    finite_mask = np.isfinite(values)
    data = np.zeros(values.size + 2, dtype=np.float64)
    data[1:-1] = np.where(finite_mask, values, 0.0)
    entries = float(np.count_nonzero(finite_mask))
    finite_values = values[finite_mask]
    if finite_values.size == 0:
        sumw = sumw2 = sumwx = sumwx2 = 0.0
    else:
        centers = x_min + (np.arange(values.size, dtype=np.float64) + 0.5) * ((x_max - x_min) / values.size)
        finite_centers = centers[finite_mask]
        sumw = float(finite_values.sum())
        sumw2 = float(np.square(finite_values).sum())
        sumwx = float((finite_values * finite_centers).sum())
        sumwx2 = float((finite_values * np.square(finite_centers)).sum())
    return to_TH1x(
        fName=None,
        fTitle=title,
        data=data,
        fEntries=entries,
        fTsumw=sumw,
        fTsumw2=sumw2,
        fTsumwx=sumwx,
        fTsumwx2=sumwx2,
        fSumw2=None,
        fXaxis=to_TAxis("xaxis", x_title, values.size, x_min, x_max),
    )


def write_folded_rates_plot(histograms: list[tuple[dict[str, object], str]], pdf_path: Path, png_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    for hist, color in histograms:
        storage = np.asarray(hist["storage"], dtype=np.float64)
        values = storage[1:-1]
        bins = int(hist["bins"])
        x_min = float(hist["min"])
        x_max = float(hist["max"])
        centers = x_min + (np.arange(bins, dtype=np.float64) + 0.5) * ((x_max - x_min) / bins)
        ax.plot(centers, values, linewidth=0.9, color=color, label=str(hist["object_name"]))
    ax.set_xlabel(str(histograms[0][0]["x_title"]))
    ax.set_ylabel("raw rising-edge counts / bin")
    ax.set_title("Folded rates vs ch0-relative time")
    ax.grid(True, linewidth=0.3, alpha=0.4)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=180)
    plt.close(fig)


def write_tot_plot(hist: dict[str, object], pdf_path: Path, png_path: Path) -> None:
    storage = np.asarray(hist["storage"], dtype=np.float64)
    counts = storage[1:-1, 1:-1]
    nonzero = counts[counts > 0]
    norm = None
    if nonzero.size > 0:
        vmax = float(np.quantile(nonzero, 0.999))
        vmax = max(vmax, 1.0)
        norm = LogNorm(vmin=1.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(10, 7))
    image = ax.imshow(
        counts,
        origin="lower",
        aspect="auto",
        extent=[float(hist["x_min"]), float(hist["x_max"]), float(hist["y_min"]), float(hist["y_max"])],
        cmap="jet",
        norm=norm,
    )
    ax.set_xlabel(str(hist["x_title"]))
    ax.set_ylabel(str(hist["y_title"]))
    ax.set_title(str(hist["object_name"]))
    fig.colorbar(image, ax=ax, label="counts / bin")
    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=180)
    plt.close(fig)


def write_asymmetry_plot(x_values: np.ndarray, series: list[tuple[np.ndarray, str, str]], pdf_path: Path, png_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 6))
    for values, label, color in series:
        ax.plot(x_values, values, linewidth=0.9, color=color, label=label)
    ax.axhline(0.5, color="#444444", linewidth=0.8, linestyle="--")
    ax.set_xlabel("ch0-referenced window number")
    ax.set_ylabel("after-ch0 fraction")
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Asymmetry coefficient vs ch0-referenced window block")
    ax.grid(True, linewidth=0.3, alpha=0.4)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=180)
    plt.close(fig)


def write_full_time_rates_plot(x_series: np.ndarray,
                               ch0_values: np.ndarray,
                               ncal_values: np.ndarray,
                               lstilbene_values: np.ndarray,
                               sstilbene_values: np.ndarray,
                               total_run_seconds: float,
                               display_bin_width_seconds: float,
                               pdf_path: Path,
                               png_path: Path) -> None:
    unit_scale, unit_label = choose_time_unit(total_run_seconds)
    fig, ax = plt.subplots(figsize=(13, 6))
    ax.plot(x_series / unit_scale, ch0_values, color="#111111", linewidth=0.35, alpha=0.85, label="ch0", rasterized=True)
    ax.plot(x_series / unit_scale, ncal_values, color="#1f77b4", linewidth=0.35, alpha=0.85, label="Ncal1", rasterized=True)
    ax.plot(x_series / unit_scale, lstilbene_values, color="#d62728", linewidth=0.35, alpha=0.8, label="Lstilbene", rasterized=True)
    ax.plot(x_series / unit_scale, sstilbene_values, color="#2ca02c", linewidth=0.45, alpha=0.9, label="Sstilbene", rasterized=True)
    ax.set_xlabel(f"Time since run start ({unit_label})")
    ax.set_ylabel("average counts / 1 us bin")
    ax.set_title(f"Full-time rates over run (display averaged over {display_bin_width_seconds:.6f} s bins)")
    ax.set_xlim(0.0, max(total_run_seconds / unit_scale, 0.0))
    ax.set_ylim(bottom=0.0)
    ax.grid(True, linewidth=0.3, alpha=0.4)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=180)
    plt.close(fig)


def write_ch0_rate_plot(x_series: np.ndarray,
                        ch0_values: np.ndarray,
                        total_run_seconds: float,
                        display_bin_width_seconds: float,
                        pdf_path: Path,
                        png_path: Path) -> None:
    unit_scale, unit_label = choose_time_unit(total_run_seconds)
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.plot(x_series / unit_scale, ch0_values, color="#111111", linewidth=0.4, alpha=0.9, label="ch0", rasterized=True)
    ax.set_xlabel(f"Time since run start ({unit_label})")
    ax.set_ylabel("average counts / 1 us bin")
    ax.set_title(f"Full-time ch0 rate over run (display averaged over {display_bin_width_seconds:.6f} s bins)")
    ax.set_xlim(0.0, max(total_run_seconds / unit_scale, 0.0))
    ax.set_ylim(bottom=0.0)
    ax.grid(True, linewidth=0.3, alpha=0.4)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    prefix = Path(args.input_prefix)

    hist_ncal1_folded = read_hist1(prefix.with_name(prefix.name + "_custom_RawRiseCountsNcal1FoldedCh0Ref_hist.txt"))
    hist_lstilbene_folded = read_hist1(prefix.with_name(prefix.name + "_custom_RawRiseCountsLstilbeneFoldedCh0Ref_hist.txt"))
    hist_sstilbene_folded = read_hist1(prefix.with_name(prefix.name + "_custom_RawRiseCountsSstilbeneFoldedCh0Ref_hist.txt"))
    hist_ncal1_tot = read_hist2(prefix.with_name(prefix.name + "_custom_ProfNcal1Ch0Ref_hist.txt"))
    hist_lstilbene_tot = read_hist2(prefix.with_name(prefix.name + "_custom_ProfLstilbeneCh0Ref_hist.txt"))
    hist_sstilbene_tot = read_hist2(prefix.with_name(prefix.name + "_custom_ProfSstilbeneCh0Ref_hist.txt"))
    asym_blocks = read_asymmetry_blocks(prefix.with_name(prefix.name + "_asymmetry_blocks.txt"))
    summary = read_summary(prefix.with_name(prefix.name + "_summary.txt"))
    rate_meta, bin_indices, ch0_counts, ncal_counts, lstilbene_counts, sstilbene_counts = read_sparse_rates(
        prefix.with_name(prefix.name + "_rates.txt")
    )

    asym_data = np.asarray(asym_blocks["data"], dtype=np.int64)
    if asym_data.size > 0:
        window_start = asym_data[:, 1].astype(np.float64)
        window_end = asym_data[:, 2].astype(np.float64)
        asym_x = 0.5 * (window_start + window_end)
        asym_x_min = float(window_start[0])
        asym_x_max = float(window_end[-1])
        ncal_asym = ratio(asym_data[:, 5].astype(np.float64), asym_data[:, 6].astype(np.float64))
        lstilbene_asym = ratio(asym_data[:, 7].astype(np.float64), asym_data[:, 8].astype(np.float64))
        sstilbene_asym = ratio(asym_data[:, 9].astype(np.float64), asym_data[:, 10].astype(np.float64))
    else:
        asym_x = np.asarray([], dtype=np.float64)
        asym_x_min = 0.0
        asym_x_max = 1.0
        ncal_asym = np.asarray([], dtype=np.float64)
        lstilbene_asym = np.asarray([], dtype=np.float64)
        sstilbene_asym = np.asarray([], dtype=np.float64)

    rate_bin_width_ns = float(rate_meta["rate_bin_width_ns"])
    time_seconds = bin_indices.astype(np.float64) * (rate_bin_width_ns * 1e-9)
    total_run_seconds = determine_total_run_seconds(summary, time_seconds, rate_bin_width_ns)
    ch0_x, ch0_display, display_bin_width_seconds = build_display_series(bin_indices, ch0_counts, total_run_seconds, rate_bin_width_ns)
    _, ncal_display, _ = build_display_series(bin_indices, ncal_counts, total_run_seconds, rate_bin_width_ns)
    _, lstilbene_display, _ = build_display_series(bin_indices, lstilbene_counts, total_run_seconds, rate_bin_width_ns)
    _, sstilbene_display, _ = build_display_series(bin_indices, sstilbene_counts, total_run_seconds, rate_bin_width_ns)

    root_output = Path(args.root_output)
    with uproot.recreate(root_output) as root_file:
        root_file["FoldedRates/hNcal1FoldedRate_ch0Ref"] = make_hist1_root_object(hist_ncal1_folded)
        root_file["FoldedRates/hLstilbeneFoldedRate_ch0Ref"] = make_hist1_root_object(hist_lstilbene_folded)
        root_file["FoldedRates/hSstilbeneFoldedRate_ch0Ref"] = make_hist1_root_object(hist_sstilbene_folded)
        root_file["FoldedTOT/hNcal1TimeTot_ch0Ref"] = make_hist2_root_object(hist_ncal1_tot)
        root_file["FoldedTOT/hLstilbeneTimeTot_ch0Ref"] = make_hist2_root_object(hist_lstilbene_tot)
        root_file["FoldedTOT/hSstilbeneTimeTot_ch0Ref"] = make_hist2_root_object(hist_sstilbene_tot)
        root_file["FoldedAsymmetry/hNcal1AfterFraction_ch0Ref"] = make_series_hist(
            ncal_asym, asym_x_min, asym_x_max, "Ncal1 after-ch0 fraction", "ch0-referenced window number"
        )
        root_file["FoldedAsymmetry/hLstilbeneAfterFraction_ch0Ref"] = make_series_hist(
            lstilbene_asym, asym_x_min, asym_x_max, "Lstilbene after-ch0 fraction", "ch0-referenced window number"
        )
        root_file["FoldedAsymmetry/hSstilbeneAfterFraction_ch0Ref"] = make_series_hist(
            sstilbene_asym, asym_x_min, asym_x_max, "Sstilbene after-ch0 fraction", "ch0-referenced window number"
        )
        root_file["FullTimeRates/hCh0_display"] = make_rate_histogram(
            f"ch0 average counts per 1 us bin (display averaged over {display_bin_width_seconds:.6f} s bins)",
            ch0_display,
            total_run_seconds,
            float(ch0_counts.sum(dtype=np.int64)),
        )
        root_file["FullTimeRates/hNcal1_display"] = make_rate_histogram(
            f"Ncal1 average counts per 1 us bin (display averaged over {display_bin_width_seconds:.6f} s bins)",
            ncal_display,
            total_run_seconds,
            float(ncal_counts.sum(dtype=np.int64)),
        )
        root_file["FullTimeRates/hLstilbene_display"] = make_rate_histogram(
            f"Lstilbene average counts per 1 us bin (display averaged over {display_bin_width_seconds:.6f} s bins)",
            lstilbene_display,
            total_run_seconds,
            float(lstilbene_counts.sum(dtype=np.int64)),
        )
        root_file["FullTimeRates/hSstilbene_display"] = make_rate_histogram(
            f"Sstilbene average counts per 1 us bin (display averaged over {display_bin_width_seconds:.6f} s bins)",
            sstilbene_display,
            total_run_seconds,
            float(sstilbene_counts.sum(dtype=np.int64)),
        )

    pdf_prefix = Path(args.pdf_prefix)
    png_prefix = Path(args.png_prefix)
    write_folded_rates_plot(
        [(hist_ncal1_folded, "#1f77b4"), (hist_lstilbene_folded, "#d62728"), (hist_sstilbene_folded, "#2ca02c")],
        pdf_prefix.with_name(pdf_prefix.name + "_FoldedRates_T.pdf"),
        png_prefix.with_name(png_prefix.name + "_FoldedRates_T.png"),
    )
    write_tot_plot(
        hist_ncal1_tot,
        pdf_prefix.with_name(pdf_prefix.name + "_Ncal1_ToT_T.pdf"),
        png_prefix.with_name(png_prefix.name + "_Ncal1_ToT_T.png"),
    )
    write_tot_plot(
        hist_lstilbene_tot,
        pdf_prefix.with_name(pdf_prefix.name + "_Lstilbene_ToT_T.pdf"),
        png_prefix.with_name(png_prefix.name + "_Lstilbene_ToT_T.png"),
    )
    write_tot_plot(
        hist_sstilbene_tot,
        pdf_prefix.with_name(pdf_prefix.name + "_Sstilbene_ToT_T.pdf"),
        png_prefix.with_name(png_prefix.name + "_Sstilbene_ToT_T.png"),
    )
    write_asymmetry_plot(
        asym_x,
        [(ncal_asym, "Ncal1", "#1f77b4"), (lstilbene_asym, "Lstilbene", "#d62728"), (sstilbene_asym, "Sstilbene", "#2ca02c")],
        pdf_prefix.with_name(pdf_prefix.name + "_FoldedAsymmetry.pdf"),
        png_prefix.with_name(png_prefix.name + "_FoldedAsymmetry.png"),
    )
    write_full_time_rates_plot(
        ch0_x,
        ch0_display,
        ncal_display,
        lstilbene_display,
        sstilbene_display,
        total_run_seconds,
        display_bin_width_seconds,
        pdf_prefix.with_name(pdf_prefix.name + "_FullTimeRates.pdf"),
        png_prefix.with_name(png_prefix.name + "_FullTimeRates.png"),
    )
    write_ch0_rate_plot(
        ch0_x,
        ch0_display,
        total_run_seconds,
        display_bin_width_seconds,
        pdf_prefix.with_name(pdf_prefix.name + "_Ch0FullTimeRate.pdf"),
        png_prefix.with_name(png_prefix.name + "_Ch0FullTimeRate.png"),
    )

    print(f"ROOT output: {root_output}")
    print(f"Ncal1 folded entries: {hist_ncal1_folded['entries']:.0f}")
    print(f"Lstilbene folded entries: {hist_lstilbene_folded['entries']:.0f}")
    print(f"Sstilbene folded entries: {hist_sstilbene_folded['entries']:.0f}")


if __name__ == "__main__":
    main()