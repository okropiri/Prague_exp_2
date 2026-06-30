from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot
from uproot.writing.identify import to_TAxis, to_TH1x


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write ROOT, PDF, and PNG outputs for trigger-window asymmetry investigation."
    )
    parser.add_argument("--input-prefix", required=True)
    parser.add_argument("--root-output", required=True)
    parser.add_argument("--pdf-output", required=True)
    parser.add_argument("--png-output", required=True)
    return parser.parse_args()


def parse_metadata_line(line: str) -> tuple[str, str] | None:
    if not line.startswith("#") or "=" not in line:
        return None
    key, value = line[1:].split("=", 1)
    return key.strip(), value.strip()


def read_blocks(path: Path) -> dict[str, object]:
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

    if not rows:
        raise ValueError(f"No block rows found in {path}")

    data = np.asarray(rows, dtype=np.int64)
    metadata["block_size_windows"] = int(metadata["block_size_windows"])
    metadata["window_min_ns"] = float(metadata["window_min_ns"])
    metadata["window_max_ns"] = float(metadata["window_max_ns"])
    metadata["data"] = data
    return metadata


def ratio(after: np.ndarray, total: np.ndarray) -> np.ndarray:
    values = np.full(after.shape, np.nan, dtype=np.float64)
    mask = total > 0
    values[mask] = after[mask] / total[mask]
    return values


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


def make_count_hist(values: np.ndarray, x_min: float, x_max: float, title: str, x_title: str):
    data = np.zeros(values.size + 2, dtype=np.float64)
    data[1:-1] = values.astype(np.float64)
    centers = x_min + (np.arange(values.size, dtype=np.float64) + 0.5) * ((x_max - x_min) / values.size)
    return to_TH1x(
        fName=None,
        fTitle=title,
        data=data,
        fEntries=float(values.size),
        fTsumw=float(values.sum()),
        fTsumw2=float(np.square(values, dtype=np.float64).sum()),
        fTsumwx=float((values * centers).sum()),
        fTsumwx2=float((values * np.square(centers)).sum()),
        fSumw2=None,
        fXaxis=to_TAxis("xaxis", x_title, values.size, x_min, x_max),
    )


def write_plot(x_values: np.ndarray,
               channel_series: list[tuple[np.ndarray, str, str]],
               lstilbene_mult_series: list[tuple[np.ndarray, str, str]],
               pdf_path: Path,
               png_path: Path,
               block_size_windows: int) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True)

    for values, label, color in channel_series:
        axes[0].plot(x_values, values, linewidth=0.9, color=color, label=label)
    axes[0].axhline(0.5, color="#444444", linewidth=0.8, linestyle="--")
    axes[0].set_ylabel("after-trigger fraction")
    axes[0].set_ylim(0.0, 1.0)
    axes[0].set_title(f"Asymmetry coefficient vs trigger window block (block size = {block_size_windows} windows)")
    axes[0].grid(True, linewidth=0.3, alpha=0.4)
    axes[0].legend(loc="best")

    for values, label, color in lstilbene_mult_series:
        axes[1].plot(x_values, values, linewidth=0.9, color=color, label=label)
    axes[1].axhline(0.5, color="#444444", linewidth=0.8, linestyle="--")
    axes[1].set_xlabel("trigger window number")
    axes[1].set_ylabel("after-trigger fraction")
    axes[1].set_ylim(0.0, 1.0)
    axes[1].set_title("Lstilbene split by per-window multiplicity")
    axes[1].grid(True, linewidth=0.3, alpha=0.4)
    axes[1].legend(loc="best")

    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    prefix = Path(args.input_prefix)
    blocks = read_blocks(prefix.with_name(prefix.name + "_blocks.txt"))
    data = np.asarray(blocks["data"], dtype=np.int64)
    block_size_windows = int(blocks["block_size_windows"])

    window_start = data[:, 1].astype(np.float64)
    window_end = data[:, 2].astype(np.float64)
    x_values = 0.5 * (window_start + window_end)
    x_min = float(window_start[0])
    x_max = float(window_end[-1])

    ncal1_values = ratio(data[:, 4].astype(np.float64), data[:, 5].astype(np.float64))
    lstilbene_values = ratio(data[:, 6].astype(np.float64), data[:, 7].astype(np.float64))
    sstilbene_values = ratio(data[:, 8].astype(np.float64), data[:, 9].astype(np.float64))
    lst_mult1_values = ratio(data[:, 10].astype(np.float64), data[:, 11].astype(np.float64))
    lst_mult2to3_values = ratio(data[:, 12].astype(np.float64), data[:, 13].astype(np.float64))
    lst_mult4to7_values = ratio(data[:, 14].astype(np.float64), data[:, 15].astype(np.float64))
    lst_mult8plus_values = ratio(data[:, 16].astype(np.float64), data[:, 17].astype(np.float64))

    channel_series = [
        (ncal1_values, "Ncal1", "#1f77b4"),
        (lstilbene_values, "Lstilbene", "#d62728"),
        (sstilbene_values, "Sstilbene", "#2ca02c"),
    ]
    lstilbene_mult_series = [
        (lstilbene_values, "Lstilbene overall", "#111111"),
        (lst_mult1_values, "Lstilbene mult 1", "#9467bd"),
        (lst_mult2to3_values, "Lstilbene mult 2-3", "#ff7f0e"),
        (lst_mult4to7_values, "Lstilbene mult 4-7", "#17becf"),
        (lst_mult8plus_values, "Lstilbene mult 8+", "#8c564b"),
    ]

    root_output = Path(args.root_output)
    with uproot.recreate(root_output) as root_file:
        root_file["Asymmetry/hNcal1AfterFraction"] = make_series_hist(
            ncal1_values, x_min, x_max, "Ncal1 after-trigger fraction", "Trigger window number"
        )
        root_file["Asymmetry/hLstilbeneAfterFraction"] = make_series_hist(
            lstilbene_values, x_min, x_max, "Lstilbene after-trigger fraction", "Trigger window number"
        )
        root_file["Asymmetry/hSstilbeneAfterFraction"] = make_series_hist(
            sstilbene_values, x_min, x_max, "Sstilbene after-trigger fraction", "Trigger window number"
        )
        root_file["Asymmetry/hLstilbeneMult1AfterFraction"] = make_series_hist(
            lst_mult1_values, x_min, x_max, "Lstilbene multiplicity 1 after-trigger fraction", "Trigger window number"
        )
        root_file["Asymmetry/hLstilbeneMult2to3AfterFraction"] = make_series_hist(
            lst_mult2to3_values, x_min, x_max, "Lstilbene multiplicity 2-3 after-trigger fraction", "Trigger window number"
        )
        root_file["Asymmetry/hLstilbeneMult4to7AfterFraction"] = make_series_hist(
            lst_mult4to7_values, x_min, x_max, "Lstilbene multiplicity 4-7 after-trigger fraction", "Trigger window number"
        )
        root_file["Asymmetry/hLstilbeneMult8plusAfterFraction"] = make_series_hist(
            lst_mult8plus_values, x_min, x_max, "Lstilbene multiplicity 8+ after-trigger fraction", "Trigger window number"
        )
        root_file["Counts/hNcal1TotalCounts"] = make_count_hist(
            data[:, 5].astype(np.float64), x_min, x_max, "Ncal1 total counts per block", "Trigger window number"
        )
        root_file["Counts/hLstilbeneTotalCounts"] = make_count_hist(
            data[:, 7].astype(np.float64), x_min, x_max, "Lstilbene total counts per block", "Trigger window number"
        )
        root_file["Counts/hSstilbeneTotalCounts"] = make_count_hist(
            data[:, 9].astype(np.float64), x_min, x_max, "Sstilbene total counts per block", "Trigger window number"
        )

    write_plot(x_values, channel_series, lstilbene_mult_series, Path(args.pdf_output), Path(args.png_output), block_size_windows)

    print(f"ROOT output: {root_output}")
    print(f"Blocks: {data.shape[0]}")
    print(f"Ncal1 overall after fraction: {np.nansum(data[:, 4]) / np.nansum(data[:, 5]):.6f}")
    print(f"Lstilbene overall after fraction: {np.nansum(data[:, 6]) / np.nansum(data[:, 7]):.6f}")
    print(f"Sstilbene overall after fraction: {np.nansum(data[:, 8]) / np.nansum(data[:, 9]):.6f}")


if __name__ == "__main__":
    main()