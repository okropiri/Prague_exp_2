from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot
from uproot.writing.identify import to_TAxis, to_TH1x


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write ROOT, PDF, and PNG outputs for Lstilbene rise/fall and rise-spacing checks."
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


def hist_centers(hist: dict[str, object]) -> np.ndarray:
    bins = int(hist["bins"])
    x_min = float(hist["min"])
    x_max = float(hist["max"])
    return x_min + (np.arange(bins, dtype=np.float64) + 0.5) * ((x_max - x_min) / bins)


def hist_values(hist: dict[str, object]) -> np.ndarray:
    return np.asarray(hist["storage"], dtype=np.float64)[1:-1]


def write_plot(hist_rise: dict[str, object], hist_fall: dict[str, object], hist_spacing: dict[str, object], pdf_path: Path, png_path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=False)

    centers_rise = hist_centers(hist_rise)
    centers_fall = hist_centers(hist_fall)
    axes[0].plot(centers_rise, hist_values(hist_rise), linewidth=0.9, color="#d62728", label="Lstilbene rises")
    axes[0].plot(centers_fall, hist_values(hist_fall), linewidth=0.9, color="#1f77b4", label="Lstilbene falls")
    axes[0].set_xlabel(str(hist_rise["x_title"]))
    axes[0].set_ylabel("edge counts / bin")
    axes[0].set_title("Lstilbene folded counts: rising vs falling edges")
    axes[0].grid(True, linewidth=0.3, alpha=0.4)
    axes[0].legend(loc="upper left")

    spacing_centers = hist_centers(hist_spacing)
    spacing_values = hist_values(hist_spacing)
    axes[1].step(spacing_centers, spacing_values, where="mid", linewidth=1.0, color="#2ca02c")
    axes[1].set_xlabel(str(hist_spacing["x_title"]))
    axes[1].set_ylabel("successive-rise pairs / bin")
    axes[1].set_title("Lstilbene rise-to-rise spacing within trigger windows")
    axes[1].set_yscale("log")
    axes[1].grid(True, linewidth=0.3, alpha=0.4)

    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    prefix = Path(args.input_prefix)

    hist_rise = read_hist1(prefix.with_name(prefix.name + "_custom_RawRiseCountsLstilbeneFolded_hist.txt"))
    hist_fall = read_hist1(prefix.with_name(prefix.name + "_custom_RawFallCountsLstilbeneFolded_hist.txt"))
    hist_spacing = read_hist1(prefix.with_name(prefix.name + "_custom_RiseToRiseSpacingLstilbene_hist.txt"))

    root_output = Path(args.root_output)
    with uproot.recreate(root_output) as root_file:
        root_file["Histograms/custom/custom_RawRiseCountsLstilbeneFolded"] = make_hist1_root_object(hist_rise)
        root_file["Histograms/custom/custom_RawFallCountsLstilbeneFolded"] = make_hist1_root_object(hist_fall)
        root_file["Histograms/custom/custom_RiseToRiseSpacingLstilbene"] = make_hist1_root_object(hist_spacing)

    write_plot(hist_rise, hist_fall, hist_spacing, Path(args.pdf_output), Path(args.png_output))

    print(f"ROOT output: {root_output}")
    print(f"Lstilbene rise entries: {hist_rise['entries']:.0f}")
    print(f"Lstilbene fall entries: {hist_fall['entries']:.0f}")
    print(f"Rise-spacing pairs: {hist_spacing['entries']:.0f}")


if __name__ == "__main__":
    main()