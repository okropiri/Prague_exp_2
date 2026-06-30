from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.colors import LogNorm
from uproot.writing.identify import to_TAxis, to_TH1x, to_TH2x


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write ROOT, PDF, and PNG outputs for colleague-style trigger-referenced DOGMA histograms."
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


def sum_hist2_over_y(hist: dict[str, object], object_name: str, title: str) -> dict[str, object]:
    storage2d = np.asarray(hist["storage"], dtype=np.float64)
    counts = storage2d[1:-1, 1:-1].sum(axis=0)
    x_bins = int(hist["x_bins"])
    x_min = float(hist["x_min"])
    x_max = float(hist["x_max"])
    storage1d = np.zeros(x_bins + 2, dtype=np.float64)
    storage1d[1:-1] = counts

    bin_width = (x_max - x_min) / x_bins
    centers = x_min + (np.arange(x_bins, dtype=np.float64) + 0.5) * bin_width
    entries = float(counts.sum(dtype=np.float64))

    return {
        "object_name": object_name,
        "title": title,
        "x_title": hist["x_title"],
        "bins": x_bins,
        "min": x_min,
        "max": x_max,
        "entries": entries,
        "sumw": entries,
        "sumw2": entries,
        "sumwx": float((counts * centers).sum(dtype=np.float64)),
        "sumwx2": float((counts * np.square(centers)).sum(dtype=np.float64)),
        "storage": storage1d,
    }


def write_folded_counts_plot(
    histograms: list[tuple[dict[str, object], str]],
    pdf_path: Path,
    png_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    for hist, color in histograms:
        storage = np.asarray(hist["storage"], dtype=np.float64)
        values = storage[1:-1]
        x_min = float(hist["min"])
        x_max = float(hist["max"])
        bins = int(hist["bins"])
        centers = x_min + (np.arange(bins, dtype=np.float64) + 0.5) * ((x_max - x_min) / bins)
        ax.plot(centers, values, linewidth=0.9, color=color, label=str(hist["object_name"]))

    ax.set_xlabel(str(histograms[0][0]["x_title"]))
    ax.set_ylabel("counts / bin")
    ax.set_title("Folded counts vs trigger-relative time")
    ax.grid(True, linewidth=0.3, alpha=0.4)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=180)
    plt.close(fig)


def write_plot(hist: dict[str, object], pdf_path: Path, png_path: Path) -> None:
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


def main() -> None:
    args = parse_args()
    prefix = Path(args.input_prefix)

    hist_ncal1 = read_hist1(prefix.with_name(prefix.name + "_custom_Ncal1_hist.txt"))
    hist_prof_ncal1 = read_hist2(prefix.with_name(prefix.name + "_custom_ProfNcal1_hist.txt"))
    hist_prof_lstilbene = read_hist2(prefix.with_name(prefix.name + "_custom_ProfLstilbene_hist.txt"))
    hist_prof_sstilbene = read_hist2(prefix.with_name(prefix.name + "_custom_ProfSstilbene_hist.txt"))
    hist_counts_ncal1 = sum_hist2_over_y(hist_prof_ncal1, "custom_CountsNcal1Folded", "Ncal1 counts vs trigger-relative time")
    hist_counts_lstilbene = sum_hist2_over_y(hist_prof_lstilbene, "custom_CountsLstilbeneFolded", "Lstilbene counts vs trigger-relative time")
    hist_counts_sstilbene = sum_hist2_over_y(hist_prof_sstilbene, "custom_CountsSstilbeneFolded", "Sstilbene counts vs trigger-relative time")

    root_output = Path(args.root_output)
    with uproot.recreate(root_output) as root_file:
        root_file["Histograms/custom/custom_Ncal1"] = make_hist1_root_object(hist_ncal1)
        root_file["Histograms/custom/custom_ProfNcal1"] = make_hist2_root_object(hist_prof_ncal1)
        root_file["Histograms/custom/custom_ProfLstilbene"] = make_hist2_root_object(hist_prof_lstilbene)
        root_file["Histograms/custom/custom_ProfSstilbene"] = make_hist2_root_object(hist_prof_sstilbene)
        root_file["Histograms/custom/custom_CountsNcal1Folded"] = make_hist1_root_object(hist_counts_ncal1)
        root_file["Histograms/custom/custom_CountsLstilbeneFolded"] = make_hist1_root_object(hist_counts_lstilbene)
        root_file["Histograms/custom/custom_CountsSstilbeneFolded"] = make_hist1_root_object(hist_counts_sstilbene)

    pdf_prefix = Path(args.pdf_prefix)
    png_prefix = Path(args.png_prefix)
    write_plot(hist_prof_ncal1, pdf_prefix.with_name(pdf_prefix.name + "_Ncal1_ToT_T.pdf"), png_prefix.with_name(png_prefix.name + "_Ncal1_ToT_T.png"))
    write_plot(hist_prof_lstilbene, pdf_prefix.with_name(pdf_prefix.name + "_Lstilbene_ToT_T.pdf"), png_prefix.with_name(png_prefix.name + "_Lstilbene_ToT_T.png"))
    write_plot(hist_prof_sstilbene, pdf_prefix.with_name(pdf_prefix.name + "_Sstilbene_ToT_T.pdf"), png_prefix.with_name(png_prefix.name + "_Sstilbene_ToT_T.png"))
    write_folded_counts_plot(
        [(hist_counts_ncal1, "#1f77b4"), (hist_counts_lstilbene, "#d62728"), (hist_counts_sstilbene, "#2ca02c")],
        pdf_prefix.with_name(pdf_prefix.name + "_FoldedCounts_T.pdf"),
        png_prefix.with_name(png_prefix.name + "_FoldedCounts_T.png"),
    )

    print(f"ROOT output: {root_output}")
    print(f"Ncal1 entries: {hist_prof_ncal1['entries']:.0f}")
    print(f"Lstilbene entries: {hist_prof_lstilbene['entries']:.0f}")
    print(f"Sstilbene entries: {hist_prof_sstilbene['entries']:.0f}")


if __name__ == "__main__":
    main()