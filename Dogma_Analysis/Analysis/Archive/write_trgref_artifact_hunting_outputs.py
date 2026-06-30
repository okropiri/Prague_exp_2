from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot
from uproot.writing.identify import to_TAxis, to_TH1x


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write ROOT, PDF, and PNG outputs for multiplicity-aware trigger-referenced artifact-hunting histograms."
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


def write_overlay_plot(histograms: list[tuple[dict[str, object], str]], pdf_path: Path, png_path: Path) -> None:
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
    ax.set_title("Multiplicity-aware folded counts vs trigger-relative time")
    ax.grid(True, linewidth=0.3, alpha=0.4)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    prefix = Path(args.input_prefix)

    hist_ncal1 = read_hist1(prefix.with_name(prefix.name + "_custom_RawRiseCountsNcal1Folded_hist.txt"))
    hist_lstilbene = read_hist1(prefix.with_name(prefix.name + "_custom_RawRiseCountsLstilbeneFolded_hist.txt"))
    hist_sstilbene = read_hist1(prefix.with_name(prefix.name + "_custom_RawRiseCountsSstilbeneFolded_hist.txt"))

    root_output = Path(args.root_output)
    with uproot.recreate(root_output) as root_file:
        root_file["Histograms/custom/custom_RawRiseCountsNcal1Folded"] = make_hist1_root_object(hist_ncal1)
        root_file["Histograms/custom/custom_RawRiseCountsLstilbeneFolded"] = make_hist1_root_object(hist_lstilbene)
        root_file["Histograms/custom/custom_RawRiseCountsSstilbeneFolded"] = make_hist1_root_object(hist_sstilbene)

    write_overlay_plot(
        [(hist_ncal1, "#1f77b4"), (hist_lstilbene, "#d62728"), (hist_sstilbene, "#2ca02c")],
        Path(args.pdf_output),
        Path(args.png_output),
    )

    print(f"ROOT output: {root_output}")
    print(f"Ncal1 raw rises: {hist_ncal1['entries']:.0f}")
    print(f"Lstilbene raw rises: {hist_lstilbene['entries']:.0f}")
    print(f"Sstilbene raw rises: {hist_sstilbene['entries']:.0f}")


if __name__ == "__main__":
    main()