from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize phi* = (phase-like rise coordinate + ToT) mod period from existing step-03 and step-04 sparse outputs."
        )
    )
    parser.add_argument("--folded3x-phase-tot-file", required=True)
    parser.add_argument("--ch0ref-tot-file", required=True)
    parser.add_argument("--channels", nargs="*", type=int, default=[2])
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--min-channel-total", type=int, default=50000)
    parser.add_argument("--phase-bin-width-ns", type=float, default=0.25)
    return parser.parse_args()


def parse_metadata_line(line: str) -> tuple[str, str] | None:
    if not line.startswith("#") or "=" not in line:
        return None
    key, value = line[1:].split("=", 1)
    return key.strip(), value.strip()


def positive_mod(value: float, period: float) -> float:
    wrapped = value % period
    return wrapped if wrapped >= 0.0 else wrapped + period


def read_sparse_metadata(path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            parsed = parse_metadata_line(line)
            if parsed is None:
                continue
            key, value = parsed
            metadata[key] = value
    return metadata


def phase_bins_for_period(period_ns: float, phase_bin_width_ns: float) -> int:
    return max(16, round(period_ns / phase_bin_width_ns))


def ensure_histogram(histograms: dict[int, list[int]], channel: int, bins: int) -> list[int]:
    histogram = histograms.get(channel)
    if histogram is None:
        histogram = [0] * bins
        histograms[channel] = histogram
    return histogram


def accumulate_folded3x_phi_star(
    path: Path,
    *,
    phase_bins: int,
    period_ns: float,
) -> tuple[dict[int, list[int]], dict[int, int]]:
    metadata = read_sparse_metadata(path)
    x_min_ns = float(metadata["x_min_ns"])
    x_max_ns = float(metadata["x_max_ns"])
    x_bins = int(metadata["x_bins"])
    y_min_ns = float(metadata["y_min_ns"])
    y_max_ns = float(metadata["y_max_ns"])
    y_bins = int(metadata["y_bins"])
    x_width_ns = (x_max_ns - x_min_ns) / x_bins
    y_width_ns = (y_max_ns - y_min_ns) / y_bins
    phase_width_ns = period_ns / phase_bins

    histograms: dict[int, list[int]] = {}
    totals: dict[int, int] = defaultdict(int)
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            channel, x_bin, y_bin, count = (int(piece) for piece in line.split())
            phase3x_ns = x_min_ns + (x_bin + 0.5) * x_width_ns
            tot_ns = y_min_ns + (y_bin + 0.5) * y_width_ns
            phi_star_ns = positive_mod(positive_mod(phase3x_ns, period_ns) + tot_ns, period_ns)
            phi_bin = min(phase_bins - 1, int(phi_star_ns / phase_width_ns))
            ensure_histogram(histograms, channel, phase_bins)[phi_bin] += count
            totals[channel] += count
    return histograms, totals


def accumulate_ch0ref_phi_star(
    path: Path,
    *,
    phase_bins: int,
    period_ns: float,
    phase_origin_ns: float,
) -> tuple[dict[int, list[int]], dict[int, int]]:
    metadata = read_sparse_metadata(path)
    x_min_ns = float(metadata["x_min_ns"])
    x_max_ns = float(metadata["x_max_ns"])
    x_bins = int(metadata["x_bins"])
    tot_min_ns = float(metadata["tot_min_ns"])
    tot_max_ns = float(metadata["tot_max_ns"])
    tot_bins = int(metadata["tot_bins"])
    x_width_ns = (x_max_ns - x_min_ns) / x_bins
    tot_width_ns = (tot_max_ns - tot_min_ns) / tot_bins
    phase_width_ns = period_ns / phase_bins

    histograms: dict[int, list[int]] = {}
    totals: dict[int, int] = defaultdict(int)
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            channel, x_bin, y_bin, count = (int(piece) for piece in line.split())
            rise_rel_ns = x_min_ns + (x_bin + 0.5) * x_width_ns
            tot_ns = tot_min_ns + (y_bin + 0.5) * tot_width_ns
            phi_star_ns = positive_mod(rise_rel_ns - phase_origin_ns + tot_ns, period_ns)
            phi_bin = min(phase_bins - 1, int(phi_star_ns / phase_width_ns))
            ensure_histogram(histograms, channel, phase_bins)[phi_bin] += count
            totals[channel] += count
    return histograms, totals


def top_peaks(histogram: list[int], total: int, period_ns: float, limit: int) -> list[tuple[float, int, float]]:
    peaks: list[tuple[float, int, float]] = []
    if total <= 0:
        return peaks
    bin_width_ns = period_ns / len(histogram)
    for index, count in sorted(enumerate(histogram), key=lambda item: item[1], reverse=True)[:limit]:
        center_ns = (index + 0.5) * bin_width_ns
        peaks.append((center_ns, count, count / total))
    return peaks


def circular_delta_ns(left_ns: float, right_ns: float, period_ns: float) -> float:
    delta = abs(left_ns - right_ns)
    return min(delta, period_ns - delta)


def main() -> None:
    args = parse_args()
    folded3x_path = Path(args.folded3x_phase_tot_file)
    ch0ref_path = Path(args.ch0ref_tot_file)

    folded3x_metadata = read_sparse_metadata(folded3x_path)
    period_ns = float(folded3x_metadata["deduced_period_ns"])
    phase_origin_ns = float(folded3x_metadata["phase_origin_ns"])
    phase_bins = phase_bins_for_period(period_ns, args.phase_bin_width_ns)

    folded3x_histograms, folded3x_totals = accumulate_folded3x_phi_star(
        folded3x_path,
        phase_bins=phase_bins,
        period_ns=period_ns,
    )
    ch0ref_histograms, ch0ref_totals = accumulate_ch0ref_phi_star(
        ch0ref_path,
        phase_bins=phase_bins,
        period_ns=period_ns,
        phase_origin_ns=phase_origin_ns,
    )

    print(f"period_ns={period_ns:.6f}")
    print(f"phase_origin_ns={phase_origin_ns:.6f}")
    print(f"phase_bins={phase_bins}")
    print("")

    for channel in args.channels:
        step03_total = ch0ref_totals.get(channel, 0)
        step04_total = folded3x_totals.get(channel, 0)
        print(f"channel=Ch{channel:02d}")
        print(f"  step03_total={step03_total}")
        if step03_total > 0:
            peaks = top_peaks(ch0ref_histograms[channel], step03_total, period_ns, args.top_n)
            print(f"  step03_top{args.top_n}_fraction={sum(frac for _, _, frac in peaks):.4f}")
            for index, (center_ns, count, fraction) in enumerate(peaks, start=1):
                print(f"  step03_peak_{index}_ns={center_ns:.3f} count={count} frac={fraction:.4f}")
        print(f"  step04_total={step04_total}")
        if step04_total > 0:
            peaks = top_peaks(folded3x_histograms[channel], step04_total, period_ns, args.top_n)
            print(f"  step04_top{args.top_n}_fraction={sum(frac for _, _, frac in peaks):.4f}")
            for index, (center_ns, count, fraction) in enumerate(peaks, start=1):
                print(f"  step04_peak_{index}_ns={center_ns:.3f} count={count} frac={fraction:.4f}")
        if step03_total > 0 and step04_total > 0:
            peak03_ns = top_peaks(ch0ref_histograms[channel], step03_total, period_ns, 1)[0][0]
            peak04_ns = top_peaks(folded3x_histograms[channel], step04_total, period_ns, 1)[0][0]
            print(f"  top_peak_circular_delta_ns={circular_delta_ns(peak03_ns, peak04_ns, period_ns):.3f}")
        print("")

    aligned_channels: list[tuple[int, float, float, float, int, int]] = []
    channel_ids = sorted(set(ch0ref_totals) | set(folded3x_totals))
    for channel in channel_ids:
        step03_total = ch0ref_totals.get(channel, 0)
        step04_total = folded3x_totals.get(channel, 0)
        if step03_total < args.min_channel_total or step04_total < args.min_channel_total:
            continue
        peak03_ns = top_peaks(ch0ref_histograms[channel], step03_total, period_ns, 1)[0][0]
        peak04_ns = top_peaks(folded3x_histograms[channel], step04_total, period_ns, 1)[0][0]
        aligned_channels.append(
            (
                channel,
                peak03_ns,
                peak04_ns,
                circular_delta_ns(peak03_ns, peak04_ns, period_ns),
                step03_total,
                step04_total,
            )
        )

    print("channels_with_sufficient_stats=")
    for channel, peak03_ns, peak04_ns, delta_ns, step03_total, step04_total in aligned_channels:
        print(
            f"  Ch{channel:02d}: step03_peak_ns={peak03_ns:.3f} step04_peak_ns={peak04_ns:.3f} "
            f"circular_delta_ns={delta_ns:.3f} totals=({step03_total},{step04_total})"
        )
    within_one_ns = sum(1 for _, _, _, delta_ns, _, _ in aligned_channels if delta_ns <= 1.0)
    print("")
    print(f"aligned_within_1ns={within_one_ns}/{len(aligned_channels)}")


if __name__ == "__main__":
    main()