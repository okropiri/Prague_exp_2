#!/usr/bin/env python3
"""Trigger-referenced RF (zebra) modulation test that does not need channel 0.

Why this exists
---------------
`step04_run_cleaned_all_channel_rf_period.py` folds hits against the channel-0
RF reference pulse. That is the right tool for the beam acquisitions, but some
runs were taken without the RF reference cabled into the TDC at all (for
example the 2026-04-27 `NCAL_Calibration_Source_Prescaler1024_*` control run,
where TDC_001049 has zero channel-0 hits). For those runs the ch0-referenced
chain cannot produce a phase at all.

The DOGMA trigger itself is the cyclotron RF clock divided by 1024, so every
trigger edge is RF-phase-locked. Hit times in the decoded `.dat` are already
measured relative to that trigger, which means an RF modulation of the hit
times is directly visible without any reference channel.

Method
------
1. Histogram trigger-referenced rising-edge times in fine bins, excluding
   channel 0 so the RF reference pulse cannot inject the signal by
   construction.
2. Divide out a wide running mean. The trigger-window acceptance varies on a
   scale of microseconds, so this removes the envelope while leaving any
   modulation at the tens-of-nanoseconds RF scale untouched.
3. Scan the trial period and take the first-harmonic amplitude of the residual.
4. Calibrate significance against the same statistic evaluated at trial periods
   far from the candidate. The amplitude is positive-definite, so the excess is
   measured against the null mean rather than against zero.

This script only reads the decoded `.dat`; it does not modify any step output.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

RAW_BIN_NS = 0.05
REBIN = 10
BIN_NS = RAW_BIN_NS * REBIN
N_RAW_BINS = 200_000
HALF_SPAN_NS = N_RAW_BINS * RAW_BIN_NS / 2.0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, required=True,
                        help="Decoded DOGMA .dat file, or a cached histogram with --from-histogram.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--label", default="dataset")
    parser.add_argument("--board", type=int, default=1,
                        help="TDC ordinal as printed in the .dat header (Prague: 0=0x1048, 1=0x1049).")
    parser.add_argument("--exclude-channels", type=int, nargs="*", default=[0],
                        help="Channels to drop; channel 0 carries the RF reference pulse.")
    parser.add_argument("--max-lines", type=int, default=0,
                        help="Stop after this many input lines (0 reads the whole file).")
    parser.add_argument("--from-histogram", action="store_true",
                        help="Treat --input as a cached histogram produced by an earlier run.")
    parser.add_argument("--time-min-ns", type=float, default=200.0)
    parser.add_argument("--time-max-ns", type=float, default=5000.0)
    parser.add_argument("--envelope-ns", type=float, default=400.0,
                        help="Running-mean width; must be far larger than the RF period.")
    parser.add_argument("--period-min-ns", type=float, default=30.0)
    parser.add_argument("--period-max-ns", type=float, default=50.0)
    parser.add_argument("--period-step-ns", type=float, default=0.0005)
    parser.add_argument("--reference-period-ns", type=float, default=None,
                        help="Also evaluate the amplitude at this fixed period, e.g. a period "
                             "measured in a beam run, to place a limit on a control run.")
    parser.add_argument("--phase-bins", type=int, default=40)
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args(argv)


def build_histogram(path: Path, board: int, excluded: set[int], max_lines: int) -> np.ndarray:
    """Histogram trigger-referenced rising-edge times from a decoded .dat."""
    counts = np.zeros(N_RAW_BINS, dtype=np.int64)
    current = -1
    seen = 0
    with path.open("r", buffering=1 << 22) as handle:
        for line in handle:
            seen += 1
            if max_lines and seen > max_lines:
                break
            head = line[0]
            if head == "T":
                parts = line.split(maxsplit=3)
                current = int(parts[1])
                continue
            if head == "G" or current != board:
                continue
            fields = line.split()
            if fields[2] != "1" or int(fields[1]) in excluded:
                continue
            index = int((float(fields[3]) + HALF_SPAN_NS) / RAW_BIN_NS)
            if 0 <= index < N_RAW_BINS:
                counts[index] += 1
    return counts


def residual_profile(counts: np.ndarray, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, int]:
    """Rebin, divide out the acceptance envelope and restrict to the analysis window."""
    rebinned = counts.reshape(-1, REBIN).sum(axis=1).astype(np.float64)
    centres = -HALF_SPAN_NS + (np.arange(rebinned.size) + 0.5) * BIN_NS
    width = int(args.envelope_ns / BIN_NS) | 1
    padded = np.pad(rebinned, width // 2, mode="edge")
    envelope = np.convolve(padded, np.ones(width) / width, mode="valid")
    keep = (centres >= args.time_min_ns) & (centres < args.time_max_ns) & (envelope > 0)
    if not keep.any():
        raise ValueError("No usable bins in the requested time window")
    return centres[keep], rebinned[keep] / envelope[keep] - 1.0, int(rebinned[keep].sum())


def harmonic_amplitude(times: np.ndarray, residual: np.ndarray, period: float) -> float:
    phase = 2.0 * np.pi * (times / period % 1.0)
    return 2.0 * np.hypot((residual * np.cos(phase)).mean(),
                          (residual * np.sin(phase)).mean())


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.from_histogram:
        counts = np.loadtxt(args.input).astype(np.int64)
        if counts.size != N_RAW_BINS:
            raise ValueError(f"Cached histogram must have {N_RAW_BINS} bins, got {counts.size}")
    else:
        counts = build_histogram(args.input, args.board, set(args.exclude_channels), args.max_lines)
        np.savetxt(args.output_dir / f"{args.label}_hist.txt", counts, fmt="%d")

    times, residual, used = residual_profile(counts, args)
    print(f"{args.label}: {used:,} hits in {times.size} bins "
          f"({args.time_min_ns:.0f}..{args.time_max_ns:.0f} ns)", flush=True)

    periods = np.arange(args.period_min_ns, args.period_max_ns, args.period_step_ns)
    amplitude = np.array([harmonic_amplitude(times, residual, p) for p in periods])
    best_index = int(np.argmax(amplitude))
    best_period = float(periods[best_index])

    # Null calibration from trial periods far enough away not to overlap the peak.
    far = np.abs(periods - best_period) > 1.0
    null_mean = float(amplitude[far].mean())
    null_sd = float(amplitude[far].std())
    excess = (amplitude[best_index] - null_mean) / null_sd if null_sd else float("nan")

    phase_index = np.clip(((times / best_period % 1.0) * args.phase_bins).astype(int),
                          0, args.phase_bins - 1)
    profile = np.array([residual[phase_index == k].mean() for k in range(args.phase_bins)])
    errors = np.array([residual[phase_index == k].std() / np.sqrt(max((phase_index == k).sum(), 1))
                       for k in range(args.phase_bins)])
    depth = float(profile.max() - profile.min())

    summary = {
        "label": args.label,
        "hits": used,
        "best_period_ns": best_period,
        "amplitude": float(amplitude[best_index]),
        "null_mean": null_mean,
        "null_sd": null_sd,
        "excess_sigma": float(excess),
        "fold_depth": depth,
    }
    if args.reference_period_ns:
        ref_amp = harmonic_amplitude(times, residual, args.reference_period_ns)
        ref_excess = (ref_amp - null_mean) / null_sd if null_sd else float("nan")
        summary["reference_period_ns"] = args.reference_period_ns
        summary["amplitude_at_reference"] = float(ref_amp)
        summary["excess_sigma_at_reference"] = float(ref_excess)
        # For a non-detection quote the null-referenced limit, not amplitude+error.
        summary["ul95_at_reference"] = float(
            ref_amp + 1.96 * null_sd if ref_excess >= 2 else null_mean + 1.96 * null_sd)

    path = args.output_dir / f"{args.label}_rf_modulation_summary.tsv"
    with path.open("w") as handle:
        handle.write("key\tvalue\n")
        for key, value in summary.items():
            handle.write(f"{key}\t{value}\n")
    for key, value in summary.items():
        print(f"  {key} = {value}")
    print("wrote", path)

    np.savetxt(args.output_dir / f"{args.label}_period_scan.tsv",
               np.column_stack([periods, amplitude]), fmt="%.6f",
               header="trial_period_ns\tharmonic_amplitude", delimiter="\t", comments="")
    np.savetxt(args.output_dir / f"{args.label}_folded_profile.tsv",
               np.column_stack([(np.arange(args.phase_bins) + 0.5) / args.phase_bins,
                                profile, errors]),
               fmt="%.6f", header="phase_turns\tresidual\terror", delimiter="\t", comments="")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots(1, 2, figsize=(13, 4.5))
        axes[0].plot(periods, amplitude * 100, lw=0.5)
        axes[0].axvline(best_period, color="tab:red", ls="--", lw=1,
                        label=f"best {best_period:.3f} ns")
        if args.reference_period_ns:
            axes[0].axvline(args.reference_period_ns, color="tab:green", ls=":", lw=1,
                            label=f"reference {args.reference_period_ns:.3f} ns")
        axes[0].set_xlabel("trial period (ns)")
        axes[0].set_ylabel("modulation amplitude (%)")
        axes[0].set_title(f"{args.label}: period scan, {used:,} hits")
        axes[0].legend(loc="upper right")

        axes[1].errorbar((np.arange(args.phase_bins) + 0.5) / args.phase_bins,
                         profile * 100, yerr=errors * 100, fmt="o-", ms=3, lw=1)
        axes[1].axhline(0, color="k", lw=0.5)
        axes[1].set_xlabel("RF phase (turns)")
        axes[1].set_ylabel("residual (%)")
        axes[1].set_title(f"folded at {best_period:.3f} ns, depth {depth*100:.2f}%")
        figure.suptitle(f"Trigger-referenced RF modulation ({args.label}, "
                        f"channels excluded: {args.exclude_channels})")
        figure.tight_layout()
        figure.savefig(args.output_dir / f"{args.label}_rf_modulation.png", dpi=150)

    return 0


if __name__ == "__main__":
    sys.exit(main())
