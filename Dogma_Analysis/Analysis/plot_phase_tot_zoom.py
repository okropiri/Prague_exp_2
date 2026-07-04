from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

from step05_estimate_gamma_neutron_background import read_channel_histogram


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a zoomed phase-ToT sparse histogram for one channel.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--channel", type=int, default=2)
    parser.add_argument("--tot-max-ns", type=float, default=40.0)
    parser.add_argument("--z-scale", choices=("linear", "log"), default="linear")
    parser.add_argument("--linear-z-quantile", type=float, default=0.995)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    histogram = read_channel_histogram(Path(args.input), args.channel)
    display_mask = histogram.y_centers_ns <= args.tot_max_ns
    display_counts = histogram.counts[display_mask, :]
    positive = display_counts[display_counts > 0.0]
    vmax = float(np.quantile(positive, args.linear_z_quantile)) if positive.size else 1.0
    masked = np.ma.masked_less_equal(display_counts, 0.0)

    figure, axis = plt.subplots(figsize=(12.5, 5.0), constrained_layout=True)
    image_kwargs = {
        "aspect": "auto",
        "origin": "lower",
        "extent": (
            float(histogram.x_centers_ns[0] - 0.5 * histogram.x_bin_width_ns),
            float(histogram.x_centers_ns[-1] + 0.5 * histogram.x_bin_width_ns),
            0.0,
            args.tot_max_ns,
        ),
        "interpolation": "nearest",
        "cmap": "magma",
    }
    if args.z_scale == "log":
        image = axis.imshow(masked, norm=colors.LogNorm(vmin=1.0, vmax=max(1.0, vmax)), **image_kwargs)
        colorbar_label = f"counts/bin, log scale, vmax=q{args.linear_z_quantile:.3f}"
    else:
        image = axis.imshow(masked, vmin=0.0, vmax=max(1.0, vmax), **image_kwargs)
        colorbar_label = f"counts/bin, linear scale, clipped at q{args.linear_z_quantile:.3f}"

    figure.colorbar(image, ax=axis, label=colorbar_label)
    run_key = histogram.metadata.get("run_key", Path(args.input).stem)
    axis.set_xlabel("RF phase, three periods (ns)")
    axis.set_ylabel("ToT (ns)")
    axis.set_title(f"{run_key} | Ch{args.channel:02d} phase-ToT | 0-{args.tot_max_ns:g} ns")
    figure.savefig(Path(args.output), dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
