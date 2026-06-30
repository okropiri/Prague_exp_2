#!/usr/bin/env python3
"""Create the analytic prompt-gamma / fast-neutron TOF reference plot."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


C_M_PER_S = 299_792_458.0
NEUTRON_REST_MEV = 939.56542052
REFERENCE_NEUTRON_MEV = 30.0
UPPER_NEUTRON_MEV = 33.0
POSITIONS_M = np.array([2.8, 3.4, 4.0, 5.0, 5.8, 6.6])


def neutron_beta(kinetic_energy_mev: float) -> float:
    gamma = 1.0 + kinetic_energy_mev / NEUTRON_REST_MEV
    return float(np.sqrt(1.0 - 1.0 / (gamma * gamma)))


def tof_ns(distance_m: np.ndarray, beta: float) -> np.ndarray:
    return distance_m / (beta * C_M_PER_S) * 1.0e9


def main() -> None:
    latex_dir = Path(__file__).resolve().parents[1]
    figure_dir = latex_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    png_path = figure_dir / "Gamma_Neutron_deltaTOF_30MeV_reference_positions.png"
    csv_path = figure_dir / "Gamma_Neutron_deltaTOF_30MeV_reference_positions.csv"

    distances_m = np.linspace(2.5, 6.8, 400)
    distances_cm = distances_m * 100.0

    gamma_tof = tof_ns(distances_m, 1.0)
    neutron_30_tof = tof_ns(distances_m, neutron_beta(REFERENCE_NEUTRON_MEV))
    neutron_33_tof = tof_ns(distances_m, neutron_beta(UPPER_NEUTRON_MEV))
    delta_30 = neutron_30_tof - gamma_tof
    delta_33 = neutron_33_tof - gamma_tof

    position_gamma_tof = tof_ns(POSITIONS_M, 1.0)
    position_neutron_30_tof = tof_ns(POSITIONS_M, neutron_beta(REFERENCE_NEUTRON_MEV))
    position_delta_30 = position_neutron_30_tof - position_gamma_tof

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "position_m",
                "position_cm",
                "prompt_gamma_tof_ns",
                "neutron_30MeV_tof_ns",
                "delta_30MeV_minus_gamma_ns",
            ]
        )
        for distance_m, gamma_ns, neutron_ns, delta_ns in zip(
            POSITIONS_M, position_gamma_tof, position_neutron_30_tof, position_delta_30
        ):
            writer.writerow(
                [
                    f"{distance_m:.3f}",
                    f"{distance_m * 100.0:.1f}",
                    f"{gamma_ns:.6f}",
                    f"{neutron_ns:.6f}",
                    f"{delta_ns:.6f}",
                ]
            )

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(9.2, 5.5), dpi=180)

    ax.plot(distances_cm, neutron_30_tof, color="#174f90", linewidth=2.4, label="30 MeV neutron TOF")
    ax.plot(distances_cm, gamma_tof, color="#d84a2b", linewidth=2.4, label="prompt gamma TOF")
    ax.plot(distances_cm, delta_30, color="#d5a100", linewidth=2.4, label="TOF separation, 30 MeV n - gamma")
    ax.fill_between(
        distances_cm,
        neutron_33_tof,
        neutron_30_tof,
        color="#174f90",
        alpha=0.12,
        label="30-33 MeV neutron TOF band",
    )
    ax.fill_between(
        distances_cm,
        delta_33,
        delta_30,
        color="#d5a100",
        alpha=0.12,
        label="30-33 MeV separation band",
    )

    ax.scatter(POSITIONS_M * 100.0, position_neutron_30_tof, color="#174f90", s=34, zorder=4)
    ax.scatter(POSITIONS_M * 100.0, position_gamma_tof, color="#d84a2b", s=34, zorder=4)
    ax.scatter(POSITIONS_M * 100.0, position_delta_30, color="#d5a100", s=34, zorder=4)

    for distance_m, delta_ns in zip(POSITIONS_M, position_delta_30):
        ax.text(distance_m * 100.0, delta_ns + 2.0, f"{distance_m:.1f} m", fontsize=8, ha="center")

    ax.set_title("Analytic prompt-gamma and fast-neutron TOF reference")
    ax.set_xlabel("Detector-face distance (cm)")
    ax.set_ylabel("Time of flight / separation (ns)")
    ax.set_xlim(240, 680)
    ax.set_ylim(0, max(neutron_30_tof) * 1.08)
    ax.legend(loc="upper left", fontsize=8, frameon=True)
    ax.text(
        0.985,
        0.03,
        "Relativistic neutron velocity; gamma at c",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="#555555",
    )

    fig.tight_layout()
    fig.savefig(png_path)
    plt.close(fig)


if __name__ == "__main__":
    main()
