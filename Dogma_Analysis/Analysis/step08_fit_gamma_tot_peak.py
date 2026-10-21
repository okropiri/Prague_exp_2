from __future__ import annotations
from analysis_provenance import record_arguments

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class PeakFitResult:
    run_key: str
    position_m: float
    duration_seconds: float
    fit_tot_min_ns: float
    fit_tot_max_ns: float
    peak_mu_ns: float
    peak_mu_err_ns: float
    peak_sigma_ns: float
    peak_sigma_err_ns: float
    peak_amplitude_hz: float
    peak_amplitude_err_hz: float
    peak_offset_hz: float
    peak_offset_err_hz: float
    chi2: float
    ndof: int
    reduced_chi2: float
    fit_success: bool
    fit_message: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit the sharp Ch02 gamma ToT peak with a Gaussian plus constant background "
            "using rate-normalized ToT distributions."
        )
    )
    parser.add_argument("--results-root", default="/data6/Dogma_analysis_by_Dachi/Results")
    parser.add_argument("--output-dir")
    parser.add_argument("--run-glob", default="NCAL_20us_Pos_*")
    parser.add_argument("--channel", type=int, default=2)
    parser.add_argument("--fit-tot-min-ns", type=float, default=12.5)
    parser.add_argument("--fit-tot-max-ns", type=float, default=15.0)
    parser.add_argument("--display-tot-max-ns", type=float, default=25.0)
    parser.add_argument("--max-iterations", type=int, default=30)
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    return record_arguments(parser.parse_args())


def parse_position_m(run_key: str) -> float:
    if "Pos_" not in run_key:
        return float("nan")
    return float(run_key.split("Pos_", 1)[1].split("m_", 1)[0])


def read_key_value_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip().strip('"')
    return values


def run_duration_seconds(run_dir: Path, run_key: str) -> float:
    summary = read_key_value_file(run_dir / f"{run_key}_cleaned_rates_summary.txt")
    return float(summary["covered_duration_seconds"])


def read_tot_distribution(path: Path, channel: int) -> tuple[str, np.ndarray, np.ndarray]:
    metadata: dict[str, str] = {}
    counts_by_bin: dict[int, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#") and "=" in line:
                key, value = line[1:].split("=", 1)
                metadata[key.strip()] = value.strip()
                continue
            if line.startswith("#"):
                continue
            channel_text, tot_bin_text, count_text = line.split()
            if int(channel_text) != channel:
                continue
            counts_by_bin[int(tot_bin_text)] = int(count_text)

    run_key = metadata.get("run_key", path.stem.replace("_TOT_distrib", ""))
    tot_min_ns = float(metadata.get("tot_min_ns", "0"))
    tot_bin_width_ns = float(metadata.get("tot_bin_width_ns", "0.1"))
    tot_bins = int(metadata.get("tot_bins", "0"))

    x = tot_min_ns + (np.arange(tot_bins, dtype=np.float64) + 0.5) * tot_bin_width_ns
    y = np.zeros(tot_bins, dtype=np.float64)
    for tot_bin, count in counts_by_bin.items():
        if 0 <= tot_bin < tot_bins:
            y[tot_bin] = count
    return run_key, x, y


def gaussian_plus_constant(x: np.ndarray, amplitude: float, mu: float, log_sigma: float, offset: float) -> np.ndarray:
    sigma = float(np.exp(log_sigma))
    z = (x - mu) / sigma
    return offset + amplitude * np.exp(-0.5 * z * z)


def jacobian(x: np.ndarray, amplitude: float, mu: float, log_sigma: float, offset: float) -> np.ndarray:
    sigma = float(np.exp(log_sigma))
    z = (x - mu) / sigma
    exp_term = np.exp(-0.5 * z * z)
    jac = np.zeros((x.size, 4), dtype=np.float64)
    jac[:, 0] = exp_term
    jac[:, 1] = amplitude * exp_term * (z / sigma)
    jac[:, 2] = amplitude * exp_term * (z * z)
    jac[:, 3] = 1.0
    return jac


def fit_peak(x: np.ndarray, rate: np.ndarray, sigma_rate: np.ndarray, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, float, int, bool, str]:
    offset0 = float(np.quantile(rate, 0.20))
    amplitude0 = max(float(np.max(rate) - offset0), max(float(np.max(rate)), 1.0e-6))
    mu0 = float(x[int(np.argmax(rate))])
    sigma0 = max((args.fit_tot_max_ns - args.fit_tot_min_ns) / 6.0, 0.15)
    params = np.asarray([amplitude0, mu0, np.log(sigma0), offset0], dtype=np.float64)

    success = False
    message = "max iterations reached"
    previous_chi2 = np.inf
    for _ in range(args.max_iterations):
        model = gaussian_plus_constant(x, *params)
        residual = rate - model
        chi2 = float(np.sum((residual / sigma_rate) ** 2))
        jac = jacobian(x, *params)
        jac_w = jac / sigma_rate[:, np.newaxis]
        residual_w = residual / sigma_rate
        delta, *_ = np.linalg.lstsq(jac_w, residual_w, rcond=None)

        step_scale = 1.0
        accepted = False
        for _ in range(12):
            candidate = params + step_scale * delta
            candidate[2] = np.clip(candidate[2], np.log(0.05), np.log(5.0))
            candidate_model = gaussian_plus_constant(x, *candidate)
            candidate_residual = rate - candidate_model
            candidate_chi2 = float(np.sum((candidate_residual / sigma_rate) ** 2))
            if candidate_chi2 <= chi2:
                params = candidate
                previous_chi2 = candidate_chi2
                accepted = True
                break
            step_scale *= 0.5

        if not accepted:
            success = True
            message = "converged (no improving step)"
            break

        scale = np.maximum(np.abs(params), np.asarray([1.0, 1.0, 1.0, 1.0], dtype=np.float64))
        if float(np.max(np.abs(step_scale * delta) / scale)) < args.tolerance:
            success = True
            message = "converged"
            break

    final_model = gaussian_plus_constant(x, *params)
    final_residual = rate - final_model
    chi2 = float(np.sum((final_residual / sigma_rate) ** 2))
    ndof = max(0, x.size - 4)
    jac = jacobian(x, *params)
    jac_w = jac / sigma_rate[:, np.newaxis]
    hessian = jac_w.T @ jac_w
    covariance = np.linalg.pinv(hessian)
    if ndof > 0:
        covariance *= max(1.0, chi2 / ndof)
    return params, covariance, chi2, ndof, success, message


def draw_fit_plot(
    run_key: str,
    position_m: float,
    x: np.ndarray,
    counts: np.ndarray,
    duration_seconds: float,
    fit_mask: np.ndarray,
    params: np.ndarray,
    result: PeakFitResult,
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    rate = counts / duration_seconds
    fit_x = x[fit_mask]
    fit_rate = rate[fit_mask]
    model_rate = gaussian_plus_constant(fit_x, *params)

    figure, axes = plt.subplots(2, 1, figsize=(10.5, 7.0), constrained_layout=True)
    display_mask = x <= args.display_tot_max_ns
    axes[0].plot(x[display_mask], rate[display_mask], color="#1d4ed8", linewidth=1.2, label="rate per ToT bin")
    axes[0].axvspan(args.fit_tot_min_ns, args.fit_tot_max_ns, color="#f59e0b", alpha=0.12, label="fit ROI")
    axes[0].plot(fit_x, model_rate, color="#dc2626", linewidth=1.5, label="Gaussian + constant fit")
    axes[0].axvline(result.peak_mu_ns, color="#059669", linestyle="--", linewidth=1.2, label=f"mu = {result.peak_mu_ns:.4f} ns")
    axes[0].set_xlim(0.0, args.display_tot_max_ns)
    axes[0].set_xlabel("ToT (ns)")
    axes[0].set_ylabel("rate per ToT bin (Hz)")
    axes[0].set_title(f"{run_key} | detector face {position_m:g} m | Ch02 gamma ToT peak fit")
    axes[0].grid(True, linewidth=0.3, alpha=0.5)
    axes[0].legend(fontsize=8)

    residual = fit_rate - model_rate
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].plot(fit_x, residual, color="#7c3aed", linewidth=1.1)
    axes[1].set_xlabel("ToT (ns)")
    axes[1].set_ylabel("fit residual (Hz)")
    axes[1].grid(True, linewidth=0.3, alpha=0.5)

    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_results(results: list[PeakFitResult], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[field.name for field in PeakFitResult.__dataclass_fields__.values()],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for result in results:
            writer.writerow({field: getattr(result, field) for field in PeakFitResult.__dataclass_fields__})


def draw_summary(results: list[PeakFitResult], output_path: Path) -> None:
    if not results:
        return
    positions = np.asarray([result.position_m for result in results], dtype=np.float64)
    mu = np.asarray([result.peak_mu_ns for result in results], dtype=np.float64)
    mu_err = np.asarray([result.peak_mu_err_ns for result in results], dtype=np.float64)
    sigma = np.asarray([result.peak_sigma_ns for result in results], dtype=np.float64)

    figure, axes = plt.subplots(2, 1, figsize=(8.2, 6.8), sharex=True, constrained_layout=True)
    axes[0].errorbar(positions, mu, yerr=mu_err, fmt="o", color="#2563eb", capsize=3)
    axes[0].set_ylabel("fitted ToT peak mu (ns)")
    axes[0].grid(True, linewidth=0.3, alpha=0.5)

    axes[1].plot(positions, sigma, "o", color="#dc2626")
    axes[1].set_xlabel("Detector-face distance (m)")
    axes[1].set_ylabel("fitted sigma (ns)")
    axes[1].grid(True, linewidth=0.3, alpha=0.5)
    figure.suptitle("Sharp gamma ToT peak fit summary | Ch02")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def write_readme(output_dir: Path, args: argparse.Namespace) -> None:
    text = f"""# Step08 Sharp Gamma ToT Peak Fits

This directory contains Gaussian fits to the sharp Ch02 gamma ToT peak.

Method:
- use the Ch02 whole-run ToT distribution already produced by step03;
- convert counts to rate per ToT bin by dividing by `covered_duration_seconds`;
- fit only the narrow ROI [{args.fit_tot_min_ns}, {args.fit_tot_max_ns}] ns;
- model = constant background + Gaussian peak.

Interpretation:
- `peak_mu_ns` is the fitted ToT peak position.
- `peak_mu_err_ns` is the statistical fit uncertainty on that mean under the chosen model and fit ROI.
- it is not the full uncertainty budget; ROI choice, model mismatch, binning, and residual time-walk effects are systematic terms.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    results_root = Path(args.results_root).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else results_root / "Ch02_Position_Comparison" / "gamma_tot_peak_fits"
    output_dir.mkdir(parents=True, exist_ok=True)

    results: list[PeakFitResult] = []
    for run_dir in sorted(results_root.glob(args.run_glob)):
        if not run_dir.is_dir():
            continue
        run_key = run_dir.name
        position_m = parse_position_m(run_key)
        tot_path = run_dir / "TOT_distrib" / f"{run_key}_TOT_distrib.tsv"
        if not tot_path.is_file():
            continue
        _, x, counts = read_tot_distribution(tot_path, args.channel)
        duration_seconds = run_duration_seconds(run_dir, run_key)
        fit_mask = (x >= args.fit_tot_min_ns) & (x <= args.fit_tot_max_ns)
        fit_x = x[fit_mask]
        fit_counts = counts[fit_mask]
        fit_rate = fit_counts / duration_seconds
        sigma_rate = np.sqrt(np.maximum(fit_counts, 1.0)) / duration_seconds

        params, covariance, chi2, ndof, success, message = fit_peak(fit_x, fit_rate, sigma_rate, args)
        amplitude, mu, log_sigma, offset = params
        sigma = float(np.exp(log_sigma))
        amplitude_err = float(np.sqrt(max(0.0, covariance[0, 0])))
        mu_err = float(np.sqrt(max(0.0, covariance[1, 1])))
        log_sigma_err = float(np.sqrt(max(0.0, covariance[2, 2])))
        sigma_err = sigma * log_sigma_err
        offset_err = float(np.sqrt(max(0.0, covariance[3, 3])))
        reduced_chi2 = chi2 / ndof if ndof > 0 else float("nan")

        result = PeakFitResult(
            run_key=run_key,
            position_m=position_m,
            duration_seconds=duration_seconds,
            fit_tot_min_ns=args.fit_tot_min_ns,
            fit_tot_max_ns=args.fit_tot_max_ns,
            peak_mu_ns=float(mu),
            peak_mu_err_ns=mu_err,
            peak_sigma_ns=sigma,
            peak_sigma_err_ns=sigma_err,
            peak_amplitude_hz=float(amplitude),
            peak_amplitude_err_hz=amplitude_err,
            peak_offset_hz=float(offset),
            peak_offset_err_hz=offset_err,
            chi2=chi2,
            ndof=ndof,
            reduced_chi2=reduced_chi2,
            fit_success=success,
            fit_message=message,
        )
        results.append(result)
        draw_fit_plot(run_key, position_m, x, counts, duration_seconds, fit_mask, params, result, output_dir / f"{run_key}_ch{args.channel:02d}_gamma_peak_fit.png", args)
        print(f"{run_key}: mu={result.peak_mu_ns:.4f} +/- {result.peak_mu_err_ns:.4f} ns, sigma={result.peak_sigma_ns:.4f} ns")

    results.sort(key=lambda result: (result.position_m, result.run_key))
    write_results(results, output_dir / "gamma_tot_peak_fit_summary.tsv")
    draw_summary(results, output_dir / "gamma_tot_peak_fit_summary.png")
    write_readme(output_dir, args)
    print(f"Wrote {output_dir}")


if __name__ == "__main__":
    main()