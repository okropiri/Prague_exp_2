from __future__ import annotations

import argparse
import csv
from collections import Counter
from decimal import Decimal
from pathlib import Path
from statistics import median


def parse_decimal_series(csv_path: Path, column_name: str) -> list[Decimal]:
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return [Decimal(row[column_name]) for row in reader if row[column_name] != ""]


def format_decimal(value: Decimal) -> str:
    return format(value, "f")


def summarize_series(label: str, values: list[Decimal]) -> str:
    counts = Counter(values)
    lines = [f"{label} summary:"]
    lines.append(f"  samples: {len(values)}")
    lines.append(f"  unique values: {len(counts)}")
    lines.append(f"  min: {format_decimal(min(values))}")
    lines.append(f"  median: {format_decimal(median(values))}")
    lines.append(f"  max: {format_decimal(max(values))}")
    lines.append("  exact value counts:")
    for value, count in sorted(counts.items()):
        lines.append(f"    {format_decimal(value)}: {count}")
    return "\n".join(lines)


def save_plot(values: list[Decimal], label: str, output_path: Path) -> bool:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    counts = Counter(values)
    ordered_values = sorted(counts)
    x_values = [float(value) for value in ordered_values]
    y_values = [counts[value] for value in ordered_values]

    figure, axis = plt.subplots(figsize=(8, 4.5))
    if len(ordered_values) <= 50:
        axis.bar(x_values, y_values, width=0.6 * min_positive_spacing(ordered_values, Decimal("1")))
    else:
        axis.hist([float(value) for value in values], bins=100)

    axis.set_title(label)
    axis.set_xlabel(label)
    axis.set_ylabel("Count")
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return True


def min_positive_spacing(values: list[Decimal], fallback: Decimal) -> float:
    if len(values) < 2:
        return float(fallback)

    spacings = [right - left for left, right in zip(values, values[1:]) if right > left]
    if not spacings:
        return float(fallback)

    return float(min(spacings))


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize Trigger_To_Trigger and Trigger_Jitter distributions.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("Results/results_normalized.csv"),
        help="CSV file to inspect.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("Results"),
        help="Directory where plots will be written if matplotlib is available.",
    )
    args = parser.parse_args()

    csv_path = args.input
    if not csv_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []

    if "Trigger_To_Trigger (us)" in fieldnames:
        trigger_column = "Trigger_To_Trigger (us)"
    elif "Trigger_To_Trigger (s)" in fieldnames:
        trigger_column = "Trigger_To_Trigger (s)"
    else:
        raise ValueError("Input CSV does not contain a Trigger_To_Trigger column.")

    if "Trigger_Jitter (ns)" not in fieldnames:
        raise ValueError("Input CSV does not contain a Trigger_Jitter (ns) column.")

    trigger_values = parse_decimal_series(csv_path, trigger_column)
    jitter_values = parse_decimal_series(csv_path, "Trigger_Jitter (ns)")

    print(summarize_series(trigger_column, trigger_values))
    print()
    print(summarize_series("Trigger_Jitter (ns)", jitter_values))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    trigger_plot_written = save_plot(
        trigger_values,
        trigger_column,
        args.output_dir / f"{csv_path.stem}_trigger_to_trigger_hist.png",
    )
    jitter_plot_written = save_plot(
        jitter_values,
        "Trigger_Jitter (ns)",
        args.output_dir / f"{csv_path.stem}_trigger_jitter_hist.png",
    )

    if trigger_plot_written and jitter_plot_written:
        print()
        print(f"Saved plots to {args.output_dir}")
    else:
        print()
        print("matplotlib is not available; printed summaries only.")


if __name__ == "__main__":
    main()