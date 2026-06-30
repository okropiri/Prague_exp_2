from __future__ import annotations

"""Digest cleaned or partially cleaned DLD event data into CSV outputs.

Run commands:
1. Activate the existing environment first:
    conda activate conda_env

2. Run with the default input/output paths:
    python Analysis/digest_dld_data.py

3. Run on a different input file:
    python Analysis/digest_dld_data.py \
         --input Dogma_test_Data/another_file.dld.dat \
         --output Results/another_results.csv \
         --normalized-output Results/another_results_normalized.csv \
         --clean-output Results/another_file_cleaned.dld.dat

What the script cleans automatically while parsing:
- Removes lines like: TDC TDC_... not found
- Ignores empty lines
- Collapses repeated TDC event headers before the trigger/data rows arrive

The raw input file is not modified in place. Cleaning happens during parsing.
"""

import argparse
import csv
import re
from contextlib import nullcontext
from decimal import Decimal
from statistics import median
from pathlib import Path


TRIGGER_PATTERN = re.compile(r"GlobalTriggerTime\s+([-+]?\d+(?:\.\d+)?)")
NOT_FOUND_PATTERN = re.compile(r"^TDC TDC_.* not found$")
MICROSECONDS_PER_SECOND = Decimal("1000000")
NANOSECONDS_PER_SECOND = Decimal("1000000000")


def format_number(value: Decimal | int | str) -> str:
    if value == "":
        return ""
    if isinstance(value, int):
        return str(value)

    decimal_value = value if isinstance(value, Decimal) else Decimal(str(value))
    return format(decimal_value, "f")


def should_skip_line(line: str) -> bool:
    return not line or bool(NOT_FOUND_PATTERN.match(line))


def iter_clean_lines(input_path: Path, clean_output_path: Path | None = None):
    if clean_output_path is not None:
        clean_output_path.parent.mkdir(parents=True, exist_ok=True)

    clean_output_context = (
        clean_output_path.open("w", encoding="utf-8") if clean_output_path else nullcontext()
    )

    with input_path.open("r", encoding="utf-8") as input_handle, clean_output_context as clean_handle:
        pending_meta: str | None = None

        def emit(clean_line: str) -> str:
            if clean_handle is not None:
                clean_handle.write(clean_line + "\n")
            return clean_line

        for raw_line in input_handle:
            line = raw_line.strip()
            if should_skip_line(line):
                continue

            if line.startswith("TDC "):
                pending_meta = line
                continue

            if pending_meta is not None:
                yield emit(pending_meta)
                pending_meta = None

            yield emit(line)


def parse_event_file(
    input_path: Path, clean_output_path: Path | None = None
) -> list[dict[str, Decimal | int | str]]:
    events: list[dict[str, Decimal | int | str]] = []
    current_meta: str | None = None
    current_trigger: Decimal | None = None
    current_leading: Decimal | None = None
    current_trailing: Decimal | None = None

    for line in iter_clean_lines(input_path, clean_output_path):
        if line.startswith("TDC "):
            # If repeated headers appear before trigger/data rows, keep only the latest one.
            current_meta = line
            current_trigger = None
            current_leading = None
            current_trailing = None
            continue

        trigger_match = TRIGGER_PATTERN.search(line)
        if trigger_match:
            if current_meta is None:
                raise ValueError("Found GlobalTriggerTime before an event metadata line.")
            current_trigger = Decimal(trigger_match.group(1))
            continue

        if current_trigger is None or not line.startswith("1 7 "):
            continue

        parts = line.split()
        if len(parts) != 4:
            raise ValueError(f"Unexpected event data row: {line}")

        edge_kind = parts[2]
        edge_value = Decimal(parts[3])

        if edge_kind == "1":
            current_leading = edge_value
        elif edge_kind == "0":
            current_trailing = edge_value
        else:
            raise ValueError(f"Unexpected edge type '{edge_kind}' in row: {line}")

        if current_leading is None or current_trailing is None:
            continue

        previous_event = events[-1] if events else None
        trigger_to_trigger = (
            current_trigger - previous_event["Global_Trigger_Time (s)"]
            if previous_event is not None
            else ""
        )
        leading_to_leading = (
            current_leading - previous_event["Leading_Edge (ns)"]
            if previous_event is not None
            else ""
        )
        trailing_to_trailing = (
            current_trailing - previous_event["Trailing_Edge (ns)"]
            if previous_event is not None
            else ""
        )

        events.append(
            {
                "Event_Index": len(events) + 1,
                "Global_Trigger_Time (s)": current_trigger,
                "Leading_Edge (ns)": current_leading,
                "Trailing_Edge (ns)": current_trailing,
                "TOT (ns)": current_trailing - current_leading,
                "Trigger_To_Trigger (s)": trigger_to_trigger,
                "Leading_To_Leading (ns)": leading_to_leading,
                "Trailing_To_Trailing (ns)": trailing_to_trailing,
                "Trigger_To_Leading (ns)": current_leading,
                "Trigger_To_Trailing (ns)": current_trailing,
            }
        )

        current_meta = None
        current_trigger = None
        current_leading = None
        current_trailing = None

    if not events:
        raise ValueError(f"No events were parsed from {input_path}")

    return events


def stringify_rows(
    rows: list[dict[str, Decimal | int | str]], fieldnames: list[str]
) -> list[dict[str, str]]:
    return [
        {fieldname: format_number(row[fieldname]) for fieldname in fieldnames}
        for row in rows
    ]


def write_results(events: list[dict[str, Decimal | int | str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Event_Index",
        "Global_Trigger_Time (s)",
        "Leading_Edge (ns)",
        "Trailing_Edge (ns)",
        "TOT (ns)",
        "Trigger_To_Trigger (s)",
        "Trigger_Jitter (ns)",
        "Leading_To_Leading (ns)",
        "Trailing_To_Trailing (ns)",
        "Trigger_To_Leading (ns)",
        "Trigger_To_Trailing (ns)",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(stringify_rows(events, fieldnames))


def determine_baseline_index(events: list[dict[str, Decimal | int | str]]) -> int:
    trigger_deltas = [
        event["Trigger_To_Trigger (s)"]
        for event in events
        if event["Trigger_To_Trigger (s)"] != "" and event["Trigger_To_Trigger (s)"] > 0
    ]
    if not trigger_deltas:
        return 0

    median_delta = median(trigger_deltas)
    stable_threshold = median_delta * 100

    for index, event in enumerate(events):
        trigger_delta = event["Trigger_To_Trigger (s)"]
        if trigger_delta == "":
            continue

        if 0 < trigger_delta <= stable_threshold:
            return index

    return 0


def trim_startup_events(
    events: list[dict[str, Decimal | int | str]]
) -> list[dict[str, Decimal | int | str]]:
    return events[determine_baseline_index(events):]


def calculate_trigger_delta_median(
    events: list[dict[str, Decimal | int | str]]
) -> Decimal | None:
    trigger_deltas = [
        event["Trigger_To_Trigger (s)"]
        for event in events
        if event["Trigger_To_Trigger (s)"] != ""
    ]
    if not trigger_deltas:
        return None

    return median(trigger_deltas)


def prepare_export_events(
    events: list[dict[str, Decimal | int | str]]
) -> list[dict[str, Decimal | int | str]]:
    prepared_events = [dict(event) for event in events]
    trigger_delta_median = calculate_trigger_delta_median(prepared_events)

    for index, event in enumerate(prepared_events):
        trigger_delta = event["Trigger_To_Trigger (s)"]
        if index == 0:
            event["Trigger_To_Trigger (s)"] = ""
            event["Trigger_Jitter (ns)"] = ""
            event["Leading_To_Leading (ns)"] = ""
            event["Trailing_To_Trailing (ns)"] = ""
            continue

        event["Trigger_Jitter (ns)"] = (
            (trigger_delta - trigger_delta_median) * NANOSECONDS_PER_SECOND
            if trigger_delta_median is not None and trigger_delta != ""
            else ""
        )

    return prepared_events


def build_normalized_events(
    events: list[dict[str, Decimal | int | str]]
) -> list[dict[str, Decimal | int | str]]:
    first_trigger = events[0]["Global_Trigger_Time (s)"]
    normalized_events: list[dict[str, Decimal | int | str]] = []

    for event in events:
        trigger_seconds = event["Global_Trigger_Time (s)"]
        trigger_to_trigger = event["Trigger_To_Trigger (s)"]

        normalized_events.append(
            {
                "Event_Index": int(event["Event_Index"]),
                "Trigger_Offset (us)": (trigger_seconds - first_trigger) * MICROSECONDS_PER_SECOND,
                "Leading_Edge (ns)": event["Leading_Edge (ns)"],
                "Trailing_Edge (ns)": event["Trailing_Edge (ns)"],
                "TOT (ns)": event["TOT (ns)"],
                "Trigger_To_Trigger (us)": (
                    trigger_to_trigger * MICROSECONDS_PER_SECOND
                    if trigger_to_trigger != ""
                    else ""
                ),
                "Trigger_Jitter (ns)": event["Trigger_Jitter (ns)"],
                "Leading_To_Leading (ns)": event["Leading_To_Leading (ns)"],
                "Trailing_To_Trailing (ns)": event["Trailing_To_Trailing (ns)"],
                "Trigger_To_Leading (ns)": event["Trigger_To_Leading (ns)"],
                "Trigger_To_Trailing (ns)": event["Trigger_To_Trailing (ns)"],
            }
        )

    return normalized_events


def write_normalized_results(
    normalized_events: list[dict[str, Decimal | int | str]], output_path: Path
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Event_Index",
        "Trigger_Offset (us)",
        "Leading_Edge (ns)",
        "Trailing_Edge (ns)",
        "TOT (ns)",
        "Trigger_To_Trigger (us)",
        "Trigger_Jitter (ns)",
        "Leading_To_Leading (ns)",
        "Trailing_To_Trailing (ns)",
        "Trigger_To_Leading (ns)",
        "Trigger_To_Trailing (ns)",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(stringify_rows(normalized_events, fieldnames))


def main() -> None:
    parser = argparse.ArgumentParser(description="Digest DLD event data into a CSV summary.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("Dogma_test_Data/FirstTestMarch23_0008.dld.dat"),
        help="Path to the raw or partially cleaned input data file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("Results/results.csv"),
        help="Path to the output CSV file.",
    )
    parser.add_argument(
        "--normalized-output",
        type=Path,
        default=Path("Results/results_normalized.csv"),
        help="Path to the normalized output CSV file.",
    )
    parser.add_argument(
        "--clean-output",
        type=Path,
        default=None,
        help="Optional path to write the cleaned intermediate text file.",
    )
    args = parser.parse_args()

    events = parse_event_file(args.input, args.clean_output)
    stable_events = prepare_export_events(trim_startup_events(events))
    write_results(stable_events, args.output)
    normalized_events = build_normalized_events(stable_events)
    write_normalized_results(normalized_events, args.normalized_output)
    print(f"Wrote {len(stable_events)} events to {args.output}")
    print(f"Wrote {len(normalized_events)} normalized events to {args.normalized_output}")
    if args.clean_output is not None:
        print(f"Wrote cleaned intermediate text to {args.clean_output}")


if __name__ == "__main__":
    main()