from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
import sys
import time

import numpy as np
import uproot


DEFAULT_TREE_NAME = "CleanedHits"
DEFAULT_METADATA_TREE_NAME = "RunMetadata"
DEFAULT_CHUNK_SIZE = 1_000_000
DEFAULT_MAX_ENTRIES_PER_FILE = 20_000_000
INVALID_RF_FOLD_N = np.iinfo(np.int32).min
INVALID_RF_SEGMENT_ID = -1

HIT_BRANCH_DTYPES: dict[str, np.dtype] = {
    "hit_id": np.dtype(np.uint64),
    "trigger_window_id": np.dtype(np.uint64),
    "global_trigger_time_s": np.dtype(np.float64),
    "tdc_id": np.dtype(np.int32),
    "channel": np.dtype(np.int32),
    "rise_time_trigger_ref_ns": np.dtype(np.float64),
    "fall_time_trigger_ref_ns": np.dtype(np.float64),
    "tot_ns": np.dtype(np.float64),
    "rise_time_ch0_ref_ns": np.dtype(np.float64),
    "fall_time_ch0_ref_ns": np.dtype(np.float64),
    "is_ch0_reference_hit": np.dtype(np.bool_),
    "rf_fold_n": np.dtype(np.int32),
    "rf_phase_ns": np.dtype(np.float64),
    "rf_segment_id": np.dtype(np.int32),
}

METADATA_BRANCH_DTYPES: dict[str, np.dtype] = {
    "rf_period_ns": np.dtype(np.float64),
    "rf_phase_origin_ns": np.dtype(np.float64),
    "rf_segment_count": np.dtype(np.int32),
    "ch0_reference_rise_min_ns": np.dtype(np.float64),
    "ch0_reference_rise_max_ns": np.dtype(np.float64),
    "ch0_reference_tot_min_ns": np.dtype(np.float64),
    "ch0_reference_tot_max_ns": np.dtype(np.float64),
}


@dataclass
class OutputPart:
    path: Path
    entries: int = 0


@dataclass(frozen=True)
class PulseRow:
    window_index: int
    tdc_ordinal: int
    has_global_trigger_seconds: bool
    global_trigger_seconds: float
    channel: int
    rise_tdc_index_one_based: int
    fall_tdc_index_one_based: int
    rise_ns: float
    fall_ns: float
    tot_ns: float


@dataclass(frozen=True)
class RfSegment:
    segment_id: int
    start_global_time_s: float
    end_global_time_s: float
    phase_origin_ns: float


@dataclass(frozen=True)
class RfAssignment:
    period_ns: float
    phase_origin_ns: float
    segment_id: int = INVALID_RF_SEGMENT_ID


@dataclass(frozen=True)
class RfSolution:
    period_ns: float
    phase_origin_ns: float
    summary_file: Path | None = None
    segments: tuple[RfSegment, ...] = ()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Write a ROOT TTree from a cleaned DOGMA pulse-table TSV."
    )
    parser.add_argument("--pulse-table", required=True)
    parser.add_argument("--output-root-file")
    parser.add_argument("--summary-file")
    parser.add_argument("--tree-name", default=DEFAULT_TREE_NAME)
    parser.add_argument("--metadata-tree-name", default=DEFAULT_METADATA_TREE_NAME)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    parser.add_argument("--max-entries-per-file", type=int, default=DEFAULT_MAX_ENTRIES_PER_FILE)
    parser.add_argument("--compression", choices=("zlib", "none"), default="zlib")
    parser.add_argument("--compression-level", type=int, default=1)
    parser.add_argument("--progress-report-interval-seconds", type=float, default=15.0)
    parser.add_argument("--ch0-valid-rise-min-ns", type=float, default=-410.0)
    parser.add_argument("--ch0-valid-rise-max-ns", type=float, default=-395.0)
    parser.add_argument("--ch0-valid-tot-min-ns", type=float, default=16.5)
    parser.add_argument("--ch0-valid-tot-max-ns", type=float, default=19.5)
    parser.add_argument("--rf-summary-file")
    parser.add_argument("--rf-period-ns", type=float)
    parser.add_argument("--rf-phase-origin-ns", type=float)
    return parser.parse_args()


def parse_metadata_line(line: str) -> tuple[str, str] | None:
    if not line.startswith("#") or "=" not in line:
        return None
    key, value = line[1:].split("=", 1)
    return key.strip(), value.strip()


def derive_output_root_file(pulse_table: Path) -> Path:
    name = pulse_table.name
    if name.endswith("_pulses.tsv"):
        return pulse_table.with_name(name[: -len("_pulses.tsv")] + "_pulses.root")
    return pulse_table.with_suffix(".root")


def derive_summary_file(pulse_table: Path) -> Path:
    name = pulse_table.name
    if name.endswith("_pulses.tsv"):
        return pulse_table.with_name(name[: -len("_pulses.tsv")] + "_summary.txt")
    return pulse_table.with_suffix(".summary.txt")


def derive_manifest_file(output_root_file: Path) -> Path:
    return output_root_file.with_name(output_root_file.stem + "_parts.txt")


def derive_part_file(output_root_file: Path, part_index: int) -> Path:
    return output_root_file.with_name(
        f"{output_root_file.stem}_part{part_index:04d}{output_root_file.suffix}"
    )


def read_key_value_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip().strip('"')
    return values


def load_rf_segments(path: Path) -> tuple[RfSegment, ...]:
    segments: list[RfSegment] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or line.startswith("segment_id"):
                continue
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"Expected at least 4 columns in RF segment file {path}, found: {line}")
            segments.append(
                RfSegment(
                    segment_id=int(parts[0]),
                    start_global_time_s=float(parts[1]),
                    end_global_time_s=float(parts[2]),
                    phase_origin_ns=float(parts[3]),
                )
            )
    return tuple(segments)


def resolve_rf_solution(args: argparse.Namespace) -> RfSolution | None:
    manual_period = args.rf_period_ns
    manual_phase_origin = args.rf_phase_origin_ns
    if (manual_period is None) != (manual_phase_origin is None):
        raise ValueError("--rf-period-ns and --rf-phase-origin-ns must be provided together")
    if manual_period is not None:
        if manual_period <= 0.0:
            raise ValueError("--rf-period-ns must be > 0")
        return RfSolution(period_ns=manual_period, phase_origin_ns=manual_phase_origin)

    if not args.rf_summary_file:
        return None

    summary_file = Path(args.rf_summary_file).expanduser().resolve()
    values = read_key_value_file(summary_file)
    try:
        period_ns = float(values["deduced_period_ns"])
        phase_origin_ns = float(values["phase_origin_ns"])
    except KeyError as error:
        raise ValueError(f"Missing RF summary key: {error.args[0]}") from error
    if period_ns <= 0.0:
        raise ValueError("deduced_period_ns from RF summary must be > 0")
    segment_file_text = values.get("rf_segment_output")
    segments: tuple[RfSegment, ...] = ()
    if segment_file_text:
        segment_file = Path(segment_file_text).expanduser().resolve()
        segments = load_rf_segments(segment_file)
    return RfSolution(
        period_ns=period_ns,
        phase_origin_ns=phase_origin_ns,
        summary_file=summary_file,
        segments=segments,
    )


def resolve_rf_assignment(rf_solution: RfSolution | None, global_time_seconds: float) -> RfAssignment | None:
    if rf_solution is None:
        return None
    if not rf_solution.segments or global_time_seconds != global_time_seconds:
        return RfAssignment(
            period_ns=rf_solution.period_ns,
            phase_origin_ns=rf_solution.phase_origin_ns,
            segment_id=0,
        )

    selected_segment = rf_solution.segments[0]
    for segment in rf_solution.segments[1:]:
        if global_time_seconds < segment.start_global_time_s:
            break
        selected_segment = segment
    return RfAssignment(
        period_ns=rf_solution.period_ns,
        phase_origin_ns=selected_segment.phase_origin_ns,
        segment_id=selected_segment.segment_id,
    )


def parse_pulse_row(line: str, pulse_table: Path) -> PulseRow:
    parts = line.split("\t")
    if len(parts) != 9:
        raise ValueError(f"Expected 9 columns in {pulse_table}, found: {line}")

    global_trigger_text = parts[2].strip()
    has_global_trigger = bool(global_trigger_text)
    global_trigger_seconds = float(global_trigger_text) if has_global_trigger else np.nan
    return PulseRow(
        window_index=int(parts[0]),
        tdc_ordinal=int(parts[1]),
        has_global_trigger_seconds=has_global_trigger,
        global_trigger_seconds=global_trigger_seconds,
        channel=int(parts[3]),
        rise_tdc_index_one_based=int(parts[4]),
        fall_tdc_index_one_based=int(parts[5]),
        rise_ns=float(parts[6]),
        fall_ns=float(parts[7]),
        tot_ns=float(parts[8]),
    )


def positive_mod(value: float, period: float) -> float:
    wrapped = math.fmod(value, period)
    return wrapped + period if wrapped < 0.0 else wrapped


def select_ch0_reference(window_rows: list[PulseRow], args: argparse.Namespace) -> PulseRow | None:
    for row in window_rows:
        if row.channel != 0:
            continue
        if row.rise_ns < args.ch0_valid_rise_min_ns or row.rise_ns > args.ch0_valid_rise_max_ns:
            continue
        if row.tot_ns < args.ch0_valid_tot_min_ns or row.tot_ns > args.ch0_valid_tot_max_ns:
            continue
        return row
    return None


def build_compression(name: str, level: int):
    if name == "none":
        return None
    return uproot.ZLIB(level)


def create_chunk_buffers(chunk_size: int) -> dict[str, np.ndarray]:
    return {name: np.empty(chunk_size, dtype=dtype) for name, dtype in HIT_BRANCH_DTYPES.items()}


def slice_chunk(buffers: dict[str, np.ndarray], size: int) -> dict[str, np.ndarray]:
    return {name: values[:size] for name, values in buffers.items()}


class SplitTreeWriter:
    def __init__(
        self,
        output_root_file: Path,
        tree_name: str,
        title: str,
        compression,
        max_entries_per_file: int,
        write_first_part_as_numbered: bool = False,
    ) -> None:
        self.output_root_file = output_root_file
        self.tree_name = tree_name
        self.title = title
        self.compression = compression
        self.max_entries_per_file = max_entries_per_file
        self.write_first_part_as_numbered = write_first_part_as_numbered
        self.parts: list[OutputPart] = []
        self._root_file = None
        self._tree = None
        self._entries_in_current_file = 0
        self._open_new_file()

    def _part_path(self, part_index: int) -> Path:
        if part_index == 1:
            if self.write_first_part_as_numbered:
                return derive_part_file(self.output_root_file, 1)
            return self.output_root_file
        return derive_part_file(self.output_root_file, part_index)

    def _open_new_file(self) -> None:
        part_path = self._part_path(len(self.parts) + 1)
        part_path.parent.mkdir(parents=True, exist_ok=True)
        self._root_file = uproot.recreate(part_path, compression=self.compression)
        self._tree = self._root_file.mktree(self.tree_name, HIT_BRANCH_DTYPES, title=self.title)
        self.parts.append(OutputPart(part_path))
        self._entries_in_current_file = 0

    def _close_current_file(self) -> None:
        if self._root_file is not None:
            self._root_file.close()
        self._root_file = None
        self._tree = None

    def extend(self, chunk: dict[str, np.ndarray]) -> None:
        if not chunk:
            return
        if self._tree is None:
            raise RuntimeError("ROOT tree is not open")
        total_size = len(next(iter(chunk.values())))
        offset = 0
        while offset < total_size:
            if self.max_entries_per_file > 0 and self._entries_in_current_file >= self.max_entries_per_file:
                self._close_current_file()
                self._open_new_file()
            remaining = total_size - offset
            if self.max_entries_per_file > 0:
                capacity = self.max_entries_per_file - self._entries_in_current_file
                take = min(remaining, capacity)
            else:
                take = remaining
            piece = {name: values[offset : offset + take] for name, values in chunk.items()}
            self._tree.extend(piece)
            self._entries_in_current_file += take
            self.parts[-1].entries += take
            offset += take

    def close(self) -> None:
        self._close_current_file()


def normalize_split_part_names(output_root_file: Path, parts: list[OutputPart]) -> Path:
    if len(parts) <= 1:
        return parts[0].path

    first_part = parts[0]
    normalized_first_path = derive_part_file(output_root_file, 1)
    if first_part.path != normalized_first_path:
        if normalized_first_path.exists():
            normalized_first_path.unlink()
        first_part.path.rename(normalized_first_path)
        first_part.path = normalized_first_path
    return normalized_first_path


def should_write_numbered_first_part(output_root_file: Path, max_entries_per_file: int) -> bool:
    if max_entries_per_file <= 0:
        return False
    return (
        derive_part_file(output_root_file, 1).exists()
        or derive_part_file(output_root_file, 2).exists()
        or derive_manifest_file(output_root_file).exists()
    )


def write_manifest_file(manifest_path: Path, tree_name: str, total_entries: int, parts: list[OutputPart]) -> None:
    with manifest_path.open("w", encoding="utf-8") as handle:
        handle.write(f"tree_name={tree_name}\n")
        handle.write(f"total_entries={total_entries}\n")
        handle.write(f"file_count={len(parts)}\n")
        handle.write("part_index\tentries\tpath\n")
        for index, part in enumerate(parts, start=1):
            handle.write(f"{index}\t{part.entries}\t{part.path}\n")


def format_progress(entries_written: int, bytes_read: int, total_bytes: int, started_at: float) -> str:
    elapsed = max(time.monotonic() - started_at, 1e-9)
    rate = entries_written / elapsed
    if total_bytes > 0:
        percent = 100.0 * float(bytes_read) / float(total_bytes)
        return f"[progress] {percent:.2f}% | entries {entries_written} | rate {rate:,.0f} rows/s"
    return f"[progress] entries {entries_written} | rate {rate:,.0f} rows/s"


def insert_summary_fields(summary_path: Path, fields: dict[str, str]) -> None:
    if not summary_path.exists():
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text("".join(f"{key}={value}\n" for key, value in fields.items()), encoding="utf-8")
        return
    content = summary_path.read_text(encoding="utf-8")
    lines = content.splitlines()
    updated_lines: list[str] = []
    updated_keys: set[str] = set()
    for line in lines:
        if line and not line.startswith("#") and "=" in line:
            key = line.split("=", 1)[0].strip()
            if key in fields:
                updated_lines.append(f"{key}={fields[key]}")
                updated_keys.add(key)
                continue
        updated_lines.append(line)
    additions = "".join(f"{key}={value}\n" for key, value in fields.items() if key not in updated_keys)
    if not additions:
        summary_path.write_text("\n".join(updated_lines) + ("\n" if content.endswith("\n") else ""), encoding="utf-8")
        return
    content = "\n".join(updated_lines)
    marker = "\n# How to read this summary\n"
    if marker in content:
        content = content.replace(marker, "\n" + additions + marker, 1)
    else:
        if content and not content.endswith("\n"):
            content += "\n"
        content += additions
    summary_path.write_text(content, encoding="utf-8")


def write_metadata_tree(
    output_root_file: Path,
    metadata_tree_name: str,
    metadata: dict[str, str],
    args: argparse.Namespace,
    rf_solution: RfSolution | None,
) -> None:
    run_key = metadata.get("run_key", output_root_file.stem)
    segment_count = 0
    if rf_solution is not None:
        segment_count = len(rf_solution.segments) if rf_solution.segments else 1
    metadata_arrays = {
        "rf_period_ns": np.array([rf_solution.period_ns if rf_solution is not None else np.nan], dtype=np.float64),
        "rf_phase_origin_ns": np.array([rf_solution.phase_origin_ns if rf_solution is not None else np.nan], dtype=np.float64),
        "rf_segment_count": np.array([segment_count], dtype=np.int32),
        "ch0_reference_rise_min_ns": np.array([args.ch0_valid_rise_min_ns], dtype=np.float64),
        "ch0_reference_rise_max_ns": np.array([args.ch0_valid_rise_max_ns], dtype=np.float64),
        "ch0_reference_tot_min_ns": np.array([args.ch0_valid_tot_min_ns], dtype=np.float64),
        "ch0_reference_tot_max_ns": np.array([args.ch0_valid_tot_max_ns], dtype=np.float64),
    }
    with uproot.update(output_root_file) as root_file:
        metadata_tree = root_file.mktree(
            metadata_tree_name,
            METADATA_BRANCH_DTYPES,
            title=f"{run_key} cleaned-hit metadata",
        )
        metadata_tree.extend(metadata_arrays)


def write_row_to_buffers(
    buffers: dict[str, np.ndarray],
    index: int,
    hit_id: int,
    row: PulseRow,
    ch0_reference_row: PulseRow | None,
    rf_assignment: RfAssignment | None,
) -> None:
    buffers["hit_id"][index] = np.uint64(hit_id)
    buffers["trigger_window_id"][index] = np.uint64(row.window_index)
    buffers["global_trigger_time_s"][index] = row.global_trigger_seconds if row.has_global_trigger_seconds else np.nan
    buffers["tdc_id"][index] = np.int32(row.tdc_ordinal + 1)
    buffers["channel"][index] = np.int32(row.channel)
    buffers["rise_time_trigger_ref_ns"][index] = row.rise_ns
    buffers["fall_time_trigger_ref_ns"][index] = row.fall_ns
    buffers["tot_ns"][index] = row.tot_ns
    if ch0_reference_row is not None:
        ch0_reference_rise_ns = ch0_reference_row.rise_ns
        rise_time_ch0_ref_ns = row.rise_ns - ch0_reference_rise_ns
        fall_time_ch0_ref_ns = row.fall_ns - ch0_reference_rise_ns
        buffers["rise_time_ch0_ref_ns"][index] = rise_time_ch0_ref_ns
        buffers["fall_time_ch0_ref_ns"][index] = fall_time_ch0_ref_ns
        buffers["is_ch0_reference_hit"][index] = row is ch0_reference_row
    else:
        rise_time_ch0_ref_ns = np.nan
        buffers["rise_time_ch0_ref_ns"][index] = np.nan
        buffers["fall_time_ch0_ref_ns"][index] = np.nan
        buffers["is_ch0_reference_hit"][index] = False

    if ch0_reference_row is not None and rf_assignment is not None:
        phase_input_ns = rise_time_ch0_ref_ns - rf_assignment.phase_origin_ns
        buffers["rf_fold_n"][index] = np.int32(math.floor(phase_input_ns / rf_assignment.period_ns))
        buffers["rf_phase_ns"][index] = positive_mod(phase_input_ns, rf_assignment.period_ns)
        buffers["rf_segment_id"][index] = np.int32(rf_assignment.segment_id)
    else:
        buffers["rf_fold_n"][index] = np.int32(INVALID_RF_FOLD_N)
        buffers["rf_phase_ns"][index] = np.nan
        buffers["rf_segment_id"][index] = np.int32(INVALID_RF_SEGMENT_ID)


def write_tree(args: argparse.Namespace) -> tuple[Path, Path, int, dict[str, str], list[OutputPart], Path | None]:
    pulse_table = Path(args.pulse_table).expanduser().resolve()
    output_root_file = (
        Path(args.output_root_file).expanduser().resolve()
        if args.output_root_file
        else derive_output_root_file(pulse_table)
    )
    summary_file = (
        Path(args.summary_file).expanduser().resolve()
        if args.summary_file
        else derive_summary_file(pulse_table)
    )

    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be > 0")
    if args.max_entries_per_file < 0:
        raise ValueError("--max-entries-per-file must be >= 0")
    if args.progress_report_interval_seconds <= 0.0:
        raise ValueError("--progress-report-interval-seconds must be > 0")
    if not args.ch0_valid_rise_min_ns < args.ch0_valid_rise_max_ns:
        raise ValueError("--ch0-valid-rise-min-ns must be < --ch0-valid-rise-max-ns")
    if not args.ch0_valid_tot_min_ns < args.ch0_valid_tot_max_ns:
        raise ValueError("--ch0-valid-tot-min-ns must be < --ch0-valid-tot-max-ns")

    output_root_file.parent.mkdir(parents=True, exist_ok=True)

    metadata: dict[str, str] = {}
    rf_solution = resolve_rf_solution(args)
    total_entries = 0
    next_hit_id = 0
    bytes_read = 0
    last_report_at = time.monotonic()
    started_at = time.monotonic()
    total_bytes = pulse_table.stat().st_size
    compression = build_compression(args.compression, args.compression_level)

    writer = SplitTreeWriter(
        output_root_file=output_root_file,
        tree_name=args.tree_name,
        title=f"{pulse_table.stem} cleaned DOGMA hits",
        compression=compression,
        max_entries_per_file=args.max_entries_per_file,
        write_first_part_as_numbered=should_write_numbered_first_part(output_root_file, args.max_entries_per_file),
    )
    try:
        with pulse_table.open("r", encoding="utf-8") as handle:
            buffers = create_chunk_buffers(args.chunk_size)
            filled = 0
            current_window_index: int | None = None
            current_window_rows: list[PulseRow] = []

            def flush_chunk() -> None:
                nonlocal buffers, filled, total_entries
                if filled == 0:
                    return
                writer.extend(slice_chunk(buffers, filled))
                total_entries += filled
                filled = 0
                buffers = create_chunk_buffers(args.chunk_size)

            def flush_window() -> None:
                nonlocal filled, current_window_rows, next_hit_id
                if not current_window_rows:
                    return
                ch0_reference_row = select_ch0_reference(current_window_rows, args)
                window_global_time_seconds = next(
                    (row.global_trigger_seconds for row in current_window_rows if row.has_global_trigger_seconds),
                    np.nan,
                )
                rf_assignment = resolve_rf_assignment(rf_solution, window_global_time_seconds)
                for row in current_window_rows:
                    write_row_to_buffers(buffers, filled, next_hit_id, row, ch0_reference_row, rf_assignment)
                    next_hit_id += 1
                    filled += 1
                    if filled == args.chunk_size:
                        flush_chunk()
                current_window_rows = []

            for raw_line in handle:
                bytes_read += len(raw_line)
                line = raw_line.rstrip("\n")
                parsed = parse_metadata_line(line)
                if parsed is not None:
                    key, value = parsed
                    metadata[key] = value
                    continue
                if not line or line.startswith("#") or line.startswith("window_index"):
                    continue

                row = parse_pulse_row(line, pulse_table)
                if current_window_index is None:
                    current_window_index = row.window_index
                elif row.window_index != current_window_index:
                    flush_window()
                    current_window_index = row.window_index
                current_window_rows.append(row)

                now = time.monotonic()
                if now - last_report_at >= args.progress_report_interval_seconds:
                    pending_entries = total_entries + filled + len(current_window_rows)
                    print(format_progress(pending_entries, bytes_read, total_bytes, started_at), file=sys.stderr)
                    last_report_at = now

            flush_window()
            flush_chunk()
    finally:
        writer.close()

    primary_output_file = normalize_split_part_names(output_root_file, writer.parts)
    write_metadata_tree(primary_output_file, args.metadata_tree_name, metadata, args, rf_solution)

    manifest_path = None
    if len(writer.parts) > 1:
        manifest_path = derive_manifest_file(output_root_file)
        write_manifest_file(manifest_path, args.tree_name, total_entries, writer.parts)

    insert_summary_fields(
        summary_file,
        {
            "root_tree_output": str(primary_output_file),
            "root_tree_name": args.tree_name,
            "root_tree_entries": str(total_entries),
            "root_tree_file_count": str(len(writer.parts)),
            "root_tree_source_pulse_table": str(pulse_table),
            "root_tree_run_key": metadata.get("run_key", pulse_table.stem),
            "root_tree_metadata_tree_name": args.metadata_tree_name,
            "root_tree_branch_schema": "cleaned_hits_v2",
            "root_tree_rf_fold_n_invalid_sentinel": str(INVALID_RF_FOLD_N),
            "root_tree_rf_segment_id_invalid_sentinel": str(INVALID_RF_SEGMENT_ID),
            **({"root_tree_rf_summary_input": str(rf_solution.summary_file)} if rf_solution and rf_solution.summary_file is not None else {}),
            **({"root_tree_rf_period_ns": str(rf_solution.period_ns)} if rf_solution is not None else {}),
            **({"root_tree_rf_phase_origin_ns": str(rf_solution.phase_origin_ns)} if rf_solution is not None else {}),
            **({"root_tree_rf_segment_count": str(len(rf_solution.segments) if rf_solution.segments else 1)} if rf_solution is not None else {}),
            **({"root_tree_rf_segment_output": read_key_value_file(rf_solution.summary_file).get("rf_segment_output", "")} if rf_solution and rf_solution.summary_file is not None else {}),
            **({"root_tree_manifest_output": str(manifest_path)} if manifest_path is not None else {}),
        },
    )
    return pulse_table, primary_output_file, total_entries, metadata, writer.parts, manifest_path


def main() -> int:
    try:
        args = parse_args()
        pulse_table, output_root_file, total_entries, metadata, parts, manifest_path = write_tree(args)
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    run_key = metadata.get("run_key", pulse_table.stem)
    print(f"Pulse table: {pulse_table}")
    print(f"ROOT tree file: {output_root_file}")
    print(f"Tree name: {args.tree_name}")
    print(f"Run key: {run_key}")
    print(f"Entries written: {total_entries}")
    print(f"ROOT files written: {len(parts)}")
    if manifest_path is not None:
        print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())