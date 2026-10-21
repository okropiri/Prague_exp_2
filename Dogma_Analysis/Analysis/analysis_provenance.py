from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import sys


def json_value(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_value(item) for item in value]
    return value


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(json_value(value), indent=2, allow_nan=False) + "\n")


def file_identity(path: Path) -> dict[str, str | int]:
    before = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError(f"File changed while hashing: {path}")
    return {"path": str(path.resolve()), "bytes": after.st_size,
            "mtime_ns": after.st_mtime_ns, "sha256": digest.hexdigest()}


def record_arguments(args: argparse.Namespace) -> argparse.Namespace:
    folder = os.environ.get("DOGMA_STEP_REPORT_DIR")
    if folder:
        source = Path(sys.argv[0]).resolve()
        output = Path(folder)
        output.mkdir(parents=True, exist_ok=True)
        write_json(output / f"{source.stem}_{os.getpid()}_settings.json", {
            "recorded_utc": datetime.now(timezone.utc).isoformat(),
            "argv": sys.argv, "cwd": str(Path.cwd()), "python": sys.executable,
            "source": file_identity(source), "effective_arguments": vars(args),
        })
    return args


def record_command(command: list[str]) -> None:
    folder = os.environ.get("DOGMA_STEP_REPORT_DIR")
    if folder:
        with (Path(folder) / "commands.jsonl").open("a") as handle:
            handle.write(json.dumps({"utc": datetime.now(timezone.utc).isoformat(),
                                     "argv": command}, allow_nan=False) + "\n")
