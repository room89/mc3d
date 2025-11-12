#!/usr/bin/env python3
"""
Helpers for working with MC3D snapshot files.

Supports both the legacy text (CSV with ';' separator) format and the new
binary format produced when `snapshots_binary=true` is enabled.
"""

from __future__ import annotations

import csv
import struct
from pathlib import Path
from typing import BinaryIO, Iterator, Sequence, Tuple

import pandas as pd

SNAPSHOT_COLUMNS: Tuple[str, ...] = (
    "x",
    "y",
    "z",
    "N",
    "ro",
    "T",
    "vx",
    "vy",
    "vz",
    "E",
)

_MAGIC = b"MC3D"
_HEADER_STRUCT = struct.Struct("<IQd")  # version, record_count, density_reference
_ROW_STRUCT = struct.Struct("<dddIdddddd")
_EXPECTED_VERSION = 1


class SnapshotFormatError(RuntimeError):
    """Raised when a snapshot file cannot be parsed."""


def _load_snapshot_binary(handle: BinaryIO, path: Path) -> pd.DataFrame:
    header = handle.read(_HEADER_STRUCT.size)
    if len(header) != _HEADER_STRUCT.size:
        raise SnapshotFormatError(f"Snapshot header truncated: {path}")
    version, record_count, density_reference = _HEADER_STRUCT.unpack(header)
    if version != _EXPECTED_VERSION:
        raise SnapshotFormatError(
            f"Unsupported snapshot version {version} in {path} "
            f"(expected {_EXPECTED_VERSION})"
        )

    columns = {name: [] for name in SNAPSHOT_COLUMNS}
    for index in range(record_count):
        chunk = handle.read(_ROW_STRUCT.size)
        if len(chunk) != _ROW_STRUCT.size:
            raise SnapshotFormatError(
                f"Snapshot payload truncated after {index} records: {path}"
            )
        (
            x,
            y,
            z,
            particle_count,
            density_value,
            temperature,
            vx,
            vy,
            vz,
            energy,
        ) = _ROW_STRUCT.unpack(chunk)
        columns["x"].append(x)
        columns["y"].append(y)
        columns["z"].append(z)
        columns["N"].append(float(particle_count))
        columns["ro"].append(density_value)
        columns["T"].append(temperature)
        columns["vx"].append(vx)
        columns["vy"].append(vy)
        columns["vz"].append(vz)
        columns["E"].append(energy)

    df = pd.DataFrame(columns)
    df.attrs["density_reference"] = density_reference
    return df


def _load_snapshot_text(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep=";", engine="python")
    except pd.errors.EmptyDataError as exc:
        raise SnapshotFormatError(f"Snapshot file is empty or malformed: {path}") from exc

    missing = [col for col in SNAPSHOT_COLUMNS if col not in df.columns]
    if missing:
        raise SnapshotFormatError(
            f"Snapshot {path} is missing expected columns: {', '.join(missing)}"
        )
    return df[SNAPSHOT_COLUMNS]


def load_snapshot(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Snapshot file not found: {path}")

    with path.open("rb") as handle:
        magic = handle.read(len(_MAGIC))
        if magic == _MAGIC:
            return _load_snapshot_binary(handle, path)

    return _load_snapshot_text(path)


def iter_snapshot_records(path: Path) -> Iterator[Sequence[float]]:
    """
    Yield snapshot rows as numeric sequences in SNAPSHOT_COLUMNS order.

    Values are converted to float (with integer counts cast to float) to match
    the behaviour of the legacy CSV reader.
    """

    with path.open("rb") as handle:
        magic = handle.read(len(_MAGIC))
        if magic == _MAGIC:
            yield from _iter_binary_records(handle, path)
            return

    yield from _iter_text_records(path)


def _iter_binary_records(handle: BinaryIO, path: Path) -> Iterator[Sequence[float]]:
    header = handle.read(_HEADER_STRUCT.size)
    if len(header) != _HEADER_STRUCT.size:
        raise SnapshotFormatError(f"Snapshot header truncated: {path}")
    version, record_count, _density_reference = _HEADER_STRUCT.unpack(header)
    if version != _EXPECTED_VERSION:
        raise SnapshotFormatError(
            f"Unsupported snapshot version {version} in {path} "
            f"(expected {_EXPECTED_VERSION})"
        )

    for index in range(record_count):
        chunk = handle.read(_ROW_STRUCT.size)
        if len(chunk) != _ROW_STRUCT.size:
            raise SnapshotFormatError(
                f"Snapshot payload truncated after {index} records: {path}"
            )
        row = _ROW_STRUCT.unpack(chunk)
        # Convert the particle count to float for compatibility with CSV stats.
        yield (
            float(row[0]),
            float(row[1]),
            float(row[2]),
            float(row[3]),
            float(row[4]),
            float(row[5]),
            float(row[6]),
            float(row[7]),
            float(row[8]),
            float(row[9]),
        )


def _iter_text_records(path: Path) -> Iterator[Sequence[float]]:
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader, None)
        if not header:
            raise SnapshotFormatError(f"No header row found in {path}")
        normalized = tuple(col.strip() for col in header)
        if normalized != SNAPSHOT_COLUMNS:
            raise SnapshotFormatError(
                f"Unexpected header in {path}. Expected "
                f"{';'.join(SNAPSHOT_COLUMNS)}, got {header}"
            )

        for row_idx, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) != len(SNAPSHOT_COLUMNS):
                raise SnapshotFormatError(
                    f"Line {row_idx} in {path} has {len(row)} columns, "
                    f"expected {len(SNAPSHOT_COLUMNS)}"
                )
            try:
                yield tuple(float(cell) for cell in row)
            except ValueError as exc:
                raise SnapshotFormatError(
                    f"Failed to parse row {row_idx} in {path}: {row}"
                ) from exc


__all__ = ["SNAPSHOT_COLUMNS", "SnapshotFormatError", "iter_snapshot_records", "load_snapshot"]

