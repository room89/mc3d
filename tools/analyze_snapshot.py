#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

from snapshot_loader import SNAPSHOT_COLUMNS, SnapshotFormatError, iter_snapshot_records


class RunningStats:
    def __init__(self, name: str):
        self.name = name
        self.count = 0
        self.valid_count = 0
        self.nan_count = 0
        self.pos_inf_count = 0
        self.neg_inf_count = 0
        self.sum = 0.0
        self.min = None
        self.max = None

    def update(self, value: float) -> None:
        self.count += 1

        if math.isnan(value):
            self.nan_count += 1
            return

        if math.isinf(value):
            if value > 0:
                self.pos_inf_count += 1
            else:
                self.neg_inf_count += 1
            return

        self.valid_count += 1
        self.sum += value
        if self.min is None or value < self.min:
            self.min = value
        if self.max is None or value > self.max:
            self.max = value

    def as_dict(self) -> dict:
        mean = None
        if self.valid_count:
            mean = self.sum / self.valid_count
        return {
            "name": self.name,
            "count": self.count,
            "valid": self.valid_count,
            "nan": self.nan_count,
            "pos_inf": self.pos_inf_count,
            "neg_inf": self.neg_inf_count,
            "min": self.min,
            "max": self.max,
            "mean": mean,
        }


def summarize_snapshot(path: Path) -> list[dict]:
    stats = [RunningStats(name) for name in SNAPSHOT_COLUMNS]

    try:
        for values in iter_snapshot_records(path):
            for stat, value in zip(stats, values):
                stat.update(float(value))
    except SnapshotFormatError as exc:
        raise RuntimeError(str(exc)) from exc

    return [stat.as_dict() for stat in stats]


def render_summary(summaries: list[dict]) -> str:
    lines = [
        "Column;count;valid;nan;pos_inf;neg_inf;min;max;mean",
    ]
    for summary in summaries:
        lines.append(
            "{name};{count};{valid};{nan};{pos_inf};{neg_inf};{min};{max};{mean}".format(
                name=summary["name"],
                count=summary["count"],
                valid=summary["valid"],
                nan=summary["nan"],
                pos_inf=summary["pos_inf"],
                neg_inf=summary["neg_inf"],
                min=summary["min"],
                max=summary["max"],
                mean=summary["mean"],
            )
        )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute basic statistics for an MC3D snapshot file."
    )
    parser.add_argument(
        "file",
        type=Path,
        help="Path to a snapshot file such as data0.0978857.dat",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional output file for the summary (defaults to stdout).",
    )
    args = parser.parse_args()

    if not args.file.exists():
        raise SystemExit(f"File not found: {args.file}")

    summary = summarize_snapshot(args.file)
    report = render_summary(summary)

    if args.output:
        args.output.write_text(report + "\n", encoding="utf-8")
    else:
        print(report)


if __name__ == "__main__":
    main()

