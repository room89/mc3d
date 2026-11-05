#!/usr/bin/env python3
"""Change a running solver's load mode through an atomically replaced file."""
import argparse
import os
from pathlib import Path
import tempfile


def set_mode(path: Path, mode: str) -> None:
    if mode not in ("full", "balanced", "quiet", "pause"):
        raise ValueError(f"Unknown load mode: {mode}")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent,
                                         prefix=f".{path.name}.", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(mode + "\n")
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Solver working/output directory")
    parser.add_argument("mode", choices=("full", "balanced", "quiet", "pause"))
    parser.add_argument("--file", default="load.control", help="Control file name or absolute path")
    parser.add_argument("--create-parent", action="store_true", help="Create the control file parent directory")
    args = parser.parse_args()
    path = args.directory / args.file
    try:
        if args.create_parent:
            path.parent.mkdir(parents=True, exist_ok=True)
        set_mode(path, args.mode)
    except OSError as error:
        parser.exit(1, f"Cannot write {path}: {error}\n")
    print(f"Requested {args.mode}: {path} (applies when solver next checks the file)")


if __name__ == "__main__":
    main()
