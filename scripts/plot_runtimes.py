#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description="Plot blossom-vi and Blossom V runtimes for all benchmark families."
    )
    parser.add_argument(
        "--blossom-vi-root",
        type=Path,
        default=repo_root / "runtimes" / "blossom-vi",
        help="Root directory containing blossom-vi family runtime CSVs.",
    )
    parser.add_argument(
        "--blossom-v-root",
        type=Path,
        default=repo_root / "runtimes" / "blossom-v",
        help="Root directory containing Blossom V family runtime CSVs.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=repo_root / "runtimes" / "runtime-comparison.png",
        help="Output image path.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show the plot interactively after saving it.",
    )
    return parser


def read_runtime_table(path: Path, runtime_column: str) -> list[tuple[int, float]]:
    rows: list[tuple[int, float]] = []
    with path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            rows.append((int(row["m"]), float(row[runtime_column])))
    rows.sort()
    return rows


def plot_family(ax, family_name: str, blossom_vi_rows: list[tuple[int, float]], blossom_v_rows: list[tuple[int, float]]) -> None:
    ax.plot(
        [m for m, _ in blossom_vi_rows],
        [runtime for _, runtime in blossom_vi_rows],
        marker="o",
        markersize=3,
        linewidth=1.5,
        label="Blossom VI",
    )
    ax.plot(
        [m for m, _ in blossom_v_rows],
        [runtime for _, runtime in blossom_v_rows],
        marker="s",
        markersize=3,
        linewidth=1.5,
        label="Blossom V",
    )
    ax.set_xlabel("m")
    ax.set_ylabel("time [s]")
    ax.set_title(family_name)
    ax.grid(True, alpha=0.3)
    ax.legend()


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required to plot the runtime comparison. "
            "Install it and rerun the script."
        ) from exc

    families = (
        "dense-random",
        "sparse-random",
        "delaunay-big-weights",
        "delaunay-small-weights",
        "geometric-big-weights",
        "geometric-small-weights",
    )
    fig, axes = plt.subplots(3, 2, figsize=(16, 15))

    for ax, family_name in zip(axes.flat, families):
        blossom_vi_path = args.blossom_vi_root / family_name / "runtimes.csv"
        blossom_v_path = args.blossom_v_root / family_name / "runtimes.csv"
        blossom_vi_rows = read_runtime_table(blossom_vi_path, "runtime_seconds")
        blossom_v_rows = read_runtime_table(blossom_v_path, "runtime_seconds")
        plot_family(ax, family_name, blossom_vi_rows, blossom_v_rows)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
