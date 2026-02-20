#!/usr/bin/env python3
"""
Publication-ready perturbation profile plotter (IEEE styling).
Author: Watosn
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot perturbation acceleration magnitude vs altitude with publication styling."
    )
    p.add_argument("input_csv", type=Path, help="CSV from perturbation_profile_cli")
    p.add_argument(
        "--output-stem",
        type=Path,
        default=Path("perturbation_vs_altitude"),
        help="Output stem (writes .pdf and .png)",
    )
    p.add_argument(
        "--column",
        choices=("single", "double"),
        default="single",
        help="Figure width: single-column (3.5 in) or double-column (7.16 in)",
    )
    return p.parse_args()


def load_csv(path: Path):
    alt_km = []
    drag = []
    third_body = []
    total = []
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            alt_km.append(float(row["altitude_km"]))
            drag.append(float(row["drag_mps2"]))
            third_body.append(float(row["third_body_mps2"]))
            total.append(float(row["total_mps2"]))
    return alt_km, drag, third_body, total


def apply_ieee_style(column: str) -> tuple[float, float]:
    width_in = 3.5 if column == "single" else 7.16
    height_in = width_in * 0.62
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "font.size": 8.0,
            "axes.labelsize": 8.0,
            "axes.titlesize": 8.0,
            "legend.fontsize": 7.0,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.2,
            "savefig.dpi": 600,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.constrained_layout.use": True,
        }
    )
    return width_in, height_in


def main() -> None:
    args = parse_args()
    alt_km, drag, third_body, total = load_csv(args.input_csv)
    width_in, height_in = apply_ieee_style(args.column)

    fig, ax = plt.subplots(figsize=(width_in, height_in))
    ax.semilogy(alt_km, drag, label="Drag", color="#1f77b4")
    ax.semilogy(alt_km, third_body, label="Third-Body", color="#d62728")
    ax.semilogy(alt_km, total, label="Total", color="#000000", linestyle="--")

    ax.set_xlabel("Altitude [km]")
    ax.set_ylabel("Acceleration Magnitude [m/s$^2$]")
    ax.grid(True, which="major", linestyle="-", linewidth=0.4, alpha=0.35)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.3, alpha=0.25)
    ax.legend(loc="best", frameon=True, framealpha=0.9, edgecolor="0.2")

    out_pdf = args.output_stem.with_suffix(".pdf")
    out_png = args.output_stem.with_suffix(".png")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight")
    print(f"wrote {out_pdf}")
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()

