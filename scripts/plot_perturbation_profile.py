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
    total = []
    component_keys = []
    series = {}
    with path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise RuntimeError("CSV has no header")
        component_keys = [
            k
            for k in reader.fieldnames
            if k.endswith("_mps2") and k not in ("total_mps2",)
        ]
        for k in component_keys:
            series[k] = []
        for row in reader:
            alt_km.append(float(row["altitude_km"]))
            for k in component_keys:
                series[k].append(float(row[k]))
            total.append(float(row["total_mps2"]))
    return alt_km, series, total


def pretty_label(key: str) -> str:
    core = key.replace("_mps2", "")
    parts = core.split("_")
    return " ".join(p.capitalize() for p in parts)


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
    alt_km, series, total = load_csv(args.input_csv)
    width_in, height_in = apply_ieee_style(args.column)

    fig, ax = plt.subplots(figsize=(width_in, height_in))
    for key, vals in series.items():
        ax.semilogy(alt_km, vals, label=pretty_label(key))
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
