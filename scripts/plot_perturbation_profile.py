#!/usr/bin/env python3
"""
Publication-ready perturbation profile plotter (IEEE styling).
Author: Watosn
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


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
    p.add_argument(
        "--cross-alt-km",
        dest="cross_alt_km",
        nargs="+",
        type=float,
        default=[],
        help="Altitude marker(s) [km] to draw as vertical lines and mark all curve intersections.",
    )
    p.add_argument(
        "--split-central-gravity",
        action="store_true",
        help="Draw gravity_central in a separate top panel and all other terms in a lower panel.",
    )
    p.add_argument(
        "--y-min",
        type=float,
        default=None,
        help="Override y-axis minimum (log scale).",
    )
    p.add_argument(
        "--y-max",
        type=float,
        default=None,
        help="Override y-axis maximum (log scale).",
    )
    p.add_argument(
        "--focus-noncentral",
        action="store_true",
        default=True,
        help="Auto-focus y-axis on non-central terms when gravity_central exists in single-panel mode.",
    )
    p.add_argument(
        "--monochrome",
        action="store_true",
        default=False,
        help="Use black lines with distinct line styles and markers for color-blind-friendly plots.",
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
    if key == "gravity_sph_tides_mps2":
        return "Gravity SPH + Tides"
    if key == "gravity_tides_mps2":
        return "Gravity Tides"
    if key == "gravity_sph_mps2":
        return "Gravity SPH"
    core = key.replace("_mps2", "")
    parts = core.split("_")
    acronyms = {"srp": "SRP", "earth_radiation": "Earth Radiation"}
    return " ".join(acronyms.get(p, p.capitalize()) for p in parts)


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


def interpolate_y(alt_km: list[float], vals: list[float], x_km: float) -> float | None:
    if not alt_km or len(alt_km) != len(vals):
        return None
    if x_km < alt_km[0] or x_km > alt_km[-1]:
        return None
    i = bisect.bisect_left(alt_km, x_km)
    if i == 0:
        y = vals[0]
        return y if (y > 0.0 and math.isfinite(y)) else None
    if i >= len(alt_km):
        y = vals[-1]
        return y if (y > 0.0 and math.isfinite(y)) else None
    if alt_km[i] == x_km:
        y = vals[i]
        return y if (y > 0.0 and math.isfinite(y)) else None
    x0, x1 = alt_km[i - 1], alt_km[i]
    y0, y1 = vals[i - 1], vals[i]
    if not (y0 > 0.0 and y1 > 0.0 and math.isfinite(y0) and math.isfinite(y1)):
        return None
    t = (x_km - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)


def main() -> None:
    args = parse_args()
    alt_km, series, total = load_csv(args.input_csv)
    series_raw = {k: list(v) for k, v in series.items()}

    # Omit central gravity from comparative perturbation plots to keep scale readable.
    if "gravity_central_mps2" in series:
        del series["gravity_central_mps2"]

    # Aggregate all gravity tide sub-terms into one "Gravity Tides" curve
    # and suppress individual tide lines to reduce clutter.
    tide_keys = [k for k in list(series.keys()) if k.startswith("gravity_tide_")]
    if tide_keys:
        n = len(alt_km)
        vals = []
        for i in range(n):
            ss = 0.0
            for key in tide_keys:
                v = series[key][i]
                if v > 0.0 and math.isfinite(v):
                    ss += v * v
            vals.append(math.sqrt(ss) if ss > 0.0 else float("nan"))
        series["gravity_tides_mps2"] = vals
        for key in tide_keys:
            del series[key]
        # Avoid double-counting: if SPH+tides is present, prefer that combined
        # representation and suppress the separate aggregated tides curve.
        if "gravity_sph_tides_mps2" in series:
            del series["gravity_tides_mps2"]

    width_in, height_in = apply_ieee_style(args.column)

    marker_cycle = ["o", "s", "^", "D", "v", "P", "X", "*", "<", ">"]

    if args.split_central_gravity and "gravity_central_mps2" in series:
        fig, (ax_top, ax) = plt.subplots(
            2,
            1,
            figsize=(width_in, height_in * 1.28),
            sharex=True,
            gridspec_kw={"height_ratios": [1, 2]},
        )
        plot_order = [("gravity_central_mps2", series["gravity_central_mps2"], ax_top)] + [
            (k, v, ax) for k, v in series.items() if k != "gravity_central_mps2"
        ]
    else:
        fig, ax = plt.subplots(figsize=(width_in, height_in))
        ax_top = None
        plot_order = [(k, v, ax) for k, v in series.items()]

    color_cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"])
    lines = []
    for i, (key, vals, axis) in enumerate(plot_order):
        mk = marker_cycle[i % len(marker_cycle)]
        mark_every = max(1, len(alt_km) // 16)
        kwargs = {
            "linestyle": "-",
            "marker": mk,
            "markevery": mark_every,
            "markersize": 3.2,
            "markeredgewidth": 0.6,
            "markerfacecolor": "black",
        }
        kwargs["color"] = "black"
        line = axis.semilogy(alt_km, vals, label=pretty_label(key), **kwargs)[0]
        lines.append((pretty_label(key), line, vals, mk, axis))

    if ax_top is not None:
        ax_top.set_ylabel("Accel [m/s$^2$]")
        ax_top.grid(True, which="major", linestyle="-", linewidth=0.4, alpha=0.35)
        ax_top.grid(True, which="minor", linestyle=":", linewidth=0.3, alpha=0.25)
        ax.set_xlabel("Altitude [km]")
    else:
        ax.set_xlabel("Altitude [km]")
    ax.set_ylabel("Acceleration Magnitude [m/s$^2$]")
    ax.grid(True, which="major", linestyle="-", linewidth=0.4, alpha=0.35)
    ax.grid(True, which="minor", linestyle=":", linewidth=0.3, alpha=0.25)

    # Single-panel scaling control.
    if ax_top is None:
        y_min = args.y_min
        y_max = args.y_max
        if y_min is None or y_max is None:
            values = []
            for key, vals in series.items():
                for v in vals:
                    if v > 0.0 and math.isfinite(v):
                        values.append(v)
            if values:
                auto_min = min(values) / 2.0
                auto_max = max(values) * 3.0
                y_min = auto_min if y_min is None else y_min
                y_max = auto_max if y_max is None else y_max
        if y_min is not None and y_max is not None and y_min > 0.0 and y_max > y_min:
            ax.set_ylim(y_min, y_max)
    color_by_label = {label: line.get_color() for label, line, _, _, _ in lines}
    readout_series = {
        k: v for k, v in series_raw.items() if k != "gravity_central_mps2"
    }
    cross_values: list[tuple[float, list[tuple[str, str, float]]]] = []
    for x_km in args.cross_alt_km:
        if x_km < alt_km[0] or x_km > alt_km[-1]:
            continue
        vals_at_cross: list[tuple[str, str, float]] = []
        # Marker overlay on plotted curves.
        for _, line, vals, mk, axis in lines:
            y = interpolate_y(alt_km, vals, x_km)
            if y is None:
                continue
            axis.plot([x_km], [y], marker=mk, markersize=3.0, color=line.get_color(), linestyle="None")
        # Numeric readout from raw components (exclude central gravity).
        for key, vals in readout_series.items():
            y = interpolate_y(alt_km, vals, x_km)
            if y is None:
                continue
            label = pretty_label(key)
            vals_at_cross.append((label, color_by_label.get(label, "black"), y))
        vals_at_cross.sort(key=lambda t: t[2], reverse=True)
        cross_values.append((x_km, vals_at_cross))

    # Sort legend entries by first (starting-altitude) acceleration magnitude.
    # Largest starting acceleration appears first in the legend.
    lines_for_legend = sorted(
        lines,
        key=lambda item: item[2][0] if (item[2] and math.isfinite(item[2][0])) else float("-inf"),
        reverse=True,
    )

    handles = []
    labels = []
    for name, line, _, mk, _ in lines_for_legend:
        handles.append(
            Line2D(
                [0],
                [0],
                color=line.get_color(),
                linestyle="-",
                marker=mk,
                markersize=4.0,
                linewidth=1.2,
                markerfacecolor="black",
                markeredgewidth=0.7,
            )
        )
        labels.append(name)
    legend = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        ncol=1,
        borderaxespad=0.0,
        frameon=True,
        framealpha=0.95,
        edgecolor="0.2",
        fontsize=5.8,
        handlelength=1.4,
        handletextpad=0.5,
        columnspacing=0.8,
    )
    for text in legend.get_texts():
        text.set_color("black")

    # Add crossing-value readout outside axes for precise numeric lookup.
    if cross_values:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        lb = legend.get_window_extent(renderer=renderer)
        p0, p1 = ax.transAxes.inverted().transform([[lb.x0, lb.y0], [lb.x1, lb.y1]])
        y0 = min(p0[1] - 0.02, 0.56)
        for x_km, vals in cross_values:
            ax.text(1.02, y0, f"@ {x_km:.0f} km", transform=ax.transAxes, fontsize=6.0, ha="left", va="top")
            y = y0 - 0.055
            for name, color, val in vals:
                ax.text(
                    1.02,
                    y,
                    f"{name}: {val:.3e}",
                    color=color,
                    transform=ax.transAxes,
                    fontsize=5.8,
                    ha="left",
                    va="top",
                )
                y -= 0.048
            y0 = y - 0.045

    out_pdf = args.output_stem.with_suffix(".pdf")
    out_png = args.output_stem.with_suffix(".png")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, bbox_inches="tight")
    print(f"wrote {out_pdf}")
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
