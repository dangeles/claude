"""Extract a JSON-serializable spec from a matplotlib Figure.

Usable two ways:
  - Import: `from figure_spec import extract_spec; extract_spec(fig)`
  - CLI:    `python3 figure_spec.py --pickle fig.pkl --out fig.json`

The lint script (style_lint.py) consumes the JSON spec via its
--figure-spec mode. This module is designed to be matplotlib-aware
but matplotlib-optional: extract_spec uses duck-typed Figure API
(.get_axes(), ax.get_xlabel(), etc.), so it works on any object
that quacks like a Figure.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys


def extract_spec(fig) -> dict:
    """Walk fig.axes and return a JSON-serializable dict.

    Captures: per-axis labels, scales, limits, line colors, legend
    presence, 3D projection flag, bar/patch count, image presence.
    Does NOT capture: tick labels (variable), font sizes (need
    rcParams), or any image pixels.
    """
    axes_spec = []
    for ax in fig.get_axes():
        xlim = tuple(ax.get_xlim())
        ylim = tuple(ax.get_ylim())
        line_colors = [_normalize_color(line.get_color()) for line in ax.get_lines()]
        axes_spec.append(
            {
                "xlabel": ax.get_xlabel() or "",
                "ylabel": ax.get_ylabel() or "",
                "title": ax.get_title() or "",
                "xscale": ax.get_xscale(),
                "yscale": ax.get_yscale(),
                "xlim": [float(xlim[0]), float(xlim[1])],
                "ylim": [float(ylim[0]), float(ylim[1])],
                "xlim_includes_zero": xlim[0] <= 0.0 <= xlim[1],
                "ylim_includes_zero": ylim[0] <= 0.0 <= ylim[1],
                "is_3d": getattr(ax, "name", "rectilinear") == "3d",
                "line_colors": line_colors,
                "has_legend": ax.get_legend() is not None,
                "n_patches": len(getattr(ax, "patches", []) or []),
                "n_images": len(getattr(ax, "images", []) or []),
            }
        )
    return {
        "n_axes": len(axes_spec),
        "axes": axes_spec,
    }


def _normalize_color(c) -> str:
    """Best-effort color → hex string. Returns '?' if unrecognized."""
    if isinstance(c, str):
        if c.startswith("#"):
            return c.upper()
        # Common named colors that matplotlib accepts
        named = {
            "k": "#000000", "black": "#000000",
            "w": "#FFFFFF", "white": "#FFFFFF",
            "r": "#FF0000", "red": "#FF0000",
            "g": "#008000", "green": "#008000",
            "b": "#0000FF", "blue": "#0000FF",
            "c": "#00FFFF", "cyan": "#00FFFF",
            "m": "#FF00FF", "magenta": "#FF00FF",
            "y": "#FFFF00", "yellow": "#FFFF00",
        }
        return named.get(c.lower(), "?")
    if isinstance(c, (tuple, list)) and len(c) >= 3:
        r, g, b = c[:3]
        # matplotlib uses 0-1 floats; convert to 0-255 hex
        if all(isinstance(v, float) and 0.0 <= v <= 1.0 for v in (r, g, b)):
            return "#{:02X}{:02X}{:02X}".format(
                int(round(r * 255)), int(round(g * 255)), int(round(b * 255))
            )
    return "?"


def cli_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Extract JSON spec from a pickled matplotlib Figure."
    )
    parser.add_argument(
        "--pickle", required=True, help="Path to a pickled matplotlib Figure"
    )
    parser.add_argument(
        "--out", required=True, help="Path to write the JSON spec"
    )
    args = parser.parse_args(argv)
    try:
        with open(args.pickle, "rb") as f:
            fig = pickle.load(f)
    except Exception as e:
        print(f"error: could not load pickle: {e}", file=sys.stderr)
        return 1
    spec = extract_spec(fig)
    with open(args.out, "w") as f:
        json.dump(spec, f, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(cli_main())
