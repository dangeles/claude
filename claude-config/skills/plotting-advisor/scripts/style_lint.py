"""Plotting Advisor lint — passive checker for figures and descriptions.

Three input modes (mutually exclusive):
  --describe TEXT       : textual description of a figure
  --figure-spec PATH    : JSON spec produced by figure_spec.extract_spec(fig)
  --image PATH          : a saved PNG/PDF (image-mode lint, optional Pillow)

Default: exit 0 (advisory). --strict exits 2 on any critical violation.

This module is built incrementally:
  Task 12: describe mode + critical rules (this task)
  Task 13: major rules
  Task 14: minor rules + figure-spec mode
  Task 15: image mode + --strict flag wiring
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field, asdict

# Sibling module path setup (palettes added when needed in T13+).
THIS_DIR = __file__.rsplit("/", 1)[0] if "/" in __file__ else "."
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

import palettes  # noqa: E402  # pyright: ignore[reportMissingImports]


SEVERITIES = ("critical", "major", "minor")


@dataclass
class Violation:
    severity: str
    rule: str
    where: str = ""
    observed: object = None
    fix: str = ""

    def to_dict(self) -> dict:
        d = asdict(self)
        # Drop None observed for clean JSON
        if d["observed"] is None:
            d.pop("observed")
        return d


@dataclass
class LintReport:
    violations: list[Violation] = field(default_factory=list)

    def add(self, v: Violation) -> None:
        self.violations.append(v)

    def critical(self) -> list[Violation]:
        return [v for v in self.violations if v.severity == "critical"]

    def major(self) -> list[Violation]:
        return [v for v in self.violations if v.severity == "major"]

    def minor(self) -> list[Violation]:
        return [v for v in self.violations if v.severity == "minor"]


# ----- Describe mode -----

_RAINBOW_PATTERN = re.compile(
    r"\b(jet|gist[\s_-]?rainbow|hsv|rainbow|gist[\s_-]?ncar|nipy[\s_-]?spectral)\b",
    re.IGNORECASE,
)
_3D_PATTERN = re.compile(r"\b3[\s-]?d\b|three[\s-]?dimensional", re.IGNORECASE)
_DUAL_AXIS_PATTERN = re.compile(
    r"dual[\s-]?(y[\s-]?)?ax(e|i)s|two[\s-]?y[\s-]?ax(e|i)s|secondary[\s-]?y",
    re.IGNORECASE,
)
_TRUNCATED_PATTERN = re.compile(
    r"(y[\s-]?axis\s+starts?\s+at\s+(?!0)\d|truncat\w+\s+(y[\s-]?)?ax|baseline\s+at\s+(?!0)\d)",
    re.IGNORECASE,
)
_NO_LABELS_PATTERN = re.compile(
    r"no\s+(axis|axes)\s+labels?|axes?\s+(are\s+)?unlabel\w+|missing\s+(axis|axes)\s+labels?",
    re.IGNORECASE,
)
_CONTINUOUS_HINT_PATTERN = re.compile(
    r"continuous|gradient|colormap|color\s*scale|heatmap|density|magnitude",
    re.IGNORECASE,
)
_PIE_SLICES_PATTERN = re.compile(
    r"pie\s+chart\s+with\s+(\d+)\s+slices?|pie\s+chart.*?(\d+)\s+slices?",
    re.IGNORECASE,
)
_DYNAMITE_PATTERN = re.compile(
    # Catches "bar chart of means with error bars", "bar chart of mean values
    # with error bars", "bar of mean expression with SD error bars", etc.
    r"bar\s+(chart|plot)?\s*of\s+means?(\s+\w+){0,3}\s+with\s+(\w+\s+)?error\s+bars?",
    re.IGNORECASE,
)
_RED_GREEN_PALETTE_PATTERN = re.compile(
    r"colored?\s+red\s+and\s+green|red[\s/]green\s+(palette|colors?|encoding)",
    re.IGNORECASE,
)
_ALPHABET_ORDER_PATTERN = re.compile(
    r"sorted\s+alphabetic\w*|alphabetic\w*\s+order(ed)?",
    re.IGNORECASE,
)
_LEGEND_FEW_SERIES_PATTERN = re.compile(
    r"(\d+)\s+series.*?legend|legend.*?(\d+)\s+series",
    re.IGNORECASE,
)


def lint_describe(text: str) -> list[dict]:
    """Run all describe-mode checks; return list of violation dicts."""
    report = LintReport()

    # --- Critical rules ---

    # Rainbow on continuous
    if _RAINBOW_PATTERN.search(text) and _CONTINUOUS_HINT_PATTERN.search(text):
        report.add(
            Violation(
                severity="critical",
                rule="anti-patterns#rainbow-on-continuous",
                where="palette",
                fix=(
                    "Use viridis / cividis (sequential) or RdBu / PuOr "
                    "(diverging). See references/palettes.md."
                ),
            )
        )

    # 3D charts
    if _3D_PATTERN.search(text):
        report.add(
            Violation(
                severity="critical",
                rule="anti-patterns#3d-charts",
                where="chart_type",
                fix="Use the 2D version. Always.",
            )
        )

    # Dual y-axes
    if _DUAL_AXIS_PATTERN.search(text):
        report.add(
            Violation(
                severity="critical",
                rule="anti-patterns#dual-axes",
                where="axes",
                fix=(
                    "Use two adjacent panels with shared x-axis, OR "
                    "normalize both series to a common unit."
                ),
            )
        )

    # Truncated baseline (look for explicit non-zero start without annotation)
    if _TRUNCATED_PATTERN.search(text) and not re.search(
        r"annotated|broken[\s-]?axis|labeled", text, re.IGNORECASE
    ):
        report.add(
            Violation(
                severity="critical",
                rule="anti-patterns#truncated-baseline",
                where="y_axis",
                fix=(
                    "Extend axis to zero, OR annotate the truncation with "
                    "a broken-axis indicator."
                ),
            )
        )

    # Missing axis labels
    if _NO_LABELS_PATTERN.search(text):
        report.add(
            Violation(
                severity="critical",
                rule="axes#labels-required",
                where="axes",
                fix="Add labels with units in parentheses, e.g. 'Time (h)'.",
            )
        )

    # --- Major rules ---

    m = _PIE_SLICES_PATTERN.search(text)
    if m:
        n_str = next((g for g in m.groups() if g), None)
        if n_str and int(n_str) > 5:
            report.add(
                Violation(
                    severity="major",
                    rule="anti-patterns#pie-misuse",
                    where="chart_type",
                    observed=int(n_str),
                    fix="Use an ordered bar or dot plot. Pie charts require ≤5 slices.",
                )
            )

    if _DYNAMITE_PATTERN.search(text):
        report.add(
            Violation(
                severity="major",
                rule="anti-patterns#dynamite-plot",
                where="chart_type",
                fix=(
                    "Show individual observations — dot plot, box plot, or violin. "
                    "See Weissgerber 2015."
                ),
            )
        )

    if _RED_GREEN_PALETTE_PATTERN.search(text):
        report.add(
            Violation(
                severity="major",
                rule="accessibility#red-green-only",
                where="palette",
                fix=(
                    "Use Okabe-Ito (e.g., #E69F00 orange and #56B4E9 sky blue). "
                    "Or add a shape/line encoding so color isn't load-bearing."
                ),
            )
        )

    if _ALPHABET_ORDER_PATTERN.search(text):
        report.add(
            Violation(
                severity="major",
                rule="anti-patterns#alphabet-order",
                where="x_axis_order",
                fix=(
                    "Sort by value (ascending or descending), unless time or "
                    "experimental order is more meaningful."
                ),
            )
        )

    # --- Minor rules ---

    m = _LEGEND_FEW_SERIES_PATTERN.search(text)
    if m:
        n_str = next((g for g in m.groups() if g), None)
        if n_str and 2 <= int(n_str) <= 5:
            report.add(
                Violation(
                    severity="minor",
                    rule="annotation#legend-when-direct-labels-fit",
                    where="annotation",
                    observed=int(n_str),
                    fix=(
                        "Direct-label each series at its terminal point "
                        "(Tufte). Legend is unnecessary at this count."
                    ),
                )
            )

    return [v.to_dict() for v in report.violations]


def lint_figure_spec(spec: dict) -> list[dict]:
    """Run all rules against a JSON spec produced by figure_spec.extract_spec."""
    report = LintReport()
    axes = spec.get("axes", [])
    n_axes = spec.get("n_axes", len(axes))

    # --- Critical rules (per axis) ---
    for i, ax in enumerate(axes):
        where = f"axes[{i}]"

        # Labels required
        if not ax.get("xlabel") or not ax.get("ylabel"):
            report.add(
                Violation(
                    severity="critical",
                    rule="axes#labels-required",
                    where=where,
                    observed={
                        "xlabel": ax.get("xlabel", ""),
                        "ylabel": ax.get("ylabel", ""),
                    },
                    fix="Add labels with units in parentheses, e.g. 'Time (h)'.",
                )
            )

        # 3D
        if ax.get("is_3d"):
            report.add(
                Violation(
                    severity="critical",
                    rule="anti-patterns#3d-charts",
                    where=where,
                    fix="Use the 2D version. Always.",
                )
            )

        # Truncated baseline (only meaningful for bar/area charts — proxy: any patches)
        if (
            ax.get("yscale") == "linear"
            and not ax.get("ylim_includes_zero")
            and ax.get("n_patches", 0) > 0
        ):
            report.add(
                Violation(
                    severity="critical",
                    rule="anti-patterns#truncated-baseline",
                    where=f"{where}.ylim",
                    observed=ax.get("ylim"),
                    fix=(
                        "Extend axis to zero, OR annotate the truncation "
                        "with a broken-axis indicator."
                    ),
                )
            )

    # Dual y-axes: two axes that share xlim but differ on ylabel/ylim
    if n_axes == 2 and len(axes) == 2:
        a, b = axes[0], axes[1]
        if (
            a.get("xlim") == b.get("xlim")
            and a.get("ylabel") != b.get("ylabel")
            and (a.get("ylabel") and b.get("ylabel"))
            and a.get("ylim") != b.get("ylim")
        ):
            report.add(
                Violation(
                    severity="critical",
                    rule="anti-patterns#dual-axes",
                    where="axes[0,1]",
                    fix=(
                        "Use two adjacent panels with shared x-axis, OR "
                        "normalize both series to a common unit."
                    ),
                )
            )

    # --- Major rules (per axis) ---
    for i, ax in enumerate(axes):
        where = f"axes[{i}]"
        colors = [c for c in ax.get("line_colors", []) if c and c != "?"]

        # Too many categorical colors
        unique_colors = list(dict.fromkeys(colors))  # preserve order, dedupe
        if len(unique_colors) > 8:
            report.add(
                Violation(
                    severity="major",
                    rule="color#too-many-categorical",
                    where=f"{where}.line_colors",
                    observed=len(unique_colors),
                    fix=(
                        "≤8 categorical levels. Collapse categories or use "
                        "small multiples."
                    ),
                )
            )

        # Palette colorblind safety
        if len(unique_colors) >= 2:
            if not palettes.is_colorblind_safe_categorical(unique_colors):
                report.add(
                    Violation(
                        severity="major",
                        rule="accessibility#colorblind-unsafe",
                        where=f"{where}.line_colors",
                        observed=unique_colors,
                        fix=(
                            "Use Okabe-Ito categorical palette, or add a "
                            "redundant encoding (shape/line style)."
                        ),
                    )
                )

    # --- Minor rules (per axis) ---
    for i, ax in enumerate(axes):
        where = f"axes[{i}]"
        colors = [c for c in ax.get("line_colors", []) if c and c != "?"]
        n_series = len(colors)

        # Legend when ≤5 series
        if ax.get("has_legend") and 2 <= n_series <= 5:
            report.add(
                Violation(
                    severity="minor",
                    rule="annotation#legend-when-direct-labels-fit",
                    where=f"{where}.legend",
                    observed=n_series,
                    fix="Direct-label each series at its terminal point (Tufte).",
                )
            )

        # Log scale with narrow range (<2 decades)
        if ax.get("yscale") == "log":
            ylo, yhi = ax.get("ylim", [None, None])
            if ylo and yhi and ylo > 0 and yhi > 0:
                import math
                decades = math.log10(yhi) - math.log10(ylo)
                if decades < 2.0:
                    report.add(
                        Violation(
                            severity="major",
                            rule="axes#log-narrow-range",
                            where=f"{where}.yaxis",
                            observed=round(decades, 2),
                            fix=(
                                "Log scale is appropriate only when range "
                                "spans ≥2 decades, or for multiplicative "
                                "processes. Use linear if not."
                            ),
                        )
                    )
        if ax.get("xscale") == "log":
            xlo, xhi = ax.get("xlim", [None, None])
            if xlo and xhi and xlo > 0 and xhi > 0:
                import math
                decades = math.log10(xhi) - math.log10(xlo)
                if decades < 2.0:
                    report.add(
                        Violation(
                            severity="major",
                            rule="axes#log-narrow-range",
                            where=f"{where}.xaxis",
                            observed=round(decades, 2),
                            fix=(
                                "Log scale is appropriate only when range "
                                "spans ≥2 decades, or for multiplicative "
                                "processes. Use linear if not."
                            ),
                        )
                    )

    return [v.to_dict() for v in report.violations]


def lint_image(path: str) -> list[dict]:
    """Inspect a saved PNG/PDF/JPG for rule violations.

    This is the coarsest mode. Without OCR, we cannot reliably read
    axis labels; we focus on what raw pixels reveal: dominant colors,
    rainbow-palette signature (continuous hue sweep), 3D perspective
    is not detectable from a raster, so we only catch palette issues.

    Returns 'error' severity violations on failure (file missing,
    Pillow unavailable, corrupt image) — never raises.
    """
    report = LintReport()

    try:
        from PIL import Image  # type: ignore
    except ImportError:
        report.add(
            Violation(
                severity="error",
                rule="lint#dependency-missing",
                where="image-mode",
                fix=(
                    "Install Pillow for image-mode lint: pip install Pillow. "
                    "Or use --figure-spec mode (no dependencies)."
                ),
            )
        )
        return [v.to_dict() for v in report.violations]

    import os
    if not os.path.exists(path):
        report.add(
            Violation(
                severity="error",
                rule="lint#file-not-found",
                where=path,
                fix="Check the path. Lint is advisory and never blocks.",
            )
        )
        return [v.to_dict() for v in report.violations]

    try:
        img = Image.open(path).convert("RGB")
    except Exception as e:
        report.add(
            Violation(
                severity="error",
                rule="lint#image-unreadable",
                where=path,
                observed=str(e),
                fix="Pillow could not open the file. Verify it's a valid PNG/JPG/etc.",
            )
        )
        return [v.to_dict() for v in report.violations]

    # Quantize to a small palette to catch dominant colors despite anti-aliasing.
    quantized = img.quantize(colors=16, method=Image.Quantize.FASTOCTREE)
    pal = quantized.getpalette() or []
    # Get color counts
    counts = quantized.getcolors(maxcolors=16) or []
    counts.sort(reverse=True)

    # Convert top non-background colors to hex
    # Heuristic: background is the most-frequent color; skip it.
    dominant_hexes: list[str] = []
    for i, (cnt, idx) in enumerate(counts):
        r, g, b = pal[idx * 3 : idx * 3 + 3]
        hex_code = "#{:02X}{:02X}{:02X}".format(r, g, b)
        # Skip near-white and near-black backgrounds for palette judgment
        is_near_white = r > 240 and g > 240 and b > 240
        is_near_black = r < 20 and g < 20 and b < 20
        if i == 0 and (is_near_white or is_near_black):
            continue
        dominant_hexes.append(hex_code)

    # Rule: rainbow signature — too many hues spread across the spectrum
    # is a soft indicator of a rainbow palette on continuous data.
    # We only flag when there are many distinct hues AND no clear primary.
    if len(dominant_hexes) >= 8:
        # Quick hue spread check
        hues = []
        for hex_code in dominant_hexes[:10]:
            r = int(hex_code[1:3], 16)
            g = int(hex_code[3:5], 16)
            b = int(hex_code[5:7], 16)
            mx = max(r, g, b)
            mn = min(r, g, b)
            if mx == mn:
                continue
            if mx == r:
                h = 60 * ((g - b) / (mx - mn)) % 360
            elif mx == g:
                h = 60 * ((b - r) / (mx - mn)) + 120
            else:
                h = 60 * ((r - g) / (mx - mn)) + 240
            hues.append(h % 360)
        if hues:
            hue_range = max(hues) - min(hues)
            if hue_range > 270:
                report.add(
                    Violation(
                        severity="critical",
                        rule="anti-patterns#rainbow-on-continuous",
                        where="image-palette",
                        observed={
                            "n_dominant_colors": len(dominant_hexes),
                            "hue_range_deg": round(hue_range, 1),
                        },
                        fix=(
                            "Dominant hues span >270° — consistent with a "
                            "rainbow palette. Use viridis / cividis (sequential) "
                            "or RdBu (diverging)."
                        ),
                    )
                )

    # Rule: colorblind-safety check on the top few dominant non-background colors
    if 2 <= len(dominant_hexes) <= 8:
        if not palettes.is_colorblind_safe_categorical(dominant_hexes[:8]):
            report.add(
                Violation(
                    severity="major",
                    rule="accessibility#colorblind-unsafe",
                    where="image-palette",
                    observed=dominant_hexes[:8],
                    fix=(
                        "Dominant palette is not colorblind-safe. Switch to "
                        "Okabe-Ito or Tol bright."
                    ),
                )
            )

    return [v.to_dict() for v in report.violations]


# ----- Render output (used by CLI) -----

def render_markdown(violations: list[dict]) -> str:
    """Render violations as a markdown report (matches SKILL.md spec)."""
    by_sev = {"critical": [], "major": [], "minor": []}
    for v in violations:
        by_sev.setdefault(v["severity"], []).append(v)
    total = len(violations)
    lines = [f"## Plotting Advisor: figure check — {total} issue{'s' if total != 1 else ''} found", ""]
    for sev in SEVERITIES:
        items = by_sev.get(sev, [])
        if not items:
            continue
        lines.append(f"### {sev.capitalize()} ({len(items)})")
        for v in items:
            head = v.get("rule", "<unknown>")
            where = v.get("where", "")
            fix = v.get("fix", "")
            line1 = f"- **{head}**"
            if where:
                line1 += f" — at `{where}`"
            lines.append(line1)
            if fix:
                lines.append(f"  - Fix: {fix}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


# ----- CLI -----

def cli_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Plotting Advisor lint — passive figure checker."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--describe", help="Textual description of the figure")
    group.add_argument(
        "--figure-spec", help="Path to JSON spec from figure_spec.extract_spec()"
    )
    group.add_argument("--image", help="Path to a saved PNG/PDF")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit 2 on any critical violation (default: always exit 0)",
    )
    parser.add_argument(
        "--json", action="store_true", help="Output JSON to stdout instead of markdown"
    )
    args = parser.parse_args(argv)

    if args.describe:
        violations = lint_describe(args.describe)
    elif args.figure_spec:
        try:
            with open(args.figure_spec) as f:
                spec = json.load(f)
        except Exception as e:
            print(f"error: could not read figure-spec JSON: {e}", file=sys.stderr)
            return 1
        violations = lint_figure_spec(spec)
    elif args.image:
        violations = lint_image(args.image)
    else:
        parser.print_help()
        return 1

    if args.json:
        print(json.dumps(violations, indent=2))
    else:
        print(render_markdown(violations), end="")

    if args.strict and any(v["severity"] == "critical" for v in violations):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(cli_main())
