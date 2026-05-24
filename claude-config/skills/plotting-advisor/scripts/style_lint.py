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

    # Major and minor rules added in later tasks
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
        # Task 13/14 will implement this.
        print("error: --figure-spec mode not yet implemented", file=sys.stderr)
        return 1
    elif args.image:
        # Task 15 will implement this.
        print("error: --image mode not yet implemented", file=sys.stderr)
        return 1
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
