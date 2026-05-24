"""Tests for style_lint.py — describe mode + critical-severity rules.

Build this file out incrementally across tasks 12-15.
"""
import json
import os
import subprocess
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, SCRIPTS_DIR)

import style_lint  # noqa: E402  # pyright: ignore[reportMissingImports]


def _run_describe(text: str) -> list[dict]:
    """Run the lint engine on a describe string, return parsed violations."""
    return style_lint.lint_describe(text)


def _severities(violations: list[dict]) -> list[str]:
    return [v["severity"] for v in violations]


def _rules(violations: list[dict]) -> list[str]:
    return [v["rule"] for v in violations]


class TestDescribeModeCritical(unittest.TestCase):
    def test_rainbow_on_continuous_flagged(self):
        v = _run_describe(
            "scatter plot of two continuous variables colored by a third "
            "continuous using jet colormap"
        )
        self.assertIn("critical", _severities(v))
        self.assertIn("anti-patterns#rainbow-on-continuous", _rules(v))

    def test_3d_chart_flagged(self):
        v = _run_describe("3D bar chart of 4 categories with z axis as count")
        self.assertIn("critical", _severities(v))
        self.assertIn("anti-patterns#3d-charts", _rules(v))

    def test_dual_axis_flagged(self):
        v = _run_describe(
            "line chart with dual y-axes, one for revenue and one for clicks"
        )
        self.assertIn("critical", _severities(v))
        self.assertIn("anti-patterns#dual-axes", _rules(v))

    def test_truncated_baseline_flagged(self):
        v = _run_describe(
            "bar chart of company revenue, y-axis starts at 42 to highlight differences"
        )
        self.assertIn("critical", _severities(v))
        self.assertIn("anti-patterns#truncated-baseline", _rules(v))

    def test_missing_axis_labels_flagged(self):
        v = _run_describe("scatter plot with no axis labels")
        self.assertIn("critical", _severities(v))
        self.assertIn("axes#labels-required", _rules(v))

    def test_clean_description_no_critical(self):
        v = _run_describe(
            "dot plot of 4 conditions using Okabe-Ito palette, "
            "y-axis labeled 'Response (counts/min)' from 0 to 100, "
            "x-axis labeled 'Condition'"
        )
        self.assertNotIn("critical", _severities(v))


if __name__ == "__main__":
    unittest.main()
