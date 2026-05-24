"""Tests for style_lint.py — describe mode + critical-severity rules.

Build this file out incrementally across tasks 12-15.
"""
import os
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


class TestDescribeModeMajor(unittest.TestCase):
    def test_pie_too_many_slices(self):
        v = _run_describe("pie chart with 8 slices showing market share")
        self.assertIn("anti-patterns#pie-misuse", _rules(v))
        # severity major
        for vv in v:
            if vv["rule"] == "anti-patterns#pie-misuse":
                self.assertEqual(vv["severity"], "major")

    def test_dynamite_plot_flagged(self):
        v = _run_describe(
            "bar chart of means with error bars showing standard error "
            "across 5 conditions"
        )
        self.assertIn("anti-patterns#dynamite-plot", _rules(v))

    def test_dynamite_plot_flagged_mean_values_phrasing(self):
        # Common alternate phrasing: "bar chart of mean values with error bars"
        v = _run_describe(
            "bar chart of mean values with error bars across 5 conditions"
        )
        self.assertIn("anti-patterns#dynamite-plot", _rules(v))

    def test_red_green_only_palette_flagged(self):
        v = _run_describe(
            "scatter plot of two groups colored red and green"
        )
        self.assertIn("accessibility#red-green-only", _rules(v))

    def test_alphabet_order_flagged(self):
        v = _run_describe(
            "bar chart of revenue per state sorted alphabetically"
        )
        self.assertIn("anti-patterns#alphabet-order", _rules(v))


class TestFigureSpecMode(unittest.TestCase):
    def _lint_spec(self, spec: dict) -> list[dict]:
        return style_lint.lint_figure_spec(spec)

    def test_missing_axis_labels_critical(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "",
                    "ylabel": "",
                    "title": "Sales",
                    "xscale": "linear",
                    "yscale": "linear",
                    "xlim": [0, 10],
                    "ylim": [0, 100],
                    "xlim_includes_zero": True,
                    "ylim_includes_zero": True,
                    "is_3d": False,
                    "line_colors": ["#1f77b4"],
                    "has_legend": False,
                    "n_patches": 0,
                    "n_images": 0,
                }
            ],
        }
        v = self._lint_spec(spec)
        self.assertIn("axes#labels-required", _rules(v))

    def test_3d_critical(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x",
                    "ylabel": "y",
                    "title": "",
                    "xscale": "linear",
                    "yscale": "linear",
                    "xlim": [0, 10],
                    "ylim": [0, 100],
                    "xlim_includes_zero": True,
                    "ylim_includes_zero": True,
                    "is_3d": True,
                    "line_colors": [],
                    "has_legend": False,
                    "n_patches": 0,
                    "n_images": 0,
                }
            ],
        }
        v = self._lint_spec(spec)
        self.assertIn("anti-patterns#3d-charts", _rules(v))

    def test_truncated_baseline_critical(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x",
                    "ylabel": "y",
                    "title": "",
                    "xscale": "linear",
                    "yscale": "linear",
                    "xlim": [0, 10],
                    "ylim": [42, 51],
                    "xlim_includes_zero": True,
                    "ylim_includes_zero": False,
                    "is_3d": False,
                    "line_colors": [],
                    "has_legend": False,
                    "n_patches": 4,  # bars present
                    "n_images": 0,
                }
            ],
        }
        v = self._lint_spec(spec)
        self.assertIn("anti-patterns#truncated-baseline", _rules(v))

    def test_dual_axes_critical(self):
        spec = {"n_axes": 2, "axes": [
            {"xlabel": "t", "ylabel": "rev", "title": "", "xscale": "linear",
             "yscale": "linear", "xlim": [0, 10], "ylim": [0, 100],
             "xlim_includes_zero": True, "ylim_includes_zero": True,
             "is_3d": False, "line_colors": ["#1f77b4"], "has_legend": False,
             "n_patches": 0, "n_images": 0},
            # Twin axis: matplotlib's ax.twinx() puts another axis with the
            # same xlim but different ylim and a different ylabel
            {"xlabel": "t", "ylabel": "clicks", "title": "", "xscale": "linear",
             "yscale": "linear", "xlim": [0, 10], "ylim": [0, 5000],
             "xlim_includes_zero": True, "ylim_includes_zero": True,
             "is_3d": False, "line_colors": ["#ff7f0e"], "has_legend": False,
             "n_patches": 0, "n_images": 0},
        ]}
        v = self._lint_spec(spec)
        self.assertIn("anti-patterns#dual-axes", _rules(v))

    def test_too_many_categorical_colors_major(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x", "ylabel": "y", "title": "",
                    "xscale": "linear", "yscale": "linear",
                    "xlim": [0, 10], "ylim": [0, 100],
                    "xlim_includes_zero": True, "ylim_includes_zero": True,
                    "is_3d": False,
                    "line_colors": [
                        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                        "#bcbd22", "#17becf",  # 10 colors > 8 limit
                    ],
                    "has_legend": True,
                    "n_patches": 0, "n_images": 0,
                }
            ],
        }
        v = self._lint_spec(spec)
        self.assertIn("color#too-many-categorical", _rules(v))


class TestRedGreenPaletteIsCaught(unittest.TestCase):
    def test_red_green_pair_flagged(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x", "ylabel": "y", "title": "",
                    "xscale": "linear", "yscale": "linear",
                    "xlim": [0, 10], "ylim": [0, 100],
                    "xlim_includes_zero": True, "ylim_includes_zero": True,
                    "is_3d": False,
                    "line_colors": ["#FF0000", "#00FF00"],
                    "has_legend": True,
                    "n_patches": 0, "n_images": 0,
                }
            ],
        }
        v = style_lint.lint_figure_spec(spec)
        self.assertIn("accessibility#colorblind-unsafe", _rules(v))


class TestMinorRules(unittest.TestCase):
    def test_legend_with_few_series_flagged_describe(self):
        v = _run_describe(
            "line plot with 3 series using a legend; "
            "x-axis 'time (s)', y-axis 'voltage (V)'"
        )
        # Minor: legend when ≤5 series — suggest direct labels
        self.assertIn("annotation#legend-when-direct-labels-fit", _rules(v))
        for vv in v:
            if vv["rule"] == "annotation#legend-when-direct-labels-fit":
                self.assertEqual(vv["severity"], "minor")

    def test_legend_with_3_series_flagged_spec(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x", "ylabel": "y", "title": "",
                    "xscale": "linear", "yscale": "linear",
                    "xlim": [0, 10], "ylim": [0, 100],
                    "xlim_includes_zero": True, "ylim_includes_zero": True,
                    "is_3d": False,
                    "line_colors": ["#E69F00", "#56B4E9", "#009E73"],
                    "has_legend": True,
                    "n_patches": 0, "n_images": 0,
                }
            ],
        }
        v = style_lint.lint_figure_spec(spec)
        self.assertIn("annotation#legend-when-direct-labels-fit", _rules(v))

    def test_log_scale_with_narrow_range_flagged(self):
        spec = {
            "n_axes": 1,
            "axes": [
                {
                    "xlabel": "x", "ylabel": "y", "title": "",
                    "xscale": "linear", "yscale": "log",
                    "xlim": [0, 10], "ylim": [10, 30],  # less than 1 decade
                    "xlim_includes_zero": True, "ylim_includes_zero": False,
                    "is_3d": False,
                    "line_colors": ["#1f77b4"],
                    "has_legend": False,
                    "n_patches": 0, "n_images": 0,
                }
            ],
        }
        v = style_lint.lint_figure_spec(spec)
        self.assertIn("axes#log-narrow-range", _rules(v))


class TestStrictFlag(unittest.TestCase):
    def test_strict_exits_2_on_critical(self):
        # Use a separate process via cli_main
        rc = style_lint.cli_main([
            "--describe",
            "3D bar chart of revenue",
            "--strict",
        ])
        self.assertEqual(rc, 2)

    def test_no_strict_exits_0_on_critical(self):
        rc = style_lint.cli_main([
            "--describe",
            "3D bar chart of revenue",
        ])
        self.assertEqual(rc, 0)

    def test_strict_exits_0_on_no_critical(self):
        rc = style_lint.cli_main([
            "--describe",
            "dot plot of 4 conditions with Okabe-Ito palette, "
            "y axis 'Response (counts/min)' from 0 to 100, "
            "x axis 'Condition'",
            "--strict",
        ])
        self.assertEqual(rc, 0)


class TestImageMode(unittest.TestCase):
    def setUp(self):
        import importlib.util
        self.have_pillow = importlib.util.find_spec("PIL") is not None

    def test_missing_image_returns_error_violation(self):
        violations = style_lint.lint_image("/nonexistent/path.png")
        self.assertTrue(any(v["severity"] == "error" for v in violations))

    def test_image_mode_skips_gracefully_without_pillow(self):
        if self.have_pillow:
            self.skipTest("Pillow is installed; cannot test no-Pillow path here")
        violations = style_lint.lint_image("/any/path.png")
        # Should return a single 'error' violation explaining Pillow is missing
        self.assertTrue(any(
            v["severity"] == "error" and "pillow" in v.get("fix", "").lower()
            for v in violations
        ))


if __name__ == "__main__":
    unittest.main()
