"""Tests for figure_spec.py — matplotlib Figure → JSON spec.

Tests use a minimal duck-typed mock so they run without matplotlib.
"""
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, SCRIPTS_DIR)

import figure_spec  # noqa: E402  # pyright: ignore[reportMissingImports]


def make_mock_axis(
    *,
    xlabel="x",
    ylabel="y",
    title="",
    xscale="linear",
    yscale="linear",
    xlim=(0.0, 10.0),
    ylim=(0.0, 100.0),
    line_colors=None,
    has_legend=False,
):
    ax = mock.MagicMock()
    ax.get_xlabel.return_value = xlabel
    ax.get_ylabel.return_value = ylabel
    ax.get_title.return_value = title
    ax.get_xscale.return_value = xscale
    ax.get_yscale.return_value = yscale
    ax.get_xlim.return_value = xlim
    ax.get_ylim.return_value = ylim
    ax.name = "rectilinear"  # not 3D
    # Lines
    lines = []
    for c in (line_colors or []):
        line = mock.MagicMock()
        line.get_color.return_value = c
        lines.append(line)
    ax.get_lines.return_value = lines
    # Legend
    ax.get_legend.return_value = mock.MagicMock() if has_legend else None
    # Patches (bars / pie wedges) — empty by default
    ax.patches = []
    # Images (heatmap content) — empty by default
    ax.images = []
    return ax


class TestExtractSpec(unittest.TestCase):
    def test_extracts_axis_labels(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [
            make_mock_axis(xlabel="Time (h)", ylabel="Response (counts/min)")
        ]
        spec = figure_spec.extract_spec(fig)
        self.assertEqual(spec["axes"][0]["xlabel"], "Time (h)")
        self.assertEqual(spec["axes"][0]["ylabel"], "Response (counts/min)")

    def test_extracts_scale(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [make_mock_axis(yscale="log")]
        spec = figure_spec.extract_spec(fig)
        self.assertEqual(spec["axes"][0]["yscale"], "log")

    def test_extracts_ylim_zero_inclusion(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [make_mock_axis(ylim=(42.0, 51.0))]
        spec = figure_spec.extract_spec(fig)
        self.assertEqual(spec["axes"][0]["ylim"], [42.0, 51.0])
        self.assertFalse(spec["axes"][0]["ylim_includes_zero"])

    def test_extracts_line_colors(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [
            make_mock_axis(line_colors=["#E69F00", "#56B4E9"])
        ]
        spec = figure_spec.extract_spec(fig)
        self.assertEqual(spec["axes"][0]["line_colors"], ["#E69F00", "#56B4E9"])

    def test_detects_legend(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [make_mock_axis(has_legend=True)]
        spec = figure_spec.extract_spec(fig)
        self.assertTrue(spec["axes"][0]["has_legend"])

    def test_n_axes(self):
        fig = mock.MagicMock()
        fig.get_axes.return_value = [make_mock_axis(), make_mock_axis()]
        spec = figure_spec.extract_spec(fig)
        self.assertEqual(spec["n_axes"], 2)

    def test_detects_3d_projection(self):
        fig = mock.MagicMock()
        ax = make_mock_axis()
        ax.name = "3d"
        fig.get_axes.return_value = [ax]
        spec = figure_spec.extract_spec(fig)
        self.assertTrue(spec["axes"][0]["is_3d"])


# --- Picklable duck-typed shims for the CLI test (MagicMock isn't picklable) ---
class _PickAxis:
    def __init__(self, xlabel="x", ylabel="y", title="",
                 xscale="linear", yscale="linear",
                 xlim=(0.0, 10.0), ylim=(0.0, 100.0)):
        self._xlabel = xlabel
        self._ylabel = ylabel
        self._title = title
        self._xscale = xscale
        self._yscale = yscale
        self._xlim = xlim
        self._ylim = ylim
        self.name = "rectilinear"
        self.patches = []
        self.images = []

    def get_xlabel(self): return self._xlabel
    def get_ylabel(self): return self._ylabel
    def get_title(self): return self._title
    def get_xscale(self): return self._xscale
    def get_yscale(self): return self._yscale
    def get_xlim(self): return self._xlim
    def get_ylim(self): return self._ylim
    def get_lines(self): return []
    def get_legend(self): return None


class _PickFig:
    def __init__(self, axes):
        self._axes = axes

    def get_axes(self):
        return self._axes


class TestCliMode(unittest.TestCase):
    def test_cli_with_pickle(self):
        # Build a picklable duck-typed Figure, pickle it, and verify
        # the CLI reads it back. MagicMock objects can't be pickled on
        # some Python versions, so we use a minimal shim instead.
        import pickle
        fig = _PickFig([_PickAxis(xlabel="t", ylabel="y")])
        with tempfile.NamedTemporaryFile(
            suffix=".pkl", delete=False
        ) as pkl:
            pickle.dump(fig, pkl)
            pkl_path = pkl.name
        with tempfile.NamedTemporaryFile(
            suffix=".json", delete=False
        ) as out:
            out_path = out.name
        try:
            rc = figure_spec.cli_main(
                ["--pickle", pkl_path, "--out", out_path]
            )
            self.assertEqual(rc, 0)
            with open(out_path) as f:
                data = json.load(f)
            self.assertEqual(data["axes"][0]["xlabel"], "t")
        finally:
            os.unlink(pkl_path)
            os.unlink(out_path)


if __name__ == "__main__":
    unittest.main()
