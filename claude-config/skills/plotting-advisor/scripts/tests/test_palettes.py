"""Tests for palettes.py — canonical hex constants and colorblind-safety helper."""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, SCRIPTS_DIR)

import palettes  # noqa: E402


class TestPaletteConstants(unittest.TestCase):
    def test_okabe_ito_has_8_colors(self):
        self.assertEqual(len(palettes.OKABE_ITO), 8)

    def test_okabe_ito_all_hex(self):
        for c in palettes.OKABE_ITO:
            self.assertRegex(c, r"^#[0-9A-Fa-f]{6}$")

    def test_tol_bright_has_7_colors(self):
        self.assertEqual(len(palettes.TOL_BRIGHT), 7)

    def test_rainbow_names_includes_jet(self):
        self.assertIn("jet", palettes.RAINBOW_NAMES)
        self.assertIn("gist_rainbow", palettes.RAINBOW_NAMES)
        self.assertIn("hsv", palettes.RAINBOW_NAMES)

    def test_viridis_names_includes_viridis_and_cividis(self):
        self.assertIn("viridis", palettes.VIRIDIS_NAMES)
        self.assertIn("cividis", palettes.VIRIDIS_NAMES)


class TestColorblindSafeHelper(unittest.TestCase):
    def test_okabe_ito_is_safe(self):
        self.assertTrue(palettes.is_colorblind_safe_categorical(palettes.OKABE_ITO))

    def test_red_green_pair_is_not_safe(self):
        # Pure red + pure green is the classic red/green confusion case
        self.assertFalse(
            palettes.is_colorblind_safe_categorical(["#FF0000", "#00FF00"])
        )

    def test_single_color_is_trivially_safe(self):
        self.assertTrue(palettes.is_colorblind_safe_categorical(["#444444"]))

    def test_empty_palette_is_trivially_safe(self):
        self.assertTrue(palettes.is_colorblind_safe_categorical([]))


if __name__ == "__main__":
    unittest.main()
