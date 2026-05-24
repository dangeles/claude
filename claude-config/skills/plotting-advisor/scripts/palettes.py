"""Canonical color palettes for plotting-advisor.

Hex constants mirror what is documented in references/palettes.md.
The colorblind-safety helper uses a whitelist of known-safe palettes
plus a coarse RGB-distance fallback. For rigorous CIELAB ΔE under
CVD simulation, callers should use the `colorspacious` library
directly — this module deliberately keeps zero non-stdlib dependencies.
"""
from __future__ import annotations

# Okabe-Ito (Okabe & Ito 2008; Wong 2011 Nature Methods)
OKABE_ITO = [
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
]

# Tol bright (Tol 2018)
TOL_BRIGHT = [
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
]

# Perceptually uniform sequential maps (matplotlib built-ins)
VIRIDIS_NAMES = frozenset({"viridis", "cividis", "plasma", "magma", "inferno"})

# ColorBrewer diverging (commonly used names available in matplotlib)
COLORBREWER_DIVERGING = frozenset(
    {"RdBu", "PuOr", "BrBG", "RdYlBu", "RdYlGn", "Spectral"}
)

# Banned for continuous data (non-monotonic luminance — Borland & Taylor 2007)
RAINBOW_NAMES = frozenset(
    {"jet", "gist_rainbow", "hsv", "rainbow", "gist_ncar", "nipy_spectral"}
)

# Palettes the lint script trusts without re-checking
_KNOWN_SAFE_CATEGORICAL = (
    tuple(OKABE_ITO),
    tuple(TOL_BRIGHT),
)


def _hex_to_rgb(h: str) -> tuple[int, int, int]:
    h = h.lstrip("#")
    if len(h) != 6:
        raise ValueError(f"not a 6-digit hex color: {h!r}")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def _rgb_distance(a: tuple[int, int, int], b: tuple[int, int, int]) -> float:
    """Euclidean RGB distance. Coarse but cheap; used only as a fallback."""
    return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5


def is_colorblind_safe_categorical(
    hex_list: list[str],
    *,
    rgb_distance_tolerance: float = 80.0,
) -> bool:
    """Return True if a categorical palette is colorblind-safe.

    Strategy (in order):
    1. Trivially safe if ≤1 color.
    2. Whitelist match against known-safe palettes (Okabe-Ito, Tol bright).
    3. Fallback: every pair must have RGB distance ≥ tolerance.
       This is a coarse approximation — it catches obviously bad pairs
       (pure red vs pure green at full saturation = distance ≈ 360,
       BUT under deuteranopia they appear similar). The tolerance is
       deliberately permissive of low-saturation pairs; the lint
       script's primary path is the whitelist, with this fallback
       only flagging clear failures.

    For rigorous checks, callers should use `colorspacious` with a
    CVD transform; this helper is for the zero-dependency advisory
    case.
    """
    if len(hex_list) <= 1:
        return True

    normalized = tuple(c.upper() for c in hex_list)
    for safe in _KNOWN_SAFE_CATEGORICAL:
        safe_upper = tuple(c.upper() for c in safe)
        # Subset check — any subset of a known-safe palette is safe
        if all(c in safe_upper for c in normalized):
            return True

    # Fallback: pairwise distance + red/green saturation check
    rgbs = [_hex_to_rgb(c) for c in hex_list]
    for i in range(len(rgbs)):
        for j in range(i + 1, len(rgbs)):
            ri, gi, bi = rgbs[i]
            rj, gj, bj = rgbs[j]
            # Flag explicit red/green-only confusion: one color
            # dominated by R channel, the other by G channel, both
            # at high saturation.
            high_sat_red_i = ri > 200 and gi < 80 and bi < 80
            high_sat_green_j = gj > 200 and rj < 80 and bj < 80
            high_sat_red_j = rj > 200 and gj < 80 and bj < 80
            high_sat_green_i = gi > 200 and ri < 80 and bi < 80
            if (high_sat_red_i and high_sat_green_j) or (
                high_sat_red_j and high_sat_green_i
            ):
                return False
            if _rgb_distance(rgbs[i], rgbs[j]) < rgb_distance_tolerance:
                return False
    return True
