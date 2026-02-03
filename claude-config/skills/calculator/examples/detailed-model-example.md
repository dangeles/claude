# Detailed Model Example: Oxygen Concentration Profile in Hepatocyte-Seeded Gel

**Version**: 1.0
**Date**: 2024-03-25
**Type**: Detailed analytical model with spatial resolution

---

## Question

What is the steady-state oxygen concentration profile in a 300 μm thick hepatocyte-seeded gel layer between two oxygenated surfaces (e.g., between two hollow fibers)? At what cell density does the center become hypoxic (<5 mmHg pO₂)?

**Context**: Refining the back-of-envelope calculation to include spatial oxygen gradients. Need to know maximum allowable cell packing density.

---

## Model Setup

**Geometry**: 1D slab geometry (gel layer of thickness $L$ between two fibers)
- Oxygen diffuses from both surfaces (x = 0 and x = L)
- Hepatocytes distributed uniformly, consuming oxygen at rate $OCR$
- Steady-state (no time dependence)

**Governing equation** (reaction-diffusion):

$$
D_{O_2} \frac{d^2 C}{dx^2} - k \cdot C_0 = 0
$$

Where:
- $D_{O_2}$ = oxygen diffusion coefficient in gel (cm²/s)
- $C(x)$ = oxygen concentration at position $x$ (mol/cm³)
- $k$ = volumetric oxygen consumption rate (mol/cm³/s)
- $C_0$ = cell density (cells/cm³) × $OCR$ (mol/cell/s)

**Boundary conditions**:
- $C(0) = C_{surface}$ (oxygen concentration at fiber surface)
- $C(L) = C_{surface}$ (symmetry; both surfaces supply oxygen)

By symmetry, we only need to solve for $0 \leq x \leq L/2$ with $dC/dx|_{x=L/2} = 0$.

---

## Parameters

| Symbol | Parameter | Value | Source |
|--------|-----------|-------|--------|
| $L$ | Gel thickness | 300 μm = 0.03 cm | Design choice |
| $D_{O_2}$ | O₂ diffusion coef. | 2.0 × 10⁻⁵ cm²/s | [1] |
| $OCR$ | Per-cell O₂ consumption | 0.8 × 10⁻¹⁵ mol/cell/s | [2] |
| $\rho_{cell}$ | Cell density | Variable (10⁷-10⁸ cells/cm³) | To determine |
| $C_{surface}$ | Surface O₂ conc. | 2 × 10⁻⁷ mol/cm³ (~160 mmHg) | [3] |
| $C_{hypoxic}$ | Hypoxia threshold | 6.6 × 10⁻⁹ mol/cm³ (~5 mmHg) | [4] |

---

## Analytical Solution

The solution to the 1D reaction-diffusion equation with symmetric BCs is:

$$
C(x) = C_{surface} - \frac{k}{2D_{O_2}} \left( \frac{L^2}{4} - \left(x - \frac{L}{2}\right)^2 \right)
$$

Where $k = \rho_{cell} \times OCR$ is the volumetric consumption rate.

**Minimum concentration occurs at center** ($x = L/2$):

$$
C_{min} = C_{surface} - \frac{k L^2}{8 D_{O_2}}
$$

Substitute $k = \rho_{cell} \times OCR$:

$$
C_{min} = C_{surface} - \frac{\rho_{cell} \cdot OCR \cdot L^2}{8 D_{O_2}}
$$

---

## Critical Cell Density Calculation

Set $C_{min} = C_{hypoxic}$ and solve for maximum allowable $\rho_{cell}$:

$$
\rho_{cell, max} = \frac{8 D_{O_2} (C_{surface} - C_{hypoxic})}{OCR \cdot L^2}
$$

Plug in values:

$$
\rho_{cell, max} = \frac{8 \times (2.0 \times 10^{-5}) \times (2 \times 10^{-7} - 6.6 \times 10^{-9})}{(0.8 \times 10^{-15}) \times (0.03)^2}
$$

$$
\rho_{cell, max} = \frac{8 \times 2.0 \times 10^{-5} \times 1.93 \times 10^{-7}}{7.2 \times 10^{-19}}
$$

$$
\rho_{cell, max} = \frac{3.1 \times 10^{-11}}{7.2 \times 10^{-19}} \approx 4.3 \times 10^7 \text{ cells/cm}^3
$$

---

## Results

**Maximum cell density to avoid hypoxia**: **~43 million cells/cm³**

**Interpretation**:
- This is **LOWER** than typical hepatocyte densities in liver tissue (~120-160 million/cm³)
- Suggests that 300 μm spacing between fibers is too large at physiological cell densities
- Need either: (1) reduce spacing to ~150-200 μm, or (2) reduce cell density to ~40-50 million/cm³

**Comparison to back-of-envelope**: The simple estimate said "keep cells within 200 μm of oxygen source." This detailed model confirms that constraint quantitatively.

---

## Sensitivity Analysis

| Parameter | 2× increase | 2× decrease | Effect on $\rho_{max}$ |
|-----------|-------------|-------------|----------------------|
| $D_{O_2}$ | +100% | -50% | HIGH |
| $L$ | -75% (L² term) | +300% | **VERY HIGH** |
| $OCR$ | -50% | +100% | HIGH |

**Conclusion**: Spacing ($L$) is the most sensitive parameter due to squared dependence. Halving the spacing (150 μm instead of 300 μm) allows 4× higher cell density.

---

## References

1. Oxygen diffusion coefficient in fibrinogen gel at 37°C
2. Kidambi et al. 2009 (see back-of-envelope example)
3. Equilibrium with 21% O₂ at 37°C, 1 atm
4. Typical hypoxia threshold for hepatocytes

---

**Status**: Complete - detailed model confirms feasibility constraints from simple estimate
