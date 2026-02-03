# Physical Constants for Bioengineering Calculations

**Last Updated**: 2026-01-28
**Purpose**: Quick reference for commonly used physical constants in bioreactor calculations

---

## Universal Constants

| Constant | Symbol | Value | Units |
|----------|--------|-------|-------|
| **Avogadro's number** | $N_A$ | 6.022 × 10²³ | molecules/mol |
| **Gas constant** | $R$ | 8.314 | J/(mol·K) |
| | | 0.08206 | L·atm/(mol·K) |
| | | 62.36 | L·mmHg/(mol·K) |
| **Boltzmann constant** | $k_B$ | 1.381 × 10⁻²³ | J/K |
| **Faraday constant** | $F$ | 96,485 | C/mol |

---

## Gas Properties (at STP: 0°C, 1 atm)

| Property | Value | Units |
|----------|-------|-------|
| **Molar volume (ideal gas)** | 22.4 | L/mol |
| **Molar volume** | 22,400 | mL/mol |
| **Oxygen density** | 1.429 | g/L |
| **Oxygen molecular weight** | 32.0 | g/mol |

---

## Oxygen Properties

| Property | Value | Conditions | Units |
|----------|-------|------------|-------|
| **Solubility in water** | 1.3 × 10⁻³ | 37°C, 1 atm air | mol/L |
| | 8.2 × 10⁻⁷ | 37°C, 1 atm air | mol/cm³ |
| | 0.23 | 37°C, 1 atm O₂ | mmol/L |
| **Henry's constant (water)** | 1.3 × 10⁻⁶ | 37°C | mol/(L·mmHg) |
| **Diffusion coefficient (water)** | 2.8 × 10⁻⁵ | 25°C | cm²/s |
| | 3.0 × 10⁻⁵ | 37°C | cm²/s |
| **Diffusion coefficient (tissue/gel)** | 1.5-2.5 × 10⁻⁵ | 37°C, hydrogel | cm²/s |

**Conversion**: pO₂ (mmHg) to concentration (mol/L) at 37°C:
$$
C_{O_2} = 1.3 \times 10^{-6} \times pO_2 \text{ (mol/L)}
$$

---

## Cell Properties

| Property | Value | Cell Type | Units |
|----------|-------|-----------|-------|
| **Hepatocyte diameter** | 20-25 | Primary, human | μm |
| **Hepatocyte volume** | 4,000-8,000 | Primary, human | μm³ |
| **Hepatocyte mass** | 3-5 | Primary | ng/cell |
| **Hepatocytes per gram liver** | 120-160 million | Human | cells/g |
| **Hepatocyte density (tissue)** | 120-160 million | In vivo | cells/cm³ |
| **Typical culture density** | 40-80 million | In vitro, 3D | cells/cm³ |

---

## Fluid Properties (37°C)

| Property | Value | Fluid | Units |
|----------|-------|-------|-------|
| **Water density** | 0.993 | Pure water | g/cm³ |
| **Media density** | ~1.00 | Cell culture media | g/cm³ |
| **Water viscosity** | 0.692 | Pure water | cP (mPa·s) |
| **Media viscosity** | 0.7-0.8 | Cell culture media | cP |
| **Blood viscosity** | 3-4 | Whole blood | cP |
| **Plasma viscosity** | 1.2-1.3 | Blood plasma | cP |

---

## Unit Conversions

### Pressure
- 1 atm = 760 mmHg = 101.325 kPa = 14.7 psi
- 1 bar = 0.987 atm = 750 mmHg
- 1 mmHg = 133.3 Pa

### Concentration
- 1 M = 1 mol/L = 1000 mmol/L = 10⁶ μmol/L
- 1 mM = 10⁻³ M = 1 mmol/L
- 1 μM = 10⁻⁶ M = 1 μmol/L

### Volume
- 1 L = 1000 mL = 1000 cm³
- 1 μL = 10⁻⁶ L = 10⁻³ mL = 1 mm³

### Flow Rate
- 1 mL/min = 16.67 μL/s
- 1 L/min = 1000 mL/min = 16.67 mL/s

---

## Useful Conversion Factors

**Oxygen consumption**:
- 1 nmol O₂/s = 60 nmol/min = 3.6 μmol/hr
- 1 μmol O₂/s = 1.344 mL O₂/min (at STP)
- 1 mL O₂/min = 0.744 μmol O₂/s (at STP)

**Cell normalization**:
- "per 10⁶ cells" = "per million cells"
- Typical: Express as nmol/(s·10⁶ cells) or μmol/(hr·10⁶ cells)

**Temperature correction (diffusion coefficients)**:
$$
D_2 = D_1 \times \frac{T_2}{T_1} \times \frac{\mu_1}{\mu_2}
$$
where $T$ is absolute temperature (K) and $\mu$ is viscosity.

---

## Quick Sanity Checks

**Typical oxygen consumption rates**:
- Hepatocyte: ~0.8 nmol/s per 10⁶ cells
- Generic mammalian cell: ~0.1-1 nmol/s per 10⁶ cells
- Whole liver: ~20-25 mL O₂/min per 100 g tissue

**Typical diffusion distances**:
- Oxygen in tissue: ~150-200 μm max before hypoxia
- Glucose in tissue: ~300-500 μm (slower consumption)
- Nutrients in general: O(100 μm) for high-metabolic tissues

**Typical flow rates (bioreactors)**:
- Small-scale (mL): 0.1-10 mL/min
- Clinical-scale (L): 50-500 mL/min
- Dialysis (comparison): 200-800 mL/min

---

**Pro tip**: When in doubt about a calculated value, check against these typical ranges. If your oxygen consumption is 100× higher than expected, check your units!
