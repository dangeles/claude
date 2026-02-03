# Unit Conversion Reference for Bioreactor Calculations

**Last Updated**: 2026-01-28
**Purpose**: Quick lookup for unit conversions commonly needed in bioengineering calculations

---

## Length

| From | To | Multiply by |
|------|----|-----------|
| meters (m) | centimeters (cm) | 100 |
| centimeters (cm) | millimeters (mm) | 10 |
| millimeters (mm) | micrometers (μm) | 1000 |
| micrometers (μm) | nanometers (nm) | 1000 |
| inches (in) | centimeters (cm) | 2.54 |
| feet (ft) | meters (m) | 0.3048 |

**Quick reference**:
- 1 μm = 10⁻⁶ m = 10⁻⁴ cm
- 1 nm = 10⁻⁹ m = 10⁻⁷ cm
- 1 mm = 1000 μm = 0.1 cm

---

## Volume

| From | To | Multiply by |
|------|----|-----------|
| liters (L) | milliliters (mL) | 1000 |
| milliliters (mL) | microliters (μL) | 1000 |
| liters (L) | cubic centimeters (cm³) | 1000 |
| cubic meters (m³) | liters (L) | 1000 |
| gallons (US) | liters (L) | 3.785 |
| fluid ounces (US) | milliliters (mL) | 29.57 |

**Key identities**:
- 1 mL = 1 cm³ (exact)
- 1 μL = 1 mm³ (exact)
- 1 L = 1 dm³

---

## Mass

| From | To | Multiply by |
|------|----|-----------|
| kilograms (kg) | grams (g) | 1000 |
| grams (g) | milligrams (mg) | 1000 |
| milligrams (mg) | micrograms (μg) | 1000 |
| pounds (lb) | kilograms (kg) | 0.4536 |
| ounces (oz) | grams (g) | 28.35 |

---

## Concentration (Molar)

| From | To | Multiply by |
|------|----|-----------|
| M (mol/L) | mM (mmol/L) | 1000 |
| mM | μM (μmol/L) | 1000 |
| μM | nM (nmol/L) | 1000 |
| M | mol/mL | 0.001 |
| M | mol/cm³ | 0.001 |

**Mass/volume conversions**:
$$
\text{Concentration (M)} = \frac{\text{mass (g/L)}}{\text{molecular weight (g/mol)}}
$$

**Example**: 5 g/L glucose (MW = 180 g/mol):
$$
C = \frac{5 \text{ g/L}}{180 \text{ g/mol}} = 0.0278 \text{ M} = 27.8 \text{ mM}
$$

---

## Pressure

| From | To | Multiply by |
|------|----|-----------|
| atmospheres (atm) | mmHg | 760 |
| atmospheres (atm) | kPa | 101.325 |
| atmospheres (atm) | psi | 14.696 |
| bar | atm | 0.987 |
| mmHg | Pa | 133.322 |
| psi | kPa | 6.895 |

**Biological context**:
- Normal atmospheric pressure = 1 atm = 760 mmHg
- Arterial O₂ (air-breathing) = ~100 mmHg
- Venous O₂ = ~40 mmHg
- "21% O₂" = 0.21 × 760 mmHg = 160 mmHg pO₂

---

## Flow Rate

| From | To | Multiply by |
|------|----|-----------|
| L/min | mL/min | 1000 |
| mL/min | μL/min | 1000 |
| mL/min | μL/s | 16.67 |
| L/min | cm³/s | 16.67 |
| gallons/min (gpm) | L/min | 3.785 |

**Conversion with velocity**:
$$
Q = v \times A
$$
where:
- $Q$ = flow rate (cm³/s or mL/s)
- $v$ = velocity (cm/s)
- $A$ = cross-sectional area (cm²)

---

## Diffusion Coefficients

**Common units**: cm²/s, m²/s, μm²/s

| From | To | Multiply by |
|------|----|-----------|
| m²/s | cm²/s | 10⁴ |
| cm²/s | mm²/s | 100 |
| cm²/s | μm²/s | 10⁸ |

**Typical values**:
- Small molecules (O₂, CO₂) in water: ~2-3 × 10⁻⁵ cm²/s
- Glucose in water: ~6 × 10⁻⁶ cm²/s
- Proteins (albumin) in water: ~6 × 10⁻⁷ cm²/s

---

## Oxygen-Specific Conversions

### pO₂ to Dissolved Concentration (37°C, aqueous)

Using Henry's Law:
$$
C_{O_2} = H \times pO_2
$$

where $H = 1.3 \times 10^{-6}$ mol/(L·mmHg) at 37°C.

**Examples**:
- **Air-equilibrated (160 mmHg)**:
  - $C = 1.3 \times 10^{-6} \times 160 = 2.08 \times 10^{-4}$ mol/L = 0.208 mM = 208 μM

- **Arterial blood (~100 mmHg)**:
  - $C = 1.3 \times 10^{-6} \times 100 = 1.3 \times 10^{-4}$ mol/L = 0.13 mM = 130 μM

- **Hypoxic (~5 mmHg)**:
  - $C = 1.3 \times 10^{-6} \times 5 = 6.5 \times 10^{-6}$ mol/L = 6.5 μM

### Oxygen Consumption Rate Units

**Common format**: nmol/(s·10⁶ cells) or μmol/(hr·10⁶ cells)

**Conversions**:
- 1 nmol/s = 60 nmol/min = 3600 nmol/hr = 3.6 μmol/hr
- 0.8 nmol/(s·10⁶ cells) = 48 nmol/(min·10⁶ cells) = 2.88 μmol/(hr·10⁶ cells)

**Volume-based (if you know cell density)**:
$$
OCR_{volumetric} = OCR_{per\ cell} \times \rho_{cells}
$$

**Example**: 0.8 nmol/(s·10⁶ cells) at 50 × 10⁶ cells/cm³:
$$
OCR = 0.8 \times 50 = 40 \text{ nmol/(s·cm}^3\text{)}
$$

---

## Temperature Conversions

| From | To | Formula |
|------|-----------|---------|
| Celsius (°C) | Kelvin (K) | K = °C + 273.15 |
| Celsius (°C) | Fahrenheit (°F) | °F = (°C × 9/5) + 32 |
| Fahrenheit (°F) | Celsius (°C) | °C = (°F - 32) × 5/9 |

**Common temperatures**:
- Room temperature: 20-25°C = 293-298 K = 68-77°F
- Body/culture temperature: 37°C = 310 K = 98.6°F
- STP (standard temperature): 0°C = 273 K = 32°F

---

## Time

| From | To | Multiply by |
|------|----|-----------|
| hours (hr) | minutes (min) | 60 |
| minutes (min) | seconds (s) | 60 |
| days | hours | 24 |
| days | seconds | 86,400 |

**Rate conversions**:
- per second (s⁻¹) → per minute (min⁻¹): multiply by 60
- per minute (min⁻¹) → per hour (hr⁻¹): multiply by 60

---

## Dimensionless Numbers (Fluid Mechanics)

### Reynolds Number
$$
Re = \frac{\rho v D}{\mu} = \frac{v D}{\nu}
$$
- $\rho$ = density (g/cm³)
- $v$ = velocity (cm/s)
- $D$ = characteristic length (cm)
- $\mu$ = dynamic viscosity (g/(cm·s) = Poise)
- $\nu$ = kinematic viscosity (cm²/s)

**Units check**: All terms dimensionless when consistent units used.

### Péclet Number (mass transfer)
$$
Pe = \frac{v L}{D}
$$
- $v$ = velocity (cm/s)
- $L$ = length scale (cm)
- $D$ = diffusion coefficient (cm²/s)

**Interpretation**:
- Pe << 1: Diffusion dominates
- Pe >> 1: Convection dominates

---

## Common Pitfalls

### Pitfall 1: Mixing mL and cm³
**Problem**: Assuming 1 mL ≠ 1 cm³
**Reality**: 1 mL = 1 cm³ exactly (by definition)

### Pitfall 2: Forgetting 10⁶ in "per million cells"
**Problem**: Using 0.8 nmol/s/cell instead of 0.8 nmol/s per 10⁶ cells
**Fix**: Always write explicitly: "nmol/(s·10⁶ cells)" or "nmol/s/million cells"

### Pitfall 3: Mishandling pO₂ vs. Percentage
**Problem**: Treating "21% O₂" as 21 mmHg
**Reality**: 21% O₂ at 1 atm = 0.21 × 760 mmHg = 160 mmHg

### Pitfall 4: Using Wrong Henry's Constant Temperature
**Problem**: Using H at 25°C for 37°C calculations
**Fix**: Always check temperature. H changes ~10% between 25°C and 37°C

### Pitfall 5: Confusing Diffusion Units (m²/s vs cm²/s)
**Problem**: Using 2.8 × 10⁻⁵ m²/s instead of cm²/s (off by 10⁴!)
**Fix**: Standard bio units are cm²/s. If you see 10⁻⁹ m²/s, that's 10⁻⁵ cm²/s.

---

## Quick Unit Checker

**Before finalizing a calculation**, verify units on both sides of equations:

**Example: Oxygen flux**
$$
J = D \frac{dC}{dx}
$$

**Left side (flux)**:
- Units: mol/(cm²·s)

**Right side**:
- $D$: cm²/s
- $dC/dx$: (mol/cm³)/cm = mol/cm⁴
- Product: cm²/s × mol/cm⁴ = mol/(cm²·s) ✅

**If units don't match, you made an error!**

---

## Useful Tools

**Online converters**:
- WolframAlpha: Great for complex unit conversions
- NIST Guide for SI Units: https://www.nist.gov/pml/weights-and-measures/metric-si/si-units

**Spreadsheet formula for common conversions**: Keep a reference sheet in Excel/Google Sheets with conversion factors for quick lookup during calculations.
