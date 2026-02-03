# INDEX.md Entry Example

**Purpose**: Show proper format for INDEX.md entries with cross-references and status tracking

---

## Example Master INDEX.md Structure

```markdown
# Bioreactor Project Documentation Index

**Last Updated**: 2026-01-28
**Documents**: 47 active, 3 in progress, 2 archived

---

## Quick Navigation

- [Literature Reviews](#literature-reviews)
- [Technical Analysis](#technical-analysis)
- [Engineering Plans](#engineering-plans)
- [Research Modules](#research-modules)
- [Meeting Materials](#meeting-materials)

---

## Literature Reviews

### Hollow Fiber Membrane Technology
- **review-hollow-fiber-oxygen-transport.md** (v2.1, Complete) - Oxygen delivery in HFM bioreactors, 15 papers synthesized
  - *Key parameters*: OCR 0.7-0.9 nmol/s/10⁶ cells, max distance 150 μm, membrane permeability values
  - *Related*: [Oxygen Transport Module](#oxygen-transport-module), [HFM Specs](#reference-documents)
  - *Status*: Fact-checked ✅, Editor-polished ✅

- **review-hepatocyte-isolation-methods.md** (v1.0, Complete) - Two-step collagenase perfusion, yield optimization
  - *Key findings*: Porcine yield >2 billion cells/liver, viability >85%, species differences
  - *Related*: [Cell Loading Analysis](#technical-analysis)

### Bioartificial Liver Clinical Studies
- **review-bal-clinical-outcomes.md** (v1.3, In Progress) - Meta-analysis of 12 clinical trials
  - *Current status*: Draft complete, awaiting Devil's Advocate review
  - *Related*: [Clinical Strategy](#engineering-plans)

---

## Technical Analysis

### Cell Loading
- **analysis-cell-loading-requirements.md** (v1.0, Complete) - 10 billion cell target justification
  - *Conclusion*: 10-20 billion hepatocytes needed for 20-30% liver function replacement
  - *Sources*: 8 clinical papers, 3 engineering analyses
  - *Related*: [Hepatocyte Isolation](#literature-reviews), [Device Sizing](#engineering-plans)

### Oxygen Transport
- **analysis-oxygen-delivery-strategies.md** (v2.0, Complete) - Comparison of 4 oxygenation approaches
  - *Options evaluated*: Membrane oxygenation, dual-lumen fibers, perfluorocarbon carriers, gas-permeable walls
  - *Recommendation*: Dual-lumen HFM architecture (best O₂ delivery per volume)
  - *Related*: [HFM Review](#hollow-fiber-membrane-technology), [Detailed O₂ Model](#research-modules)

---

## Engineering Plans

### Vision Documents
- **2025-01-25-exo-organ-bioreactor-vision.md** (v1.0, Complete) - Overarching project strategy
  - *Scope*: Portable BAL for acute liver failure, 3-7 day bridge-to-transplant
  - *Key decisions*: Target 10-20 billion cells, focus on hollow fiber architecture

### Design Specifications
- **2026-01-15-device-architecture-specification.md** (v0.5, In Progress) - Detailed engineering spec
  - *Status*: Awaiting oxygen transport calculations before finalizing dimensions

---

## Research Modules

### Oxygen Transport Module
**Location**: `modules/oxygen-transport/`

**Documents**:
- `oxygen-diffusion-model.md` (Complete) - 1D reaction-diffusion model, critical cell density calculation
- `membrane-mass-transfer.md` (Complete) - Boundary layer effects, mass transfer coefficients
- `system-level-oxygen-budget.md` (In Progress) - Full-device O₂ balance

**Key Results**: Maximum 150 μm spacing between oxygenated surfaces at 50 million cells/cm³

---

## Reference Documents

### Equipment Specifications
- **reference-hollow-fiber-membrane-specs.md** (v1.1, Complete)
  - *Vendors*: Membrana/3M, Spectrum Labs, Repligen
  - *Key specs*: Inner diameter 300-500 μm, wall thickness 30-80 μm, O₂ permeability values
  - *Commercial modules*: 8,000-15,000 fiber bundles available off-shelf

### Biological Constants
- **reference-hepatocyte-parameters.md** (v1.0, Complete)
  - *Compilation*: OCR, albumin secretion, urea synthesis, CYP450 activities across species
  - *Sources*: 20 papers, values cross-validated

---

## Meeting Materials

- **2026-01-27-research-synthesis-presentation.md** - Quarterly progress review deck

---

## Document Status Legend

- **Complete**: Content finalized, reviewed, cited, ready for use
- **In Progress**: Actively being written/revised
- **Placeholder**: Title reserved, not yet started
- **Archived**: Superseded by newer version or no longer relevant

---

## Cross-Reference Guidelines

When adding new documents:
1. Add entry under appropriate section
2. Include version number and status
3. Write 1-sentence description of content/findings
4. Link to related documents in other sections
5. Update document count at top of INDEX.md
6. If creating new subsection, add to Quick Navigation

---

**Maintenance**: Archivist reviews INDEX.md weekly, updates status, adds new documents, removes archived entries
```

---

## Key Formatting Principles

### 1. Hierarchical Structure
- Use `##` for major sections (Literature Reviews, Technical Analysis, etc.)
- Use `###` for subsections (topic areas)
- Use `####` for individual documents (rarely needed)

### 2. Document Entry Format
```
- **filename.md** (vX.Y, Status) - One-sentence description
  - *Key findings/parameters*: Brief bullet points if quantitative
  - *Related*: [Links to related docs with anchor tags]
  - *Status notes*: Workflow status (reviewed, polished, etc.)
```

### 3. Status Indicators
- Use emojis sparingly: ✅ (complete), ⚠️ (needs review), 🚧 (in progress)
- Or use text: (Complete), (In Progress), (Draft), (Archived)

### 4. Cross-References
- Use `[Display Text](#anchor-tag)` for internal links
- Anchor tags auto-generated from headings: `## Hollow Fiber` → `#hollow-fiber`
- Link bidirectionally: if Doc A mentions Doc B, ensure Doc B entry mentions Doc A

### 5. Metadata at Top
- **Last Updated**: ISO date (YYYY-MM-DD)
- **Documents**: Count of active, in-progress, archived files
- Helps track index freshness

---

## Common Mistakes to Avoid

### Mistake 1: Vague Descriptions
❌ `review-oxygen.md (Complete) - About oxygen`
✅ `review-hollow-fiber-oxygen-transport.md (v2.1, Complete) - Oxygen delivery in HFM bioreactors, 15 papers synthesized. Key: OCR 0.7-0.9 nmol/s/10⁶ cells, 150 μm max distance.`

### Mistake 2: Missing Cross-References
❌ Listing document in isolation with no "Related" links
✅ Always add "Related: [Doc A], [Doc B]" when documents inform each other

### Mistake 3: Stale Status
❌ Marking document "In Progress" months after completion
✅ Update status as documents move through workflow (Draft → Reviewed → Complete)

### Mistake 4: Inconsistent Naming
❌ `oxygen_transport_2026.md`, `OxygenTransportReview.md`, `o2-transport.md`
✅ `review-oxygen-transport-hollow-fiber.md` (follows CLAUDE.md convention)

---

**Use this example as a template when creating or updating INDEX.md files.**
