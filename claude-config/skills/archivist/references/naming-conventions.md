# File Naming Conventions Reference

**Last Updated**: 2026-01-28
**Source**: CLAUDE.md project guidelines
**Purpose**: Quick reference for consistent file naming across bioreactor project

---

## Core Principles

1. **Use kebab-case**: Words separated by hyphens, all lowercase
2. **Type prefix**: Start with type identifier (review-, analysis-, etc.)
3. **Descriptive**: Name should indicate content without opening file
4. **Consistent**: Follow patterns within each category

---

## Document Type Prefixes

| Type | Prefix | Example | When to Use |
|------|--------|---------|-------------|
| **Literature review** | `review-` | `review-hepatocyte-oxygen-consumption.md` | Multi-paper synthesis on a specific topic |
| **Technical analysis** | `analysis-` | `analysis-oxygen-delivery-feasibility.md` | Focused deep-dive on one design question |
| **Reference document** | `reference-` | `reference-commercial-hollow-fibers.md` | Spec sheets, constants, vendor info |
| **Paper notes** | `<author>-<year>-` | `wang-2023-hollow-fiber-oxygen.md` | Single-paper summary |
| **Plan/vision** | `<YYYY-MM-DD>-` | `2025-01-25-device-architecture-plan.md` | Strategic planning docs (date-stamped) |
| **Meeting materials** | `<YYYY-MM-DD>-` | `2026-01-27-quarterly-review-deck.md` | Presentations, agendas (date-stamped) |

---

## Detailed Naming Patterns

### Literature Reviews

**Pattern**: `review-<topic>-<subtopic>.md`

**Examples**:
- `review-oxygen-transport-hollow-fiber.md`
- `review-hepatocyte-isolation-methods.md`
- `review-bal-clinical-outcomes-2015-2025.md`

**Guidelines**:
- Start broad, then narrow: `review-[general-topic]-[specific-aspect]`
- Avoid version numbers in filename (use git history or internal version metadata)
- Include date range if review is time-bounded: `review-topic-2015-2025.md`

### Technical Analysis

**Pattern**: `analysis-<question-or-topic>.md`

**Examples**:
- `analysis-cell-loading-requirements.md`
- `analysis-pressure-drop-fiber-bundle.md`
- `analysis-cost-comparison-oxygenation-methods.md`

**Guidelines**:
- Frame as question being analyzed: "What are the cell loading requirements?"
- Or name the system/component: "Pressure drop in fiber bundle"
- Keep focused: analysis should address ONE design question

### Reference Documents

**Pattern**: `reference-<category>-<type>.md`

**Examples**:
- `reference-hollow-fiber-membrane-specs.md`
- `reference-hepatocyte-parameters-human.md`
- `reference-vendor-contact-list.md`

**Guidelines**:
- Organized by what you're looking up, not the source
- Include species/system if multiple versions exist: `reference-ocr-values-human.md`, `reference-ocr-values-porcine.md`

### Paper Notes

**Pattern**: `<first-author-lastname>-<year>-<brief-topic>.md`

**Examples**:
- `wang-2023-oxygen-hollow-fiber.md`
- `kidambi-2009-hepatocyte-oxygen-metabolism.md`
- `gerlach-2003-bal-clinical-trial.md`

**Guidelines**:
- Use first author's last name only (lowercase)
- 4-digit year
- Topic should distinguish from other papers by same author in same year
- If same author published multiple papers same year: `smith-2024-oxygen-transport-a.md`, `smith-2024-oxygen-transport-b.md`

### Plans and Vision Documents

**Pattern**: `<YYYY-MM-DD>-<descriptive-title>.md`

**Examples**:
- `2025-01-25-exo-organ-bioreactor-vision.md`
- `2026-01-15-device-architecture-specification-v2.md`
- `2026-02-01-q1-research-priorities.md`

**Guidelines**:
- ISO 8601 date format (YYYY-MM-DD) for sortability
- Date represents when plan was created, not when it's executed
- Use present/future tense in title: "device architecture" not "architected device"

### Meeting Materials

**Pattern**: `<YYYY-MM-DD>-<meeting-type>-<topic>.md`

**Examples**:
- `2026-01-27-quarterly-review-research-synthesis.md`
- `2026-02-10-design-meeting-oxygen-system.md`
- `2026-03-01-presentation-conference-abstract.md`

**Guidelines**:
- Date = meeting date (not prep date)
- Include meeting type (review, design-meeting, presentation, etc.)
- Topic should identify content

---

## Directory Naming

**Pattern**: `kebab-case`, descriptive

**Examples**:
- `docs/literature/hollow-fiber-membranes/`
- `docs/literature/hepatocyte-biology/`
- `modules/oxygen-transport/`
- `models/fluid-dynamics/`

**Guidelines**:
- Plural for collections: `membranes`, not `membrane`
- Specific over generic: `hollow-fiber-membranes`, not `fibers`
- Maximum 2-3 words if possible

---

## File Extension Standards

| Type | Extension | Notes |
|------|-----------|-------|
| Documentation | `.md` | Markdown (project standard) |
| Python code | `.py` | Scripts, models, analysis |
| Jupyter notebooks | `.ipynb` | Interactive analyses |
| Data (structured) | `.csv`, `.json` | Prefer CSV for tabular, JSON for nested |
| Images | `.png`, `.jpg`, `.svg` | PNG for diagrams, SVG for schematics when possible |
| PDFs (literature) | `.pdf` | Stored in `docs/literature/<topic>/pdfs/` |

---

## Anti-Patterns (DON'T Use These)

### ❌ Version Numbers in Filename
- ❌ `review-oxygen-v2.1.md`
- ✅ `review-oxygen.md` (version tracked in git or metadata block)

**Why**: Leads to proliferation of files; git history is the version system

### ❌ Dates in Middle of Filename
- ❌ `review-2024-oxygen-transport.md`
- ✅ `review-oxygen-transport.md` (date in metadata)

**Why**: Inconsistent sorting; date not the primary classifier

### ❌ Underscores or CamelCase
- ❌ `review_oxygen_transport.md`, `ReviewOxygenTransport.md`
- ✅ `review-oxygen-transport.md`

**Why**: Inconsistent with project standard (kebab-case)

### ❌ Ambiguous Abbreviations
- ❌ `review-o2-hfm-bal.md`
- ✅ `review-oxygen-hollow-fiber-bioreactor.md`

**Why**: Not immediately understandable; searchability reduced

### ❌ Personal Identifiers
- ❌ `review-oxygen-david-v3-final.md`
- ✅ `review-oxygen.md`

**Why**: Documents are project assets, not personal files; git tracks authorship

### ❌ "Final" or "Latest"
- ❌ `analysis-oxygen-final-FINAL-v2.md`
- ✅ `analysis-oxygen.md`

**Why**: Nothing is ever truly final; causes filename clutter

---

## Renaming Existing Files

**When to rename**:
- File doesn't follow conventions (found during archivist audit)
- Scope has changed (paper notes expanded into full review)
- Original name was ambiguous

**How to rename safely**:
1. Check if file is referenced in INDEX.md or other docs
2. Rename file (use `git mv` to preserve history)
3. Update all references in INDEX.md
4. Update cross-references in related documents
5. Commit with message: "Rename [old] to [new] for convention compliance"

**Example**:
```bash
git mv docs/literature/oxygen-review.md docs/literature/review-oxygen-transport-hollow-fiber.md
# Update INDEX.md references
git add docs/INDEX.md
git commit -m "Rename oxygen-review.md to follow naming convention"
```

---

## Quick Decision Tree

**Need to name a new file?**

```
Is it a multi-paper synthesis?
├─ YES → Use `review-<topic>.md`
└─ NO ↓

Is it analyzing one design question?
├─ YES → Use `analysis-<question>.md`
└─ NO ↓

Is it notes on a single paper?
├─ YES → Use `<author>-<year>-<topic>.md`
└─ NO ↓

Is it a specification or lookup table?
├─ YES → Use `reference-<category>.md`
└─ NO ↓

Is it a strategic plan or meeting material?
└─ YES → Use `<YYYY-MM-DD>-<description>.md`
```

---

## References

- **CLAUDE.md**: Section "Repository Organization" → "Naming conventions" table
- **docs/INDEX.md**: Examples of properly named documents in navigation structure

---

**Remember**: Consistent naming is a courtesy to your future self and collaborators. If you can't tell what a file contains from its name alone, the name needs improvement.
