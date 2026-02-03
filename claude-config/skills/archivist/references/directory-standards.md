# Directory Structure Standards

**Last Updated**: 2026-01-28
**Source**: CLAUDE.md project guidelines
**Purpose**: Repository organization principles and directory hierarchy

---

## Repository Root Structure

```
bioreactor/
├── docs/                   # All documentation
│   ├── INDEX.md           # Master navigation
│   ├── WORK-LOG.md        # Agent task tracking
│   ├── literature/        # Research papers and reviews
│   ├── reports/           # Technical reports
│   ├── plans/             # Vision and strategy documents
│   ├── meeting-materials/ # Presentations, agendas
│   └── presentations/     # Conference/external talks
├── modules/               # Engineering subsystems
│   ├── oxygen-transport/  # Subsystem documentation + models
│   ├── cell-biology/      # Cell-related engineering work
│   └── fluid-dynamics/    # Flow modeling
├── models/                # Mathematical models (Python, notebooks)
├── src/                   # Python code (reusable functions)
├── experiments/           # Lab protocols, data, analysis
├── scratchpad/            # Exploratory work (not git-tracked)
├── .claude/               # User-level agent config (not part of bioreactor project)
├── CLAUDE.md              # Project guidelines
└── README.md              # Repository overview
```

---

## Core Directories

### `docs/` - All Documentation

**Purpose**: Human-readable documents in Markdown format

**Subdirectories**:

**`docs/literature/`**
- Topic-based subdirectories: `hollow-fiber-membranes/`, `hepatocyte-biology/`, etc.
- Each topic directory contains:
  - `README.md` - Overview of topic area
  - `review-*.md` - Literature reviews
  - `<author>-<year>-*.md` - Paper notes
  - `pdfs/` - Acquired PDF papers
- **Naming**: Use biological/engineering concepts, not project components: `hepatocyte-metabolism`, not `device-module-3`

**`docs/reports/`**
- Technical reports synthesizing multiple analyses
- Format: `<YYYY-MM-DD>-<report-title>.md`
- Example: `2026-01-27-technical-research-synthesis.md`

**`docs/plans/`**
- Vision documents, specifications, strategic plans
- Format: `<YYYY-MM-DD>-<plan-title>.md`
- Example: `2025-01-25-exo-organ-bioreactor-vision.md`

**`docs/meeting-materials/`**
- Meeting agendas, slide decks, prep materials
- Format: `<YYYY-MM-DD>-<meeting-type>-<topic>.md`

**`docs/presentations/`**
- External presentations (conferences, seminars)
- Subdirectories by date: `2026-01-27-research-synthesis/`
- Contains: slides (markdown, PDF), assets (images), speaker notes

---

### `modules/` - Engineering Subsystems

**Purpose**: Self-contained engineering modules with documentation + models

**Structure per module**:
```
modules/<module-name>/
├── README.md                # Module overview
├── <analysis-files>.md      # Focused analyses
├── models/                  # Module-specific models
├── calculations/            # Back-of-envelope, detailed calculations
└── data/                    # Module-specific data files
```

**Examples**:
- `modules/oxygen-transport/` - All oxygen delivery work (diffusion models, membrane specs, system budgets)
- `modules/cell-biology/` - Cell loading, viability, metabolic function analyses
- `modules/fluid-dynamics/` - Pressure drop, flow distribution, shear stress

**When to create a module**:
- Subsystem is complex enough to warrant multiple documents
- Analyses build on each other (prerequisite relationships)
- Topic has dedicated models or calculations

---

### `models/` - Mathematical Models

**Purpose**: Reusable Python scripts and Jupyter notebooks

**Organization**:
```
models/
├── <topic>/                # Group by physical phenomenon
│   ├── *.py                # Python scripts
│   ├── *.ipynb             # Jupyter notebooks
│   └── README.md           # Model documentation
```

**Examples**:
- `models/diffusion/` - Reaction-diffusion solvers
- `models/fluid-flow/` - CFD scripts, pressure drop calculators
- `models/mass-transfer/` - Mass balance, metabolic network models

**Guidelines**:
- Include docstrings in all Python functions
- Notebooks should have markdown cells explaining approach
- README.md explains: purpose, inputs, outputs, dependencies

---

### `src/` - Python Source Code

**Purpose**: Reusable functions and utilities

**Organization**:
```
src/
├── utils/                  # General utilities (unit conversions, etc.)
├── biology/                # Cell biology calculations (OCR, yields, etc.)
├── physics/                # Physics models (diffusion, flow, etc.)
└── visualization/          # Plotting functions
```

**Guidelines**:
- Production-quality code (not exploratory)
- Tested and documented
- Importable from models/ and experiments/

---

### `experiments/` - Lab Work

**Purpose**: Protocols, raw data, analysis notebooks

**Organization**:
```
experiments/
├── <YYYY-MM-DD>-<experiment-name>/
│   ├── protocol.md         # Step-by-step procedure
│   ├── data/               # Raw data files
│   ├── analysis.ipynb      # Data analysis notebook
│   └── results.md          # Summary and conclusions
```

**Guidelines**:
- Date-stamped experiment directories for chronological tracking
- Keep raw data immutable (don't edit original files)
- Analysis notebooks should be self-contained (reproducible)

---

### `scratchpad/` - Exploratory Work

**Purpose**: Temporary workspace for rough drafts and exploration

**Characteristics**:
- **NOT git-tracked** (in `.gitignore`)
- Agents use `scratchpad/<task-name>/` subdirectories to avoid collisions
- Delete subdirectories when outputs are finalized and moved to `docs/` or `modules/`

**Typical contents**:
- Rough outlines before formal writing
- Exploratory calculations
- Comparison tables while weighing options
- Raw reading notes before structuring

**Cleanup protocol**:
- When task completes, move useful content to appropriate `docs/` or `modules/` location
- Delete `scratchpad/<task-name>/` directory
- Scratchpad should be <1 GB total (clean out old work)

---

## Directory Naming Principles

### Use Plural for Collections
- ✅ `docs/literature/hollow-fiber-membranes/`
- ❌ `docs/literature/hollow-fiber-membrane/`

### Specific Over Generic
- ✅ `modules/oxygen-transport/`
- ❌ `modules/gas-exchange/` (too vague)

### Biological/Engineering Concepts, Not Project Jargon
- ✅ `docs/literature/hepatocyte-metabolism/`
- ❌ `docs/literature/module-A-stuff/`

### Maximum 2-3 Words
- ✅ `fluid-dynamics`, `oxygen-transport`
- ❌ `fluid-flow-modeling-and-pressure-drop-analysis`

### Lowercase Kebab-Case
- ✅ `hollow-fiber-membranes`
- ❌ `HollowFiberMembranes`, `hollow_fiber_membranes`

---

## Creating New Directories

**Before creating a new directory**, check:

1. **Does a similar directory already exist?**
   - Search `ls -R docs/` or use file explorer
   - If similar topic exists, add to existing directory rather than creating duplicate

2. **Is this a long-term category or one-off?**
   - If one-off analysis → put in relevant existing directory
   - If new research area → create new directory

3. **Will this directory contain 3+ files?**
   - If only 1-2 files expected → don't create directory yet
   - Create directory when it becomes unwieldy to keep files in parent dir

**When creating a new directory**:

```bash
# Example: Creating new literature topic area
mkdir -p docs/literature/plasma-separation/pdfs
touch docs/literature/plasma-separation/README.md
```

**Then populate `README.md`**:
```markdown
# Plasma Separation Technology

**Scope**: Methods for separating plasma from whole blood for ex vivo processing

**Key Questions**:
- What separation methods exist? (centrifugation, membrane, etc.)
- What are typical efficiency and throughput values?
- How does separation affect hepatocyte function?

**Documents**:
- (none yet)

**Status**: Initial exploration phase
```

**Finally, update parent INDEX.md** to include new topic area.

---

## Moving Files Between Directories

**When to move**:
- File was created in wrong location during rapid drafting
- Scope of work expanded (paper notes → full review)
- Repository was reorganized

**How to move safely**:
```bash
# Use git mv to preserve history
git mv docs/literature/oxygen-transport.md docs/literature/oxygen-transport/review-oxygen-transport.md

# Update INDEX.md with new path
# Update cross-references in related documents

git commit -m "Move oxygen-transport review to dedicated subdirectory"
```

**Checklist**:
- [ ] Use `git mv` (not `mv`) to preserve git history
- [ ] Update INDEX.md references
- [ ] Update cross-references in other documents that link to this file
- [ ] Verify relative paths (e.g., image links) still work after move
- [ ] Commit with clear message explaining move

---

## Special Directories

### `.claude/skills/` - Agent Configuration

**Important**: This is a **user-level directory**, NOT part of the bioreactor project repository.

**Location**: `~/.claude/skills/` (home directory, not project root)

**Contents**: Skill definitions for agents (researcher, calculator, etc.)

**Do NOT**:
- Commit skill files to bioreactor repo
- Confuse with project documentation structure
- Move project files into `.claude/`

**Project guideline**: Changes to skills are user workflow improvements, not project deliverables. No need to commit skill changes to bioreactor repo.

---

## Maintenance Tasks

**Archivist responsibilities**:

**Weekly**:
- Check for misnamed files (run naming convention audit)
- Update INDEX.md with new documents
- Verify directories have README.md files

**Monthly**:
- Clean out `scratchpad/` (delete completed tasks)
- Archive outdated documents (move to `docs/archive/` with explanation)
- Check for orphaned files (not referenced in INDEX.md)

**Quarterly**:
- Review directory structure for emerging patterns
  - If multiple related files in one directory, consider splitting
  - If directory has <3 files, consider consolidating
- Update CLAUDE.md if conventions evolve

---

## Anti-Patterns

### ❌ Deeply Nested Directories (>4 levels)
- ❌ `docs/literature/membranes/hollow-fiber/oxygen-transport/high-density/`
- ✅ `docs/literature/hollow-fiber-membranes/` (flat, searchable)

**Why**: Hard to navigate, easy to lose files

### ❌ Date-Based Organization
- ❌ `docs/2024/Q1/january/week3/`
- ✅ `docs/literature/<topic>/` (organized by topic, dates in filenames if needed)

**Why**: Time-based structure doesn't reflect how you search ("What did we learn about oxygen transport?")

### ❌ Personal Directories
- ❌ `docs/david-notes/`, `docs/researchers/alice/`
- ✅ `docs/literature/<topic>/` (shared space, git tracks authorship)

**Why**: Documents are project assets; personal directories fragment knowledge

### ❌ "Misc" or "Other" Directories
- ❌ `docs/misc/`, `docs/other-stuff/`
- ✅ Create specific topic directories or leave in parent until category emerges

**Why**: "Misc" becomes a dumping ground; indicates poor categorization

---

## References

- **CLAUDE.md**: Section "Repository Organization"
- **docs/INDEX.md**: Working example of directory structure
- **Git best practices**: Use `git mv` for file moves, commit with clear messages

---

**Principle**: A well-organized directory structure makes knowledge findable. If you can't locate a document in <30 seconds, the structure needs improvement.
