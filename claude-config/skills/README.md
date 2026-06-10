# Claude Skills Ecosystem

**Version**: 3.0.0
**Last Updated**: 2026-06-09
**Total Skills**: 55 (53 active + 2 deprecated tombstones)
**Platform**: Claude Code on Fable 5 / Opus 4.8

## Overview

This directory holds the local skill library that drives `~/.claude/skills/` via
`sync-config.py` (the repository is the source of truth — edit here, never `~/.claude/`
directly). Skills are reusable reference guides for proven techniques and workflows:
literature research, bioinformatics, software development, multi-agent orchestration,
quantitative modeling, scientific writing, and document tooling.

Each skill is a directory containing `SKILL.md` (required) plus optional
`references/`, `examples/`, `assets/`, and `scripts/` subdirectories for progressive
disclosure. Claude loads a skill's `description` to decide relevance, then reads the
body on demand.

> To **add or change** any skill in this repo, use the `skill-editor` skill — it is the
> mandated entry point and wraps the 7-step sync workflow (planning entry → edit repo →
> validate → `push --dry-run` → `push` → test → commit). See `../../CLAUDE.md` and
> `../../CONFIG_MANAGEMENT.md`.

---

## Finding the Right Skill

Every active skill's `description` follows the repo convention: it states **when to use**
the skill and, where overlap exists, an explicit **NOT for… (use X instead)** boundary.
Read the description to disambiguate between neighboring skills. The tables below group
skills by domain; a skill may relate to more than one domain but is listed once in its
primary home.

---

## Skills by Domain

### 🗂 Orchestration & Project Management

| Skill | Use when… |
|-------|-----------|
| **technical-pm** | Coordinating GENERIC multi-agent work with custom skill sequences/dependencies. Default orchestrator when no domain-specific PM applies. |
| **programming-pm** | Coordinating Python software development across multi-specialist pipelines with quality gates (architecture, pre-mortem, review, testing, version control). |
| **lit-pm** | Coordinating comprehensive adaptive literature reviews (4–24h, 9 stages, complexity-calibrated checkpoints). |
| **research-pipeline** | A quick automated research→polish chain (fixed 5 stages, no checkpoints, 2–8h). |
| **program-officer** | Adaptive multi-source research coordination with checkpoints and confidence-rated recommendations. |
| **principal-investigator** | Directing a research project end-to-end — team feedback, final scientific decisions, publication prose, delegation. |
| **brainstorming-pm** | Orchestrating multi-perspective brainstorming via 5 parallel agents with confidence-weighted synthesis. |
| **scientific-analysis-architect** | Planning multi-chapter scientific analyses (RNA-seq, proteomics) with pseudocode and expert consultation. |

### 🔗 Coordination Infrastructure

| Skill | Use when… |
|-------|-----------|
| **workflow-coordinator** | Coordinating handoffs between workflows — universal handoff schema v3.0, validation, distributed tracing, discovery. |
| **perspective-swarm** | Protocol specification for the 5-agent parallel brainstorming swarm (executed by brainstorming-pm). |

### 🔍 Research & Literature

| Skill | Use when… |
|-------|-----------|
| **researcher** | Comprehensive literature research sourcing quantitative parameters from primary literature with citations and context. |
| **literature-researcher** | Deep review-discovery and targeted section research (15–30 papers/section) with convergence tracking (dispatched by lit-pm). |
| **lit-synthesizer** | Synthesizing a literature review into publication-quality prose with authority to restructure (dispatched by lit-pm). |
| **synthesizer** | Integrating multiple reviews/notes, identifying cross-cutting themes and project implications. |
| **fact-checker** | Verifying quantitative claims with inline superscript citations before publication. |
| **consistency-auditor** | Verifying parameter-value consistency across multiple documents. |
| **devils-advocate** | Adversarially reviewing substantive documents to strengthen arguments before editorial polish. |
| **editor** | Prose improvement, CLAUDE.md style enforcement, final polish before archival. |

### ✍️ Scientific Writing & Essays

| Skill | Use when… |
|-------|-----------|
| **essay-pipeline** | Collaboratively writing a science blog essay (interactive thesis → structure → paragraphs) with fact-checking + voice matching. |
| **essay-fact-checker** | Verifying factual claims in a science blog essay with source URLs. |
| **essay-voice-matcher** | Evaluating text against a user's writing style profile for voice consistency. |

### 🧬 Bioinformatics, Notebooks & Data

| Skill | Use when… |
|-------|-----------|
| **bioinformatician** | Implementing -omics data analysis pipelines, statistical tests, or workflows in Python/R. |
| **biologist-commentator** | Evaluating biological relevance, methodological appropriateness, or choosing between analysis methods. |
| **data-pipeline-manager** | Designing/troubleshooting robust data pipelines — stage validation, error handling, monitoring (ETL, FASTQ→BAM→counts). |
| **experimental-planner** | Designing experiments/protocols with hypotheses, controls, success criteria, and resource estimates. |
| **notebook-writer** | Creating reproducible Jupyter notebooks (narrative, environment/session info, git-friendly Jupytext). |
| **notebook-debugger** | Debugging Jupyter errors — kernel crashes, environment conflicts, import/memory issues. |

### 💻 Software Development & Code Quality

| Skill | Use when… |
|-------|-----------|
| **systems-architect** | Designing software architecture, data structures, scalability for complex Python systems. |
| **senior-developer** | Implementing production-quality Python with component-level decisions and formal review of junior outputs. |
| **junior-developer** | Implementing well-scoped Python tasks with unit tests for senior-developer review. |
| **copilot** | Inline code review during development — real-time, adversarial-but-collaborative second opinion. |
| **completion-verifier** | Before marking a task complete or handing off — verifying requirements, edge cases, tests, no regressions. |
| **edge-case-analyst** | BEFORE implementation — identifying failure scenarios while planning features. |
| **systematic-troubleshooter** | Systematic debugging of errors/bugs, with extended thinking for multi-layer issues. |
| **git-strategy-advisor** | Deciding integration strategy — branch vs. commit, when to push, PR vs. merge. |

### 📐 Quantitative, Modeling & Visualization

| Skill | Use when… |
|-------|-----------|
| **calculator** | Physics/engineering feasibility — order-of-magnitude estimates, mass-balance/diffusion/flow models. |
| **mathematician** | Designing deterministic algorithms, complexity analysis, numerical method selection, correctness proofs. |
| **statistician** | Selecting statistical methods, power analysis, uncertainty quantification, MCMC/Monte Carlo validation. |
| **economist** | Order-of-magnitude cost estimates and cost-effectiveness comparisons. |
| **procurement** | Matching equipment specs to vendors and mapping sourcing/lead-time landscape. |
| **plotting-advisor** | BEFORE writing plotting code — chart type, palette, axis, accessibility (Tufte/Cleveland/Wong/Wilke). |

### 🌊 CFD / Bioreactor Simulation

| Skill | Use when… |
|-------|-----------|
| **cfd-bioreactor** | Simulating fluid flow / O2 transport through bioprocess cartridges and bioreactor geometries (FEniCSx). |
| **cfd-mathematician** | (Invoked by cfd-bioreactor) Rigorous FEM variational analysis, stability, convergence, dimensionless numbers. |
| **cfd-reviewer** | (Invoked by cfd-bioreactor) Adversarial engineering review of CFD plans and generated code. |

### 🧭 Strategy, Planning & Requirements

| Skill | Use when… |
|-------|-----------|
| **strategist** | Assessing research direction, identifying knowledge gaps, recommending priorities. |
| **ai-strategist** | Evaluating AI tools and agentic workflows against workflow gaps; quarterly landscape scans. |
| **requirements-analyst** | Refining vague requirements, scope, and success criteria before implementation. |

### 📄 Documents & Domain Tools

| Skill | Use when… |
|-------|-----------|
| **latex-document-manager** | Examining, editing, proofreading, or compiling LaTeX documents on macOS. |
| **patent-review** | Scientist-led review of a lawyer-drafted patent (technical accuracy, claim gaps, revised claim language). |
| **archive-workflow** | Organizing projects — clutter detection, naming conventions, structure, gitignore, expandability. |
| **web-presence-manager** | Reviewing/improving professional web presence across GitHub-hosted sites. |
| **pov-expansion** | Seeking analogous solutions from other domains; cross-domain perspective generation. |

### 🛠 Meta / Skill Tooling

| Skill | Use when… |
|-------|-----------|
| **skill-editor** | Creating, modifying, or refactoring any skill in this repository (the mandated entry point). |

### 🪦 Deprecated (tombstones, retained for redirects)

| Skill | Replaced by |
|-------|-------------|
| **software-developer** | `senior-developer` (deprecated 2026-05-12) |
| **parallel-coordinator** | `superpowers:dispatching-parallel-agents` (deprecated 2026-06-09) |

---

## Frontmatter Conventions

Skills use YAML frontmatter delimited by `---`. The schema is minimal and intentionally
loose; only two fields are enforced.

### Required (enforced at `push` time)

| Field | Notes |
|-------|-------|
| `name` | kebab-case identifier; letters/numbers/hyphens only. Must match the directory name. |
| `description` | Third-person, starts with **"Use when…"**, states triggering conditions (not workflow), with an explicit **NOT for… (use X)** boundary where neighbors overlap. Keep under ~500 chars. |

`sync-config.py push` blocks if any `SKILL.md` is missing `name` or `description`.

### Optional (use consistently when present)

| Field | Notes |
|-------|-------|
| `last_updated` | ISO date (`YYYY-MM-DD`). Update **only when the skill's content actually changes** — it records provenance, so do not bump it on untouched skills. |
| `version` | Optional semver. Bump only on a notable rewrite; the CHANGELOG is the authority for version history. Not required — and intentionally absent on skills that don't track one, rather than backfilled. |
| `success_criteria` | List of measurable completion benchmarks. |
| `deprecated` / `deprecated_date` / `replaced_by` | Tombstone markers (see deprecated table above). |

### Handoff metadata (custom extension)

Skills that participate in cross-workflow handoffs carry a `handoff:` block plus
top-level `categories`, following the **Universal Handoff Schema v3.0**:

```yaml
handoff:
  accepts_from:
    - "*"
  provides_to:
    - "*"            # or a list of specific downstream skill names
  schema_version: "3.0"
  schema_type: universal
categories:
  - research
  - analysis
```

These fields are consumed by `workflow-coordinator` for discovery and are ignored by
Claude Code's built-in skill loader. The legacy v2.0 form (`protocol_version`,
`accepts_handoff`, `handoff_categories`, `requires`, `optional_consumes`) is retired.
Full field definitions and the controlled category vocabulary live in
`workflow-coordinator/references/frontmatter-metadata-standard.md`.

---

## References

- `CHANGELOG.md` — skill-level change history
- `../../CLAUDE.md` — repo invariants and standard workflow
- `../../CONFIG_MANAGEMENT.md` — full 7-step sync workflow with rollback procedures
- `workflow-coordinator/references/` — handoff schema, metadata standard, tracing

**Maintainer**: David Angeles Albores
