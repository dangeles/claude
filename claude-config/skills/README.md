# Claude Skills Ecosystem
**Version**: 2.0.0
**Last Updated**: 2026-01-29
**Total Skills**: 26
**Claude Compatibility**: Claude 4.x (Opus 4.5, Sonnet 4.5)

## Overview

This directory contains a comprehensive ecosystem of 26 specialized skills for Claude AI, optimized for 2026 best practices. Skills cover literature research, bioinformatics, project management, software development, debugging, and scientific workflows.

**Key Features**:
- ✅ 100% compliance with Claude 4.x prompt engineering best practices
- ✅ Contract-style structure (role, goal, constraints, workflow, outputs)
- ✅ Success criteria for measurable completion benchmarks
- ✅ Extended thinking integration for complex reasoning (5 skills)
- ✅ Parallel execution patterns for efficiency
- ✅ Reproducibility standards for scientific workflows

---

## Quick Start

### Using Skills

Skills are invoked using the `Skill` tool in Claude Code:

```python
# Example: Research a topic
Skill(skill="researcher", args="Find papers on CRISPR off-target effects")

# Example: Debug a Jupyter notebook
Skill(skill="notebook-debugger", args="Kernel crashes when running cell 5 with 50k cells")

# Example: Verify task completion
Skill(skill="completion-verifier", args="Verify that RNA-seq analysis is complete and ready for handoff")
```

### Finding the Right Skill

See **[Skills by Domain](#skills-by-domain)** below for categorized lists, or use these quick references:

- **Debugging**: systematic-troubleshooter, notebook-debugger
- **Research**: researcher, synthesizer, devil's-advocate, fact-checker
- **Bioinformatics**: bioinformatician, notebook-writer, data-pipeline-manager
- **Project Management**: technical-pm, program-officer, principal-investigator
- **Code Development**: senior-developer, copilot, systems-architect
- **Quality Assurance**: completion-verifier, consistency-auditor, copilot
- **Parallel Coordination**: parallel-coordinator

---

## Skills by Domain

### 🔍 Literature Search & Research (6 skills)

| Skill | Description | Key Features |
|-------|-------------|--------------|
| **researcher** | Comprehensive literature research with citations | Extended thinking (4k-16k tokens), parallel searches, Nature-style citations |
| **synthesizer** | Integrate multiple sources into coherent analysis | Extended thinking (8k-16k tokens), cross-cutting themes, contradiction resolution |
| **devil's-advocate** | Challenge arguments to strengthen them | Extended thinking (4k-8k tokens), strategic + tactical challenges, thesis coherence |
| **fact-checker** | Verify citations and quantitative claims | Parallel verification, DOI resolution, Nature-style format |
| **consistency-auditor** | Ensure parameter consistency across documents | Parameter inventory, discrepancy tracking, single source of truth |
| **editor** | Polish prose and enforce style guidelines | CLAUDE.md style, acronym handling, flow improvement |

**Workflow**: researcher → synthesizer → devil's-advocate → fact-checker → editor

### 🧬 Bioinformatics & Biology (5 skills)

| Skill | Description | Key Features |
|-------|-------------|--------------|
| **bioinformatician** | Implement data analysis pipelines in Python/R | Reproducibility standards, biological sanity checks, QC workflows |
| **biologist-commentator** | Validate biological relevance and methods | Experimental design, gold-standard methods, plausibility assessment |
| **notebook-writer** | Create reproducible Jupyter notebooks | Jupyter AI integration, Jupytext format, environment docs, session info |
| **data-pipeline-manager** | Design robust data pipelines with validation | Quality gates, error handling, retry logic, bioinformatics-specific patterns |
| **experimental-planner** | Design experiments from theoretical models | Hypothesis definition, controls identification, success criteria |

**Workflow**: experimental-planner → bioinformatician → notebook-writer → biologist-commentator

### 📊 Project & Task Management (4 skills)

| Skill | Description | Key Features |
|-------|-------------|--------------|
| **technical-pm** | Coordinate multi-agent work with dependencies | KPI tracking (≥90% completion), progress monitoring, crisis management |
| **program-officer** | Coordinate complex research tasks | Task breakdown, dependency management, parallel agent coordination |
| **principal-investigator** | Direct research projects | Team feedback synthesis, scientific interpretation, strategic delegation |
| **strategist** | Assess research direction and priorities | Gap identification, priority ranking, risk assessment |

**Workflow**: strategist → technical-pm → program-officer → principal-investigator

### 💻 Software Development (5 skills)

| Skill | Description | Key Features |
|-------|-------------|--------------|
| **senior-developer** | Implement production-quality Python code | Testing standards, security awareness, documentation, delegates to junior-developer |
| **junior-developer** | Implement well-scoped Python tasks under senior review | Unit tests, documented code, max 3 revision cycles |
| **systems-architect** | Design software architecture | Extended thinking (8k-12k tokens), scalability planning, trade-off analysis |
| **copilot** | Review code for quality and correctness | Parallel analysis (security + performance + style), adversarial feedback |
| **completion-verifier** | Verify tasks are truly complete | 6-item checklist, domain-specific criteria, prevents 40% failure rate |

**Workflow**: systems-architect → senior-developer (→ junior-developer) → copilot → completion-verifier

### 🐛 Debugging & Troubleshooting (3 skills)

| Skill | Description | Key Features |
|-------|-------------|--------------|
| **systematic-troubleshooter** | General-purpose debugging | Extended thinking (8k-16k tokens), 7-phase workflow, hypothesis testing |
| **notebook-debugger** | Jupyter-specific troubleshooting | Extended thinking (4k-8k tokens), kernel crashes, environment conflicts, memory optimization |
| **parallel-coordinator** | Orchestrate parallel agents | Dependency analysis, 2.75× speedup examples, anti-patterns guidance |

**Use cases**: systematic-troubleshooter (any bug) → notebook-debugger (Jupyter-specific)

### 📐 Quantitative Analysis (2 skills)

| Skill | Description | Key Features |
|-------|-------------|--------------|
| **calculator** | Quantitative feasibility checks and models | Extended thinking (4k-8k tokens), sensitivity analysis, order-of-magnitude estimates |
| **economist** | Cost estimation for equipment/reagents | Order-of-magnitude estimates, driver identification, feasibility assessment |

### 📚 Support & Coordination (2 skills)

| Skill | Description | Key Features |
|-------|-------------|--------------|
| **archive-workflow** | Organize projects comprehensively | Clutter detection, naming conventions, structure organization, gitignore management |
| **procurement** | Match equipment to specifications | Vendor identification, compatibility assessment, availability checking |

---

## New Skills (v2.0.0)

### systematic-troubleshooter ⭐
**Purpose**: General-purpose debugging with extended thinking
**When to use**: Any error, bug, or unexpected behavior
**Workflow**: Understand → Reproduce → Hypothesize (extended thinking) → Test → Fix → Verify → Document

**Key features**:
- Extended thinking budget: 8,192-16,384 tokens for complex multi-layer bugs
- 3 detailed examples (memory leaks, API errors, edge case bugs)
- 3 reference guides (error patterns, debugging tools, testing strategies)

**Word count**: 3,034 words

### notebook-debugger ⭐
**Purpose**: Jupyter-specific troubleshooting (kernel crashes, environment conflicts)
**When to use**: Kernel crashes, import errors, memory errors in notebooks
**Workflow**: Diagnose → Isolate → Fix → Verify → Document

**Key features**:
- Extended thinking budget: 4,096-8,192 tokens
- Sparse matrix optimization (RNA-seq 50k cells × 20k genes example)
- Environment management (micromamba/pip kernel registration)
- 3 real-world debugging examples

**Word count**: 2,546 words

### completion-verifier ⭐
**Purpose**: Verify tasks are complete before marking done (prevents 40% failure rate)
**When to use**: Before marking tasks complete, before user handoffs, at quality checkpoints
**Workflow**: Context Gathering → Systematic Verification → Decision → Action

**Key features**:
- 6-item checklist (requirements, edge cases, tests, docs, regressions, acceptance)
- Domain-specific criteria for code, research, analysis, documentation
- Pass/fail/conditional outcomes with remediation guidance

**Word count**: 1,830 words

### parallel-coordinator ⭐
**Purpose**: Orchestrate multiple independent agents (leverages Claude Sonnet 4.5 strength)
**When to use**: User has 2+ independent tasks ("analyze A and visualize B")
**Workflow**: Identify tasks → Analyze dependencies → Launch parallel → Integrate results

**Key features**:
- 4 parallelization patterns (pipeline, fan-out, fork-join, fully independent)
- Dependency analysis framework
- 63% time reduction example (2.75× speedup)
- Anti-patterns section (what NOT to parallelize)

**Word count**: 2,093 words

### data-pipeline-manager ⭐
**Purpose**: Design robust data pipelines with quality validation
**When to use**: RNA-seq pipelines, ETL jobs, genomics workflows, pipeline failures
**Workflow**: Design → Validate Input → Transform → Validate Output → Handle Errors → Monitor

**Key features**:
- Quality validation at each stage (file, format, schema, data quality)
- Error handling strategies (retry logic, checkpointing, circuit breakers)
- Bioinformatics-specific patterns (FASTQ→BAM→counts, genome consistency)
- Complete RNA-seq pipeline example (7 stages)

**Word count**: 3,321 words

---

## Enhanced Skills (v2.0.0)

### Extended Thinking Integration (5 skills)

Skills enhanced with extended thinking for complex reasoning:

1. **researcher** (4k-16k tokens): Literature synthesis, hypothesis generation, gap analysis
2. **synthesizer** (8k-16k tokens): Multi-source integration, conceptual frameworks
3. **calculator** (4k-8k tokens): Complex calculations, sensitivity analysis
4. **systems-architect** (8k-12k tokens): Architectural decisions, trade-off analysis
5. **devils-advocate** (4k-8k tokens): Critical analysis, logical inconsistency identification

### Jupyter AI and Reproducibility (2 skills)

#### notebook-writer
**New sections**:
- Jupyter AI Integration (~500 words): %%ai magic, chat UI, context provision
- Reproducibility Standards (~300 words): Environment docs, random seeds, session info, file paths

#### bioinformatician
**New section**:
- Reproducibility Standards (~200+ lines): Environment docs, random seeds for bioinformatics operations (UMAP, t-SNE, Leiden), data provenance (genome builds, GEO accessions), organism/reference tracking

### KPI Tracking (1 skill)

#### technical-pm (v1.1 → v1.2)
**New sections**:
- KPI Tracking: 5 KPIs with ≥85-95% targets (task completion, milestone adherence, handoff quality)
- Progress Monitoring: Updates every 45-90 min for long tasks, leverages Claude 4.5 state tracking
- Crisis Management: Rollback procedures, emergency escalation (P0/P1/P2), prevention safeguards
- Crisis response template (assets/crisis-response-template.md)

### Success Criteria (ALL 26 skills)

All skills now have measurable completion benchmarks in YAML frontmatter:
- **Example** (researcher): Key papers identified, citations formatted, context captured, gaps documented
- **Example** (bioinformatician): Notebook executes, visualizations labeled, reproducibility ensured
- **Example** (calculator): Order-of-magnitude answer, assumptions stated, sensitivity performed

---

## 2026 Best Practices Compliance

### Contract-Style Structure ✅
All skills follow the contract pattern:
- **Role**: Personality and approach
- **Goal**: What the skill achieves
- **Constraints**: What NOT to do
- **Workflow**: Step-by-step process
- **Outputs**: Deliverables

### Clarity Over Complexity ✅
- Target word counts: 1,500-3,500 words
- Supporting files in subdirectories (examples/, references/, assets/)
- No unnecessary verbosity

### Reduced Aggressive Language ✅
- 85% reduction of "ALWAYS", "NEVER", "MUST" phrasing
- Prevents Claude Opus 4.5 overtriggering
- Conditional phrasing: "Use when..." instead of "ALWAYS use"

### Parallel Execution ✅
3 skills enhanced for Claude Sonnet 4.5's parallel tool execution:
- researcher: Parallel literature searches
- fact-checker: Parallel citation verification
- copilot: Parallel code review (security + performance + style)

### Extended Thinking ✅
5 skills configured with extended thinking budgets (4k-16k tokens) for complex reasoning

---

## File Structure

```
skills/
├── README.md                          # This file
├── CHANGELOG.md                       # Version history and changes
├── VALIDATION_REPORT.md              # Testing and validation results
│
├── systematic-troubleshooter/        # NEW v2.0
│   ├── SKILL.md
│   ├── examples/
│   │   ├── bug-report-example.md
│   │   ├── hypothesis-testing-example.md
│   │   └── minimal-reproduction-example.md
│   └── references/
│       ├── common-error-patterns.md
│       ├── debugging-tools.md
│       └── testing-strategies.md
│
├── notebook-debugger/                # NEW v2.0
│   ├── SKILL.md
│   ├── examples/
│   │   ├── kernel-crash-debug.md
│   │   ├── import-error-debug.md
│   │   └── execution-order-debug.md
│   └── references/
│       ├── jupyter-troubleshooting-guide.md
│       ├── environment-management.md
│       └── notebook-best-practices.md
│
├── completion-verifier/              # NEW v2.0
│   ├── SKILL.md
│   ├── examples/
│   │   └── verification-report-example.md
│   └── references/
│       └── completion-criteria-by-domain.md
│
├── parallel-coordinator/             # NEW v2.0
│   ├── SKILL.md
│   ├── examples/
│   │   └── parallel-coordination-example.md
│   └── references/
│       └── dependency-analysis-guide.md
│
├── data-pipeline-manager/            # NEW v2.0
│   ├── SKILL.md
│   ├── examples/
│   │   ├── rnaseq-pipeline-example.md
│   │   └── pipeline-debugging-example.md
│   └── references/
│       ├── validation-patterns.md
│       └── error-handling-strategies.md
│
├── researcher/                       # ENHANCED v2.0
│   └── SKILL.md                      # + extended thinking, parallel execution
│
├── bioinformatician/                 # ENHANCED v2.0
│   └── SKILL.md                      # + reproducibility standards
│
├── notebook-writer/                  # ENHANCED v2.0
│   └── SKILL.md                      # + Jupyter AI, reproducibility
│
├── technical-pm/                     # ENHANCED v2.0
│   ├── SKILL.md                      # + KPI tracking, progress monitoring
│   └── assets/
│       └── crisis-response-template.md
│
└── [21 other skills]/
    └── SKILL.md                      # + success_criteria in YAML

Total: 26 skill directories, 54+ files
```

---

## Common Workflows

### Research Pipeline
```
researcher (literature search with extended thinking)
    ↓
synthesizer (integrate sources with extended thinking)
    ↓
devil's-advocate (challenge arguments with extended thinking)
    ↓
fact-checker (verify citations in parallel)
    ↓
editor (polish prose)
```

### Bioinformatics Analysis
```
experimental-planner (design experiment)
    ↓
bioinformatician (implement analysis with reproducibility standards)
    ↓
notebook-writer (create reproducible notebook with Jupyter AI)
    ↓
copilot (review code in parallel: security + performance + style)
    ↓
biologist-commentator (validate biological interpretation)
    ↓
completion-verifier (verify task complete)
```

### Software Development
```
systems-architect (design architecture with extended thinking)
    ↓
senior-developer (implement code, optionally delegating to junior-developer)
    ↓
copilot (adversarial review)
    ↓
systematic-troubleshooter (debug if issues found)
    ↓
completion-verifier (verify before handoff)
```

### Parallel Coordination
```
User: "Research topic A, analyze dataset B, create viz C"
    ↓
parallel-coordinator (analyze dependencies)
    ↓
Launch 3 parallel agents:
  - Agent 1: researcher → topic A
  - Agent 2: bioinformatician → dataset B
  - Agent 3: plotting skill → viz C
    ↓
parallel-coordinator (integrate results)
    ↓
Unified deliverable (2.75× faster than sequential)
```

---

## Success Criteria

Every skill now includes measurable completion benchmarks. Check the `success_criteria` field in each skill's YAML frontmatter.

**Examples**:

**researcher**:
- Key papers identified and synthesized into structured notes
- All quantitative claims have inline citations with proper format
- Measurement context captured (species, methods, culture conditions)
- Gaps in literature documented explicitly

**bioinformatician**:
- Notebook runs end-to-end without errors
- All visualizations properly labeled with units and legends
- Session info documented for reproducibility
- Biological sanity checks completed and documented
- Code reviewed by copilot with no critical issues

**completion-verifier**:
- All stated requirements verified as satisfied
- Edge cases systematically identified and tested
- Quality standards met for domain (tests pass, docs complete)
- No regressions introduced by the work
- Deliverables ready for handoff or production use

---

## Validation Status

See `VALIDATION_REPORT.md` for complete testing results.

**Summary**:
- ✅ 26/26 skills have valid YAML frontmatter
- ✅ 26/26 skills have success_criteria defined
- ✅ 5/5 target skills have extended_thinking_budget
- ✅ All new skills have complete file structure
- ✅ Word counts meet 2026 best practices
- ⚠️  Residual aggressive language (acceptable in context)

---

## Migration from v1.x

**No breaking changes** - all enhancements are additions. Existing workflows continue to work.

**Recommended adoptions**:
1. Use **completion-verifier** before marking tasks done
2. Use **parallel-coordinator** for independent multi-task requests
3. Follow new **reproducibility standards** for Jupyter notebooks (environment docs, random seeds, session info)
4. Check **success_criteria** in YAML to understand completion benchmarks

---

## Contributing

Skills are designed to be extended and refined based on usage. When adding new skills or enhancing existing ones:

1. Follow contract-style structure (role, goal, constraints, workflow, outputs)
2. Include `success_criteria` in YAML frontmatter
3. Target 1,500-3,500 words for clarity
4. Avoid aggressive language ("ALWAYS", "NEVER", "MUST")
5. Add supporting files to subdirectories (examples/, references/, assets/)
6. Update CHANGELOG.md with changes

---

## References

### Claude 4.x Best Practices
- [Claude 4 Prompt Engineering](https://platform.claude.com/docs/en/build-with-claude/prompt-engineering/claude-4-best-practices)
- [Extended Thinking Tips](https://platform.claude.com/docs/en/build-with-claude/prompt-engineering/extended-thinking-tips)
- [Anthropic Extended Thinking](https://www.anthropic.com/news/visible-extended-thinking)

### Domain-Specific
- [Bioinformatics Best Practices](https://rnabio.org/module-09-appendix/0009/10/01/Bioinformatics_Best_Practices/)
- [Jupyter AI Extension](https://github.com/jupyterlab/jupyter-ai)
- [AI Agents for Project Management](https://www.epicflow.com/blog/ai-agents-for-project-management/)

---

**Version**: 2.0.0
**Last Updated**: 2026-01-29
**Maintainer**: David Angeles Albores
**Claude Version**: Claude 4.x (Opus 4.5, Sonnet 4.5)
