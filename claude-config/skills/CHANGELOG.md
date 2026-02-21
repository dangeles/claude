# Skills Changelog
All notable changes to the Claude skills ecosystem.

## [3.0.0] - 2026-02-21

### skill-editor v3.0 — Major Refactor

- Consolidated from 9 agents to 5 (hybrid swarm + specialist model)
- Integrated brainstorming-pm swarm delegation for Phase 2 analysis
- Reduced SKILL.md from 2,206 to ~645 lines via bash extraction + architectural simplification
- Simplified from 3 modes (SIMPLE/STANDARD/EXPERIMENTAL) to 2 (QUICK/FULL)
- Consolidated from 5 quality gates to 3
- Eliminated conditional Phase 2.5 (absorbed into swarm)
- Updated all model references to opus (resolves to latest Opus)
- Added swarm challenge template (references/swarm-challenge-templates.md)

### Agents Removed
- skill-editor-best-practices-reviewer (absorbed into swarm)
- skill-editor-knowledge-engineer (absorbed into swarm)
- skill-editor-external-researcher (absorbed into swarm)
- skill-editor-strategy-consultant (absorbed into swarm)
- skill-editor-decision-synthesizer (inline synthesis by orchestrator)

## [2.1.0] - 2026-02-04

### Summary
Added comprehensive project organization skill (archive-workflow) with 1 PM orchestrator and 5 specialist agents. Replaced the single-purpose archivist skill with a full-featured multi-agent workflow.

### Added - New Skills (1)

#### archive-workflow (v1.0)
- **Purpose**: Comprehensive project organization with clutter detection, naming conventions, structure organization, and gitignore management
- **Architecture**: 1 PM (library-pm) + 5 specialist agents (clutter-analyst, nomenclature-enforcer, structure-organizer, expandability-reviewer, decision-integrator)
- **Features**:
  - 4-wave workflow with quality gates between each wave
  - Project type detection (code/research/data/mixed)
  - Pre-flight validation (git status, permissions, disk space)
  - Conflict resolution rules (7 rules for handling analyst disagreements)
  - Rollback procedures (pre-commit and post-commit)
  - Editor skill integration for documentation
- **Use cases**: Organizing new projects, cleaning up existing projects, enforcing naming conventions, managing gitignore
- **References**: 8 reference documents (naming conventions, structure templates, gitignore patterns, clutter detection rules)
- **Examples**: 3 worked examples (code project, research project, mixed project)

### Removed

#### archivist (REPLACED by archive-workflow)
- **Reason**: Limited to document filing for specific project type; archive-workflow provides comprehensive project organization across all project types
- **Migration**: Use `/archive-workflow` instead of `/archivist`
- **Preserved**: Naming conventions and directory standards absorbed into archive-workflow references

### Changed

#### Updated References (7 files)
- editor/SKILL.md: Handoff to archive-workflow instead of Archivist
- research-pipeline/SKILL.md: Updated handoff documentation
- technical-pm/SKILL.md: Updated agent assignment guide
- technical-pm/references/coordination-patterns.md: Updated literature pipeline flow
- technical-pm/references/parallel-execution.md: Updated skill list
- technical-pm/references/dependency-detection.md: Updated estimation heuristics
- skills/README.md: Updated skill catalog

## [2.0.0] - 2026-01-29

### Summary
Major update implementing Claude 4.x (2026) best practices across all 26 skills. Added 5 new skills, enhanced 15 existing skills, and achieved 100% compliance with success criteria standards.

### Added - New Skills (5)

#### systematic-troubleshooter (v1.0)
- **Purpose**: General-purpose debugging with extended thinking for complex issues
- **Features**:
  - 7-phase debugging workflow (Understand → Reproduce → Hypothesize → Test → Fix → Verify → Document)
  - Extended thinking budget: 8,192-16,384 tokens for multi-layer bugs
  - 3 detailed examples (bug reports, hypothesis testing, minimal reproduction)
  - 3 reference guides (error patterns, debugging tools, testing strategies)
- **Use cases**: Any error, bug, unexpected behavior across all domains
- **Word count**: 3,034 words

#### notebook-debugger (v1.0)
- **Purpose**: Jupyter-specific troubleshooting (kernel crashes, environment conflicts, data pipeline failures)
- **Features**:
  - 5-phase workflow (Diagnose → Isolate → Fix → Verify → Document)
  - Extended thinking budget: 4,096-8,192 tokens
  - Sparse matrix memory optimization patterns
  - Environment management with micromamba/pip kernel registration
- **Use cases**: Kernel crashes, import errors, memory errors, execution order issues
- **Examples**: RNA-seq analysis kernel crash (50k cells × 20k genes), environment conflicts
- **Word count**: 2,546 words

#### completion-verifier (v1.0)
- **Purpose**: Verify tasks are truly complete before marking done, preventing 40% failure rate from premature completion
- **Features**:
  - 6-item verification checklist (requirements, edge cases, tests, docs, regressions, acceptance)
  - 4-phase workflow (Context → Verification → Decision → Action)
  - Domain-specific completion criteria (code, research, analysis, docs)
- **Use cases**: Before marking tasks complete, before user handoffs, at quality checkpoints
- **Integration**: Called by technical-pm, program-officer, execution skills
- **Word count**: 1,830 words

#### parallel-coordinator (v1.0)
- **Purpose**: Orchestrate multiple independent agents running simultaneously (leverages Claude Sonnet 4.5 parallel execution strength)
- **Features**:
  - Dependency analysis framework
  - 4 parallelization patterns (pipeline, fan-out, fork-join, fully independent)
  - Anti-patterns section (what NOT to parallelize)
- **Use cases**: "Analyze dataset A and create visualization B" → 2 parallel agents
- **Performance**: Example shows 63% time reduction (2.75× speedup) for 3 parallel research tasks
- **Word count**: 2,093 words

#### data-pipeline-manager (v1.0)
- **Purpose**: Design and troubleshoot data pipelines with quality validation (addresses prevalent AI agent error cause)
- **Features**:
  - 6-stage workflow (Design → Input Validation → Transform → Output Validation → Error Handling → Monitoring)
  - Quality validation at each stage (file system, format, schema, data quality)
  - Error handling strategies (retry logic, checkpointing, circuit breakers)
  - Bioinformatics-specific patterns (FASTQ→BAM→counts, genome consistency, metadata tracking)
- **Use cases**: RNA-seq pipelines, ETL jobs, genomics workflows
- **Examples**: Complete RNA-seq pipeline (7 stages), pipeline debugging scenarios
- **Word count**: 3,321 words

### Changed - Enhanced Skills (15)

#### Aggressive Language Reduction (8 skills)
Reduced "ALWAYS", "NEVER", "MUST" phrasing by 85% to prevent Claude Opus 4.5 overtriggering:
- **researcher** (v1.2 → v1.3): 7 instances removed
- **calculator** (v1.1 → v1.2): 3 instances removed
- **copilot**: "MUST fix" → "Fix before proceeding"
- **program-officer**: "Always" → "Escalate"
- **technical-pm**: "you must monitor" → "monitor"
- **principal-investigator**: "Always frame" → "Frame"
- **archivist**: "Always ask" → "Ask"
- **fact-checker**: "Each must have" → "Each should have"

#### Parallel Execution Patterns (3 skills)
Added guidance for Claude Sonnet 4.5's parallel tool execution:
- **researcher**: Parallel literature searches across PubMed, bioRxiv, OpenAlex
- **fact-checker**: Parallel citation verification for independent claims
- **copilot**: Parallel code review (security + performance + style simultaneously)

#### Extended Thinking Integration (5 skills)
Added extended thinking budgets with metadata:
- **researcher** (v1.2 → v1.3): 4,096-16,384 tokens
  - Use cases: Literature synthesis, hypothesis generation, gap analysis, contradiction resolution
- **synthesizer** (v1.0 → v1.1): 8,192-16,384 tokens
  - Use cases: Multi-source integration, cross-cutting patterns, conceptual frameworks
- **calculator** (v1.1 → v1.2): 4,096-8,192 tokens
  - Use cases: Complex calculations, uncertainty propagation, sensitivity analysis, trade-offs
- **systems-architect** (v1.0 → v1.1): 8,192-12,288 tokens
  - Use cases: Architectural decisions, scalability planning, technology selection
- **devils-advocate** (v1.1 → v1.2): 4,096-8,192 tokens
  - Use cases: Deep critical analysis, logical inconsistency identification, alternative interpretations

#### Jupyter AI and Reproducibility (2 skills)

**notebook-writer** (enhanced):
- Added "Jupyter AI Integration" section (~500 words):
  - %%ai magic command usage and best practices
  - Chat UI guidance
  - Context provision (API docs, dataset descriptions)
  - When to use AI vs manual coding
  - JetBrains AI Assistant features
- Added "Reproducibility Standards" section (~300 words):
  - Environment documentation requirements
  - Random seed setting for stochastic processes
  - Session info output mandates
  - File path best practices (Path variables, no hardcoded paths)
  - Reproducibility checklist

**bioinformatician** (enhanced):
- Added comprehensive "Reproducibility Standards" section (~200+ lines):
  - Environment documentation (micromamba/pip, kernel specification)
  - Random seed setting (numpy, scanpy, PyTorch, TensorFlow)
  - Bioinformatics-specific operations requiring seeds (UMAP, t-SNE, Leiden clustering, permutation tests)
  - Session info output (session_info.show(), scanpy logging)
  - Data provenance documentation (GEO accessions, genome builds, download dates)
  - Pre-flight and handoff checklists
  - Organism and reference version tracking (genome builds, annotations)
  - QC parameter documentation

#### KPI Tracking and Progress Monitoring (1 skill)

**technical-pm** (v1.1 → v1.2):
- Added "KPI Tracking and Success Metrics" section:
  - 5 KPIs with specific targets:
    * Task completion rate: ≥90%
    * Milestone adherence: ≥85%
    * Handoff quality (no rework): ≥95%
    * User intervention rate: ≤15%
    * Bug discovery rate in review: ≤10%
  - Tracking methodology (weekly reviews, per-project retrospectives)
  - Improvement actions when targets missed
- Added "Progress Monitoring for Long-Running Tasks":
  - Update schedules (45-90 min intervals based on task length)
  - Leverages Claude 4.5's extended context for multi-turn workflows
  - Incremental milestones with checkpoints
- Added "Crisis Management" section:
  - Crisis definitions (data loss, security breaches, critical failures, system issues)
  - 3-step response protocol (STOP → Notify → Execute)
  - Rollback procedures with evidence preservation
  - Emergency escalation path (P0/P1/P2 severity levels)
  - Prevention safeguards
- Created crisis response template (assets/crisis-response-template.md)

#### Success Criteria (26 skills)
Added `success_criteria` field to YAML frontmatter of all 26 skills:
- 2-7 specific, measurable completion benchmarks per skill
- Domain-appropriate criteria (code testing, research citations, analysis reproducibility)
- Examples:
  - **researcher**: Key papers identified, citations formatted, context captured, gaps documented
  - **bioinformatician**: Notebook executes, visualizations labeled, reproducibility ensured, sanity checks done
  - **calculator**: Order-of-magnitude answer, assumptions stated, sensitivity analysis, implications interpreted
  - **software-developer**: Tests pass, style compliance, documentation complete, security verified

### Fixed
- Task completion verification gaps (addressed by completion-verifier skill)
- Missing reproducibility standards for Jupyter notebooks (notebook-writer, bioinformatician)
- No guidance for parallel task coordination (parallel-coordinator skill)
- Data pipeline failure handling (data-pipeline-manager skill)
- Inconsistent success criteria across skills (now standardized)

### Validation Results
See `VALIDATION_REPORT.md` for complete validation details.

**Summary**:
- ✅ 26/26 skills have valid YAML frontmatter
- ✅ 26/26 skills have success_criteria defined
- ✅ 5/5 target skills have extended_thinking_budget
- ✅ All new skills have complete file structure (examples + references)
- ✅ Word counts meet 2026 best practices (1,800-3,500 words)
- ⚠️  Residual aggressive language (14 skills, 1-6 instances) - acceptable in context

### Metrics Achieved

| Metric | Target | Achieved | Status |
|--------|--------|----------|--------|
| Skills with success_criteria | 100% | 100% (26/26) | ✅ |
| Skills with extended thinking | 5 | 5 | ✅ |
| New priority skills | 2 minimum | 5 | ✅ (250%) |
| Aggressive language reduction | 80% | ~85% | ✅ |
| Parallel execution patterns | 3 | 3 | ✅ |
| Reproducibility enhancements | 2 | 2 | ✅ |
| Word count compliance | 100% | 100% | ✅ |

---

## Implementation Timeline

**Month 1 (Tasks 1-8)** - Completed 2026-01-29:
1. Aggressive language audit (8 skills)
2. Create systematic-troubleshooter skill
3. Create notebook-debugger skill
4. Add parallel execution patterns (3 skills)
5. Verify skill lengths (all compliant, no work needed)
6. Integrate extended thinking (5 skills)
7. Expand notebook-writer (Jupyter AI + reproducibility)
8. Enhance bioinformatician (reproducibility)

**Month 2 (Tasks 9-13)** - Completed 2026-01-29:
9. Create completion-verifier skill
10. Create parallel-coordinator skill
11. Create data-pipeline-manager skill
12. Add KPI tracking to technical-pm
13. Add success_criteria to all 26 skills

**Validation & Documentation (Tasks 14-15)** - Completed 2026-01-29:
14. Testing and validation (VALIDATION_REPORT.md)
15. Comprehensive documentation (this CHANGELOG, README)

---

## Migration Guide

### For Existing Workflows

**No breaking changes** - all enhancements are additions or refinements. Existing workflows continue to work.

**Recommended adoptions**:

1. **Use completion-verifier before marking tasks done**:
   ```
   [Your work] → completion-verifier → mark complete
   ```

2. **Use parallel-coordinator for independent tasks**:
   ```
   User: "Analyze dataset A and create viz B"
   You: parallel-coordinator → launch 2 agents → integrate results
   ```

3. **For Jupyter notebooks, follow new reproducibility standards**:
   - Document environment (requirements.txt or environment.yml)
   - Set random seeds for stochastic processes
   - Add session info cell at end
   - Use Path variables for file paths

4. **For bioinformatics analyses, include**:
   - Genome build and annotation version
   - Organism specification
   - Data provenance (GEO accessions, download dates)

5. **Check success_criteria in YAML frontmatter** to understand completion benchmarks for each skill

### For New Projects

Start with these skills for different use cases:

**Debugging**: systematic-troubleshooter (general) or notebook-debugger (Jupyter-specific)

**Task completion**: completion-verifier (verify before marking done)

**Parallel work**: parallel-coordinator (orchestrate 2+ independent tasks)

**Data pipelines**: data-pipeline-manager (design robust pipelines with validation)

**Research**: researcher → synthesizer → devil's-advocate → fact-checker (with extended thinking)

**Bioinformatics**: bioinformatician → notebook-writer (with reproducibility standards)

**Project management**: technical-pm (with KPI tracking) → completion-verifier

---

## File Statistics

**Total files**:
- 26 skill directories
- 26 SKILL.md files
- 28+ supporting files (examples, references, assets)
- 54+ total files

**Content added**:
- New skills: ~12,800 words
- Enhanced skills: ~7,000 words
- Supporting files: ~15,000 words
- **Total**: ~35,000 words

---

## Known Issues

### Minor
- Residual aggressive language in 14 skills (1-6 instances each)
  - **Impact**: Minimal - most are in examples or legitimate technical contexts
  - **Status**: Acceptable, no action required

### Future Enhancements
- Monitor actual KPI metrics over 2-4 weeks to validate targets
- Collect user feedback on new skills for refinements
- Calibrate extended thinking token budgets based on usage
- Integration testing for skill handoff workflows

---

## Contributors
- Implementation: Claude Sonnet 4.5 (claude-code)
- Validation: Claude Sonnet 4.5 (claude-code)
- Research: Based on Anthropic 2026 best practices, academic literature, industry implementations
- User: David Angeles Albores

---

## References

### 2026 Best Practices Sources
- [Claude 4 Prompt Engineering Best Practices](https://platform.claude.com/docs/en/build-with-claude/prompt-engineering/claude-4-best-practices)
- [Extended Thinking Tips](https://platform.claude.com/docs/en/build-with-claude/prompt-engineering/extended-thinking-tips)
- [Anthropic's Extended Thinking Announcement](https://www.anthropic.com/news/visible-extended-thinking)
- [Best Practices for AI Agent Implementations 2026](https://onereach.ai/blog/best-practices-for-ai-agent-implementations/)

### Domain-Specific References
- Bioinformatics: [Griffith Lab Best Practices](https://rnabio.org/module-09-appendix/0009/10/01/Bioinformatics_Best_Practices/)
- Jupyter AI: [jupyter-ai extension](https://github.com/jupyterlab/jupyter-ai), [JetBrains AI Assistant](https://www.jetbrains.com/help/ai-assistant/ai-in-jupyter-notebooks.html)
- Project Management: [AI Agents for Project Management 2026](https://www.epicflow.com/blog/ai-agents-for-project-management/)

---

For detailed validation results, see `VALIDATION_REPORT.md`.
For skills overview and usage guide, see `README.md`.
