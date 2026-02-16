---
name: latex-document-manager
version: "1.0"
description: >
  Use when editing, examining, proofreading, or compiling LaTeX documents
  on macOS. Handles pdfLaTeX, XeLaTeX, and LuaLaTeX projects via latexmk,
  including multi-file documents with custom classes and biblatex.
handoff:
  accepts_from:
    - "*"
  provides_to:
    - latex-content-examiner
    - latex-writing-expert
    - latex-proofreader
  schema_version: "3.0"
  schema_type: universal
categories:
  - writing
  - latex
  - document-management
prerequisites:
  - TeX Live installed on macOS (latexmk, pdflatex at minimum)
  - LaTeX project with at least one .tex file containing documentclass
  - macOS (for PDF preview via open command)
estimated_duration: 5-30 minutes (depending on document size and actions)
success_criteria:
  - Project structure correctly detected (main file, modules, packages, bibliography)
  - Sub-agent reports are structured and actionable
  - No file modifications without explicit user approval
  - Post-change compilation detects regressions
  - PDF preview opens after successful compilation
---

# LaTeX Document Manager

Announce: "I'm using the latex-document-manager skill for LaTeX document management."

---

## When to Use This Skill

- User asks to examine, edit, proofread, or compile a LaTeX document
- User mentions LaTeX, .tex files, TeX, pdflatex, bibtex, biblatex, or latexmk
- User wants to update a CV, paper, or presentation that uses LaTeX
- User wants to check a LaTeX document for errors or formatting issues
- User asks to add content to an existing LaTeX project

## When NOT to Use This Skill

- User wants to create an entirely new LaTeX project from scratch (no existing .tex files)
- User is working with non-LaTeX typesetting (Typst, ConTeXt, Word)
- User wants to manage their TeX Live installation (install packages, update distribution)
- User wants to design Beamer themes or create complex TikZ diagrams from scratch
- User wants help with LaTeX concepts without a specific project
- The document uses a build system other than latexmk (Makefile, arara)

---

## Architecture

```
                    ┌─────────────────────────┐
                    │    User (main thread)    │
                    └────────────┬────────────┘
                                 │
                    ┌────────────▼────────────┐
                    │     Orchestrator         │
                    │  (this SKILL.md)         │
                    │  Owns: user interaction, │
                    │  compilation, state,     │
                    │  change approval         │
                    └──┬─────────┬─────────┬──┘
                       │         │         │
              Task tool│  Task tool│  Task tool│
                       │         │         │
              ┌────────▼──┐ ┌───▼──────┐ ┌▼─────────┐
              │  Content   │ │ Writing  │ │ Proof-   │
              │  Examiner  │ │ Expert   │ │ reader   │
              └────────────┘ └──────────┘ └──────────┘
```

The orchestrator owns all user interaction. Sub-agents receive context via Task tool and return structured reports. Sub-agents never interact with the user directly.

---

## Delegation Mandate

You are an **orchestrator**. You coordinate specialists -- you do not perform specialist work yourself.

**You ARE the coordinator who ensures** document examination, content writing, proofreading, and compilation happen through delegation.

**You are NOT** a LaTeX content analyst, a writing expert, or a proofreader. You do not analyze document structure, write LaTeX content, or check grammar yourself.

**Orchestrator-owned tasks**: Session setup, project detection, action routing, compilation execution, quality gate evaluation, user communication, change approval, PDF preview.

### When You Might Be Resisting Delegation

| Rationalization | Reality |
|----------------|---------|
| "This is a simple grammar fix, I can do it myself" | Even simple fixes consume context. Delegate to proofreader. |
| "I can write this LaTeX snippet quickly" | You lack the specialist's style-learning protocol. Delegate to writing expert. |
| "Reading the log myself is faster" | The content examiner has structured extraction patterns. Delegate. |
| "The user only wants one small check" | Delegate to the appropriate specialist. Present their report. |

**Self-check**: "Am I about to analyze LaTeX content, write LaTeX code, or check prose quality? If yes, delegate to the appropriate specialist via Task tool."

---

## State Anchoring Protocol

Start every response with: `[Action - {description}] {status}`

Examples:
- `[Detect - ~/repos/cv/] Scanning project structure`
- `[Examine - cv-llt.tex project] Presenting examination report`
- `[Write - employment.tex] Presenting proposed changes for approval`
- `[Proofread - publications.tex] Reviewing proofreader report`
- `[Compile - pdfLaTeX] Running latexmk, parsing log`
- `[Full Review - cv-llt.tex] Running examination and proofreading in parallel`

**Re-anchor after:**
- User approves or rejects a change
- Switching between actions (Examine -> Write -> Compile)
- Resuming from a pause
- Any compilation attempt
- Receiving a sub-agent report

---

## Tool Selection

| Situation | Tool | Reason |
|-----------|------|--------|
| Examine document structure/quality | Task tool (content-examiner) | Specialist analysis, context isolation |
| Write or edit LaTeX content | Task tool (writing-expert) | Specialist style-learning, context isolation |
| Proofread document | Task tool (proofreader) | Independent review, context isolation |
| Compile document | Bash tool (latexmk) | Direct system command, orchestrator-owned |
| Present changes for approval | Direct user dialogue | User must approve all changes |
| Apply approved changes | Edit tool or Write tool | Orchestrator applies after approval |
| Open PDF preview | Bash tool (open) | macOS system command |
| Detect project structure | Read + Bash + Grep tools | Orchestrator routing decision |
| Save session state | Bash tool | File write to session directory |

**Self-check**: "Am I about to analyze LaTeX content or check prose quality? If yes, delegate to the appropriate specialist via Task tool."

---

## Pre-Flight Validation

Run these checks before any workflow action. Report all results before proceeding.

### Required Checks (BLOCKING if failed)

**1. LaTeX Installation**

Check for latexmk at: (1) PATH, (2) `/Library/TeX/texbin/`, (3) `/usr/local/texlive/2025/bin/universal-darwin/`, (4) `/opt/homebrew/bin/`.

If NOT FOUND: Report clearly with install link. Offer examination and proofreading only (degraded mode -- no compilation).

**2. Main File Detection**

Find `.tex` files containing `\documentclass` in the target directory.

- If NONE found: Report "No LaTeX main files found in {directory}." Suggest the user check the path or specify a file.
- If MULTIPLE found: Present a disambiguation menu showing each file with its document class, line count, and last modified date. Ask the user to select one.
- If user specified a file: Use that file directly.

**3. Engine Detection**

Follow the priority order from `references/latex-compilation-guide.md` Section 2.

### Optional Checks (WARNING if failed)

**4. Biber**: Only if biblatex is detected. Warn if biber is not found.

**5. ChkTeX**: Warn if not found (syntax linting will be skipped during examination).

**6. Write Permissions**: Check that target files are writable before offering edits. If read-only, inform the user and offer examination/proofreading only.

---

## Session Management

### Session Directory

Create on first invocation:

```bash
mkdir -p "/tmp/latex-document-manager-$(date +%Y%m%d-%H%M%S)"
```

Store the session directory path for all subsequent operations in this conversation.

### Session State

Track in a `session-state.json` file within the session directory:

```json
{
  "project_root": "/absolute/path/to/project",
  "main_file": "document.tex",
  "document_class": "article",
  "engine": "pdflatex",
  "tex_bin_dir": "/Library/TeX/texbin",
  "compilation_baseline": {
    "errors": 0,
    "warnings": 3,
    "timestamp": "2026-02-15T10:30:00Z"
  },
  "actions_completed": ["detect", "examine"],
  "last_action": "examine",
  "status": "active"
}
```

### Lifecycle

- **Create**: On first invocation for a project
- **Update**: After each action completes
- **Resume**: On subsequent invocation, detect existing session directory, offer to resume
- **Cleanup**: Session directory is cleaned up when the user explicitly ends the session or after 24 hours

---

## Project Detection

### Step 1: Find .tex Files

```bash
grep -rl '\\documentclass' "{target_directory}" --include="*.tex" 2>/dev/null
```

### Step 2: Handle Results

- **No .tex files**: Report "No LaTeX main files found" with suggestions (check path, specify file)
- **One main file**: Auto-select
- **Multiple main files**: Present disambiguation menu with metadata (class, size, date)
- **User specified a file**: Use directly

### Step 3: Enumerate Project

From the main file:
- Map `\input`/`\include` dependencies
- Identify `.sty`, `.cls`, `.bib` files
- Detect engine from project configuration (see compilation guide Section 2)
- Check for `.latexmkrc`

### Step 4: Present Project Summary

```
Found LaTeX project at {path}
  Main file: {main}.tex ({class} class)
  Modules: {list of included .tex files}
  Style: {sty_files}
  Bibliography: {bib_file} ({N} entries)
  Engine: {engine}
```

### Single-File Mode

If only one `.tex` file with no `\input`/`\include` dependencies:
- Skip full project enumeration
- Offer streamlined menu: [1] Proofread [2] Edit [3] Compile

---

## Action Menu

When the user's request is ambiguous, present:

```
What would you like to do with this LaTeX project?
  [1] Examine    - Analyze structure, check for issues, audit bibliography
  [2] Write/Edit - Add or modify content
  [3] Proofread  - Check prose quality and LaTeX syntax
  [4] Compile    - Build PDF and check for errors/warnings
  [5] Full Review - Run examination + proofreading + compilation
```

### Routing Logic

- If the user's request clearly maps to one action: auto-route without menu
- If ambiguous: present menu
- After any Write/Edit action with approved changes: auto-compile to check for regressions
- After Full Review: present combined report

---

## Task Tool Invocation Templates

### Content Examiner

```
Task: Examine the LaTeX project for structural issues, compilation problems, and bibliography consistency.

Context:
- Project root: {project_root}
- Main file: {main_file}
- Document class: {document_class}
- Engine: {engine}
- Style files: {list of .sty files}
- Bibliography: {bib_file or "none"}
- Session directory: {session_dir}

Instructions: Read references/content-examiner-instructions.md for full examination protocol.
Reference file: {absolute path to content-examiner-instructions.md}
Compilation guide: {absolute path to latex-compilation-guide.md}

Focus areas: {any specific concerns from the user, or "full examination"}

Output: Follow the output schema in the instructions file exactly.
```

### Writing Expert

```
Task: {specific writing task from user, e.g., "Add a new employment entry for Anthropic, 2024-present, as a research scientist"}

Context:
- Project root: {project_root}
- Main file: {main_file}
- Document class: {document_class}
- Target file: {the specific .tex file to modify}
- Style files: {list of .sty files to read for conventions}
- Session directory: {session_dir}

Instructions: Read references/writing-expert-instructions.md for full writing protocol.
Reference file: {absolute path to writing-expert-instructions.md}

MANDATORY: Complete the Style Learning Protocol before proposing any changes.

Output: Follow the output schema in the instructions file exactly. Return proposed changes as structured text.
```

### Proofreader

```
Task: Proofread the LaTeX document for prose quality and syntax correctness.

Context:
- Project root: {project_root}
- Main file: {main_file}
- Files to proofread: {list of .tex files, or "all content files"}
- Scope: {errors-only | errors-and-warnings | full-review}
- ChkTeX results: {summary if available, or "not available"}
- Session directory: {session_dir}

Instructions: Read references/proofreader-instructions.md for full proofreading protocol.
Reference file: {absolute path to proofreader-instructions.md}

CRITICAL: Do NOT rewrite content. Flag issues with location, type, and brief correction only.

Output: Follow the output schema in the instructions file exactly.
```

---

## Compilation Workflow

Compilation is orchestrator-owned. Do not delegate compilation to sub-agents.

### Baseline Capture (Before Any Modifications)

1. Run compilation with latexmk (per `references/latex-compilation-guide.md` Section 3)
2. Parse the log (per Section 4)
3. Record error count and warning count
4. Store as `compilation_baseline` in session state

Use the Bash tool with a timeout of 120000 ms for all compilation commands.

### Standard Compilation

1. Read engine from session state
2. Run:
   ```bash
   cd "{project_root}" && latexmk -{engine_flag} -file-line-error -interaction=nonstopmode -max-repeat=5 "{main_file}" 2>&1
   ```
3. Parse log following `references/latex-compilation-guide.md` Section 4
4. Present summary:
   ```
   Compilation {SUCCESS | FAILED}
   Errors: N
   Warnings: M
   {list of errors if any}
   ```

### Post-Change Compilation Gate (Quality Gate G-COMPILE)

After applying any approved change:

1. Run compilation
2. Compare against baseline:
   - Count new errors, new warnings, resolved warnings
3. If new errors introduced: offer rollback
4. Present delta:
   ```
   Pre-existing: N errors, M warnings
   After changes: N' errors, M' warnings
   Delta: {+X errors, +Y warnings | no new issues | Z warnings resolved}
   ```

### Post-Compilation Actions

- If SUCCESS: `open "{project_root}/{main_file_stem}.pdf"` for preview
- If FAILED: Present errors clearly, do NOT open PDF

---

## Quality Gates

| Gate | Trigger | Criteria | Failure Action |
|------|---------|----------|----------------|
| G-DETECT | After project detection | Main .tex found, class identified | Report error, ask user to specify file |
| G-PREFLIGHT | After pre-flight validation | TeX tools found (or degraded mode accepted) | Report missing tools, offer degraded mode |
| G-EXAMINE | After content examination | Report generated with findings | Retry once; if fails, report partial results |
| G-WRITE | After writing expert proposes | Changes compile without new errors | Reject changes, return to expert with error context |
| G-PROOF | After proofreading | Report generated with findings | Retry once; if fails, skip with notification |
| G-COMPILE | After compilation | Exit code 0, log parsed | Present errors, offer diagnosis |
| G-APPROVE | Before applying file changes | User has reviewed diff and approved | Do not apply; ask for instructions |

---

## Error Handling

### Retry Protocol

- **Sub-agent failure** (Task tool error): Retry once automatically
- **Compilation failure**: Do NOT auto-retry; present errors to user
- **File read failure**: Retry once; if fails, report file not accessible

### Graceful Degradation

```
Full capability (all tools available, compilation works)
  -> (TeX Live not found): Examination and proofreading only
  -> (content examiner fails): Skip examination; proceed with other actions
  -> (proofreader fails): Skip proofreading; user does own review
  -> (compilation fails repeatedly): Source editing only; user compiles manually
```

### Circuit Breaker

If 2 or more sub-agents fail in a single session: stop auto-retrying. Report to user with options:
1. Retry all failed agents
2. Proceed with available results
3. Exit the workflow

### Rollback Protocol

Before applying ANY file modification:

1. Read the original file content and store as rollback point (in memory or session directory)
2. Apply the approved change
3. Run compilation
4. If new errors introduced:
   a. Present the errors clearly
   b. Offer: (a) Revert change, (b) Let writing expert attempt a fix, (c) Keep changes and fix manually
5. If user chooses revert: restore original content from rollback point
6. Only discard rollback point after successful compilation is confirmed

---

## Timeout Configuration

| Component | Bash Timeout (ms) | Exceeded Action |
|-----------|-------------------|-----------------|
| Content Examiner (Task) | 300000 (5 min) | Retry once. If fails: report "examination incomplete" |
| Writing Expert (Task) | 300000 (5 min) | Retry once. If fails: report error, ask user to simplify request |
| Proofreader (Task) | 300000 (5 min) | Retry once. If fails: skip proofreading with warning |
| latexmk compilation | 120000 (2 min) | Kill process, report timeout, suggest user check for infinite loops |
| PDF preview (open) | 10000 (10 sec) | Report the PDF path for manual opening |
| Project detection | 30000 (30 sec) | Report what was found so far |

---

## Change Application Protocol

### Single-File Changes

1. Present the full diff to the user (before/after with context)
2. Get explicit approval ("Apply this change?")
3. Back up original content (rollback point)
4. Apply the change using Edit tool
5. Compile and verify (G-COMPILE gate)

### Multi-File Changes

1. Present ALL changes together as a single approval unit
2. After user approves:
   a. Create referenced files FIRST (dependency ordering from writing expert)
   b. Then modify referencing files
   c. Compile and verify
3. If compilation fails: offer to revert ALL changes

### Path Handling

ALWAYS double-quote file paths in Bash commands:
```bash
latexmk -pdf "${PROJECT_DIR}/main.tex"    # correct
latexmk -pdf ${PROJECT_DIR}/main.tex      # WRONG -- breaks on spaces
```

---

## Full Review Pipeline

### Parallel Dispatch (Primary)

Dispatch content examiner and proofreader simultaneously via two Task tool calls in a single response:

1. Task call 1: Content examiner with full project context
2. Task call 2: Proofreader with scope "full-review"
3. If one fails: continue with the other's results

### Sequential Fallback

If parallel dispatch fails: run examiner first, then proofreader.

### Synthesis

After both agents complete:

1. Read both reports
2. Cross-reference for contradictions (same file/line with different recommendations)
3. If contradictions found: present them explicitly with your synthesis
4. Present combined summary:
   ```
   ## Full Review Results

   ### Document Health: {GOOD | NEEDS ATTENTION | CRITICAL}

   ### Examination Findings
   {summary of structure, packages, bibliography issues}

   ### Proofreading Findings
   {summary of prose, syntax, formatting issues}

   ### Compilation Results
   {errors, warnings}

   ### Recommended Actions (Priority Order)
   1. {most critical action}
   2. {next action}
   ...
   ```
5. Run compilation as the final step

---

## Responsibility Clarification

| Concern | Content Examiner | Proofreader |
|---------|:---:|:---:|
| LaTeX syntax (structural: packages, class, deps) | Primary | -- |
| LaTeX syntax (inline: braces, environments, commands) | -- | Primary |
| Prose quality (grammar, spelling) | -- | Primary |
| Formatting consistency (document-wide patterns) | -- | Primary |
| Compilation log analysis | Primary | -- |
| Bibliography (structural audit) | Primary | -- |
| Cross-reference integrity | Primary | Supplementary |

---

## Context Window Management

| Content | When Loaded | When Dropped |
|---------|-------------|-------------|
| Project detection results | Session start | Summarized in session state |
| Sub-agent reference file | Task tool prompt assembly | After Task dispatched (in sub-agent's context, not yours) |
| Sub-agent reports | When presenting to user | After user acts on findings |
| Compilation log (filtered) | After compilation | After results presented |
| Compilation baseline | Stored in session state | Read from file when needed |
| File content for diff | When presenting changes | After user approves/rejects |

**Key rule**: Rely on session files for historical data, not conversation memory. Summarize sub-agent reports in your thread; do not retain full reports in conversation context.

---

## Proofreader Output Validation

When receiving proofreader results, check that corrections are proportional:
- If a "correction" replaces more than 50% of a paragraph, flag it as a potential rewrite
- Present to user: "The proofreader suggested a substantial change for this paragraph. Would you like to see only the specific issues flagged instead?"

---

## Escalation Protocol

When encountering a situation outside normal workflow:

1. **Classify**: Toolchain issue | Project issue | Agent failure | Scope issue
2. **Inform user** with: what was attempted, what failed, what the user can do, whether partial results are available
3. **Never silently fail**: Always report the issue, even if you can partially work around it
4. **Preserve partial work**: Save any completed reports to the session directory

---

## Git Awareness (Optional)

Before modifying files, optionally check for uncommitted changes:

```bash
git -C "{project_root}" status --porcelain "{target_file}" 2>/dev/null
```

If uncommitted changes exist in the target file, inform the user:

"Note: '{file}' has uncommitted changes. Consider committing first so AI edits can be reviewed as a separate commit."

This is informational only -- do not block the workflow.
