---
name: latex-document-manager
last_updated: 2026-05-24
description: >
  Use when user asks to examine, edit, proofread, or compile LaTeX documents.
  Triggers on mentions of .tex files, pdflatex, xelatex, lualatex, latexmk,
  bibtex, biblatex, CV updates, or LaTeX paper editing on macOS.

# Handoff metadata (custom extension -- see workflow-coordinator/references/frontmatter-metadata-standard.md)
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

## Prerequisites

- TeX Live installed on macOS (latexmk, pdflatex at minimum)
- LaTeX project with at least one .tex file containing `\documentclass`
- macOS (for PDF preview via `open` command)

## Success Criteria

- Project structure correctly detected (main file, modules, packages, bibliography)
- Sub-agent reports are structured and actionable
- No file modifications without explicit user approval
- Post-change compilation detects regressions
- PDF preview opens after successful compilation

---

## Architecture

You are the orchestrator. You own user interaction, compilation, session state, and change approval. Three specialists receive context via the Task tool and return structured reports: the content examiner (structure, logs, bibliography), the writing expert (authoring and editing), and the proofreader (prose and inline syntax). Sub-agents never interact with the user directly.

---

## Workflow

1. **Pre-Flight Validation**: Check TeX installation, find main file, detect engine (blocking checks)
2. **Project Detection**: Enumerate .tex files, map dependencies, identify style and bibliography files
3. **Action Routing**: Present action menu or auto-route based on user request
4. **Delegate to Specialist(s)**: Dispatch appropriate sub-agent(s) via Task tool
5. **Present Results**: Show report to user, handle follow-up actions
6. **Compilation** (if applicable): Build PDF, compare against baseline, present delta

Read `references/latex-compilation-guide.md` (resolved to absolute path) before any compilation step to obtain engine detection rules and log parsing protocol.

---

## Delegation

You are an orchestrator. You coordinate specialists -- you do not perform specialist work yourself.

Delegate specialist work when the subtask is substantial and independent -- a full document examination, a content-writing task that needs the style-learning protocol, a prose review across several files. Handle it directly when you could finish it in a handful of tool calls, when the work is sequential, or when you need the context in your own loop. See `../references/delegation-and-scope.md`.

**Orchestrator-owned tasks**: Session setup, project detection, action routing, compilation execution, user communication, change approval, PDF preview.

---

## State Anchoring

Open responses with an action anchor -- `[Compile - pdfLaTeX] Running latexmk, parsing log`, `[Write - employment.tex] Presenting proposed changes for approval` -- and re-anchor when the action changes, after a compilation attempt, or after a sub-agent report comes back. It keeps a long editing session legible.

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

---

## Pre-Flight Validation

Run these checks before any workflow action. Report all results before proceeding.

### Required Checks (blocking if failed)

**1. LaTeX Installation**

Check for latexmk at: (1) PATH, (2) `/Library/TeX/texbin/`, (3) `/usr/local/texlive/2025/bin/universal-darwin/`, (4) `/opt/homebrew/bin/`.

If NOT FOUND: Report clearly with install link. Offer examination and proofreading only (degraded mode -- no compilation).

If latexmk IS found, verify the engine binary also exists:

```bash
which pdflatex 2>/dev/null || "${TEX_BIN_DIR}/pdflatex" --version 2>/dev/null
```

If the engine is not found: Report "latexmk was found but pdflatex is not available. You may have latexmk installed without a full TeX distribution." Offer degraded mode.

**2. Main File Detection**

Find `.tex` files containing `\documentclass` in the target directory.

- If NONE found: Report "No LaTeX main files found in {directory}." Suggest the user check the path or specify a file.
- If MULTIPLE found: Present a disambiguation menu showing each file with its document class, line count, and last modified date. Ask the user to select one.
- If user specified a file: Use that file directly.

**3. Engine Detection**

Follow the priority order from `references/latex-compilation-guide.md` Section 2.

### Optional Checks (warning if failed)

**4. Biber**: Only if biblatex is detected. Warn if biber is not found.

**5. ChkTeX**: Warn if not found (syntax linting will be skipped during examination).

**6. Write Permissions**: Check that target files are writable before offering edits. If read-only, inform the user and offer examination/proofreading only.

---

## Session Management

Create a session directory on first invocation and store its path for the rest of the conversation:

```bash
mkdir -p "/tmp/latex-document-manager-$(date +%Y%m%d-%H%M%S)-$$"
```

Track detected facts in `session-state.json` within it: `project_root`, `main_file`, `document_class`, `engine`, `tex_bin_dir`, `compilation_baseline` (errors, warnings, timestamp), `actions_completed`, `last_action`, `status`. Update it after each action. On a later invocation, detect an existing session directory and offer to resume. The directory is cleaned up when the user ends the session or after 24 hours.

Rely on these session files for historical data rather than conversation memory: summarize sub-agent reports in your thread instead of retaining them in full, and read the compilation baseline back from the file when you need it.

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

If only one `.tex` file with no `\input`/`\include` dependencies, skip full project enumeration and offer a streamlined menu: [1] Proofread [2] Edit [3] Compile.

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

If the request clearly maps to one action, auto-route without the menu. After any Write/Edit action with approved changes, compile to check for regressions. After a Full Review, present a combined report.

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

Instructions: Read the content examiner instructions and the compilation guide. Resolve these to absolute paths before dispatching:
- Content examiner instructions: {skill_dir}/references/content-examiner-instructions.md
- Compilation guide: {skill_dir}/references/latex-compilation-guide.md

where {skill_dir} is determined by reading this SKILL.md's own absolute path and extracting its directory.

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

Instructions: Read the writing expert instructions, resolved to an absolute path before dispatching:
- Writing expert instructions: {skill_dir}/references/writing-expert-instructions.md

Complete the Style Learning Protocol before proposing any changes.

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

Instructions: Read the proofreader instructions, resolved to an absolute path before dispatching:
- Proofreader instructions: {skill_dir}/references/proofreader-instructions.md

Do NOT rewrite content. Flag issues with location, type, and brief correction only.

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
2. Compare against baseline: count new errors, new warnings, resolved warnings
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
| G-PREFLIGHT | After pre-flight validation | TeX tools found (or degraded mode accepted) | Report missing tools, offer degraded mode |
| G-WRITE | After writing expert proposes | Changes compile without new errors | Reject changes, return to expert with error context |
| G-COMPILE | After compilation | Exit code 0, log parsed | Present errors, offer diagnosis |
| G-APPROVE | Before applying file changes | User has reviewed diff and approved | Do not apply; ask for instructions |

G-APPROVE is not optional: no file is modified without the user seeing the diff first.

### G-WRITE Failure Recovery

If the writing expert's proposed changes introduce compilation errors:

1. Extract relevant error messages from the compilation log (use grep extraction Steps 1-2)
2. Dispatch a new Task to the writing expert with the original task description, the proposed changes that failed, the relevant compilation-error excerpt, and the instruction "Revise proposed changes to resolve these compilation errors while preserving the original intent"
3. After 2 retry cycles, report to the user: "The writing expert was unable to produce changes that compile cleanly after 2 attempts. Here are the latest proposed changes and the remaining errors."

---

## Error Handling

Retry a failed sub-agent Task once. Do not auto-retry a failed compilation -- present the errors to the user instead. Retry a failed file read once, then report the file as inaccessible. If two or more sub-agents fail in one session, stop auto-retrying and ask the user whether to retry, proceed with partial results, or exit.

Degrade rather than stop: without TeX Live, offer examination and proofreading; if the content examiner fails, proceed with the other actions; if the proofreader fails, say so and let the user review; if compilation keeps failing, edit the source and let the user compile manually.

### Rollback Protocol

Before applying any file modification:

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

### CREATE Action Safety Check

Before executing any CREATE action from the writing expert:

```bash
test -f "{filepath}" && echo "EXISTS" || echo "NEW"
```

If the file already exists:
1. Warn the user: "The file '{filepath}' already exists. The writing expert proposed creating it as a new file."
2. Offer options: (a) View the existing content first, (b) Overwrite it (backs up original first), (c) Abort this change
3. If overwriting, store the original content as a rollback point before proceeding

### APPEND Action Handling

APPEND actions from the writing expert are handled as MODIFY operations where the insertion point is at the end of the file or at the line specified in the Line range field. No Before block is required for APPEND -- the After block contains the content to append.

### Path Handling

ALWAYS double-quote file paths in Bash commands:
```bash
latexmk -pdf "${PROJECT_DIR}/main.tex"    # correct
latexmk -pdf ${PROJECT_DIR}/main.tex      # WRONG -- breaks on spaces
```

---

## Full Review Pipeline

Dispatch the content examiner and the proofreader simultaneously, as two Task tool calls in a single response (examiner with full project context, proofreader with scope "full-review"). If one fails, continue with the other's results; if parallel dispatch itself fails, run them sequentially. See `../references/delegation-and-scope.md` for scoping fan-out.

Then read both reports, cross-reference for contradictions (same file and line with different recommendations, presented explicitly with your synthesis), and give the user a combined summary: overall document health, examination findings (structure, packages, bibliography), proofreading findings (prose, syntax, formatting), compilation results, and recommended actions in priority order. Run compilation as the final step.

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

## Proofreader Output Validation

When receiving proofreader results, check that corrections are proportional. If a "correction" replaces more than 50% of a paragraph, flag it as a potential rewrite and ask the user: "The proofreader suggested a substantial change for this paragraph. Would you like to see only the specific issues flagged instead?"

---

## Escalation Protocol

For situations outside the normal workflow, classify the problem (toolchain, project, agent failure, scope), then tell the user what was attempted, what failed, what they can do, and whether partial results exist. Report the issue even when you can partially work around it, and save completed reports to the session directory.

---

## Git Awareness (Optional)

Before modifying files, optionally check for uncommitted changes:

```bash
git -C "{project_root}" status --porcelain "{target_file}" 2>/dev/null
```

If uncommitted changes exist in the target file, inform the user:

"Note: '{file}' has uncommitted changes. Consider committing first so AI edits can be reviewed as a separate commit."

This is informational only -- do not block the workflow.
