# Example: CV Editing Workflow

This example demonstrates the latex-document-manager skill editing a CV built with the `curve` document class.

---

## User Request

> "Use latex-document-manager to add a new position to my CV at ~/repos/cv"

---

## Step 1: Project Detection

The orchestrator scans `~/repos/cv/` for `.tex` files containing `\documentclass`.

```
[Detect - ~/repos/cv/] Scanning project structure

Found LaTeX project at ~/repos/cv
  Main file: cv-llt.tex (curve class)
  Modules: employment.tex, education.tex, publications.tex, talks.tex, service.tex, skills.tex, leadership.tex
  Style: settings.sty
  Bibliography: own-bib.bib (15 entries)
  Engine: pdfLaTeX (detected from .latexmkrc)
```

The orchestrator auto-routes to Write/Edit because the user's request ("add a new position") clearly maps to content authoring.

---

## Step 2: Baseline Compilation

Before any changes, the orchestrator captures a compilation baseline:

```bash
cd ~/repos/cv && latexmk -pdf -file-line-error -interaction=nonstopmode -max-repeat=5 cv-llt.tex 2>&1
```

Log parsing results:
```
Compilation SUCCESS
Errors: 0
Warnings: 2
  - Underfull \hbox (badness 10000) in paragraph at lines 12--15 (employment.tex)
  - Underfull \hbox (badness 10000) in paragraph at lines 8--11 (education.tex)
```

Baseline stored: 0 errors, 2 warnings.

---

## Step 3: Writing Expert Dispatch

The orchestrator dispatches the writing expert via Task tool:

```
Task: Add a new employment entry for the position "Research Scientist at Anthropic, 2024-present, working on AI safety research and Claude model development."

Context:
- Project root: /Users/david/repos/cv
- Main file: cv-llt.tex
- Document class: curve
- Target file: /Users/david/repos/cv/employment.tex
- Style files: /Users/david/repos/cv/settings.sty

Instructions: Read references/writing-expert-instructions.md for full writing protocol.

MANDATORY: Complete the Style Learning Protocol before proposing any changes.
```

---

## Step 4: Writing Expert Response

The writing expert reads `settings.sty` and `employment.tex`, completes the Style Learning Protocol, and returns:

```
### 1. Context
- Files read: settings.sty (145 lines), employment.tex (89 lines), cv-llt.tex (52 lines)
- Custom commands available:
  - \entry{dates}{title}{organization}{location}{description} (from curve class)
  - \makefield (contact info), \circled (icon), \prefixmarker (bullet styling)
- Style conventions detected:
  - Entries use \entry command with 5 arguments
  - Dates in "YYYY -- Present" or "YYYY -- YYYY" format
  - Descriptions are 2-3 sentences, no trailing period
  - Indentation: 2 spaces
- Document class: curve (local class not found; standard CTAN version)

### 2. Proposed Changes

#### Change 1: Add Anthropic employment entry
- **File**: /Users/david/repos/cv/employment.tex
- **Action**: MODIFY
- **Line range**: 1-3 (insert before first existing entry)
- **Before**:
  ```latex
  \begin{rubric}{Employment}
  \entry*[2022 -- 2024]
  ```
- **After**:
  ```latex
  \begin{rubric}{Employment}
  \entry*[2024 -- Present]
    {Research Scientist}
    {Anthropic}
    {San Francisco, CA}
    {AI safety research and Claude model development, focusing on
    alignment techniques and capability evaluation}
  \entry*[2022 -- 2024]
  ```
- **Rationale**: Follows existing \entry format with 5 arguments. Date format matches "YYYY -- Present" convention. Description is 2 lines without trailing period, matching existing entries.

### 3. Verification Notes
- Compilation should be run to verify the new entry does not cause spacing issues
- The curve class may produce an Underfull \hbox warning for short descriptions

### 4. Dependencies
- No dependencies; single-file change
```

---

## Step 5: User Approval

The orchestrator presents the diff to the user:

```
[Write - employment.tex] Presenting proposed changes for approval

The writing expert proposes adding a new entry at the top of employment.tex:

  + \entry*[2024 -- Present]
  +   {Research Scientist}
  +   {Anthropic}
  +   {San Francisco, CA}
  +   {AI safety research and Claude model development, focusing on
  +   alignment techniques and capability evaluation}

This follows the existing \entry format from the curve class.

Apply this change? (yes/no/modify)
```

User responds: "yes"

---

## Step 6: Apply and Compile

The orchestrator:
1. Backs up original `employment.tex` content (rollback point)
2. Applies the change using Edit tool
3. Runs compilation:

```
[Compile - pdfLaTeX] Running latexmk, parsing log

Compilation SUCCESS
Errors: 0
Warnings: 3
  - Underfull \hbox (badness 10000) in paragraph at lines 4--8 (employment.tex) [NEW]
  - Underfull \hbox (badness 10000) in paragraph at lines 18--21 (employment.tex)
  - Underfull \hbox (badness 10000) in paragraph at lines 8--11 (education.tex)

Pre-existing: 0 errors, 2 warnings
After changes: 0 errors, 3 warnings
Delta: +1 warning (Underfull \hbox in new entry -- typical for curve class short descriptions)

No new errors. Change is safe.
```

---

## Step 7: PDF Preview

```bash
open ~/repos/cv/cv-llt.pdf
```

The orchestrator opens the PDF for visual inspection:

"PDF opened in your default viewer. Please verify the new entry appears correctly. The new Underfull hbox warning is typical for the curve class with shorter description text and does not affect visual output."

---

## Edge Cases Demonstrated

1. **Baseline comparison**: The pre-existing 2 warnings are tracked so they are not confused with new issues after the edit.
2. **Style Learning Protocol**: The writing expert reads `settings.sty` and existing entries before proposing content, ensuring the `\entry` command is used correctly with the right date format.
3. **Underfull hbox awareness**: The new warning is identified as typical for the curve class and reported as non-critical.
