# Example: Academic Paper Full Review

This example demonstrates the latex-document-manager skill performing a Full Review on an academic paper, including parallel dispatch of the content examiner and proofreader.

---

## User Request

> "Use latex-document-manager to do a full review of my paper at ~/repos/papers/myproject/"

---

## Step 1: Project Detection

```
[Detect - ~/repos/papers/myproject/] Scanning project structure

Found LaTeX project at ~/repos/papers/myproject
  Main file: paper.tex (article class)
  Modules: introduction.tex, methods.tex, results.tex, discussion.tex, appendix.tex
  Style: none (standard article class)
  Bibliography: references.bib (47 entries)
  Engine: pdfLaTeX (default, no .latexmkrc or magic comment found)
```

The user requested "full review," so the orchestrator routes to the Full Review pipeline.

---

## Step 2: Baseline Compilation

```
[Compile - pdfLaTeX] Capturing compilation baseline

Compilation SUCCESS
Errors: 0
Warnings: 5
  - Citation 'martinez2024' undefined (paper.tex)
  - Overfull \hbox (badness 7142) at line 45 (results.tex)
  - LaTeX Warning: Float too large for page at line 112 (results.tex)
  - Package hyperref Warning: Token not allowed in PDF string
  - Underfull \hbox (badness 10000) at line 23 (discussion.tex)

Baseline stored: 0 errors, 5 warnings (1 undefined citation).
```

---

## Step 3: Parallel Dispatch

The orchestrator dispatches the content examiner and proofreader simultaneously via two Task tool calls in a single response.

**Task 1: Content Examiner**

```
Task: Examine the LaTeX project for structural issues, compilation problems, and bibliography consistency.

Context:
- Project root: /Users/david/repos/papers/myproject
- Main file: paper.tex
- Document class: article
- Engine: pdflatex
- Bibliography: references.bib

Instructions: Read references/content-examiner-instructions.md
Focus areas: full examination
```

**Task 2: Proofreader**

```
Task: Proofread the LaTeX document for prose quality and syntax correctness.

Context:
- Project root: /Users/david/repos/papers/myproject
- Main file: paper.tex
- Files to proofread: introduction.tex, methods.tex, results.tex, discussion.tex
- Scope: full-review
- ChkTeX results: not yet available

Instructions: Read references/proofreader-instructions.md
CRITICAL: Do NOT rewrite content.
```

---

## Step 4: Content Examiner Report (Summary)

```
*Note: This example shows an abbreviated report. A full report includes all 8 sections per the output schema in content-examiner-instructions.md.*

### 1. Project Summary
- Main file: paper.tex
- Document class: article
- Engine: pdfLaTeX
- Modules: 5 (introduction.tex, methods.tex, results.tex, discussion.tex, appendix.tex)
- Packages: 14
- Bibliography: references.bib (47 entries)

### 3. ChkTeX Results
| File | Line | Col | Warning# | Severity | Message |
|------|------|-----|----------|----------|---------|
| methods.tex | 23 | 15 | 2 | Warning | Non-breaking space (~) should have been used |
| results.tex | 67 | 8 | 2 | Warning | Non-breaking space (~) should have been used |

### 5. Bibliography Audit
biber --validate-datamodel output:
- WARNING: Entry 'martinez2024' not found in references.bib
- WARNING: Entry 'chen2023' has missing 'doi' field
- INFO: 12 entries in .bib are not cited in the document

Cross-reference check:
- \cite{martinez2024} on line 34 of introduction.tex -- key not in .bib (ERROR)
- 12 uncited .bib entries (INFO)

### 7. Statistics
- Errors: 1
- Warnings: 4
- Info: 13

### 8. Recommendations
1. [ERROR] Add 'martinez2024' entry to references.bib or correct the citation key
2. [WARNING] Add non-breaking spaces before \ref commands in methods.tex:23 and results.tex:67
3. [WARNING] Add 'doi' field to chen2023 entry in references.bib
4. [INFO] Consider removing 12 uncited bibliography entries to keep .bib clean
```

---

## Step 5: Proofreader Report (Summary)

```
### 1. Summary
- Files reviewed: introduction.tex, methods.tex, results.tex, discussion.tex
- Total findings: 8
- Errors: 2
- Warnings: 3
- Suggestions: 3
- Scope used: full-review

### 2. Findings

#### Finding 1
- **Severity**: ERROR
- **Category**: prose
- **File**: methods.tex
- **Line**: 45
- **Current**: `The data were collected and it shows`
- **Suggested**: `it shows` -> `they show`
- **Explanation**: Subject-verb disagreement; "data" is plural, matching "were collected"

#### Finding 2
- **Severity**: ERROR
- **Category**: syntax
- **File**: results.tex
- **Line**: 89
- **Current**: `\begin{figure}...\end{fgure}`
- **Suggested**: `\end{fgure}` -> `\end{figure}`
- **Explanation**: Misspelled environment closing tag

#### Finding 3
- **Severity**: WARNING
- **Category**: typography
- **File**: introduction.tex
- **Line**: 12
- **Current**: `results from 2020-2024`
- **Suggested**: `2020-2024` -> `2020--2024`
- **Explanation**: Number ranges should use en-dash (--)

...additional findings omitted for brevity...

### 3. Style Observations
- Date ranges inconsistently use hyphens and en-dashes across sections
- Some \cite commands use ~ prefix and others do not
```

---

## Step 6: Combined Report

The orchestrator synthesizes both reports:

```
[Full Review - paper.tex] Presenting combined review results

## Full Review Results

### Document Health: NEEDS ATTENTION

### Examination Findings
- 1 undefined citation: martinez2024 (not in references.bib)
- 2 missing non-breaking spaces before \ref
- 1 missing DOI field in bibliography
- 12 uncited bibliography entries

### Proofreading Findings
- 2 errors: subject-verb disagreement (methods.tex:45), misspelled \end{figure} (results.tex:89)
- 3 warnings: en-dash usage, non-breaking space inconsistency, capitalization inconsistency
- 3 suggestions: minor clarity and style improvements

### Compilation Results
- 0 errors, 5 warnings (1 undefined citation, 1 overfull hbox, 1 float warning, 1 hyperref warning, 1 underfull hbox)

### Recommended Actions (Priority Order)
1. Fix misspelled \end{fgure} -> \end{figure} in results.tex:89 (syntax error)
2. Add 'martinez2024' to references.bib or correct the citation key in introduction.tex:34
3. Fix subject-verb disagreement in methods.tex:45
4. Standardize en-dash usage for number ranges across all sections
5. Add non-breaking spaces before \ref commands
6. Add DOI field to chen2023

Would you like me to fix any of these issues? I can route specific fixes to the writing expert.
```

---

## Step 7: User Follow-Up

> "Fix the misspelled end{figure} and the subject-verb issue."

The orchestrator routes to the writing expert:

```
[Write - results.tex, methods.tex] Dispatching writing expert for targeted fixes

Task: Fix two specific issues:
1. In results.tex line 89: change \end{fgure} to \end{figure}
2. In methods.tex line 45: change "it shows" to "they show"
```

The writing expert returns proposed changes. The orchestrator presents diffs for approval.

After user approves, the orchestrator applies changes and recompiles:

```
[Compile - pdfLaTeX] Verifying changes

Compilation SUCCESS
Errors: 0
Warnings: 4 (was 5)

Pre-existing: 0 errors, 5 warnings
After changes: 0 errors, 4 warnings
Delta: -1 warning (Float too large resolved by fixing the figure environment)

No new errors. Changes verified.
```

```bash
open ~/repos/papers/myproject/paper.pdf
```

---

## Edge Cases Demonstrated

1. **Parallel dispatch**: Content examiner and proofreader run simultaneously, reducing total review time.
2. **Undefined citation handling**: The `martinez2024` citation key is detected by both the compilation log (undefined citation warning) and the content examiner's cross-reference check (key not in .bib). The combined report deduplicates this finding.
3. **Proofreader respects scope**: With `full-review` scope, the proofreader reports suggestions in addition to errors and warnings. With `errors-only`, findings 3-8 would have been omitted.
4. **Post-fix recompilation**: After the targeted fixes, recompilation shows the warning count decreased, confirming no regressions.
5. **Synthesis of reports**: The orchestrator cross-references both reports and presents a unified priority list rather than two separate reports.
