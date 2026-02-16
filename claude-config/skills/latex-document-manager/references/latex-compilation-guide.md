# LaTeX Compilation Guide

Shared reference for compilation commands, log parsing, and macOS TeX integration. Used by both the orchestrator (compilation execution) and the content examiner (log analysis).

---

## 1. TeX Live Detection on macOS

Claude Code's Bash tool may not inherit the user's full PATH. Always detect TeX binaries explicitly.

### Detection Priority Order

Check these paths in order. Use the first one that succeeds.

```bash
# Step 1: Check if latexmk is already in PATH
which latexmk 2>/dev/null

# Step 2: Standard MacTeX location
/Library/TeX/texbin/latexmk --version 2>/dev/null

# Step 3: Direct TeX Live 2025 path
/usr/local/texlive/2025/bin/universal-darwin/latexmk --version 2>/dev/null

# Step 4: Homebrew location
/opt/homebrew/bin/latexmk --version 2>/dev/null
```

### Engine Verification

After finding latexmk, verify that the detected engine binary also exists:

```bash
# After latexmk is found, verify the engine is available
which pdflatex 2>/dev/null || "${TEX_BIN_DIR}/pdflatex" --version 2>/dev/null
```

If latexmk is found but the engine binary is not, report:

```
latexmk was found but pdflatex is not available. You may have latexmk installed
without a full TeX distribution.

Install TeX Live: https://www.tug.org/mactex/
Homebrew: brew install --cask mactex
```

The orchestrator should enter degraded mode (examination and proofreading only).

### PATH Setup

Once the correct directory is found, prepend it to PATH for the session:

```bash
export PATH="/Library/TeX/texbin:$PATH"
```

Store the discovered directory as `TEX_BIN_DIR` in session state. All subsequent commands use this PATH.

### If Not Found

If none of the paths above yield a working `latexmk`, report:

```
TeX Live does not appear to be installed or is not in any standard macOS location.

Install options:
- MacTeX: https://www.tug.org/mactex/
- Homebrew: brew install --cask mactex

After installation, restart Claude Code so PATH changes take effect.
```

The orchestrator can still run examination and proofreading in degraded mode (no compilation).

---

## 2. Engine Detection Priority Order

Determine which LaTeX engine to use. Check in this order; use the first match.

### Priority 1: .latexmkrc Configuration

```bash
# Check project root for latexmkrc
if [ -f ".latexmkrc" ] || [ -f "latexmkrc" ]; then
  grep -E '^\$pdf_mode|^\$pdflatex|^\$xelatex|^\$lualatex' .latexmkrc latexmkrc 2>/dev/null
fi
```

- `$pdf_mode = 1` or `$pdflatex` defined -> pdfLaTeX
- `$pdf_mode = 5` or `$xelatex` defined -> XeLaTeX
- `$pdf_mode = 4` or `$lualatex` defined -> LuaLaTeX

### Priority 2: Magic Comment

```bash
grep -m1 "^%!TEX program" "{main_file}" 2>/dev/null
```

- `%!TEX program = pdflatex` -> pdfLaTeX
- `%!TEX program = xelatex` -> XeLaTeX
- `%!TEX program = lualatex` -> LuaLaTeX

### Priority 3: fontspec Package Usage

```bash
grep -l '\\usepackage.*{fontspec}' "{main_file}" *.tex 2>/dev/null
```

If fontspec is loaded, the document requires XeLaTeX or LuaLaTeX. Default to XeLaTeX unless LuaLaTeX indicators are also present.

### Priority 4: Engine Conditionals

```bash
grep -l '\\ifxetex\|\\ifluatex\|\\ifxetexorluatex' "{main_file}" *.sty 2>/dev/null
```

If conditionals are present, the document supports multiple engines. Use pdfLaTeX as the default since it is the broadest compatibility choice.

### Priority 5: Default

If none of the above match, default to pdfLaTeX.

---

## 3. Compilation Commands

### Standard Invocations

Use the Bash tool's timeout parameter (set to 120000 ms) instead of a shell `timeout` command, which is not available by default on macOS.

| Engine | Command |
|--------|---------|
| pdfLaTeX | `latexmk -pdf -file-line-error -interaction=nonstopmode -max-repeat=5 "{main_file}"` |
| XeLaTeX | `latexmk -pdfxe -file-line-error -interaction=nonstopmode -max-repeat=5 "{main_file}"` |
| LuaLaTeX | `latexmk -pdflua -file-line-error -interaction=nonstopmode -max-repeat=5 "{main_file}"` |

### Flag Explanation

| Flag | Purpose |
|------|---------|
| `-pdf` / `-pdfxe` / `-pdflua` | Select engine and produce PDF output |
| `-file-line-error` | Format errors as `file:line: message` for easier parsing |
| `-interaction=nonstopmode` | Do not pause for user input on errors |
| `-max-repeat=5` | Prevent infinite recompilation loops |

### Working Directory

Always `cd` to the project root before running latexmk:

```bash
cd "{project_root}" && latexmk -pdf -file-line-error -interaction=nonstopmode -max-repeat=5 "{main_file}" 2>&1
```

### Clean Build (When Needed)

If compilation produces stale artifacts or unexplained errors:

```bash
cd "{project_root}" && latexmk -C "{main_file}" && latexmk -pdf -file-line-error -interaction=nonstopmode -max-repeat=5 "{main_file}" 2>&1
```

---

## 4. Log Parsing Protocol

### Critical Rule

**NEVER read the raw .log file directly with the Read tool.** LaTeX log files can be tens of thousands of lines and will consume the entire context window.

### Size Guard

```bash
wc -l "{logfile}"
```

If the log exceeds 5000 lines, use ONLY the grep extraction steps below. Never attempt to read the full file.

### Extraction Steps

Run these steps in order. Each step targets a specific category of log output.

**Step 1 -- Errors (highest priority)**:
```bash
grep -n "^!" "{logfile}" | head -50
```

**Step 2 -- File-line-error format** (from `-file-line-error` flag):
```bash
grep -n "^.*:[0-9]*:.*" "{logfile}" | head -50
```

**Step 3 -- LaTeX warnings**:
```bash
grep -n "^LaTeX Warning\|^LaTeX Font Warning\|^LaTeX3 Warning" "{logfile}" | head -100
```

**Step 4 -- Package and class warnings**:
```bash
grep -n "^Package\|^Class\|^Module" "{logfile}" | grep "Warning" | head -50
```

**Step 5 -- Box warnings** (overfull/underfull):
```bash
grep -n "Overfull\|Underfull" "{logfile}" | head -50
```

**Step 6 -- Missing references and citations**:
```bash
grep -n "Citation.*undefined\|Reference.*undefined\|multiply defined" "{logfile}" | head -50
```

**Step 6b -- Emergency stop / Fatal error**:
```bash
grep -n "Emergency stop\|Fatal error\|==> Fatal error" "{logfile}" | head -10
```

**Step 6c -- Missing characters** (font encoding issues):
```bash
grep -n "Missing character" "{logfile}" | head -20
```

**Step 6d -- Backend warnings** (biber/bibtex):
```bash
grep -n "^Biber\|^BibTeX" "{logfile}" | head -20
```

**Step 7 -- Summary** (last 20 lines contain error/warning counts):
```bash
tail -20 "{logfile}"
```

### Severity Classification

| Pattern | Severity | Action |
|---------|----------|--------|
| `^!` (error) | ERROR | Must fix before PDF is usable |
| `undefined` reference/citation | ERROR | Fix source or run biber/bibtex |
| `multiply defined` | WARNING | Review labels for duplicates |
| `Overfull \\hbox` (badness > 5000) | WARNING | Review line breaking |
| `Overfull \\hbox` (badness <= 5000) | INFO | Usually acceptable |
| `Underfull \\hbox` | INFO | Usually acceptable |
| Package/class warning | WARNING | Review, may be benign |
| Font warning | INFO | Usually substitution, acceptable |

---

## 5. Post-Compilation Verification

### Exit Code Check

- Exit code 0: Compilation succeeded (but check for warnings)
- Non-zero exit code: Compilation failed (check errors)

### Undefined Reference/Citation Check

Even with exit code 0, check for unresolved references:

```bash
grep -c "Citation.*undefined\|Reference.*undefined" "{logfile}"
```

If citations are undefined AND a `.bib` file exists in the project:

```bash
# Check if biber needs to run
grep -l "\\\\usepackage.*{biblatex}" "{main_file}" && echo "biblatex detected -- biber may need to run"
```

Suggest: "Undefined citations detected. This may resolve after running biber. Shall I recompile?"

### Multi-Pass Check

If `latexmk` reports "Latexmk: Run number N of rule 'pdflatex'" where N equals `max-repeat`, compilation may not have converged. Report this to the user.

---

## 6. Missing Package Detection

### Detection Pattern

```bash
grep -n "File '.*\.sty' not found\|File '.*\.cls' not found\|LaTeX Error: File .* not found" "{logfile}"
```

### Response

When a missing package is detected:

```
Missing LaTeX package: {package_name}

To install:
  sudo tlmgr install {package_name}

Then recompile. Do NOT remove the \usepackage line from the document.
```

**Never** suggest removing a `\usepackage` line as a fix for a missing package. The package is intentionally used by the document.

---

## 7. Common Error Diagnostic Reference

| Error | Likely Cause | Fix |
|-------|-------------|-----|
| `Undefined control sequence` | Typo in command name, or missing package | Check spelling; verify package is loaded |
| `Missing $ inserted` | Math-mode character outside math environment | Wrap in `$...$` or escape the character |
| `Mismatched braces` / `Extra }` | Unbalanced `{` and `}` | Count braces in the affected region |
| `Runaway argument` | Missing closing brace, often in `\begin{}` | Check for unclosed environments or arguments |
| `Missing \begin{document}` | Preamble error prevents document body | Fix the preamble error reported above this one |
| `File not found` | Missing .tex, .sty, .cls, or image file | Check file exists at the referenced path |
| `Environment undefined` | Missing package that defines the environment | Identify which package provides the environment |
| `Too many unprocessed floats` | Too many figures/tables without text between them | Add `\clearpage` or use `[H]` float specifier |
| `Dimension too large` | TikZ or geometry calculation overflow | Review coordinate values and page dimensions |
| `TeX capacity exceeded` | Infinite loop in macros or enormous table | Check for recursive macro definitions |

---

## 8. PDF Preview on macOS

### Open PDF

After successful compilation:

```bash
open "{main_file_stem}.pdf"
```

This opens the PDF in the default viewer (Preview.app, Skim, or Adobe Acrobat).

### Adobe Acrobat File Locking

Adobe Acrobat locks the PDF file while it is open. This prevents `latexmk` from overwriting it on recompilation.

If recompilation fails with a permission error on the PDF:
- Close the PDF in Acrobat before recompiling, OR
- Switch to Preview.app or Skim, which do not lock PDF files
- Skim auto-refreshes on file change, making it ideal for LaTeX workflows

### Path Quoting

Always double-quote the PDF path in case it contains spaces:

```bash
open "{project_root}/{main_file_stem}.pdf"
```
