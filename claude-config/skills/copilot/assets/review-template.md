# Code Review Template

## Standard Review Format

```
CODE REVIEW: [File/Notebook Name]
Reviewer: Copilot
Date: [Date]

=== CRITICAL ISSUES (🔴) ===
[List all critical issues that MUST be fixed]

Issue 1:
  Location: [Line number or cell]
  Problem: [What's wrong]
  Impact: [What will happen]
  Fix: [Specific solution]

---

=== MAJOR ISSUES (🟠) ===
[List issues that should be fixed]

Issue 1:
  Problem: [Description]
  Suggestion: [How to improve]

---

=== MINOR SUGGESTIONS (🟡) ===
[List nice-to-have improvements]

---

=== GOOD PRACTICES (✅) ===
[Acknowledge what's done well]

✅ [Specific good practice observed]

---

VERDICT: [APPROVED | NEEDS REVISION]
[If NEEDS REVISION, state what must be fixed before approval]
```

## Example Review

```
CODE REVIEW: Gene_Expression_Analysis.ipynb
Reviewer: Copilot
Date: 2026-01-28

=== CRITICAL ISSUES (🔴) ===

🔴 Division by zero in normalization (Cell 5, Line 3)
  Location: Cell 5, Line 3
  Problem: cpm = (counts / counts.sum()) * 1e6
  Impact: Will fail if all counts are zero (returns inf)
  Fix: Add check:
    total = counts.sum()
    cpm = (counts / total * 1e6) if total > 0 else np.zeros_like(counts)

🔴 No multiple testing correction (Cell 8)
  Location: Cell 8
  Problem: sig_genes = genes[genes['p_value'] < 0.05]
  Impact: ~1000 false positives expected with 20,000 genes
  Fix: Apply FDR correction:
    from statsmodels.stats.multitest import multipletests
    _, genes['p_adj'], _, _ = multipletests(genes['p_value'], method='fdr_bh')
    sig_genes = genes[genes['p_adj'] < 0.05]

---

=== MAJOR ISSUES (🟠) ===

🟠 Missing random seed (Cell 10)
  Problem: PCA applied without random_state parameter
  Suggestion: Add: pca = PCA(n_components=10, random_state=42)

🟠 Unclear variable name (Cell 6)
  Problem: Variable "x" used for filtered gene list
  Suggestion: Rename to "filtered_genes" for clarity

---

=== MINOR SUGGESTIONS (🟡) ===

🟡 Add docstring to custom function (Cell 4)
🟡 Consider using seaborn for more polished plots
🟡 Export figures in PDF format for publication

---

=== GOOD PRACTICES (✅) ===

✅ Random seed set for reproducibility (Cell 2)
✅ Clear comments explaining biological context
✅ Proper axis labels on all figures
✅ Session info included at end
✅ Data validation before analysis (Cell 3)

---

VERDICT: NEEDS REVISION

Must fix:
1. Division by zero in normalization (CRITICAL)
2. Multiple testing correction (CRITICAL)

Once these are addressed, re-submit for approval.
```
