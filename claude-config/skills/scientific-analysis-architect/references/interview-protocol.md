# Interview Protocol

Specification for Phase 6 statistical fact-checking interview mode.

## Contents

- Initial Presentation
- Concern Template
- User Response Handling
- Batch Options
- Rejection Handling
- Correction Application
- Progress and Resume

---

## Initial Presentation

Concern count sets the entry point: zero concerns means showing the justification and finishing; 5 or fewer goes straight to one-at-a-time interview; more than 5 shows the summary and batch options first.

### Zero concerns

When no statistical concerns are found, the fact-checker states why:

```
Statistical Review Complete

No statistical concerns identified.

Justification:
- {N} analysis documents reviewed
- {M} statistical methods examined
- All tests match data types correctly
- Multiple testing corrections present where needed
- Assumptions documented appropriately

Confidence: {high/medium/low}

Proceed to completion? [yes/request-second-review]
```

If confidence is low or the user asks, re-run the fact-checker with stricter criteria; otherwise proceed with the warning logged.

### Summary presentation (>5 concerns)

```
Statistical Review Found {N} Concerns

Summary by Severity:
- Critical: {X} (incorrect conclusions possible)
- Standard: {Y} (best practice violations)
- Minor: {Z} (improvement opportunities)

Summary by Chapter:
- Chapter 1: {A} concerns
- Chapter 2: {B} concerns
- Chapter 3: {C} concerns

How would you like to proceed?

(A) Review all {N} concerns one-by-one
(B) Review only Critical ({X}) concerns
(C) Accept all Critical, review Standard ({Y})
(D) Accept all recommended corrections
(E) Skip statistical review (not recommended)

Enter choice [A/B/C/D/E]:
```

---

## Concern Template

```
Statistical Concern {current} of {total}

Document: {document_path}
Section: {section_path}
Code Block: {code_block_index}
Severity: {Critical|Standard|Minor}

Issue: {issue_description}

Current:
{current_content}

Concern: {detailed_explanation}

Recommendation:
{recommended_fix}

Accept? [yes/no/skip/explain]
```

Section paths use hierarchical notation to disambiguate duplicate headings: `"Analysis Steps > Step 4: Differential Expression"`.

### Example concern

```
Statistical Concern 1 of 7

Document: chapter1_data-atlas/analysis1_2_differential-expression.md
Section: Analysis Steps > Step 4: Differential Expression
Code Block: 0
Severity: Critical

Issue: Multiple testing correction missing

Current:
# Run t-test for each gene
# p_values = [ttest(group1, group2) for gene in genes]
# significant = p_values < 0.05

Concern: Testing 20,000+ genes without correction will produce
~1,000 false positives at alpha=0.05. This could lead to
incorrect biological conclusions.

Recommendation:
# Run t-test for each gene
# p_values = [ttest(group1, group2) for gene in genes]
# Apply Benjamini-Hochberg correction
# adjusted_p = multipletests(p_values, method='fdr_bh')[1]
# significant = adjusted_p < 0.05

Accept? [yes/no/skip/explain]
```

---

## User Response Handling

Each response records `user_decision` and `decision_timestamp` on the concern, appends it to the matching list in `session_state["corrections"]`, and moves to the next concern with the remaining count shown.

- **yes** — decision "accepted".
- **no** — decision "rejected". Ask for a brief reason (optional, encouraged) and store it as `rejection_reason`.
- **skip** — decision "skipped", deferred rather than decided.
- **explain** — set `explanation_requested`, then present statistical background, why it matters, alternative approaches, and references, and re-prompt with [yes/no/skip] only.

---

## Batch Options

Once 5 concerns have been reviewed and more remain:

```
You've reviewed 5 of {total} concerns.

Remaining by severity:
- Critical: {X}
- Standard: {Y}
- Minor: {Z}

Continue options:
(A) Continue one-by-one ({remaining} remaining)
(B) Accept all remaining Critical, skip others
(C) Accept all remaining
(D) Reject all remaining
(E) Accept Critical+Standard, skip Minor

Enter choice [A/B/C/D/E]:
```

Batch decisions apply to every remaining concern: B accepts critical and skips the rest, C accepts all, D rejects all, E accepts critical and standard while skipping minor.

---

## Rejection Handling

If the user rejects 80% or more of concerns:

```
High Rejection Rate Detected

You've rejected {X} of {Y} concerns ({Z}%).

This may indicate:
- Domain knowledge the fact-checker lacks
- Different statistical philosophy
- Concerns not applicable to your context

Would you like to:
(A) Continue reviewing (your expertise is respected)
(B) Explain your reasoning (helps improve future reviews)
(C) Skip remaining statistical review

Enter choice [A/B/C]:
```

Group rejections by category (multiple testing, test selection, assumption checking, effect size, other) with the user's stated reason, and store them in `session_state["rejection_patterns"]`.

---

## Correction Application

```
Interview Complete

Summary:
- {X} corrections accepted
- {Y} corrections rejected
- {Z} corrections skipped

Accepted corrections:
1. [{severity}] {brief_description} (document {path})
2. [{severity}] {brief_description} (document {path})
...

Documents affected: {list}

Apply all accepted corrections now? [yes/no]
```

Apply corrections grouped by document. For each one, locate the section by its hierarchical path and confirm `current_content` still appears in the file before replacing it — this guards against stale references from earlier edits. Replace the first occurrence only and prepend an HTML comment recording the issue: `<!-- CORRECTED: {issue} -->`. Report the count applied per document.

### Correction manifest

Save all decisions to `corrections-manifest.json`:

```json
{
  "session_id": "session-20260204-143022-12345",
  "interview_completed": "2026-02-04T15:45:00Z",
  "total_concerns": 7,
  "summary": {
    "accepted": 4,
    "rejected": 2,
    "skipped": 1
  },
  "corrections": {
    "accepted": [
      {
        "id": 1,
        "document": "chapter1/analysis1_2_differential-expression.md",
        "section_path": "Analysis Steps > Step 4: Differential Expression",
        "code_block_index": 0,
        "severity": "critical",
        "issue": "Multiple testing correction missing",
        "recommendation": "# Apply BH correction...",
        "applied": true,
        "applied_at": "2026-02-04T15:46:00Z"
      }
    ],
    "rejected": [
      {
        "id": 3,
        "document": "chapter2/analysis2_1_differential-expression.md",
        "section_path": "Analysis Steps > Step 3: Run DE",
        "code_block_index": 0,
        "severity": "standard",
        "issue": "Consider Mann-Whitney instead of t-test",
        "rejection_reason": "Data confirmed to be normally distributed"
      }
    ],
    "skipped": [
      {
        "id": 7,
        "document": "chapter3/analysis3_2_gene-regulatory-networks.md",
        "section_path": "Analysis Steps > Step 5: Network Inference",
        "code_block_index": 0,
        "severity": "minor",
        "issue": "Effect size not reported"
      }
    ]
  }
}
```

---

## Progress and Resume

Show progress after each response:

```
Correction accepted. (9 remaining)
[======>             ] 4/12 concerns reviewed
Accepted: 2 | Rejected: 1 | Skipped: 0
```

If the interview is interrupted, save `interview_progress` into session state with `current_index`, `total_concerns`, `decisions_made`, and `last_saved`, and checkpoint Phase 6. On resume:

```
Resuming Statistical Review

Progress: {X} of {Y} concerns reviewed
Last concern reviewed: {last_concern_summary}

Continue from concern {X+1}? [yes/restart]
```
