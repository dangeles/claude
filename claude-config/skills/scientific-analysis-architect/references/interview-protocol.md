# Interview Protocol

Specification for Phase 6 statistical fact-checking interview mode.

## Initial Presentation

### Decision Logic

```python
def present_concerns(concerns: list):
    """Determine how to present statistical concerns."""

    total = len(concerns)

    if total == 0:
        # No concerns - show justification
        return present_zero_concerns()
    elif total <= 5:
        # Few concerns - go straight to interview
        return start_interview_mode(concerns)
    else:
        # Many concerns - show summary first
        return present_summary_then_interview(concerns)
```

### Zero Concerns Handling

When no statistical concerns are found, the fact-checker must provide justification:

```
Statistical Review Complete

No statistical concerns identified.

Justification:
- {N} notebooks reviewed
- {M} statistical methods examined
- All tests match data types correctly
- Multiple testing corrections present where needed
- Assumptions documented appropriately

Confidence: {high/medium/low}

Proceed to completion? [yes/request-second-review]
```

If confidence is low or user requests:
- Re-run fact-checker with stricter criteria
- Or proceed with warning logged

### Summary Presentation (>5 concerns)

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

## Concern Template

### Standard Format

```
Statistical Concern {current} of {total}

Notebook: {notebook_path}
Cell: {cell_number}
Severity: {Critical|Standard|Minor}

Issue: {issue_description}

Current pseudocode:
```
{current_code}
```

Concern: {detailed_explanation}

Recommendation:
```
{recommended_fix}
```

Accept? [yes/no/skip/explain]
```

### Example Concern

```
Statistical Concern 1 of 7

Notebook: chapter1_data-atlas/notebook1_2_differential-expression.ipynb
Cell: 5
Severity: Critical

Issue: Multiple testing correction missing

Current pseudocode:
```python
# Run t-test for each gene
# p_values = [ttest(group1, group2) for gene in genes]
# significant = p_values < 0.05
```

Concern: Testing 20,000+ genes without correction will produce
~1,000 false positives at alpha=0.05. This could lead to
incorrect biological conclusions.

Recommendation:
```python
# Run t-test for each gene
# p_values = [ttest(group1, group2) for gene in genes]
# Apply Benjamini-Hochberg correction
# adjusted_p = multipletests(p_values, method='fdr_bh')[1]
# significant = adjusted_p < 0.05
```

Accept? [yes/no/skip/explain]
```

## User Response Handling

### Response: "yes"

```python
def handle_yes(concern: dict, session_state: dict):
    """User accepts correction."""

    concern["user_decision"] = "accepted"
    concern["decision_timestamp"] = datetime.now().isoformat()

    session_state["corrections"]["accepted"].append(concern)

    print(f"Correction accepted. ({remaining} remaining)")
    return next_concern()
```

### Response: "no"

```python
def handle_no(concern: dict, session_state: dict):
    """User rejects correction."""

    # Ask for reason (optional but encouraged)
    reason = AskUserQuestion(
        "Briefly explain why (or press Enter to skip): "
    )

    concern["user_decision"] = "rejected"
    concern["rejection_reason"] = reason if reason else None
    concern["decision_timestamp"] = datetime.now().isoformat()

    session_state["corrections"]["rejected"].append(concern)

    print(f"Noted. ({remaining} remaining)")
    return next_concern()
```

### Response: "skip"

```python
def handle_skip(concern: dict, session_state: dict):
    """User defers decision."""

    concern["user_decision"] = "skipped"
    concern["decision_timestamp"] = datetime.now().isoformat()

    session_state["corrections"]["skipped"].append(concern)

    print(f"Skipped. ({remaining} remaining)")
    return next_concern()
```

### Response: "explain"

```python
def handle_explain(concern: dict, session_state: dict):
    """User requests more detail."""

    concern["explanation_requested"] = True

    # Generate expanded explanation
    explanation = generate_explanation(concern)

    print(f"""
Extended Explanation:

{explanation["statistical_background"]}

Why this matters:
{explanation["impact"]}

Alternative approaches:
{explanation["alternatives"]}

References:
{explanation["references"]}

Accept? [yes/no/skip]
""")

    # Re-prompt (no second explain option)
    return handle_response(get_user_input(["yes", "no", "skip"]))
```

## Batch Options

### After 5 Concerns Reviewed

If total concerns > 5 and user has reviewed 5, offer batch options:

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

### Batch Processing

```python
def apply_batch_decision(decision: str, remaining: list):
    """Apply batch decision to remaining concerns."""

    if decision == "B":  # Accept critical only
        for c in remaining:
            if c["severity"] == "critical":
                c["user_decision"] = "accepted"
            else:
                c["user_decision"] = "skipped"

    elif decision == "C":  # Accept all
        for c in remaining:
            c["user_decision"] = "accepted"

    elif decision == "D":  # Reject all
        for c in remaining:
            c["user_decision"] = "rejected"

    elif decision == "E":  # Accept critical+standard
        for c in remaining:
            if c["severity"] in ["critical", "standard"]:
                c["user_decision"] = "accepted"
            else:
                c["user_decision"] = "skipped"
```

## Rejection Handling

### High Rejection Rate Detection

If user rejects >= 80% of concerns:

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

### Recording Rejection Patterns

```python
def log_rejection_pattern(session_state: dict):
    """Log rejection patterns for future improvement."""

    rejected = session_state["corrections"]["rejected"]

    patterns = {
        "multiple_testing": [],
        "test_selection": [],
        "assumption_checking": [],
        "effect_size": [],
        "other": []
    }

    for concern in rejected:
        category = classify_concern(concern)
        patterns[category].append({
            "concern": concern["issue"],
            "reason": concern.get("rejection_reason")
        })

    session_state["rejection_patterns"] = patterns
```

## Correction Application

### Summary Before Application

```
Interview Complete

Summary:
- {X} corrections accepted
- {Y} corrections rejected
- {Z} corrections skipped

Accepted corrections:
1. [{severity}] {brief_description} (notebook {path})
2. [{severity}] {brief_description} (notebook {path})
...

Notebooks affected: {list}

Apply all accepted corrections now? [yes/no]
```

### Application Process

```python
def apply_corrections(session_state: dict):
    """Apply accepted corrections to notebooks."""

    accepted = session_state["corrections"]["accepted"]

    # Group by notebook
    by_notebook = group_corrections_by_notebook(accepted)

    for notebook_path, corrections in by_notebook.items():
        print(f"Updating {notebook_path}...")

        # Load notebook
        with open(notebook_path) as f:
            nb = nbformat.read(f, as_version=4)

        # Apply each correction
        for correction in corrections:
            cell_index = correction["cell_number"]
            new_code = correction["recommendation"]

            # Update cell
            nb.cells[cell_index].source = new_code

            # Add correction comment
            nb.cells[cell_index].source = (
                f"# CORRECTED: {correction['issue']}\n" +
                nb.cells[cell_index].source
            )

        # Validate
        nbformat.validate(nb)

        # Save
        with open(notebook_path, "w") as f:
            nbformat.write(nb, f)

        print(f"  Applied {len(corrections)} corrections")

    print(f"\nTotal: {len(accepted)} corrections applied to {len(by_notebook)} notebooks")
```

### Correction Manifest

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
        "notebook": "chapter1/notebook1_2.ipynb",
        "cell": 5,
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
        "notebook": "chapter2/notebook2_1.ipynb",
        "cell": 3,
        "severity": "standard",
        "issue": "Consider Mann-Whitney instead of t-test",
        "rejection_reason": "Data confirmed to be normally distributed"
      }
    ],
    "skipped": [
      {
        "id": 7,
        "notebook": "chapter3/notebook3_2.ipynb",
        "cell": 8,
        "severity": "minor",
        "issue": "Effect size not reported"
      }
    ]
  }
}
```

## Progress Indicators

### During Interview

```
Statistical Review Progress
[=====>              ] 3/12 concerns reviewed
Accepted: 2 | Rejected: 1 | Skipped: 0
Estimated time remaining: ~8 minutes
```

### After Each Response

```
Correction accepted. (9 remaining)
[======>             ] 4/12 concerns reviewed
```

## Interview State Persistence

If interview is interrupted, save state:

```python
def save_interview_progress(session_state: dict, current_index: int):
    """Save interview progress for resume."""

    session_state["interview_progress"] = {
        "current_index": current_index,
        "total_concerns": len(session_state["concerns"]),
        "decisions_made": current_index,
        "last_saved": datetime.now().isoformat()
    }

    checkpoint_phase(6, session_state, {})
```

Resume message:

```
Resuming Statistical Review

Progress: {X} of {Y} concerns reviewed
Last concern reviewed: {last_concern_summary}

Continue from concern {X+1}? [yes/restart]
```
