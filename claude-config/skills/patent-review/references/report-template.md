# Patent Review Report Template

Use this exact structure when compiling the Phase 4 report. Scale each section to the findings
from the session — don't include empty sections.

---

```markdown
# Patent Review Report: [Patent Title]
**Date**: YYYY-MM-DD
**Reviewer role**: Scientific review (technical accuracy and claim coverage)

## Executive Summary
[2–4 sentences on the invention and overall assessment]
**Critical issues**: N | **Minor issues**: N

## Section-by-Section Findings
[Grouped by section. Each finding on its own line:]
- [SEVERITY][TAG] Section X — Issue description → Resolution (accepted/modified/dismissed)

## Claims Analysis
| Claim # | Type | Status | Notes |
|---------|------|--------|-------|
[Every claim listed. Status: accepted / revised / flagged / new (proposed)]

**Proposed new claims:**
[Any new claim language drafted during the gap pass — use full formal claim text]

## Recommendations
[Prioritized list — critical issues first, then minor. Each item actionable for the attorney.]

## Literature
[Populated only if Phase 5 was triggered. Each result with a note on whether it materially
affects a specific claim. Omit this section entirely if Phase 5 was not triggered.]
```

---

## Notes on Each Section

**Executive Summary**: Lead with what the invention is in plain language, then give the overall
assessment. State the critical/minor counts as a quick triage signal for the attorney.

**Section-by-Section Findings**: Group by document section (Background, Detailed Description,
Claims, Abstract). Use the `[CRITICAL][ACCURACY]`, `[MINOR][GAP]`, `[CRITICAL][LANGUAGE]` labeling
from Phase 3. Every logged issue — including dismissed ones — should appear, with its resolution.

**Claims Analysis table**: Every original claim must have a row. Don't omit claims that were
accepted unchanged — list them with `Status: Accepted`. Proposed new claims from the gap pass go in
the "Proposed new claims" subsection below the table.

**Recommendations**: This is what the attorney acts on. Write each item as a directive, not a
question. "Add a composition-of-matter claim covering [X]" not "Should there be a claim for [X]?"
Critical items first; minor items follow.

**Literature**: Only populate if Phase 5 was triggered. For each result, note the specific claim it
bears on and whether it strengthens, weakens, or is neutral to patentability. If no result is
materially relevant to a claim, say so explicitly — don't just list papers.
