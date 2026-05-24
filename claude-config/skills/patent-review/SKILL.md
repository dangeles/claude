---
name: patent-review
last_updated: 2026-05-24
description: >
  This skill should be used when the user asks to "review this patent", "check the claims",
  "patent draft review", "is this claim accurate", "what's missing from this patent", or
  when a user pastes a patent draft or provides a .docx file path with patent content.
  Use this skill whenever a scientist needs to review a patent draft written by a lawyer —
  to check technical accuracy, identify gaps in claim coverage, or propose revised or new
  claim language. This is a multi-phase interactive workflow: ingest → triage → section-by-section
  review → structured report. Do not use this skill for legal advice, filing assistance, or
  non-scientific patent tasks.
---

# Patent Review

Guide a scientist through structured review of a lawyer-drafted patent. The goal is to catch
technical errors, identify coverage gaps, and produce better claims — not to give legal advice.
The scientist is the domain expert; this workflow provides the structure.

Work through five phases in order. Surface one issue at a time, get a decision, log it, move on.

---

## Phase 1 — Ingest

Ask the user to paste the patent text or provide a `.docx` file path.

- **Text paste**: parse directly
- **`.docx` file**: extract text using `markitdown` (preferred) or `python-docx` as fallback;
  if neither is available, ask the user to paste the text manually
- **PDF**: not supported — ask the user to convert to `.docx` or paste the text

Once text is in hand, identify and label these sections (use only what's present — not all sections
appear in every patent):

- Title
- Abstract
- Background
- Summary of the Invention
- Brief Description of Drawings
- Detailed Description / Embodiments
- Claims (enumerate all; note which are independent vs. dependent)
- Examples / Sequence Listings

Record the total claim count (independent + dependent separately).

---

## Phase 2 — Triage

Give the scientist a quick orientation before diving in:

1. **Invention summary** — 2–3 sentences: what is this patent about?
2. **Inferred domain** — e.g., "small-molecule therapeutics", "gene editing", "semiconductor fabrication"
3. **Structural checklist** — which standard sections are present, which are missing
4. **Initial flags** — anything immediately anomalous: claims that contradict the description,
   undefined terms, missing definitions, obvious internal inconsistencies

Ask the scientist to confirm or correct before proceeding. A wrong domain inference can skew the
entire review.

---

## Phase 3 — Deep Review (Interactive)

Work through the document section by section. For each issue:

1. Surface the issue, explain why it's a concern
2. Label it **[ACCURACY]**, **[GAP]**, or **[LANGUAGE]**; assign severity:
   - **Critical**: affects claim validity, patentability, or scientific accuracy — attorney must address
   - **Minor**: stylistic, terminology, or clarity improvement — recommended but not essential
3. Wait for the scientist to choose:
   > **(a) accept** — use the suggested fix as-is
   > **(b) modify** — scientist provides the corrected version
   > **(c) dismiss** — skip this issue (log as dismissed with reason if given)

Log every decision, including dismissals. Be explicit: "This is a critical issue because…" /
"This is a minor issue — worth fixing but not legally significant."

### 3.1 Background & Summary
- Is the problem statement technically accurate?
- Does the summary fairly and completely represent what the invention actually does?

### 3.2 Detailed Description & Embodiments
- Are the described mechanisms scientifically sound?
- Are any embodiments vague, ambiguous, or technically incorrect?
- Are there known embodiments or variations that exist but aren't described?

### 3.3 Claims — three sequential passes

Run in order. Don't conflate them — each has a distinct goal.

**Accuracy pass**: Go through each claim. Flag technical errors: wrong terminology, impossible
limitations, overly broad scientific assertions unsupported by the description.

**Gap pass**: Ask: what aspects of this invention are NOT claimed? Think about alternative
embodiments, different methods of making or using it, different materials, compositions, use cases.
Prompt the scientist with questions if the space of variations is unclear.

**Language pass**: For each claim flagged in the accuracy or gap passes, propose specific revised
claim language or draft new claims. Be concrete — provide actual claim text, not just "fix the wording."

### 3.4 Abstract (reviewed last)
Once the claims are settled, verify the abstract accurately reflects the agreed-upon scope. A
patent's abstract should match its claims, not the other way around.

### Checkpoints

At the end of each sub-section (3.1, 3.2, 3.3 accuracy, 3.3 gap, 3.3 language, 3.4), output a
checkpoint block:

```
--- CHECKPOINT START ---
[All logged findings so far, with labels, severities, and accepted/modified/dismissed status]
--- CHECKPOINT END ---
```

Then prompt: "Checkpoint above — copy this block if you want a record. Type `continue` to proceed."

If the user doesn't respond (session ends), the session is simply over. There is no auto-save. The
final report is only written in Phase 4.

---

## Phase 4 — Synthesis

Before compiling, present this prompt:

> "Ready to generate the report. Where should I save it? [default: `./`] — or type `literature` to
> run a literature search first."

After the user answers, Phase 5 is no longer available. If they type `literature`, run Phase 5 first,
then compile.

Compile the report using the structure in `references/report-template.md`.

**Filename**: `YYYY-MM-DD-<slug>-review.md`
- Slug: lowercase patent title, spaces/specials → hyphens, collapse consecutive hyphens, truncate to 60 characters
- Example: `"Compositions and Methods for Treatment of Inflammatory Disease"` →
  `2026-04-14-compositions-and-methods-for-treatment-of-inflammato-review.md`

**Save path behavior**: Use the user's specified path, or default `./`. If the path is unwritable,
ask once for an alternative. If that also fails, print the full report as a fenced Markdown block
so the user can copy it.

---

## Phase 5 — Literature Review (Optional, On-Demand)

**When to trigger**: During Phase 3 (user can say "check the literature on X" at any point), or at
the Phase 4 entry prompt (user types `literature`). Once the Phase 4 prompt is answered with
anything else, Phase 5 is no longer available.

**On trigger**, ask two quick questions:
1. Scope — scientific literature, prior art patents, or both?
2. Any constraints — date range, jurisdiction, specific databases?

Run the search, summarize relevant hits, and for each result note whether it materially affects a
specific claim under review. Append findings to the Literature section of the report.

---

## Scope Boundaries

This skill does NOT:
- Give legal advice or opinions on patentability
- Perform a formal freedom-to-operate analysis
- File anything or format for USPTO/EPO submission
- Make the final call on claim language — the scientist decides, you propose

---

## Additional Resources

- **`references/report-template.md`** — exact Markdown template for the Phase 4 report, including
  the claims analysis table format and section structure
- **`evals/evals.json`** — benchmark test cases (DMD gene therapy, KRAS inhibitor, SCD base editing)
  for evaluating skill behavior across domains
