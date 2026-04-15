# Patent Review Workflow — Design Spec

**Date:** 2026-04-14
**Status:** Approved

---

## Overview

A Claude Code skill (`/patent-review`) that guides a scientist through a structured, interactive review of a lawyer-drafted patent. The scientist validates technical accuracy, identifies coverage gaps, and proposes improved or new claim language. A final structured report is saved at the end of each session.

**Primary user:** A scientist (domain expert) reviewing a patent drafted by a patent attorney.
**Domain:** Domain-agnostic — the skill makes no assumptions about the field of science.
**Input formats:** Plain text pasted into the conversation, or a `.docx` file path.
**Output:** Interactive review session + a saved Markdown report.

---

## Workflow Phases

The skill progresses through five phases. Phases 1–4 run sequentially and automatically. Phase 5 is optional and user-triggered at any point during Phases 3 or 4.

### Phase 1 — Ingest

The skill prompts the user to either paste the patent text or provide a `.docx` file path. For `.docx` files, the skill extracts text using available tooling. The document is then parsed into named sections:

- Title
- Abstract
- Background
- Summary of the Invention
- Brief Description of Drawings (if present)
- Detailed Description / Embodiments
- Claims (independent and dependent, enumerated)
- Examples / Sequence Listings (if present)

The skill identifies and records the total number of claims, how many are independent, and how many are dependent.

### Phase 2 — Triage

Before deep review, the skill produces a brief orientation snapshot:

- **Invention summary** (2–3 sentences) — what is being claimed at a high level
- **Inferred technology domain**
- **Structural checklist** — which standard sections are present or missing
- **Initial flags** — anything immediately anomalous (e.g., claims that contradict the description, missing definitions, obvious internal inconsistencies)

The scientist confirms or corrects the triage before the skill proceeds to Phase 3.

### Phase 3 — Deep Review (Interactive)

The core phase. The skill works through the document section by section in this order:

1. **Background & Summary**
   - Is the problem statement technically accurate?
   - Does the summary fairly and completely represent the invention?

2. **Detailed Description & Embodiments**
   - Are the described mechanisms scientifically sound?
   - Are any embodiments vague, ambiguous, or technically incorrect?
   - Are there known embodiments or variations not described?

3. **Claims** (three sequential passes)
   - **Accuracy pass** — each claim is checked for technical correctness: wrong terminology, impossible limitations, overly broad scientific assertions
   - **Gap pass** — what aspects of the invention are unclaimed? What variations, methods, compositions, materials, or uses are left unprotected?
   - **Language pass** — for each flagged claim, the skill proposes specific revised claim language or drafts new claims to address identified gaps

4. **Abstract** (reviewed last)
   - Once the claims are finalized, the abstract is checked to ensure it accurately reflects the agreed-upon scope

**Interaction model:** At each step, the skill surfaces one issue at a time, explains the concern, and asks the scientist to: (a) accept the suggested fix, (b) modify it, or (c) dismiss it. All decisions — accepted, modified, and dismissed — are logged for the final report.

### Phase 4 — Synthesis

The skill compiles all logged findings into a structured Markdown report saved to the working directory as:

```
YYYY-MM-DD-<patent-title>-review.md
```

**Report structure:**

| Section | Content |
|---|---|
| Executive Summary | Invention overview, overall assessment, count of critical vs. minor issues |
| Section-by-Section Findings | Grouped by patent section; each finding tagged `[ACCURACY]`, `[GAP]`, or `[LANGUAGE]` |
| Claims Analysis | Table of all claims: status (accepted / revised / flagged), proposed new claims |
| Recommendations | Prioritized action list for the attorney |
| Literature | Populated only if Phase 5 was triggered |

### Phase 5 — Literature Review (Optional, On-Demand)

**Trigger:** The scientist requests a literature search at any time during Phase 3 or 4 using natural language (e.g., "check the literature on X", "find prior art on Y", "are there papers about Z?").

**On trigger, the skill asks two quick questions:**
1. Scope — scientific literature, prior art patents, or both?
2. Any constraints — date range, jurisdiction, specific databases?

The skill then runs the search, summarizes relevant hits, and flags any results that materially affect the claims under review. Findings are appended to the Literature section of the final report.

---

## Design Decisions

| Decision | Rationale |
|---|---|
| Domain-agnostic | The scientist brings domain expertise; the workflow should not constrain the field |
| Literature review is opt-in | Adding it by default would slow every session; most reviews don't require it |
| Claims reviewed in three passes | Separating accuracy, gaps, and language keeps each pass focused and prevents issues from being conflated |
| Abstract reviewed last | The abstract should reflect finalized claims, not drive them |
| Issues surfaced one at a time | Prevents overwhelm; keeps the scientist in a clear decision-making posture |
| All decisions logged | Dismissed issues are as important as accepted ones — the attorney needs the full picture |
| Report saved as Markdown | Portable, version-controllable, easy to share |

---

## Out of Scope

- Automated filing or submission to patent offices
- Legal advice or legal review (the skill is a scientific review tool)
- Claim novelty determination (that requires formal prior art search, which is Phase 5 if requested)
- Formatting to USPTO or EPO template standards

---

## Success Criteria

- The scientist can complete a full patent review in a single session
- Every claim in the draft is explicitly assessed (accepted, revised, or flagged)
- The final report is complete enough for the attorney to act on without follow-up questions
- Literature review, when triggered, integrates cleanly into the session without disrupting the review flow
