# Patent Review Workflow — Design Spec

**Date:** 2026-04-14
**Status:** Approved

---

## Overview

A Claude Code skill (`/patent-review`) that guides a scientist through a structured, interactive review of a lawyer-drafted patent. The scientist validates technical accuracy, identifies coverage gaps, and proposes improved or new claim language. A final structured report is saved at the end of each session.

**Primary user:** A scientist (domain expert) reviewing a patent drafted by a patent attorney.
**Domain:** Domain-agnostic — the skill makes no assumptions about the field of science.
**Input formats:** Plain text pasted into the conversation, or a `.docx` file path. PDF input is out of scope (patent PDFs vary widely in extraction quality; users should convert to `.docx` or paste text).
**Output:** Interactive review session + a saved Markdown report.

---

## Workflow Phases

The skill progresses through five phases. Phases 1–4 run sequentially and automatically. Phase 5 is optional and user-triggered at any point during Phase 3, or by typing `literature` at the Phase 4 entry prompt — before the skill begins compiling the report (see Phase 4). After the entry prompt is answered with anything other than `literature`, Phase 5 is no longer available for that session.

## Invocation

The skill is invoked as `/patent-review` with no required arguments. Optionally, the user may provide a `.docx` file path as the first argument:

```
/patent-review
/patent-review path/to/draft.docx
```

If a file path is provided on invocation, Phase 1 skips the input prompt and proceeds directly to extraction. If no argument is given, Phase 1 prompts the user to paste text or enter a file path interactively.

### Phase 1 — Ingest

The skill prompts the user to either paste the patent text or provide a `.docx` file path. For `.docx` files, text is extracted using `markitdown` (preferred) or `python-docx` as a fallback. If neither is available, the skill instructs the user to paste the text manually. The document is then parsed into named sections:

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

**Interaction model:** This pattern applies uniformly across all Phase 3 sub-sections (Background, Detailed Description, and Claims). At each step, the skill surfaces one issue at a time, explains the concern, and asks the scientist to: (a) accept the suggested fix, (b) modify it, or (c) dismiss it. All decisions — accepted, modified, and dismissed — are logged for the final report.

**Critical vs. minor classification:** When logging a finding, the skill assigns a severity:
- **Critical** — the finding would affect claim validity, patentability, or scientific accuracy in a way the attorney must address (e.g., a technically impossible claim, a missing claim for a core embodiment)
- **Minor** — the finding is a stylistic, terminology, or clarity improvement that is recommended but not legally or scientifically essential

The severity is stated explicitly when the skill surfaces each issue ("This is a critical issue because…" / "This is a minor issue…"). The Executive Summary counts both categories separately.

### Phase 4 — Synthesis

**Phase 4 entry prompt:** Before compiling, the skill presents a single prompt:

> "Ready to generate the report. Where should I save it? [default: ./] — or type 'literature' to run a literature search first."

This is the last opportunity to trigger Phase 5. After the user responds, Phase 5 is no longer available for this session.

The skill compiles all logged findings into a structured Markdown report. Default save location is the current working directory. If the user specifies a custom path, that path is used. If the path is not writable, the skill alerts the user and asks once for an alternative. If the second path is also unwritable, the skill prints the full report to the conversation as a fenced Markdown block so the user can copy it manually.

Report filename:
```
YYYY-MM-DD-<patent-title>-review.md
```

The `<patent-title>` slug is derived from the document title by: lowercasing, replacing spaces and special characters with hyphens, collapsing consecutive hyphens, and truncating to 60 characters. Example: `"Compositions and Methods for Treatment of Inflammatory Disease"` → `compositions-and-methods-for-treatment-of-inflammatory-disea`.

**Report structure:**

| Section | Content |
|---|---|
| Executive Summary | Invention overview, overall assessment, count of critical vs. minor issues |
| Section-by-Section Findings | Grouped by patent section; each finding tagged `[ACCURACY]`, `[GAP]`, or `[LANGUAGE]` |
| Claims Analysis | Table of all claims: status (accepted / revised / flagged), proposed new claims |
| Recommendations | Prioritized action list for the attorney |
| Literature | Populated only if Phase 5 was triggered |

### Phase 5 — Literature Review (Optional, On-Demand)

**Trigger:** The scientist requests a literature search at any time during Phase 3, or before the report is finalized at the start of Phase 4 (i.e., before the skill asks about the save path). Once the report is being written, Phase 5 cannot be triggered in that session. If triggered during Phase 4, the literature results are incorporated before the report is saved.

Trigger examples: "check the literature on X", "find prior art on Y", "are there papers about Z?"

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

## Session Interruption

Patent reviews can be long. If the session is interrupted (terminal closed, context reset), partial work is not automatically recoverable — Claude Code does not persist conversation state across sessions.

**Defined behavior:**
- The skill does not attempt to resume a prior session
- If a session is interrupted mid-review, the user must restart by reinvoking `/patent-review` and re-providing the document
- To mitigate loss, the skill offers a **checkpoint summary** at the end of each Phase 3 sub-section. The checkpoint is a fenced Markdown block, delimited with `--- CHECKPOINT START ---` and `--- CHECKPOINT END ---`, listing all logged findings to that point. The skill prompts: "Checkpoint above. Copy this block if you want a record in case the session is interrupted. Type `continue` to proceed." If the user does not respond (e.g., closes the terminal), the session simply ends — no further action is taken and the in-progress findings are lost. This is the defined and acceptable behavior; no background saving or retry is attempted
- The final report is only written at Phase 4; there is no auto-save of partial findings

This is acceptable given the typical session length (one patent per sitting) and the checkpoint mitigation.

---

## Out of Scope

- Automated filing or submission to patent offices
- Legal advice or legal review (the skill is a scientific review tool)
- Claim novelty determination (that requires formal prior art search, which is Phase 5 if requested)
- Formatting to USPTO or EPO template standards

---

## Success Criteria

- Every claim in the draft receives an explicit status: accepted, revised, or flagged — no claim is skipped
- Every tagged finding in the report (`[ACCURACY]`, `[GAP]`, `[LANGUAGE]`) maps to a specific section and claim number, with a proposed resolution
- The report contains a prioritized Recommendations section with at least one actionable item per critical finding
- When Phase 5 is triggered, the Literature section of the report is populated with results, each annotated with whether it materially affects a specific claim. If Phase 5 is never triggered, the Literature section is omitted from the report entirely
- The session completes and the report is saved without requiring the user to manually collate findings
