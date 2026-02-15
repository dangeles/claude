# Stage 4: Paragraph Writing (Per-Paragraph)

## Stage Overview

- **Duration**: 5-15 minutes per paragraph (user-paced)
- **Pushback type**: Minimal (execution quality only)
- **Goal**: Draft ONE paragraph at a time in the user's voice, following the argument map
- **State anchor prefix**: `[Stage 4/4 - Paragraph Writing - Section N, Paragraph M]`

## Circumscribed Responsibility

Your ONLY job in this stage is to draft paragraphs. You do NOT:
- Restructure arguments (that was Stage 3; if restructuring is needed, flag it and offer to go back)
- Restructure the essay (that was Stage 2)
- Refine the thesis (that was Stage 1)
- Verify facts yourself (delegate to the fact-checker for NEW claims only)
- Evaluate voice yourself (delegate to the voice-matcher)

You propose text. The user approves, edits, or rejects it. The user's edits are always final. You never argue with the user's prose preferences.

## Process (Per Paragraph)

### 1. Identify What This Paragraph Covers

From the argument map for the current section, determine which point(s) this paragraph handles. Each paragraph typically covers 1-2 argument map points.

State clearly: "This paragraph covers Point [N] from the argument map: [brief description]."

### 2. Draft Paragraph in User's Voice

Using the style profile as the primary voice reference:

- Match sentence length patterns from the profile
- Use the user's vocabulary preferences (include characteristic words, avoid forbidden words)
- Follow the user's rhetorical patterns (analogies, evidence presentation, hedging)
- Embed factual claims with their verified sources
- Match the user's paragraph length range

### 3. Self-Audit: Argument Map Compliance

Before presenting the draft, check:

- **MAPPED content**: Tag each sentence or claim as covering a specific argument map point
- **NOT IN MAP content**: Flag any content that was NOT in the argument map

If NOT IN MAP content exists, explicitly notify the user:
"Note: The sentence '[sentence]' introduces content not in your argument map for this section. Options: (a) keep it and add to map retroactively, (b) remove it, (c) go back to Stage 3 to revise the argument map."

### 4. Present Draft with Annotations

Present the paragraph followed by:

```markdown
---

**Argument map coverage**: Points [N, M] from Section [X]
**Sources used**:
- [Claim]: [Source title] ([URL])
**Voice notes**: [Observations about how the draft matches the style profile]
**Flags**: [Any NOT IN MAP content, uncertain voice matches, or other issues]

---
```

### 5. User Reviews

Accept one of four responses:

| Response | Action |
|----------|--------|
| **Approve** | Save paragraph; move to next |
| **Edit** | User provides modified text; accept modifications as-is without debate |
| **Reject** | User provides direction; rewrite from scratch following guidance |
| **Rewrite with feedback** | User gives specific notes ("make it more conversational", "cut the analogy", "lead with the data"); revise accordingly |

**Never argue with user edits.** If the user changes factual content, note if the new version needs fact-checking, but do not resist the edit.

## Voice-Matcher Integration

The orchestrator invokes the voice-matcher:

- **After each section is complete**: Section-level voice consistency check
- **After the full essay is complete**: Essay-level voice consistency check
- **On user request at any time**: User asks "does this sound like me?"

### Handling Voice-Matcher Results

| Score | Action |
|-------|--------|
| 5/5 | Proceed; no changes needed |
| 4/5 | Proceed with notes; mention specific observations to user |
| 3/5 | Present assessment to user; ask if they want revisions |
| 2/5 | Flag for revision; present specific suggestions |
| 1/5 | Pause; suggest reviewing the style profile for clarity |

### Calibration

On the first voice-matcher invocation, present the assessment to the user and ask: "How accurate is this voice assessment? (1-5)" Save their rating in session state for subsequent invocations.

## Fact-Checker Integration (Stage 4)

In Stage 4, invoke the fact-checker ONLY for NEW claims -- claims introduced during paragraph writing that were NOT already verified in Stage 3.

Do NOT re-verify claims that were verified during Stage 3 (they are already in the argument maps with sources).

**Identifying new claims:**
- Compare each factual statement in the draft against the argument map
- If a claim appears in the argument map with `Verified: YES`, it does not need re-verification
- If a claim is new (not in the argument map), prepare it for fact-checking

## Paragraph Presentation Format

```markdown
**Section [N], Paragraph [M] -- Draft**

[Draft paragraph text here]

---

**Covers**: Argument map points [list]
**Sources**:
- "[claim]" -- [Source] ([URL])
- "[claim]" -- [Source] ([URL])
**Voice**: [Brief note on voice match, e.g., "Matched your pattern of concrete-number-then-implication. Closing sentence uses your characteristic short-punchy style."]
**Flags**: [None / "Sentence 3 introduces content not in argument map" / etc.]

Approve / Edit / Reject / Rewrite with feedback?
```

## Section Completion Protocol

After all paragraphs in a section are approved:

1. **Present full section**: Show all approved paragraphs assembled in order
2. **Bird's-eye review**: Ask user to review the section as a whole
   - "Does the section read well as a unit?"
   - "Are there any transitions that feel rough?"
   - "Does the section accomplish its stated purpose from the outline?"
3. **Voice check**: Orchestrator invokes voice-matcher on the complete section
4. **User options**: Approve section, request specific paragraph revisions, reorder paragraphs
5. **Save**: Write completed section to `stage-4-draft.md`

## Essay Completion Protocol

After all sections are complete:

1. **Present full essay**: Show all sections assembled in order
2. **Final voice check**: Orchestrator invokes voice-matcher on the entire essay
3. **Final fact-check sweep**: Orchestrator invokes fact-checker on all claims in the complete essay
   - This catches any claims that slipped through or were modified during user edits
4. **Resolve deferred verifications**: Any claims in the deferred queue must be resolved now
5. **Present completion summary**:
   ```markdown
   ## Pipeline Completion Summary

   - **Word count**: [N] words
   - **Sections**: [N]
   - **Sources verified**: [N] claims with [M] unique sources
   - **Voice consistency**: [score]/5
   - **Deferred claims resolved**: [N] of [N]
   - **User overrides**: [N] (logged in fact-check-log.md)
   ```
6. **User final approval**: This is Quality Gate G5
7. **Write final essay**: Save to `{session_dir}/final-essay.md`
8. **Offer export**: Ask if user wants the essay written to a specific path

## Exit Criteria

### Per Paragraph
- [ ] Paragraph covers intended argument map points
- [ ] All factual claims have verified sources (or are explicitly deferred)
- [ ] User has explicitly approved the paragraph text

### Per Section
- [ ] All paragraphs approved
- [ ] Section reads well as a unit (user confirmed)
- [ ] Voice consistency checked (voice-matcher invoked)
- [ ] Section accomplishes its stated purpose from the outline

### Final (Full Essay)
- [ ] All sections approved
- [ ] Full essay voice check completed
- [ ] Full essay fact-check sweep completed
- [ ] All deferred verifications resolved
- [ ] User has given final approval (Quality Gate G5)

## Output Format

### stage-4-draft.md (Incremental)

Append each approved paragraph to this file as it is completed. After each section, add a section separator.

```markdown
# [Essay Working Title]

## [Section 1 Title]

[Approved paragraph 1]

[Approved paragraph 2]

[Approved paragraph 3]

---

## [Section 2 Title]

[Approved paragraph 1]

...
```

### stage-4-paragraph-log.md

Log each paragraph's journey:

```markdown
# Paragraph Writing Log

## Section 1, Paragraph 1
- **Draft presented**: [timestamp]
- **User action**: Approved with edit
- **Edit summary**: Changed "fifty-fold" to "roughly fifty-fold"
- **Argument map points covered**: 1.1, 1.2
- **Sources**: [list]
- **Voice notes**: [assessment summary]

## Section 1, Paragraph 2
...
```

### final-essay.md

```markdown
# [Essay Title]

[Complete essay text with all approved paragraphs]

---

## Sources

1. [Source title]. [URL]. Accessed [date].
2. [Source title]. [URL]. Accessed [date].
...

---

## Metadata

- **Pipeline session**: [session_id]
- **Thesis**: [thesis statement]
- **Word count**: [N]
- **Voice consistency score**: [N]/5
- **Verified claims**: [N]
- **User overrides**: [N]
- **Completed**: [date/time]
```
