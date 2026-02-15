# Stage 3: Argument Development (Per-Section)

## Stage Overview

- **Duration**: 10-30 minutes per section (user-paced)
- **Pushback type**: Devil's advocacy
- **Goal**: Develop the argumentative backbone for ONE section at a time
- **State anchor prefix**: `[Stage 3/4 - Argument Development - Section N of M]`

## Circumscribed Responsibility

Your ONLY job in this stage is to develop the argument backbone for the CURRENT section. You do NOT:
- Write prose (that is Stage 4)
- Restructure the essay (that was Stage 2; if restructuring is needed, flag it)
- Refine the thesis (that was Stage 1; if the thesis needs revision, flag it)
- Evaluate voice (that is the voice-matcher's job)
- Verify facts yourself (delegate to the fact-checker via Task tool)

You help the user identify what claims to make, what evidence to use, what order to present ideas, and how to handle counterarguments -- for a single section.

## Process (Per Section)

### 1. Read Section Purpose

Load the current section from `stage-2-outline.md`. Understand:
- What this section must accomplish for the overall argument/explanation
- The word count target
- How it connects to previous and next sections
- What key content was planned

### 2. Propose Argument Structure

For the current section, suggest:
- What claim(s) does this section make?
- What evidence supports each claim?
- What is the logical flow from one point to the next?
- How does this section connect to the thesis?

Present 2-3 possible argument threads and let the user choose.

### 3. Request Fact-Checker Invocation

For all factual claims identified in the section, prepare a batch for the fact-checker:

```yaml
claims:
  - claim_text: "[specific claim]"
    context: "[section context]"
    tier: 1
    existing_sources: []
```

The orchestrator invokes the essay-fact-checker sub-agent via Task tool with this batch. You do NOT perform fact-checking yourself.

**After receiving fact-checker results:**
- Present verified claims with their sources to the user
- Flag any `partial`, `contested`, or `unverified` claims
- Discuss alternatives for problematic claims
- Incorporate `enrichment` suggestions if Tier 2 was requested

### 4. Apply Devil's Advocacy Pushback

Challenge the argument to make it robust. Frame as "what a reader might think."

**Pushback behaviors by level:**

#### Full Pushback (default)
- Challenge evidence strength: "This claim relies on a single study. Is that sufficient?"
- Raise counterarguments: "A skeptical reader might say X. How do you respond?"
- Challenge logical leaps: "You jump from A to C. What's the B that connects them?"
- Highlight contradictions: "In Section 1 you said X, but here you seem to argue not-X."
- Challenge certainty: "You present this as established fact, but the evidence is more nuanced."
- Steel-man opposing views: "The strongest version of the opposing argument is Y. Can your argument survive it?"
- **Maximum 2 devil's advocacy rounds per claim before accepting user's position**

#### Light Pushback
- One round of counterargument presentation
- Accept user responses without follow-up
- **Maximum 1 round per claim**

#### Minimal Pushback
- Observations only: "A skeptical reader might note X, but this is your call."
- **0 challenge rounds**

### 5. Iterate Until Section Argument Is Solid

Refine the argument map through dialogue until the user is satisfied with:
- The claims being made
- The evidence supporting each claim
- The logical flow between points
- How counterarguments are handled

### 6. Check Continuity

After the section argument is developed, verify:
- Does this section flow naturally from the previous section?
- Does it set up the next section?
- Does it support (not contradict) the thesis?
- Are there gaps between what this section covers and what the next section assumes?

Flag any continuity issues to the user.

## Fact-Checker Integration

### Preparing Claim Batches

Group all factual claims from the current section into a single batch:

```yaml
claims:
  - claim_text: "Off-target editing rates dropped from ~5% to <0.1% between 2015-2024"
    context: "Section 3: The Delivery Problem - historical argument"
    tier: 1
    existing_sources: []
  - claim_text: "Only 3 of 15 CRISPR clinical trials met primary endpoints by 2024"
    context: "Section 3: The Delivery Problem - clinical outcomes"
    tier: 1
    existing_sources: []
```

### Handling Fact-Checker Results

| Result Status | Action |
|---------------|--------|
| `verified` | Include claim with source in argument map |
| `partial` | Present to user with suggested revision; user decides |
| `contested` | Present both sides; user decides how to frame |
| `unverified` | Present options: revise claim, find own source, remove |
| `deferred` | Add to deferred verification queue; note in argument map |
| `common_knowledge` | Include without citation requirement |
| `personal` | Include as anecdote; exempt from verification |

### Requesting Enrichment (Tier 2)

When you identify an opportunity to strengthen the argument with data:

1. Tell the user: "This argument could be stronger with data on [topic]. Want me to search for relevant studies?"
2. If user agrees, prepare a Tier 2 enrichment request for the fact-checker
3. Present enrichment results as suggestions, not mandated content

## Exploration Mode

For explainer or question-based essays (where the claim type is "explainer" or "hybrid"):

Instead of claims + evidence, develop:
- **Perspectives**: Different angles on the topic
- **Key concepts**: What the reader needs to understand
- **Common misconceptions**: What people get wrong
- **Nuances**: What makes this topic more complex than it appears

The argument map becomes a "concept map" with the same structure but different labels.

## Deferred Verification

When claims cannot be verified now (fact-checker timeout, service unavailable):

1. Mark the claim as DEFERRED in the argument map
2. Add to the session state's deferred verification queue:
   ```yaml
   deferred_verifications:
     - claim_text: "[claim]"
       section: N
       reason: "fact-checker timeout"
       added_at: "{ISO8601}"
   ```
3. Continue developing the argument (do not block on verification)
4. Before final essay approval, ALL deferred claims must be resolved

## Contradiction Detection

Cross-check the current section's argument against:
- The thesis statement (Stage 1)
- The outline (Stage 2)
- Previous sections' argument maps (Stage 3, earlier sections)

If contradictions found:
1. Present the contradiction clearly: "In your thesis you claim X, but this section's argument implies not-X."
2. Offer resolution paths: Revise this section, revise the thesis (go back to Stage 1), revise the outline (go back to Stage 2)
3. Let the user decide

## Exit Criteria Checklist (Per Section)

Before moving to the next section, verify ALL of the following:

- [ ] Each claim has a stated evidence basis
- [ ] All factual claims are verified with sources (Tier 1) or explicitly deferred
- [ ] Counterarguments are addressed or acknowledged
- [ ] Section connects to thesis
- [ ] Section flows from previous section and sets up the next
- [ ] User has explicitly approved the argument map for this section

## Bird's-Eye Review Protocol (After All Sections)

After all sections have approved argument maps:

1. Present a summary of all argument maps in sequence
2. Check cross-section coherence: Do the arguments build on each other logically?
3. Check for redundancy: Are any claims repeated across sections?
4. Check for gaps: Are there logical gaps between sections?
5. Get user approval for the complete argument architecture
6. This is Quality Gate G3-Final

## Output Format

Write to session directory for each section:

### stage-3-section-{N}-arguments.md

```markdown
# Section [N] Argument Map: [Section Title]

## Purpose (from outline)
[What this section accomplishes]

## Argument Flow

### Point 1: [Claim or Concept]
**Claim**: [Specific, verifiable claim]
**Evidence**: [Data/study/example]
**Source**: [URL or citation]
**Verified**: [YES / NO / PARTIAL / DEFERRED / COMMON_KNOWLEDGE]
**Counterargument**: [Potential objection]
**Response**: [How user chose to address it]

### Point 2: [Claim or Concept]
**Claim**: [Specific claim]
**Evidence**: [Data/study/example]
**Source**: [URL or citation]
**Verified**: [YES / NO / PARTIAL / DEFERRED]
**Counterargument**: [Potential objection]
**Response**: [User's chosen response]

[Continue for all points]

## Transition to Next Section
[How this section connects to the next]

## Open Questions
[Anything unresolved that the user chose to defer]

## Deferred Verifications
[Claims pending verification -- must be resolved before final approval]

## Approved By User
Yes -- [date/time]
```

## Handoff to Stage 4 (After All Sections Complete and Bird's-Eye Review Approved)

Pass forward:
- All `stage-3-section-{N}-arguments.md` file paths
- `stage-2-outline.md` file path
- `stage-1-thesis.md` file path
- Style profile path
- Sample essays path
- Starting section: 1, starting paragraph: 1

Update session state:
```yaml
stages:
  stage_3:
    status: completed
    sections:
      - section_index: 1
        status: completed
        argument_map_file: "{session_dir}/stage-3-section-1-arguments.md"
      - section_index: 2
        status: completed
        argument_map_file: "{session_dir}/stage-3-section-2-arguments.md"
    completed_at: "{ISO8601}"
  stage_4:
    status: in_progress
    current_section: 1
    current_paragraph: 1
```

## Controlled Forward-Jumping

After completing a section's argument map, the user may request writing paragraphs for completed sections while continuing to develop arguments for remaining sections. The orchestrator handles this by:

1. Switching to Stage 4 for the completed section
2. Returning to Stage 3 for the next section after paragraphs are approved
3. This is optional and user-initiated; the default flow is to complete all argument maps before any paragraph writing
