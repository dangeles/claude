# Stage Specifications

Reference document for detailed 8-stage pipeline process specifications.

---

## Stage 1: Scope Refinement

**Owner**: requirements-analyst (existing skill)
**Checkpoint**: ALWAYS (required)
**Duration**: 15-30 minutes
**Timeout**: 45 minutes

### Input
- User's initial prompt

### Process
1. Clarify research question, target audience, intended use
2. Define success criteria (decision support? landscape map? comprehensive review?)
3. Set boundaries (in-scope topics, out-of-scope tangents)
4. **Complexity detection** -> Propose checkpoint plan
5. User approves scope + checkpoint plan

### Output
```yaml
scope_document:
  research_question: string  # Specific, testable
  success_criteria: list     # Measurable outcomes
  in_scope: list            # Topics to cover
  out_of_scope: list        # Explicit exclusions
  complexity_tier: enum     # Simple | Medium | Complex | High-Stakes
  checkpoint_plan: object   # Which stages have checkpoints
```

### Quality Gate Checklist
- [ ] Research question is specific (not "what is X?" but "what governs X?")
- [ ] Success criteria are measurable
- [ ] In-scope and out-of-scope boundaries clear
- [ ] User approves scope and checkpoint plan

### Failure Handling
- User rejects scope: requirements-analyst revises based on feedback
- Timeout: Escalate to user with partial scope for approval

---

## Stage 2: Parallel Review Discovery

**Owner**: lit-pm orchestrates 2-3 literature-researcher agents in parallel
**Checkpoint**: Only if High-Stakes complexity
**Duration**: 45-90 minutes (parallel execution)
**Timeout**: 120 minutes total

### Input
- Refined scope document

### Search Strategies

lit-pm generates 2-3 diverse search strategies:
- Strategy A: Broad keywords (e.g., "CAR-T manufacturing review")
- Strategy B: Specific technical terms (e.g., "transduction efficiency optimization")
- Strategy C: Application focus (e.g., "clinical CAR-T production challenges")

### Process
1. Launch 2-3 literature-researcher agents in parallel (using Task tool)
2. Each agent finds 2-3 high-quality reviews using their assigned strategy
3. Each writes brief annotation per review (key findings, coverage, quality)
4. lit-pm collects results and analyzes convergence

### Convergence Tracking Algorithm

1. Collect reviews from all agents with source attribution
2. Match duplicates using (in order):
   - DOI exact match (highest confidence)
   - PubMed ID exact match
   - Title fuzzy match (Levenshtein distance < 0.2)
   - First author + year + journal match
3. Score by convergence:
   - Found by 3/3 agents: HIGH priority (must-read)
   - Found by 2/3 agents: MEDIUM priority
   - Found by 1/3 agents: Unique perspective (coverage breadth)

### Output
```yaml
review_collection:
  reviews:
    - title: string
      doi: string | null
      authors: string
      year: integer
      source_agent: string
      convergence_count: integer
      annotation: string
  convergence_analysis:
    high_priority: list
    medium_priority: list
    unique_perspectives: list
    themes_covered: list
    gaps_identified: list
```

### Quality Gate Checklist
- [ ] 6-9 reviews collected total (2-3 per agent)
- [ ] At least 2 reviews show convergence (found by multiple agents)
- [ ] Reviews published in last 10 years (unless older is definitive)
- [ ] Coverage of all major themes in scope

### Failure Handling
- If convergence < 2 reviews: Adjust search strategies, re-run
- If total reviews < 4: Broaden search terms, extend to 15 years
- If quality low: Prioritize high-citation reviews
- Timeout: Proceed with available reviews, flag incomplete coverage

---

## Stage 3: Layered Outline Synthesis

**Owner**: lit-pm (high-level structure) + literature-researcher agents (section details)
**Checkpoint**: Medium/High complexity
**Duration**: 30-60 minutes
**Timeout**: 90 minutes

### Input
- Review collection with convergence data

### Process
1. **lit-pm reads all reviews, creates high-level structure**:
   - 3-5 sections (not too broad, not too narrow)
   - Each section has: Title, Core thesis (testable claim or question), Key question to address

2. **For each section: Assign detail-work to literature-researcher agent**:
   - Agent proposes 2-4 subsections with specific questions
   - Agent identifies which reviews/papers to cite per subsection
   - Agent notes potential gaps or areas needing targeted research

3. **lit-pm presents complete outline to user** (if checkpoint enabled)

### Output
```yaml
outline:
  introduction:
    thesis: string
    assigned_to: lit-synthesizer
  sections:
    - title: string
      thesis: string
      subsections: list
      assigned_to: string  # literature-researcher agent
      key_reviews: list
      research_mandate: string
  balance_check:
    largest_section_pct: float
    smallest_section_pct: float
```

### Quality Gate Checklist
- [ ] 3-5 sections (balanced scope)
- [ ] Each section has specific thesis (testable claim or question)
- [ ] Sections cover all major themes from reviews
- [ ] Section balance (~15-30% of total per section)
- [ ] User approves (if checkpoint enabled)

### Failure Handling
- If user rejects: lit-pm revises outline based on feedback
- If user rejects 2x: Compensation - return to Stage 2 with adjusted scope
- If unbalanced: Merge small sections, split large ones

---

## Stage 4: Introduction Writing

**Owner**: lit-synthesizer + editor (quick polish)
**Checkpoint**: Never (always runs automatically)
**Duration**: 30-45 minutes
**Timeout**: 60 minutes

### Input
- Approved outline

### Process
1. **lit-synthesizer writes introduction** that:
   - Frames the research question and its importance
   - Previews the document structure (sections and their theses)
   - Sets expectations for depth/scope
   - Provides context that section writers will reference

2. **editor does quick polish**:
   - Clarity check, no deep fact-checking yet
   - Ensure introduction aligns with outline

### Output
- Polished introduction (2-4 paragraphs)
- Available to Stage 5 section writers as context

### Quality Gate Checklist
- [ ] Introduction frames research question clearly
- [ ] Preview of structure matches approved outline
- [ ] Length: 2-4 paragraphs (not too long)
- [ ] Editor quick-polish applied

### Failure Handling
- If inconsistent with outline: Return to lit-synthesizer with specific issues
- Timeout: Escalate to user

---

## Stage 5: Parallel Section Research & Writing

**Owner**: literature-researcher agents (specialist continuity from Stage 3)
**Checkpoint**: Never (runs automatically, gated by per-section fact-checking)
**Duration**: 3-5 hours per section (parallel execution)
**Timeout**: 6 hours per section

### Input
- Introduction (context)
- Outline section assignment

### Research Depth Requirements

**CRITICAL**: Section writers are active researchers, not review summarizers.

Stage 2 reviews provide:
- High-level landscape (what's known, what's controversial)
- Key research groups and seminal papers
- Entry points for deeper investigation
- Identification of gaps

Section writers conduct:
- **Targeted primary literature search** for their specific section thesis
- **Forward/backward citation tracking** from review papers
- **Parameter-specific searches** if quantitative data needed
- **Recent papers** (last 2-3 years) that reviews might have missed

### Paper Requirements per Section
- **15-30 papers total**:
  - 10-15 foundational papers (established findings)
  - 5-10 recent papers (last 2-3 years)
  - 3-5 very recent papers (last 6-12 months) for recency survey

### Required: Recency Survey Subsection

Each section MUST include a brief (1-2 paragraph) subsection highlighting most recent literature:

```markdown
### X.X.X Recent Developments (2025-2026)

Three recent studies have begun addressing [specific challenge]. [Author] et al.
(2025) demonstrated that [finding] [citation]. Concurrently, [another group]
reported [finding] [citation]. Most recently, [third study] identified [finding],
suggesting [implication] [citation].

These developments suggest [emerging trend], though [remaining uncertainty].
```

### Writer Autonomy (Moderate Flexibility)
- **Must address assigned thesis** - Core question cannot be changed
- **Can add subsections** - If research reveals need
- **Can adjust structure** - Reorder subsections if logical flow improves
- **Cannot expand scope** - Stay within section boundaries

### Output
- Drafted section with 15-30 papers cited
- Recency survey subsection included
- Ready for Stage 6a validation

### Quality Gate Checklist
- [ ] 15-30 primary papers cited (not just reviews)
- [ ] Recency survey present (3-5 papers from last 6-12 months)
- [ ] Citations include publication dates
- [ ] Section addresses assigned thesis
- [ ] No placeholder text ("TODO", "[CITE]", "[INSERT]")
- [ ] Length reasonable (2000-3000 words)

### Failure Handling
- Timeout: Proceed without section, flag gap for synthesis stage
- Low paper count: Return to writer with specific research directions

---

## Stage 6a: Per-Section Quick Validation (BLOCKING)

**Owner**: fact-checker
**Checkpoint**: Blocking per section
**Duration**: 5-10 minutes per section
**Timeout**: 15 minutes per section

### Input
- Completed section draft

### Quick Checks
- [ ] 15-30 primary papers cited (not just reviews)
- [ ] Recency survey present (3-5 papers from last 6-12 months)
- [ ] Citations include publication dates
- [ ] Quantitative data with units and measurement context
- [ ] Section addresses assigned thesis
- [ ] No contradictions with introduction
- [ ] No placeholder text ("TODO", "[CITE]", "[INSERT]")
- [ ] Length reasonable (2000-3000 words)

### Output
- PASS or REVISION-NEEDED with specific issues

### Blocking Behavior
Section cannot proceed to synthesis until PASS.

### Escape Hatch Protocol (Prevents Infinite Loops)

```yaml
blocking_gate:
  max_revision_cycles: 3

  on_max_reached:
    notify_user: |
      Section {name} has failed fact-check 3 times.

      Issues: [specific issues]

      Options:
      1. Accept section as-is (waive requirement)
      2. Adjust requirements for this section
      3. Assign different researcher to section
      4. Remove section from outline

    record_decision: workflow_state.quality_overrides[]
```

### Failure Handling
- Minor issues: fact-checker provides revision list, section writer fixes (15-30 min)
- Major issues: Section returned to writer for substantial revision
- Max cycles reached: Escalate to user with options
- Timeout: Pass with warning

---

## Stage 6b: Comprehensive Final Fact-Check (NON-BLOCKING)

**Owner**: fact-checker
**Checkpoint**: Never (runs automatically, non-blocking)
**Duration**: 45-90 minutes (after synthesis)
**Timeout**: 120 minutes

### Input
- All sections + introduction

### Deep Checks
- [ ] Cross-section consistency - No contradictory claims between sections
- [ ] Citation accuracy - Spot-check 10 random citations (retrieve papers, verify claims)
- [ ] Quantitative claim verification - Values match sources, units correct, context preserved
- [ ] Gap analysis - Are there obvious missing topics given the scope?
- [ ] Methodological context - Are measurement methods and conditions noted for key data?

### Output
Revision list with priority levels:
- **P0 (critical)**: Factual errors, contradictions - Must fix before delivery
- **P1 (important)**: Missing citations, unclear claims - Should fix
- **P2 (nice-to-have)**: Minor formatting, style consistency - Editor handles

### Non-Blocking
Revision list goes to Stage 8 (editor incorporates during polish).

### Failure Handling
- Timeout: Skip deep check, proceed with P0 items only from quick checks
- Major P0 issues: May require return to Stage 7 for correction

---

## Stage 7: Active Synthesis & Augmentation

**Owner**: lit-synthesizer (senior author role)
**Checkpoint**: High-stakes only
**Duration**: 2-4 hours
**Timeout**: 5 hours

### Input
- Introduction + all validated section drafts

### Process (Active Curation, Not Passive Assembly)

1. **Read all sections as a scientist**:
   - Do the sections collectively answer the research question?
   - What's the narrative arc these sections are telling?
   - Where do sections contradict, overlap, or leave gaps?

2. **Identify cross-cutting themes**:
   - Patterns mentioned in multiple sections
   - Implicit connections that section writers hinted at but didn't state explicitly
   - Emergent insights visible only when viewing all sections together

3. **Authority to restructure** (senior author powers):
   - Reorder sections if logical flow improves
   - Merge similar subsections across different sections
   - Add new subsections to fill gaps (with targeted research, ~30 min per gap)
   - Rewrite transitions and connecting paragraphs for narrative continuity
   - Restructure section endings/beginnings for better handoffs

4. **Write conclusion section**:
   - Synthesize key takeaways across all sections
   - Highlight remaining uncertainties and contradictions
   - Draw project-specific implications (tie back to original research question)
   - Identify future research directions

5. **Flag major additions**:
   - If >20% new content added: Notify user of substantial augmentation
   - Provide rationale for major additions

### Output
- Cohesive research document with conclusion
- Cross-cutting themes identified
- Major additions flagged (if any)

### Quality Gate Checklist
- [ ] Narrative flow is logical (each section builds on previous)
- [ ] Cross-cutting themes identified and woven through document
- [ ] Gaps filled (if lit-synthesizer added content)
- [ ] Conclusion synthesizes findings across all sections
- [ ] Major additions flagged (if >20% new content)

### Failure Handling
- If user rejects at checkpoint: lit-synthesizer revises based on feedback
- If internal quality check fails: lit-synthesizer self-revises before presenting
- If cross-section inconsistency detected: Return to sections with specific feedback
- Timeout: Escalate to user

---

## Stage 8: Editorial Polish

**Owner**: editor
**Checkpoint**: Never (always runs automatically)
**Duration**: 30-60 minutes
**Timeout**: 90 minutes

### Input
- Synthesized document + fact-check revision list from Stage 6b

### Process
1. **Incorporate fact-check revisions** (P0 and P1 priority items from Stage 6b)
2. **Polish for clarity and conciseness**:
   - Remove redundancy
   - Clarify ambiguous phrasing
   - Ensure consistent terminology
3. **Formatting consistency**:
   - Citation format (Nature-style inline)
   - Section numbering
   - Figure/table references (if applicable)
4. **Voice consistency** across sections (smooth over stylistic differences)
5. **Final read-through** for flow and coherence

### Output
- Final polished document ready for delivery

### Quality Gate Checklist
- [ ] Fact-check revision list (P0/P1) incorporated
- [ ] Consistent voice across sections
- [ ] Clear, concise, properly formatted
- [ ] Final read-through complete

### Failure Handling
- Timeout: Deliver as-is with note about incomplete polish
- Major issues found during polish: Return to Stage 7 for correction

---

## Stage Summary Table

| Stage | Owner | Checkpoint | Duration | Timeout |
|-------|-------|------------|----------|---------|
| 1: Scope | requirements-analyst | ALWAYS | 15-30 min | 45 min |
| 2: Reviews | lit-pm + researchers | High-Stakes | 45-90 min | 120 min |
| 3: Outline | lit-pm + researchers | Medium+ | 30-60 min | 90 min |
| 4: Intro | lit-synthesizer + editor | Never | 30-45 min | 60 min |
| 5: Sections | researchers (parallel) | Never | 3-5 hr/section | 6 hr/section |
| 6a: Quick FC | fact-checker | Blocking | 5-10 min/section | 15 min/section |
| 6b: Deep FC | fact-checker | Never | 45-90 min | 120 min |
| 7: Synthesis | lit-synthesizer | High-Stakes | 2-4 hr | 5 hr |
| 8: Polish | editor | Never | 30-60 min | 90 min |
