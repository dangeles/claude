# Literature Review Engine Design

**Date**: 2026-02-03
**Purpose**: Transform technical-pm/researcher/synthesizer into specialized literature curation/compilation/review engine
**Target Use Cases**: Internal research synthesis (decision-focused) and literature surveys (landscape mapping)
**Orchestration Style**: Adaptive with smart checkpoints based on complexity, stakes, and user hints

---

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [The 8-Stage Pipeline](#the-8-stage-pipeline)
3. [Adaptive Orchestration Logic](#adaptive-orchestration-logic)
4. [lit-synthesizer: The Senior Scientific Author](#lit-synthesizer-the-senior-scientific-author)
5. [Quality Gates & Handoff Protocols](#quality-gates--handoff-protocols)
6. [Implementation Roadmap](#implementation-roadmap)

---

## Architecture Overview

### Three-Tier Skill Architecture

**Tier 1: Orchestrator (NEW)**
- **`lit-pm`** - Literature Pipeline Manager
  - Coordinates 8-stage literature review pipeline
  - Implements adaptive orchestration (complexity detection → checkpoint plan)
  - Manages parallel review discovery with convergence tracking
  - Handles workflow state, handoffs, and quality gates
  - Built on Phase 3 parallel execution capabilities from technical-pm

**Tier 2: Specialized Literature Skills (NEW + Enhanced)**
- **`literature-researcher`** (enhanced from `researcher`)
  - Retains existing research methodology and personality
  - Adds literature-specific search strategies (review discovery, convergence tracking)
  - Can draft section outlines from discovered literature
  - Conducts deep targeted research (15-30 papers per section)
  - Includes recency surveys (last 6-12 months)

- **`lit-synthesizer`** (NEW - distinct from `synthesizer`)
  - **Personality: Senior Scientific Author**
  - Logical yet creative - finds non-obvious connections
  - Shapes narrative coherence across disparate sections
  - Authority to restructure, rewrite, add analysis
  - Treats sections as "material" not "final text"
  - Writes introductions/conclusions that frame the whole document

**Tier 3: Supporting Skills (Existing)**
- `fact-checker` - Quick validation + comprehensive final review
- `editor` - Final polish (clarity, conciseness, formatting)
- `requirements-analyst` - Prompt refinement (Stage 1)

### Key Architectural Decisions

1. **lit-pm orchestrates, doesn't execute** - Delegates to specialists
2. **lit-synthesizer is elevated role** - Not just stitching, but authoring with senior authority
3. **Existing skills preserved** - General researcher/synthesizer still available for non-literature work
4. **Parallel execution ready** - Leverages Phase 3 technical-pm parallel capabilities
5. **Adaptive checkpoints** - System proposes checkpoint plan, user can override

---

## The 8-Stage Pipeline

### Stage 1: Scope Refinement
**Owner**: `requirements-analyst` (existing skill)
**Checkpoint**: Always (required)
**Duration**: 15-30 minutes

**Input**: User's initial prompt
**Output**: Refined research scope document

**Process**:
1. Clarify research question, target audience, intended use
2. Define success criteria (decision support? landscape map? comprehensive review?)
3. Set boundaries (in-scope topics, out-of-scope tangents)
4. **Complexity detection** → Propose checkpoint plan
5. User approves scope + checkpoint plan

**Quality Gate**:
- [ ] Research question is specific (not "what is X?" but "what governs X?")
- [ ] Success criteria are measurable
- [ ] In-scope and out-of-scope boundaries clear
- [ ] User approves scope and checkpoint plan

**Handoff**: Scope document → lit-pm

---

### Stage 2: Parallel Review Discovery
**Owner**: `lit-pm` orchestrates 2-3 `literature-researcher` agents in parallel
**Checkpoint**: Only if high-stakes
**Duration**: 45-90 minutes (parallel execution)

**Input**: Refined scope
**Output**: Convergence-tracked review collection (6-9 reviews total)

**Process**:
1. **lit-pm generates 2-3 diverse search strategies**:
   - Strategy A: Broad keywords (e.g., "CAR-T manufacturing review")
   - Strategy B: Specific technical terms (e.g., "transduction efficiency optimization")
   - Strategy C: Application focus (e.g., "clinical CAR-T production challenges")

2. **Launch 3 literature-researcher agents in parallel** (using Phase 3 parallel execution):
   - Each finds 2-3 high-quality reviews using their strategy
   - Each writes brief annotation per review (key findings, coverage, quality assessment)

3. **lit-pm analyzes convergence**:
   - Which reviews appeared in multiple searches? → High signal (must-read)
   - Which unique perspectives from divergent searches? → Coverage breadth

**Quality Gate**:
- [ ] 6-9 reviews collected total (2-3 per agent)
- [ ] At least 2 reviews show convergence (found by multiple agents)
- [ ] Reviews published in last 10 years (unless older review is definitive)
- [ ] Coverage of all major themes in scope

**Failure Handling**:
- If convergence low (<2 reviews in common): Adjust search strategies, re-run
- If coverage gaps: Add focused search for missing theme
- If low quality reviews: Expand search, prioritize high-citation reviews

**Handoff**: Annotated review collection + convergence analysis → Stage 3

---

### Stage 3: Layered Outline Synthesis
**Owner**: `lit-pm` (high-level structure) + `literature-researcher` agents (section details)
**Checkpoint**: Medium/High complexity
**Duration**: 30-60 minutes

**Input**: Review collection with convergence data
**Output**: Structured outline with section theses

**Process**:
1. **lit-pm reads all reviews, creates high-level structure**:
   - 3-5 sections (not too broad, not too narrow)
   - Each section has: Title, Core thesis (testable claim or question), Key question to address

2. **For each section: Assign detail-work to literature-researcher agent**:
   - Agent proposes 2-4 subsections with specific questions
   - Agent identifies which reviews/papers to cite per subsection
   - Agent notes potential gaps or areas needing targeted research

3. **lit-pm presents complete outline to user** (if checkpoint):
   - User can: Approve, request changes, reprioritize sections

**Quality Gate**:
- [ ] 3-5 sections (balanced scope)
- [ ] Each section has specific thesis (testable claim or question)
- [ ] Sections cover all major themes from reviews
- [ ] Section balance (~15-30% of total per section)
- [ ] User approves (if checkpoint enabled)

**Failure Handling**:
- If user rejects: lit-pm revises outline based on feedback
- If unbalanced: Merge small sections, split large ones

**Handoff**: Approved outline with section assignments → Stage 4

---

### Stage 4: Introduction Writing
**Owner**: `lit-synthesizer`
**Checkpoint**: Never (always runs automatically)
**Duration**: 30-45 minutes

**Input**: Approved outline
**Output**: Introduction section (2-4 paragraphs)

**Process**:
1. **lit-synthesizer writes introduction** that:
   - Frames the research question and its importance
   - Previews the document structure (sections and their theses)
   - Sets expectations for depth/scope
   - Provides context that section writers will reference

2. **`editor` does quick polish**:
   - Clarity check, no deep fact-checking yet
   - Ensure introduction aligns with outline

**Quality Gate**:
- [ ] Introduction frames research question clearly
- [ ] Preview of structure matches approved outline
- [ ] Length: 2-4 paragraphs (not too long)
- [ ] Editor quick-polish applied

**Handoff**: Polished introduction → Stage 5 (section writers can reference during writing)

---

### Stage 5: Parallel Section Research & Writing
**Owner**: `literature-researcher` agents (specialist continuity from Stage 3)
**Checkpoint**: Never (runs automatically, gated by per-section fact-checking)
**Duration**: 3-5 hours per section (parallel execution)

**Input**: Introduction + Outline section assignment
**Output**: Drafted sections with moderate flexibility

**Process**:

**CRITICAL: Section writers are active researchers, not review summarizers**

The reviews from Stage 2 provide **foundational landscape mapping**, not final research. Each section writer must conduct **substantial new targeted research** for their specific thesis.

**Stage 2 reviews provide**:
- High-level landscape (what's known, what's controversial)
- Key research groups and seminal papers
- Entry points for deeper investigation
- Identification of gaps

**Stage 5 section writers conduct**:
- **Targeted primary literature search** for their specific section thesis
- **Forward/backward citation tracking** from review papers to find detailed studies
- **Parameter-specific searches** if quantitative data needed
- **Recent papers** (last 2-3 years) that reviews might have missed
- **15-30 primary papers per section** (not just 2-3 reviews)

**Research depth per section**:
- **15-30 papers total**:
  - 10-15 foundational papers (established findings)
  - 5-10 recent papers (last 2-3 years)
  - 3-5 very recent papers (last 6-12 months) for recency survey

**Required: Recency Survey Subsection**

Each section must include a brief (1-2 paragraph) subsection highlighting most recent literature:

**Example**:
```markdown
### 3.2.4 Recent Developments (2025-2026)

Three recent studies have begun addressing the MOI optimization challenge
at clinical scale. Wang et al. (2025) demonstrated that adaptive MOI
titration based on real-time cell density monitoring improved transduction
efficiency by 23% in suspension culture [47]. Concurrently, the Berlin
group reported that pulsed vector addition (rather than single bolus)
reduced vector waste while maintaining >80% transduction [48]. Most
recently, Lee et al. (2026) identified a previously overlooked pH
sensitivity in the transduction process, suggesting that bioreactor pH
control strategies require revision [49].

These developments suggest the field is moving toward closed-loop process
control for vector addition, though clinical validation remains pending.
```

**Purpose of recency survey**:
- Captures cutting-edge work that reviews might have missed
- Shows reader the document is current (not relying on 5-year-old reviews)
- Identifies emerging trends and future directions
- Prevents document from being outdated on publication

**Writer autonomy (moderate flexibility)**:
- **Must address assigned thesis** - Core question cannot be changed
- **Can add subsections** - If research reveals need (e.g., discovered important sub-mechanism)
- **Can adjust structure** - Reorder subsections if logical flow improves
- **Cannot expand scope** - Stay within section boundaries (avoid scope creep)

**Section structure guidance**:
- Core content: 3-4 subsections addressing thesis
- Final subsection: Recent Developments (2025-2026)
- Length: 4-6 pages per section

**As each section completes** → Quick fact-check gate (Stage 6a)

**Handoff**: All drafted sections (after fact-checking) → Stage 7

---

### Stage 6a: Per-Section Quick Validation (Blocking)
**Owner**: `fact-checker`
**Checkpoint**: Blocking per section
**Duration**: 5-10 minutes per section

**Input**: Completed section draft
**Output**: PASS or REVISION-NEEDED with specific issues

**Quick checks** (fast validation):
- [ ] **15-30 primary papers cited** (not just reviews)
- [ ] **Recency survey present** (3-5 papers from last 6-12 months)
- [ ] Citations include publication dates (to verify recency)
- [ ] Quantitative data with units and measurement context
- [ ] Section addresses assigned thesis
- [ ] No contradictions with introduction
- [ ] No placeholder text ("TODO", "[CITE]", "[INSERT]")
- [ ] Length reasonable (4-6 pages)

**Blocking behavior**: Section cannot proceed to synthesis until PASS

**Failure Handling**:
- **Minor issues**: fact-checker provides revision list → section writer fixes (15-30 min)
- **Major issues** (off-topic, insufficient citations, missing recency survey): Section returned to writer for substantial revision

---

### Stage 6b: Comprehensive Final Fact-Check (Non-Blocking)
**Owner**: `fact-checker`
**Checkpoint**: Never (runs automatically, non-blocking)
**Duration**: 45-90 minutes (after synthesis)

**Input**: All sections + introduction
**Output**: Revision list (flagged issues, not blocking)

**Deep checks** (comprehensive validation):
- [ ] **Cross-section consistency** - No contradictory claims between sections
- [ ] **Citation accuracy** - Spot-check 10 random citations (retrieve papers, verify claims)
- [ ] **Quantitative claim verification** - Values match sources, units correct, context preserved
- [ ] **Gap analysis** - Are there obvious missing topics given the scope?
- [ ] **Methodological context** - Are measurement methods and conditions noted for key data?

**Output**: Revision list with priority levels
- **P0 (critical)**: Factual errors, contradictions → Must fix before delivery
- **P1 (important)**: Missing citations, unclear claims → Should fix
- **P2 (nice-to-have)**: Minor formatting, style consistency → Editor handles

**Non-blocking**: Revision list goes to Stage 8 (editor incorporates during polish)

---

### Stage 7: Active Synthesis & Augmentation
**Owner**: `lit-synthesizer` (senior author role)
**Checkpoint**: High-stakes only
**Duration**: 2-4 hours

**Input**: Introduction + all validated section drafts
**Output**: Cohesive research document

**Process** (active curation, not passive assembly):

1. **Read all sections as a scientist**:
   - Do the sections collectively answer the research question?
   - What's the narrative arc these sections are telling?
   - Where do sections contradict, overlap, or leave gaps?

2. **Identify cross-cutting themes**:
   - Patterns mentioned in multiple sections (e.g., "three sections mentioned oxygen gradients")
   - Implicit connections that section writers hinted at but didn't state explicitly
   - Emergent insights visible only when viewing all sections together

3. **Authority to restructure** (senior author powers):
   - **Reorder sections** if logical flow improves
   - **Merge similar subsections** across different sections
   - **Add new subsections** to fill gaps (with targeted research, ~30 min per gap)
   - **Rewrite transitions** and connecting paragraphs for narrative continuity
   - **Restructure section endings/beginnings** for better handoffs

4. **Write conclusion section**:
   - Synthesize key takeaways across all sections
   - Highlight remaining uncertainties and contradictions
   - Draw project-specific implications (tie back to original research question)
   - Identify future research directions

5. **Flag major additions**:
   - If >20% new content added: Notify user of substantial augmentation
   - Provide rationale for major additions

**Quality Gate**:
- [ ] Narrative flow is logical (each section builds on previous)
- [ ] Cross-cutting themes identified and woven through document
- [ ] Gaps filled (if lit-synthesizer added content)
- [ ] Conclusion synthesizes findings across all sections
- [ ] Major additions flagged (if >20% new content)

**Failure Handling**:
- If user rejects at checkpoint: lit-synthesizer revises based on feedback
- If internal quality check fails: lit-synthesizer self-revises before presenting

**Handoff**: Synthesized document + revision notes → Stage 8

---

### Stage 8: Editorial Polish
**Owner**: `editor`
**Checkpoint**: Never (always runs automatically)
**Duration**: 30-60 minutes

**Input**: Synthesized document + fact-check revision list from Stage 6b
**Output**: Final polished document

**Process**:
1. **Incorporate fact-check revisions** (P0 and P1 priority items from Stage 6b)
2. **Polish for clarity and conciseness**:
   - Remove redundancy
   - Clarify ambiguous phrasing
   - Ensure consistent terminology
3. **Formatting consistency**:
   - Citation format (Nature-style inline)
   - Section numbering
   - Figure/table references (if applicable)
4. **Voice consistency** across sections (smooth over stylistic differences between section writers)
5. **Final read-through** for flow and coherence

**Quality Gate**:
- [ ] Fact-check revision list (P0/P1) incorporated
- [ ] Consistent voice across sections
- [ ] Clear, concise, properly formatted
- [ ] Final read-through complete

**Handoff**: Deliver final document to user

---

## Adaptive Orchestration Logic

### Complexity Detection (Stage 1)

When `lit-pm` receives a refined scope from `requirements-analyst`, it analyzes three dimensions:

**Dimension 1: Scope Indicators**
- **Paper count estimate**: <5 (simple), 5-15 (medium), 15+ (complex)
- **Topic breadth**: Single topic, 2-3 themes, or cross-domain
- **Literature maturity**: Well-established field vs. emerging/contradictory

**Dimension 2: Stakes Indicators**
- **Use case**: Exploration, decision support, or high-stakes deliverable
- **Keyword signals**:
  - Low stakes: "quick survey", "what's known about", "explore"
  - Medium stakes: "inform decision", "compare approaches", "assess feasibility"
  - High stakes: "grant background", "comprehensive review", "publication"

**Dimension 3: User Hints**
- **Explicit flags**: `--review-outline`, `--review-drafts`, `--full-auto`
- **Time constraints**: "quick", "thorough", "comprehensive"
- **Prior context**: Has user iterated on similar documents? (Trust built)

### Checkpoint Plan Generation

Based on detection, `lit-pm` proposes one of four checkpoint plans:

| Complexity Level | Stage 1 (Scope) | Stage 2 (Reviews) | Stage 3 (Outline) | Stage 7 (Synthesis) | Rationale |
|------------------|-----------------|-------------------|-------------------|---------------------|-----------|
| **Simple** | ✓ Always | Auto | Auto | Auto | Scope approval sufficient |
| **Medium** | ✓ Always | Auto | ✓ Checkpoint | Auto | Direction check before heavy lifting |
| **Complex** | ✓ Always | Auto | ✓ Checkpoint | ✓ Checkpoint | Multiple approval points |
| **High-Stakes** | ✓ Always | ✓ Checkpoint | ✓ Checkpoint | ✓ Checkpoint | Maximum oversight |

### Example Detection Scenarios

**Example 1: Simple**
```
User: "Quick survey of CAR-T manufacturing challenges"

Detection:
- Scope: <10 papers, single domain → Simple
- Stakes: "quick survey" → Low
- No explicit flags

Proposed Plan: SIMPLE (only Stage 1 checkpoint)

System message:
"I'll run the full pipeline automatically after scope approval.
You can add --review-outline if you want to check structure first."
```

**Example 2: High-Stakes**
```
User: "Comprehensive review of mRNA delivery for tissue engineering -
       this will inform our $2M grant proposal"

Detection:
- Scope: 15+ papers, cross-domain (mRNA + tissue eng) → Complex
- Stakes: "grant proposal", "$2M" → High stakes
- "Comprehensive" → High thoroughness

Proposed Plan: HIGH-STAKES (all checkpoints)

System message:
"Given the grant context, I'll pause for approval at:
- Review discovery (ensure we have the right foundational papers)
- Outline structure (get section framing right before writing)
- Final synthesis (review before editorial polish)"
```

### User Override

After proposing checkpoint plan, user can:
- **Accept**: "Looks good, proceed"
- **Reduce checkpoints**: "Skip review discovery checkpoint, trust you on that"
- **Add checkpoints**: "Also pause after section drafts"
- **Full auto**: `--full-auto` (skip all optional checkpoints, only Stage 1 remains)

### Mid-Pipeline Adaptation

If during execution `lit-pm` detects:
- **Unexpected complexity**: Literature more contradictory than expected → Offer checkpoint before synthesis
- **Faster than expected**: Simple topic, high-quality reviews → Suggest removing remaining checkpoints
- **Blocker**: Missing critical information → Mandatory checkpoint for user guidance

---

## lit-synthesizer: The Senior Scientific Author

### Personality & Philosophy

**You are a senior scientific author** - the kind who reads a stack of disparate papers and sees the narrative thread others miss. You think in **logical structures and creative connections simultaneously**.

You're the researcher who:
- Reads five papers on different topics and notices they're all dancing around the same underlying mechanism
- Recognizes when two "contradictory" findings actually illuminate different aspects of the same phenomenon
- Knows when to preserve technical precision and when to step back for the big picture
- Can write the paragraph that makes a reader say "Oh, *that's* why this matters"

**You are not a compiler**. You don't just stitch sections together with transitions. You're the **senior author** who shapes the scientific narrative. Section drafts are your **material**, not your final text.

### Core Responsibilities

#### 1. Active Curation (Not Passive Assembly)

When you receive section drafts:
- **Read as a scientist first**: Do the sections collectively answer the research question?
- **Identify the narrative arc**: What's the story these sections are trying to tell?
- **Find the tensions**: Where do sections contradict, overlap, or leave gaps?

#### 2. Logical Structure

You ensure:
- **Argument flow**: Each section builds on previous sections logically
- **Evidence hierarchy**: Key claims are supported before dependent claims are made
- **Scope boundaries**: Document stays focused, tangents are pruned
- **Progressive disclosure**: Complex ideas are introduced before they're built upon

#### 3. Creative Synthesis

You add value by:
- **Cross-cutting insights**: "Three sections mentioned oxygen gradients—there's a unifying principle here"
- **Implicit connections**: Making explicit what section writers hinted at but didn't state
- **Emergent themes**: Patterns visible only when viewing all sections together
- **Non-obvious implications**: Drawing conclusions that follow from the evidence but weren't explicitly stated

#### 4. Gap Filling Authority

You can:
- **Add subsections**: If sections reveal "wait, no one covered manufacturing costs"
- **Conduct targeted research**: Quick literature searches to fill specific gaps (~30 min per gap)
- **Restructure sections**: Reorder, split, or merge sections for better flow
- **Rewrite transitions**: Not just add connecting paragraphs, but rewrite section endings/beginnings for continuity

#### 5. Introduction & Conclusion Ownership

These are **your sections**, not delegated:
- **Introduction**: Frames the entire document, sets expectations, provides context others will reference
- **Conclusion**: Synthesizes findings, highlights uncertainties, draws project-specific implications

### What You DON'T Do

- **Don't preserve section text verbatim if it hurts the narrative** - You're senior author, not copy editor
- **Don't add major new research directions without flagging** - If you want to add >20% new content, tell the user
- **Don't remove technical precision for "readability"** - This is scientific writing, precision matters
- **Don't resolve genuine scientific uncertainty by picking a side** - Acknowledge and explain tensions

### Writing Style: Logical Yet Creative

- Use **active voice** and **clear causal language** ("X causes Y" not "X is associated with Y" when causal relationship is established)
- **Signpost transitions**: "This raises the question..." "Building on this mechanism..." "However, three studies contradict..."
- **Integrate without flattening**: Preserve important nuance while creating coherent flow
- **Write for the skeptical expert**: Assume reader is smart but unfamiliar with this specific topic

### Example: Authority in Action

**Section Draft** (from section writer):
```markdown
### 3.2 Oxygen Transport Limitations

Hepatocyte oxygen consumption is 0.4-0.8 nmol/s per 10^6 cells [3].
Diffusion limits transport beyond 100-150 μm [5]. Several studies
have examined bioreactor oxygenation strategies [7,8,9].
```

**lit-synthesizer Active Curation**:
```markdown
### 3.2 Why Oxygen Becomes the Bottleneck

The mismatch between hepatocyte oxygen demand and delivery defines the
core engineering challenge. Hepatocytes consume oxygen at 0.4-0.8 nmol/s
per 10^6 cells [3]—among the highest rates of any mammalian cell type.
Passive diffusion can only supply this rate within 100-150 μm of an
oxygen source [5], creating a fundamental spatial constraint.

This explains the contradictory scaling results in Section 2.3:
bioreactors that worked at 10^6 cells failed at 10^8 cells not because
the design was wrong, but because the oxygen diffusion distance grew
beyond the critical threshold. Three engineering strategies attempt to
address this [7,8,9], each with distinct trade-offs we examine next.
```

**Changes made**:
- Reframed heading to emphasize the "why" (logical structure)
- Connected to earlier section (Section 2.3) that section writer didn't reference (creative synthesis)
- Explained the mechanism, not just cited it (scientific authorship)
- Signposted the transition to next content (narrative flow)
- Added context ("among the highest rates") that section writer knew but didn't state

---

## Quality Gates & Handoff Protocols

### Handoff Document Format

Every stage produces a **structured handoff document** (building on Phase 1 `handoff-format.md`):

```yaml
stage: 3  # Outline Synthesis
status: complete
producer: lit-pm
consumer: literature-researcher (section writers)
timestamp: 2026-02-03T14:30:00Z

# Core content
outline:
  introduction:
    thesis: "CAR-T manufacturing faces scalability challenges"
    assigned_to: lit-synthesizer

  section_1:
    title: "Vector Production Bottlenecks"
    thesis: "Lentiviral vector yield limits clinical-scale production"
    subsections:
      - "Transduction efficiency trade-offs"
      - "MOI optimization in suspension culture"
    assigned_to: literature-researcher-agent-1
    key_reviews: [review_id_3, review_id_7]
    research_mandate: |
      15-30 primary papers for this section.
      Include recency survey (last 6-12 months).
      Focus on quantitative yield data and MOI optimization.

  section_2:
    title: "T-Cell Expansion Kinetics"
    thesis: "Ex vivo expansion time directly impacts product cost"
    subsections:
      - "Growth factor optimization"
      - "Bioreactor geometry effects"
    assigned_to: literature-researcher-agent-2
    key_reviews: [review_id_2, review_id_5]
    research_mandate: |
      15-30 primary papers for this section.
      Include recency survey (last 6-12 months).
      Focus on expansion kinetics data and cost modeling.

# Quality metadata
quality_checks_passed:
  - user_approved: true
  - sections_balanced: true  # No single section >40% of total
  - theses_specific: true    # Each thesis is testable/addressable

# Context for next stage
notes:
  - "Section 1 writer: Review #3 has excellent vector production kinetics data"
  - "Section 2 writer: Be mindful of overlap with Section 3 (cell viability)"
  - "All: Introduction will be available before you start writing"
```

### State Management for Resumption

If workflow is interrupted (user Ctrl+C, system issue), state is preserved:

```yaml
workflow_state:
  workflow_id: "lit-review-cart-manufacturing-20260203"
  stage_current: 5  # Section Writing
  stage_completed: [1, 2, 3, 4]
  checkpoints_remaining: [7]  # Synthesis checkpoint still pending

  artifacts:
    scope: "/path/to/refined-scope.yaml"
    reviews: "/path/to/review-collection.yaml"
    outline: "/path/to/approved-outline.yaml"
    introduction: "/path/to/intro-draft.md"
    sections_complete: [1, 2]  # Sections 1 and 2 finished
    sections_in_progress: [3]  # Section 3 was in-progress when interrupted

  resume_options:
    - "Continue from Section 3 (in-progress)"
    - "Restart Section 3 from outline"
    - "Skip to synthesis with Sections 1-2 only (not recommended)"
```

User can resume with: `lit-pm --resume workflow-id`

---

## Implementation Roadmap

### Phase 1: Create lit-pm Orchestrator
- New skill: `lit-pm`
- Implements 8-stage pipeline workflow
- Adaptive complexity detection and checkpoint planning
- Parallel review discovery (leverages Phase 3 parallel execution)
- Layered outline synthesis
- Workflow state management and resumption

### Phase 2: Create lit-synthesizer
- New skill: `lit-synthesizer`
- Senior scientific author personality
- Active curation capabilities
- Introduction/conclusion writing
- Gap filling with targeted research authority
- Cross-cutting theme identification

### Phase 3: Enhance literature-researcher
- Enhanced from existing `researcher` skill
- Add literature-specific modes
- Review discovery with convergence tracking
- Deep section research (15-30 papers per section)
- Recency survey requirements (last 6-12 months)
- Section outline detailing capabilities

### Phase 4: Integration & Testing
- Integrate with existing fact-checker and editor
- Test full pipeline on representative use cases:
  - Simple: Quick survey (<5 papers)
  - Medium: Decision support (5-15 papers)
  - Complex: Comprehensive review (15+ papers, cross-domain)
- Validate adaptive checkpoint logic
- Validate workflow resumption after interruption

### Phase 5: Documentation & Examples
- Create example workflows in `examples/`
- Document checkpoint override flags
- Create troubleshooting guide
- Add pipeline visualization diagram

---

## Success Criteria

1. **Full pipeline execution**: User can invoke `lit-pm` with prompt, receive polished document
2. **Adaptive checkpoints work**: System correctly detects complexity and proposes appropriate checkpoint plan
3. **Parallel review discovery**: 2-3 agents find reviews with convergence tracking
4. **Deep section research**: Sections contain 15-30 primary papers + recency survey
5. **Active synthesis**: lit-synthesizer demonstrably restructures/augments sections (not just stitching)
6. **Quality gates pass**: Fact-checking catches issues, blocking gates prevent bad sections from proceeding
7. **Workflow resumption**: Interrupted workflows can be resumed from last checkpoint
8. **User satisfaction**: Documents are decision-useful and current (not outdated)

---

## Next Steps

1. **Design approval**: Review this document with user
2. **Write to repository**: Commit design to `docs/plans/`
3. **Implementation**: Use `/skill-editor` workflow to create:
   - `lit-pm` skill (new)
   - `lit-synthesizer` skill (new)
   - `literature-researcher` skill (enhanced from `researcher`)
4. **Testing**: Validate on representative use cases
5. **Iteration**: Refine based on real-world usage

---

**Design Status**: Ready for implementation
**Estimated Implementation Time**: 2-3 weeks (3 new/enhanced skills + integration)
**Risk Level**: Medium (new orchestrator pattern, but builds on proven Phase 3 parallel execution)
