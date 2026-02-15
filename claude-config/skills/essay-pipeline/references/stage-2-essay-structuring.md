# Stage 2: Essay Structuring

## Stage Overview

- **Duration**: 10-20 minutes (user-paced)
- **Pushback type**: Audience awareness
- **Goal**: Negotiate and agree upon the essay's architecture -- section order, length, audience, tone
- **State anchor prefix**: `[Stage 2/4 - Essay Structuring]`

## Circumscribed Responsibility

Your ONLY job in this stage is to negotiate structure. You do NOT:
- Develop arguments (that is Stage 3)
- Write prose (that is Stage 4)
- Refine the thesis (that was Stage 1)
- Do research (that is the fact-checker's job)
- Evaluate voice (that is the voice-matcher's job)

You produce a structural blueprint that all subsequent stages follow.

## Process

### 1. Read Thesis and Claim Type

Load `stage-1-thesis.md` from the session directory. Understand:
- What needs to be argued/explained/predicted
- The claim type (argument, explainer, prediction, hybrid)
- The domain and preliminary audience

### 2. Propose Initial Structure

Based on claim type, suggest a section breakdown:

**Argument essays:**
- Setup / Hook
- The Common Narrative (what people think)
- The Core Claim (your argument)
- Evidence and Analysis
- Counterarguments and Responses
- Implications / Takeaway

**Explainer essays:**
- Hook / Why This Matters
- Background (what you need to know)
- Core Explanation (the main content)
- Nuance and Caveats
- Takeaway / What's Next

**Prediction essays:**
- Current State / Context
- The Trend (what's changing)
- The Projection (your prediction)
- Caveats and Uncertainty
- Implications / If This Happens

**Hybrid essays:**
- Combine elements as appropriate
- Note which sections are argument-driven vs. explanatory

Present the proposal with word count estimates per section.

### 3. Negotiate Length

Propose a word count based on complexity and section count:

| Length Category | Word Range | When Appropriate |
|----------------|-----------|-----------------|
| Short | 800-1,500 | Single-point arguments, brief explainers, focused predictions |
| Medium | 1,500-3,500 | Standard blog essays, multi-section arguments |
| Long | 3,500-7,000+ | Deep dives, multi-part arguments, comprehensive explainers |

**Check proportionality:** Section word targets should sum to approximately the total. Flag if any section is too thin (<150 words) to develop its purpose properly.

### 4. Apply Audience Awareness Pushback

Challenge the structure based on how the target audience will experience the essay.

**Pushback behaviors by level:**

#### Full Pushback (default)
- Challenge reader assumptions: "Your structure assumes the reader understands [concept]. Is that safe for a general audience?"
- Challenge section necessity: "Section N doesn't advance your argument. What purpose does it serve for the reader?"
- Challenge length allocation: "You're giving 200 words to your strongest argument and 500 to background. Consider flipping those proportions."
- Challenge transition logic: "A reader going from Section 2 to Section 3 will experience a jarring shift. How should we bridge?"
- Challenge pacing: "Three technical sections in a row may lose general readers. Consider alternating with a narrative section."
- **Maximum 3 pushback rounds per structural point**

#### Light Pushback
- One round of audience-focused observations
- Accept user's structural decisions without follow-up
- **Maximum 1 round per point**

#### Minimal Pushback
- Observations only: "A general reader might find [X] challenging here."
- **0 challenge rounds**

### 5. Iterate Until Agreement

The user and you converge on a structure. Accept unusual structures if the user articulates a reason.

## Contradiction Detection

Cross-check user instructions against the inherited thesis:
- Does the proposed structure actually support the thesis?
- Does the audience match what was discussed in Stage 1?
- Are there sections that contradict the thesis direction?

If contradictions found: "I notice your proposed structure for Section 3 [describes X], but your thesis argues [Y]. These seem to pull in different directions. Which should we adjust?"

## Science Blog Conventions

Offer these as defaults (user can override):

- **Lead with the conclusion**: Science blog readers want the payoff early, not at the end
- **One idea per post**: If the thesis has two independent claims, consider two posts
- **Define jargon early**: Technical terms should be defined on first use, not assumed
- **Use analogies**: One good analogy per complex concept helps general readers
- **End with forward-looking implications**: Blog readers want to know "so what" and "what's next"

## Exit Criteria Checklist

Before moving to Stage 3, verify ALL of the following:

- [ ] Outline has defined sections with word count targets
- [ ] Total word count is agreed upon
- [ ] Audience is specified (general, semi-technical, expert)
- [ ] Each section has a stated purpose (not just a topic)
- [ ] User has explicitly approved the outline
- [ ] Section word targets sum to approximately the total target (+/- 10%)
- [ ] No contradictions between outline and thesis

## Output Format

Write to session directory:

### stage-2-outline.md

```markdown
# Essay Outline: [Working Title]

## Thesis
[Approved thesis from Stage 1]

## Parameters
- **Target length**: [N] words
- **Audience**: [general / semi-technical / expert] -- [description]
- **Tone**: [description, referencing style profile]
- **Sections**: [N]
- **Claim type**: [argument / explainer / prediction / hybrid]

## Structure

### Section 1: [Title] (~[N] words)
**Purpose**: [What this section accomplishes for the argument/explanation]
**Key content**:
- [Bullet point 1]
- [Bullet point 2]
- [Bullet point 3]
**Transition to Section 2**: [How this connects to the next section]

### Section 2: [Title] (~[N] words)
**Purpose**: [What this section accomplishes]
**Key content**:
- [Bullet point 1]
- [Bullet point 2]
**Transition to Section 3**: [Connection]

[Continue for all sections]

## Structural Notes
- [Any special considerations discussed during negotiation]
- [Audience assumptions noted]
- [Conventions adopted or rejected]

## Approved By User
Yes -- [date/time]
```

## Handoff to Stage 3

Pass forward:
- `stage-2-outline.md` file path
- `stage-1-thesis.md` file path
- Style profile path
- Total section count
- Starting section index: 1

Update session state:
```yaml
stages:
  stage_2:
    status: completed
    outline_file: "{session_dir}/stage-2-outline.md"
    completed_at: "{ISO8601}"
  stage_3:
    status: in_progress
    total_sections: N
    current_section: 1
```
