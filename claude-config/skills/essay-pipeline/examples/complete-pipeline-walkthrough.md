# Complete Pipeline Walkthrough

This example demonstrates the essay pipeline from start to finish using a CRISPR delivery bottleneck topic. Exchanges are condensed -- a real session would have more back-and-forth.

---

## Session Initialization

```
User: /essay-pipeline

Orchestrator: I'm using the essay-pipeline skill for interactive essay writing.

[Pre-Flight] Checking prerequisites...

Style profile: Found at essays/style-profile.md
Sample essays: Found 4 essays in essays/samples/
WebSearch: Available via fact-checker sub-agent

Session created: /tmp/essay-pipeline-20260214-143022/

Before we begin, how much intellectual pushback would you like?

- Full (default): I'll challenge your ideas vigorously to strengthen them.
  3 rounds of pushback in thesis development, 2 in argument development.
- Light: I'll raise concerns once and accept your response.
  1 round per point.
- Minimal: I'll offer observations but won't challenge.
  No pushback rounds.

You can change this at any time during the session.

User: Full.

Orchestrator: Pushback level set to Full. Let's begin with thesis development.
What topic do you want to write about?
```

---

## Stage 1: Thesis Development

```
[Stage 1/4 - Thesis Development] Beginning Socratic dialogue

User: I want to write about how CRISPR is changing medicine.

Orchestrator: That's a broad topic. Before I ask questions, let me
make sure I understand: you're interested in the intersection of
CRISPR technology and clinical medicine -- how it's being applied
to treat diseases. Is that right?

User: Yes, specifically that people focus on the wrong things about
CRISPR.

Orchestrator: Interesting. "People focus on the wrong things" -- that
suggests you have an argument, not just an explanation. Let me probe:

1. What specifically do people focus on that you think is misplaced?
2. What should they be focusing on instead?
3. Who are "people" here -- the media, researchers, investors, the
   general public?

User: People -- mainly the media and general public -- focus on
editing precision. But the real bottleneck is delivery: getting the
editor into the right cells. We've solved the precision problem but
can't get the therapy where it needs to go.

Orchestrator: Now we're getting somewhere. Let me push on this:

You're claiming the precision problem is "solved." That's a strong
claim. Is off-target editing truly a solved problem, or is it
more accurate to say it's been dramatically improved while delivery
remains comparatively neglected?

User: Fair point. It's dramatically improved, not solved. The
contrast is what matters -- precision has improved enormously while
delivery hasn't kept pace.

Orchestrator: Much more defensible. Here's a draft thesis:

"The biggest obstacle to CRISPR-based medicine is not editing
precision -- which has improved roughly fifty-fold in a decade --
but the unsolved problem of delivering gene editors to the right
cells in a living patient."

This is an argument thesis. The claim type is "argument." You're
arguing that the bottleneck is misidentified by the public, and
that delivery deserves more attention.

Is this the thesis you want to build your essay around?

User: Yes, that captures it.

Orchestrator: [Stage 1/4 - Thesis Development] Quality Gate G1 check:
- Thesis is a clear, specific claim: YES
- Defensible with evidence: YES (precision improvement data, clinical
  trial outcomes, delivery challenges)
- Interesting ("so what?"): YES (reframes public understanding of CRISPR)
- Claim type: argument
- User approved: YES

Thesis approved. Saving to session directory.

Ready to move to Stage 2: Essay Structuring?
```

---

## Stage 2: Essay Structuring

```
[Stage 2/4 - Essay Structuring] Negotiating structure

Orchestrator: Based on your argument thesis, here's a proposed
structure for a medium-length essay (~2,000 words):

1. Hook (~200 words): Open with a specific CRISPR clinical trial
   that succeeded at editing but struggled with delivery
2. The Precision Narrative (~350 words): What the media and public
   focus on -- the editing accuracy story
3. The Delivery Problem (~500 words): Core argument -- why delivery
   is the real bottleneck
4. Current Delivery Approaches (~400 words): What's being tried
   (lipid nanoparticles, AAV, etc.)
5. What Would Change the Game (~300 words): Forward-looking
6. Takeaway (~250 words): Restate thesis with implications

Two concerns:

First, Section 4 risks becoming a technology survey. For a general
audience, you might want to focus on 2-3 approaches and explain
them well rather than listing everything.

Second, at 2,000 words, you're committing to roughly Quanta Magazine
feature depth. Is that right, or would a tighter 1,200-word piece
serve the argument better?

User: 2,000 words is right. And I'll focus on just lipid
nanoparticles and AAV in Section 4.

Orchestrator: Good choices. One more thing: your audience. You said
"general public" during thesis development. Your outline includes
terms like "lipid nanoparticles" and "AAV." A truly general audience
won't know these. Options:

(a) Define them in context (adds ~50 words per term)
(b) Narrow your audience to "readers with basic science literacy"
(c) Use analogies instead of technical terms

Which approach?

User: B -- readers with basic science literacy. But I'll still
define LNP and AAV briefly on first use.

Orchestrator: [Stage 2/4 - Essay Structuring] Quality Gate G2 check:
- Sections defined with word targets: YES (6 sections, 2000 words)
- Total word count agreed: YES (2,000)
- Audience specified: YES (general readers with basic science literacy)
- Each section has a purpose: YES
- User approved: YES
- Word targets sum correctly: 200+350+500+400+300+250 = 2,000

Outline approved. Saving to session directory.

Ready to move to Stage 3: Argument Development?
```

---

## Stage 3: Argument Development (Section 3 Excerpt)

```
[Stage 3/4 - Argument Development - Section 3 of 6] Developing argument map

Orchestrator: Let's develop Section 3: "The Delivery Problem." This
is your core argument section (~500 words). It needs to convince the
reader that delivery, not editing precision, is the real bottleneck.

I see three possible argument threads:

A) Historical: Show precision has improved dramatically but clinical
   outcomes haven't improved proportionally -- pointing to delivery.
B) Mechanistic: Explain what "delivery" means concretely and why
   each step is hard.
C) Economic: Show delivery R&D dominates CRISPR therapy budgets.

Which approach resonates? Or a combination?

User: A and B combined. The history shows the gap, then the
mechanism explains why.

Orchestrator: Good combination. Let me prepare claims for fact-
checking. I'm identifying these factual claims for verification:

1. "Off-target editing rates dropped from ~5% to <0.1% between
   2015-2024" (precision improvement data)
2. "Only a minority of CRISPR clinical trials met primary endpoints
   by 2024" (clinical outcomes data)
3. "Delivery requires reaching the right cells, in the right organ,
   at the right dose" (mechanistic claim)

Claim 3 is a conceptual/definitional claim, not a specific data
point. Claims 1 and 2 are factual and need verification.

Invoking fact-checker for claims 1 and 2...

[Fact-checker results received]

Results:
- Claim 1: PARTIAL. The review reports off-target rates dropped from
  4.8% to 0.08% (roughly 60-fold, not 50-fold). Rates vary by
  editor type. Source: Nature Medicine 2024 Review (doi:10.1038/...)
- Claim 2: PARTIAL. The source reports 4 of 17 trials met primary
  endpoints, not "only a minority of 15." Source: same review.

For Claim 1: The data supports your general point but the specific
numbers are slightly different. Want to use the actual figures
(4.8% to 0.08%) or round (roughly ~5% to <0.1%)?

User: Use the actual figures. More credible.

Orchestrator: Agreed. For Claim 2: The source says 4 of 17, not
3 of 15 as originally estimated. Want to update?

User: Yes, use 4 of 17.

Orchestrator: Updated. Now let me play devil's advocate:

A skeptical reader could say: "The trials that failed had many
problems beyond delivery -- patient selection, dosing, disease
complexity. You can't attribute all failures to delivery." How do
you want to address that?

User: Good point. I'll acknowledge that delivery isn't the only
factor but argue it was a contributing factor in most cases.

Orchestrator: That's a more defensible position. Here's the argument
map for Section 3:

Point 1: Precision has improved dramatically (4.8% to 0.08%)
  Source: Nature Medicine 2024 Review
  Verified: YES

Point 2: Clinical outcomes haven't kept pace (4 of 17 trials)
  Source: Same review
  Verified: YES

Point 3: The gap points to delivery as a major bottleneck
  Type: Analytical claim (not separately verifiable)

Point 4: Delivery means reaching right cells, right organ, right dose
  Type: Definitional/conceptual

Counterargument addressed: Other factors contribute to trial
failures; delivery is one major factor among several.

Approve this argument map for Section 3?

User: Approved.
```

---

## Stage 4: Paragraph Writing (Two Paragraph Excerpts)

```
[Stage 4/4 - Paragraph Writing - Section 3, Paragraph 1] Drafting

Orchestrator: This paragraph covers Points 1 and 2 from the
argument map: precision improvement data and clinical outcomes gap.

---

In 2015, the first wave of CRISPR experiments in human cells hit
their targets with roughly 95% accuracy -- impressive for a brand-
new tool, but that remaining 4.8% meant thousands of unintended
edits scattered across the genome. A decade of engineering later,
base editors and prime editors have pushed off-target rates below
0.1%, a sixty-fold improvement that by any measure counts as a
technical triumph. Yet the clinical scorecard tells a different
story: of the seventeen CRISPR-based trials that had reported
results by early 2024, only four met their primary endpoints. The
editing works. Getting it where it needs to go is another matter.

---

Sources:
- Off-target rates (4.8% to 0.08%): Nature Medicine 2024 Review
  (doi:10.1038/s41591-024-xxxxx)
- Clinical trial outcomes (4/17): Same review

Voice: Matched your pattern of "concrete number followed by
implication." The closing uses your characteristic short-punchy
sentence style. Paragraph is 5 sentences (within your 3-5 range,
at the upper end).

Flags: None. All content mapped to argument map points 1 and 2.

Approve / Edit / Reject / Rewrite with feedback?

User: The "sixty-fold improvement" is too precise given the
uncertainty. Say "roughly sixty-fold." Otherwise approved.

Orchestrator: Updated: "...a roughly sixty-fold improvement that
by any measure counts as a technical triumph."

Saving approved paragraph.

[Stage 4/4 - Paragraph Writing - Section 3, Paragraph 2] Drafting

This paragraph covers Points 3 and 4: delivery as bottleneck and
what delivery actually means.

---

The problem isn't the scissors -- it's the envelope. CRISPR-based
therapies face a delivery challenge that no amount of editing
precision can solve: the molecular machinery has to reach the
right cells, in the right organ, at a dose high enough to matter
but low enough to be safe. For a liver disease, that means getting
billions of editor molecules past the immune system, through the
bloodstream, across cell membranes, and into hepatocyte nuclei --
all without triggering a dangerous inflammatory response. Even the
most precise gene editor is useless if it never arrives.

---

Sources: No specific factual claims requiring citation (this
paragraph is conceptual/mechanistic).

Voice: Used your analogy pattern ("the scissors and the envelope")
-- one analogy per complex concept, drawn from everyday life. The
liver example makes the abstract concrete, matching your "specific
before general" pattern.

Flags: None. Content maps to argument points 3 and 4.

Approve / Edit / Reject / Rewrite with feedback?

User: Approved.
```

---

## Navigation: Show Full State

```
User: Show full state

Orchestrator:
[Stage 4/4 - Paragraph Writing - Section 3, Paragraph 3]

Pipeline State:
- Session: /tmp/essay-pipeline-20260214-143022/
- Pushback level: Full

Stage 1 (Thesis): COMPLETED
  Thesis: "The biggest obstacle to CRISPR-based medicine is not
  editing precision...but the unsolved problem of delivering gene
  editors to the right cells."

Stage 2 (Outline): COMPLETED
  6 sections, 2,000 words target
  Audience: General readers with basic science literacy

Stage 3 (Arguments): COMPLETED
  6 of 6 sections have approved argument maps
  12 verified claims, 0 deferred, 1 user override

Stage 4 (Paragraphs): IN PROGRESS
  Section 1: 2/2 paragraphs approved
  Section 2: 3/3 paragraphs approved
  Section 3: 2/4 paragraphs approved  <-- YOU ARE HERE
  Section 4: Not started
  Section 5: Not started
  Section 6: Not started

Voice checks completed: Sections 1, 2 (scores: 4/5, 5/5)
Word count so far: ~1,100 of ~2,000 target

Continuing with Section 3, Paragraph 3...
```

---

## Final Essay Excerpt

After all paragraphs are approved, the pipeline presents the completion summary:

```
[Stage 4/4 - Pipeline Completion] Final checks

Final voice check: 4/5 (strong match; minor note about slightly
more formal register in Section 5)

Final fact-check sweep: All 14 factual claims verified. 0 new
issues found.

Deferred claims: 0 (all resolved)

Pipeline Complete!

Word count: 2,047 words
Sections: 6
Sources verified: 14 claims with 9 unique sources
Voice consistency: 4/5
User overrides: 1 (logged in fact-check-log.md)
Deferred claims resolved: 0 of 0
Time spent: 3 hours 22 minutes

Essay saved to: /tmp/essay-pipeline-20260214-143022/final-essay.md

Would you like me to write the essay to a specific file path?

User: Yes, save it to essays/drafts/crispr-delivery.md

Orchestrator: Essay written to essays/drafts/crispr-delivery.md.
Session marked as completed.
```

---

## Key Observations From This Walkthrough

1. **Every stage was interactive**: The user was a participant at every step, not a reviewer of autonomous output.
2. **Pushback was stage-appropriate**: Logical rigor in Stage 1 (challenging the thesis), audience awareness in Stage 2 (challenging assumptions about readers), devil's advocacy in Stage 3 (raising counterarguments), minimal in Stage 4 (execution focus).
3. **Fact-checker was delegated**: The orchestrator never verified facts itself. It prepared claim batches and delegated to the fact-checker sub-agent.
4. **Voice notes accompanied every paragraph**: The paragraph writer explained its voice matching decisions.
5. **Argument map compliance was tracked**: Each paragraph explicitly noted which argument map points it covered.
6. **Navigation was available**: The user could check full state at any time.
7. **User edits were accepted without debate**: When the user changed "sixty-fold" to "roughly sixty-fold," it was accepted immediately.
