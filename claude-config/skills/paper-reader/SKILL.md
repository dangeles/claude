---
name: paper-reader
version: 1.0
last_updated: 2026-07-30
description: >
  Use when reading, summarizing, or discussing one or a few scientific papers — what the
  authors found, whether it holds up, and what it means biologically. Triggers on "read this
  paper", "summarize this preprint", "what does this paper say", a pasted DOI/PMID/PDF path.
  NOT for sourcing quantitative parameters with citations (use researcher), orchestrated
  multi-stage literature reviews (use lit-pm), or judging an analysis method's biological
  validity (use biologist-commentator).
categories:
  - research
  - reading
prerequisites:
  - A paper — DOI, PMID, URL, or local PDF
success_criteria:
  - The reader learns what the paper means, not just what it measured
  - Claims are separated from evidence, and the gap between them is named
  - Full texts are read out-of-context; the main thread receives notes, not papers
---

# Paper Reader

Announce: "I'm using the paper-reader skill."

## The point

A paper is an argument, not a data dump. Your job is to understand what the authors are
claiming about how biology works, whether their evidence supports it, and what would change
if they're right. Someone who reads your notes should understand the biology — not receive
a table of numbers stripped of meaning.

**Do not default to extracting statistics.** Effect sizes, n's, p-values, and parameter
values matter when they bear on whether the claim holds, and are noise otherwise. A
p-value copied out of a paper you didn't understand is worse than nothing: it looks like
knowledge. If you find yourself building a table of every number in the results section,
stop — you're doing the wrong job. That job is `researcher`, and it's the right job only
when the user needs parameters for a model or calculation.

## Read the paper out of context

Full texts are 8,000–15,000 tokens each. Reading three of them into the main thread costs
more context than the entire conversation that follows, and most of those tokens are
methods boilerplate and reference lists you will never refer to again.

**Default: dispatch a subagent per paper.** Give it the paper and the reading brief below;
it returns notes. The main thread holds notes, never full texts.

Read directly in the main thread only when the paper is short (a two-page comment or
preprint abstract), when the user explicitly wants to work through the text together, or
when you're returning to a specific passage you already know you need.

When several papers are in play, dispatch them concurrently — one subagent each, in a single
message.

## The reading brief

For each paper, work out:

**What question were they actually asking?** Often narrower than the title suggests. The
title is marketing; the question is in the last paragraph of the introduction.

**What did they actually do?** Not the method name — the thing they physically did. "Knocked
out gene X in mouse liver and looked at what happened to bile acid transport" beats
"conditional knockout with LC-MS metabolomics." System, model organism, and conditions
matter enormously for what the result means: an effect in immortalized cells at 21% O₂ is a
different claim about biology than the same effect in primary tissue.

**What did they find, and what do they claim it means?** These are two different things and
the gap between them is where papers go wrong. Note both, separately.

**Does the evidence carry the claim?** The questions that usually matter:
- Is the control the right control, or a convenient one?
- How many independent replicates — biological, or technical dressed up as biological?
- Correlation being asked to do causation's work?
- Is the effect size biologically meaningful, or just statistically detectable? A 4% change
  with p < 0.001 and n = 10,000 may mean nothing at all.
- Does the abstract's confidence match the data in figure 4?
- What experiment would have settled this that they didn't run?

**What does it mean biologically?** The question everything else serves. What does this
imply about the mechanism? What does it predict? What does it rule out? If it's true, what
else has to be true, and does that match what the field believes? Say this in mechanistic
prose, not bullet fragments.

**What would change your mind?** Name the result that would overturn it. A paper you can't
imagine being wrong is a paper you haven't understood.

## Numbers, when they earn their place

Include a number when it carries the argument — the effect that motivates the conclusion,
the parameter someone would reuse, the discrepancy with prior work. Always with its
context: species, system, conditions, and what was actually measured versus inferred. A
number without its measurement context is not reusable, it's just a digit.

Skip the rest. If the user later needs a full parameter extraction, that's `researcher`, and
this note tells them whether the paper is worth that effort.

## What to return

Prose that a colleague could read and understand the paper from. Roughly:

- The claim, in one or two sentences.
- What they did, enough to judge it.
- What they found, and where finding and claim diverge.
- The weaknesses that actually matter — not a symmetric list of caveats. If the paper is
  solid, say it's solid and name the one thing you'd still want.
- What it means, and what it changes.

Length follows the paper's substance. A landmark paper deserves more than a modest one; a
paper that turns out to be thin deserves two sentences saying so. Don't pad to a template
and don't fill sections that have nothing in them.

## Several papers at once

When the user hands you a stack, the value is in what they say *together*: where they agree,
where they contradict each other, and what the disagreement is actually about — different
systems, different measurements, or a genuine dispute about mechanism. Read each one, then
write the synthesis. Don't concatenate summaries; that's a reading list, not an
understanding.

For a real literature review across many papers, hand off to `lit-pm`.

## Adjacent skills

- `researcher` — you need parameters with citations for a model or calculation.
- `lit-pm` / `literature-researcher` — a comprehensive review across many papers.
- `biologist-commentator` — judging whether an analysis approach is biologically sound.
- `synthesizer` — integrating findings across notes you've already written.
