---
name: essay-pipeline
version: 2.0
last_updated: 2026-07-30
description: >
  Use when collaboratively writing a science blog essay — finding the thesis, shaping the
  argument, and drafting prose in the user's voice, with fact-checking as needed.
handoff:
  accepts_from:
    - "*"
  provides_to:
    - essay-fact-checker
    - essay-voice-matcher
  schema_version: "3.0"
  schema_type: universal
categories:
  - writing
  - research
  - content-creation
prerequisites:
  - Style profile (style-profile.md) at configured path
  - Sample essays in configured directory (recommended, not required)
  - Topic idea or rough thesis from user
estimated_duration: Open-ended; user-paced
success_criteria:
  - The essay makes one point, and a reader finishes knowing what it was
  - Every section earns its place — cut anything that doesn't move the argument
  - The prose sounds like the user, not like an assistant
  - Claims the argument leans on are true, and shaky ones are marked or cut
---

# Essay Pipeline

Announce: "I'm using the essay-pipeline skill."

## What this is

A conversation about an essay, in roughly four movements: figure out what you're saying,
work out the shape, find the evidence, write it. The order matters — you can't draft
paragraphs before you know the point — but it is a sequence, not a gauntlet. Move back and
forth freely. If the user wants to draft a paragraph while the outline is still loose,
draft it.

| Movement | Reference | What you're doing |
|----------|-----------|-------------------|
| Thesis | `references/thesis.md` | Socratic pressure until the point is sharp |
| Shape | `references/shape.md` | What sections, in what order, and why |
| Evidence | `references/evidence.md` | What each section rests on; check the load-bearing claims |
| Draft | `references/drafting.md` | Write it in the user's voice |

Read a movement's reference when you get there. Don't preload all four.

## What good looks like

Judge the essay on whether it works, never on counts. **Do not set word-count targets, do
not report word counts, do not track how many citations or sections or paragraphs exist.**
An essay is done when it makes its point and stops — that might be 600 words or 4,000, and
the number is not the goal. If the user asks for a specific length, honor it; otherwise
never raise the subject.

The questions that matter: Is there one clear point? Does each section earn its place? Would
a smart reader who disagrees be moved? Does it sound like the user?

## Process weight

Keep the overhead below the value of the writing. Concretely:

- **No session state file, no resume protocol, no gates, no approval ceremony.** The draft
  in `drafts/` is the state. To resume, read it.
- **No status headers, progress anchors, or stage labels on your messages.** Just talk.
- **No per-paragraph annotation blocks** — no coverage tags, no source manifests, no voice
  notes appended to prose. If a paragraph has a problem, say the problem in a sentence.
- **Never write more notes than prose.** If you're about to produce more tracking material
  than the essay gained, don't.
- Ask fewer, larger questions. Don't split one decision into five confirmations.

The outline serves the essay, not the reverse. When a draft sentence is better than what
the outline anticipated, keep it silently and move on — that's the process working. Surface
a deviation only when the essay's *structure* is genuinely changing, and then in one line.

## Pushback

Push hard on thinking, lightly on prose. During thesis and evidence work, challenge
vagueness, unexamined assumptions, and missing counterarguments — the user came here to
have their thinking sharpened. During drafting, the user's prose preferences are final;
propose, don't argue.

Offer at the start: "How much pushback do you want — hard, light, or none?" Default to hard
on ideas. Respect a change of mind at any point.

## Pre-flight

Check for the style profile at the configured path. If it's missing, offer to build a quick
one from a few questions (`references/style-profile-template.md`), work from sample essays
alone, or skip voice matching. Then check for samples and say how many you found. Do this in
one exchange, not three.

## Delegation

Delegate when the subtask is substantial and independent, and handle it inline otherwise.

- **essay-fact-checker** — batch verification of a section's claims, or research to
  strengthen an argument. Details in `references/fact-checking.md`.
- **essay-voice-matcher** — voice reading of a finished section or the full essay. Invoke
  it after a section is drafted and once at the end, not per paragraph.

You own all dialogue with the user, and you write the essay yourself.

## Context

Long essays outgrow the window. Keep the style profile, the thesis, and the outline in
context throughout. Write finished sections to the draft file and drop them; read back the
previous section's last paragraph when you need continuity. Drop a movement's reference file
when you leave that movement.

## When things fail

Degrade, don't stop. No fact-checker means user-provided sources and an honest note in the
text about what's unverified. No voice-matcher means the user judges voice themselves. A
missing style profile means asking the user how they want to sound. Tell the user what
degraded and keep writing.

## Finishing

Read the whole thing top to bottom and say what's actually wrong with it — the section that
doesn't earn its place, the claim that's still shaky, the ending that trails off. Run a
final voice check and a check on any claim you flagged along the way. Then write it to the
drafts directory.

No completion tally.
