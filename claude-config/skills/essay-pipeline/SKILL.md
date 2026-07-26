---
name: essay-pipeline
version: 1.0
last_updated: 2026-05-24
description: >
  Use when collaboratively writing a science blog essay through interactive
  thesis development, essay structuring, argument development, and paragraph
  writing with tiered fact-checking and voice matching.
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
  - WebSearch tool available (for fact-checking)
  - Topic idea or rough thesis from user
estimated_duration: 2-8 hours (depending on essay length and user pace)
success_criteria:
  - User completes all 4 stages from topic idea to finished essay
  - Every stage involves interactive dialogue with user approval
  - Stage-appropriate pushback applied at configured level
  - Final essay matches user voice (voice-matcher confirms)
  - All factual claims verified with source URLs
  - Session state preserved across pauses and resumes
---

# Essay Pipeline

Announce: "I'm using the essay-pipeline skill for interactive essay writing."

## Architecture

Orchestrator-as-conductor. You run all four interactive stages in your own thread, in direct conversation with the user; only non-interactive specialists are delegated as sub-agents.

```
User
  |
  v
essay-pipeline orchestrator (YOU -- main thread)
  |
  |-- Stage 1: Thesis Development (interactive)
  |-- Stage 2: Essay Structuring (interactive)
  |-- Stage 3: Argument Development (interactive, per-section)
  |     Support: essay-fact-checker (via Task tool)
  |-- Stage 4: Paragraph Writing (interactive, per-paragraph)
  |     Support: essay-fact-checker, essay-voice-matcher (via Task tool)
  v
Final Essay Output
```

**Why this architecture**: Claude Code sub-agents cannot use AskUserQuestion. Since every stage requires interactive dialogue, the orchestrator runs all stages directly. Sub-agents handle only non-interactive work (fact verification and voice evaluation).

| Stage | Reference file | Pushback type | Support agents | Gate |
|-------|---------------|---------------|----------------|------|
| 1. Thesis Development | `references/stage-1-thesis-development.md` | Logical rigor | none | G1 |
| 2. Essay Structuring | `references/stage-2-essay-structuring.md` | Audience awareness | none | G2 |
| 3. Argument Development | `references/stage-3-argument-development.md` | Devil's advocacy | essay-fact-checker | G3, G3-Final |
| 4. Paragraph Writing | `references/stage-4-paragraph-writing.md` | Minimal | essay-fact-checker, essay-voice-matcher | G4, G4-Section, G5 |

Fact-check tiers, deferred-verification queue, source-quality hierarchy, and the user-override protocol live in `references/fact-check-tiers.md`. Session state schema and change-impact assessment live in `references/session-state-schema.md`.

## Delegation

Delegate specialist work when the subtask is substantial and independent — batch verification of a section's factual claims (essay-fact-checker), voice evaluation of a completed section or the full essay (essay-voice-matcher), proactive research enrichment for an argument. Handle it directly when you could finish it in a handful of tool calls, when the work is sequential, or when you need the context in your own loop. See `../references/delegation-and-scope.md`.

You own everything interactive and stateful: all dialogue via AskUserQuestion, session create/save/resume, stage transitions and approvals, context assembly, and final essay assembly.

## State anchoring

Start every response with a state anchor: `[Stage N/4 - {stage_name}] {status}`, for example `[Stage 3/4 - Argument Development - Section 2 of 5] Developing argument map`. Re-anchor after navigation commands, stage transitions, and resumes.

## Pre-flight

Check for the style profile:

```bash
test -f "{configured_path}/style-profile.md" && echo "FOUND" || echo "NOT_FOUND"
```

If found, read it and keep it in context for the whole session. If not, offer three paths via AskUserQuestion: (a) a quick 5-question express profile (questions in `references/style-profile-template.md`), synthesized, saved, and used; (b) degraded mode, where the voice-matcher works from sample essays only, or voice matching is skipped if there are none; (c) exit so the user can fill out the full template.

Then check for samples with `ls "{configured_path}/samples/"*.md 2>/dev/null | head -5` and tell the user how many the voice-matcher will have. WebSearch reaches the pipeline only through essay-fact-checker; if that invocation fails, fall back to user-provided sources.

## Session initialization

Scan for existing sessions with `ls -d /tmp/essay-pipeline-*/ 2>/dev/null`. If any exist, present them with status summaries and ask whether to resume or start fresh.

For a new session, create `/tmp/essay-pipeline-$(date +%Y%m%d-%H%M%S)`, initialize state per `references/session-state-schema.md`, then offer a pushback level via AskUserQuestion:

```
Before we begin, how much intellectual pushback would you like?

- Full (default): I'll challenge your ideas vigorously to strengthen them.
  3 rounds of pushback in thesis development, 2 in argument development.
- Light: I'll raise concerns once and accept your response.
  1 round per point.
- Minimal: I'll offer observations but won't challenge.
  No pushback rounds.

You can change this at any time during the session.
```

Save the initial state with the atomic write protocol, then begin Stage 1.

## Running a stage

Read the stage's reference file, anchor state, and follow its instructions, applying pushback at the configured level. Save incremental progress to session files as you go — after each completed argument point (Stage 3), each approved paragraph (Stage 4), each stage transition, and each navigation command.

At the end of a stage, present a completion summary and get the user's approval at the gate. Then drop that stage's reference file from active consideration, keep its output files (thesis, outline, argument maps) as context, and load the next one.

Stage 3 loops over sections: develop the argument map, batch-invoke the fact-checker for that section's claims, work through the results with the user, get approval, save. After the last section, present a bird's-eye review of all argument maps (G3-Final). The user may ask to jump forward and write paragraphs for sections already complete.

Stage 4 loops over paragraphs within sections: draft against the argument map and style profile, present for approval, save. Invoke the voice-matcher once a section is complete, then present the section for bird's-eye review. When a draft paragraph contains material the argument map doesn't cover, say so plainly — "This paragraph includes content not in your argument map: [description]" — and offer to keep it and extend the map, remove it, or go back to Stage 3.

### Gates

| Gate | Where | User approves |
|------|-------|---------------|
| G1 | End of Stage 1 | Thesis is clear and defensible, claim type identified |
| G2 | End of Stage 2 | Outline complete, length agreed, sections purposeful |
| G3 | End of each Stage 3 section | Claims evidenced, counterarguments addressed, sources verified |
| G3-Final | After all Stage 3 sections | Argument maps cohere |
| G4 | End of each Stage 4 paragraph | Voice matches, facts verified |
| G4-Section | After a section's paragraphs | Section reads well, voice consistent |
| G5 | End of Stage 4 | Full essay, all facts verified, all deferred claims resolved, all overrides confirmed |

A gate that doesn't pass means iterating on that stage, section, or paragraph — not proceeding.

## User navigation

| Command | Action |
|---------|--------|
| "Show full state" / "Where am I?" | Display current stage, section, paragraph; completed work; pending work |
| "Show essay so far" | Read and display all approved paragraphs from `stage-4-draft.md` |
| "Go back to thesis" / "Stage 1" | Assess impact, archive downstream if needed, return to Stage 1 |
| "Go back to outline" / "Stage 2" | Assess impact, archive downstream if needed, return to Stage 2 |
| "Go back to section N arguments" | Return to Stage 3 for section N |
| "Pause" / "Save and exit" | Save state, display resume instructions |
| "Resume" | Load state, present summary, continue |
| "Change pushback level" | Update `pushback_level` in session state |

Going back cascades. Assess the change impact (see `references/session-state-schema.md`): minor flags downstream work for review, moderate archives the affected artifacts, major archives all downstream artifacts. Tell the user what will happen — "Going back to [stage] will [impact]. Your current work will be [preserved/archived/flagged]. Continue?" — and on confirmation version the old artifacts and restart from there. **Never silently discard work.**

## Fact-checker invocation

Tiers, deferred queue, claim classification, source hierarchy, and override handling: `references/fact-check-tiers.md`.

In Stage 3, batch all of a section's factual claims into one invocation after the argument map is developed, then update the map with verification status and sources. In Stage 4, send only claims that are new — those already verified in Stage 3 don't need rechecking. After all paragraphs are approved and before G5, run one final sweep over the complete essay to catch claims that shifted during editing and to resolve anything deferred.

```
Task: essay-fact-checker

Verify the following claims from Section [N] of a science blog essay about [topic].

[YAML claim batch]

Return results as YAML with verification status, source URLs, and notes for each claim.
```

Track unverified claims in session state under `deferred_verifications`; all of them are resolved before G5, by retrying, by a user-supplied source, by revising the claim, or by removing it.

## Voice-matcher invocation

Invoke after each completed section, after the full essay, and whenever the user asks about voice.

On the first result, present the assessment and ask "How accurate is this voice assessment? Rate 1-5, where 5 is very accurate." Save the rating under `voice_calibration` and include it in later invocations.

| Score | Action |
|-------|--------|
| 5/5 | Proceed; tell the user it's a strong match |
| 4/5 | Proceed; mention specific observations |
| 3/5 | Present the assessment; ask whether they want revisions |
| 2/5 | Flag for revision with specific suggestions |
| 1/5 | Pause and suggest reviewing the style profile |

```
Task: essay-voice-matcher

Evaluate the following text for voice consistency against the user's style profile.

Text to evaluate:
[text]

Style profile path: [path]
Sample essays path: [path]
Context: Stage [N], Section [N], [content type]
Previous assessments: [YAML list of prior assessments]
Calibration: [calibration data if available]
Voice reference mode: [profile | sample | hybrid]

Return assessment as YAML with score (1-5), observations, and suggestions.
```

## Error handling

| Failure | Detection | Action |
|---------|-----------|--------|
| Fact-checker sub-agent fails | Task tool returns error | Retry once; then mark claims DEFERRED, inform user, proceed |
| Voice-matcher sub-agent fails | Task tool returns error | Skip voice check, note in metadata, inform user, proceed |
| WebSearch unavailable | Fact-checker reports service failure | User-provided sources only; inform user |
| Session state corrupt | YAML parse fails on Read | Recover from `.bak`; else reconstruct from artifacts; else ask user |
| Style profile missing | Pre-flight finds no file | Offer the three paths above |
| Style profile malformed | Read returns unparseable content | Warn user; offer best-effort interpretation or exit to fix |
| literature-researcher unavailable | Pre-flight or runtime | Inform user; Tier 3 escalation unavailable, expanded Tier 2 as substitute |
| Stage reference file missing | Read tool returns error | Fatal — broken installation. Report and stop |

Degrade rather than stop: losing the fact-checker means user-provided sources and DEFERRED claims; losing the voice-matcher means the user self-evaluates voice; a corrupt session state means reconstruction from artifacts, and only if the artifacts are gone too does recovery become a question for the user.

## Session state and context

Session state updates are atomic: write to a `.tmp` file, copy the current state to `.bak`, then rename `.tmp` over the primary. On load, try primary state, then backup, then artifact reconstruction (`references/session-state-schema.md`).

Long essays (5,000+ words) put pressure on the context window, so rely on files for historical content rather than conversation memory. The style profile, thesis, and outline stay in context throughout; a stage's reference file and a section's argument map and paragraphs are dropped when that stage or section finishes. Once a section is complete in Stage 4, write its full text to `stage-4-draft.md` and drop it; when continuity requires it, read back only the previous section's last paragraph.

## Completion

At G5, after the final voice check, the final fact-check sweep, and resolution of every deferred verification, present the summary:

```
Pipeline Complete!

Word count: [N] words
Sections: [N]
Sources verified: [N] claims with [M] unique sources
Voice consistency: [score]/5
User overrides: [N] (see fact-check-log.md)
Deferred claims resolved: [N] of [N]
```

Write the essay to `{session_dir}/final-essay.md`, offer to export it to a path of the user's choosing, and mark the session "completed" in state.
