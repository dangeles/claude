---
name: genesis-design
last_updated: 2026-08-02
description: >
  Use when building or changing any Genesis surface — a white paper, the paper
  library, a KOL record, a deck, an explorable, the bioreactor screens, or
  genesis.bio — or when a design question needs the house answer: colour,
  contrast, type scale, layout, focus rings, chart palettes, alt text, voice.
  Triggers on CSS work in ~/repos/white_papers, ~/repos/papers or ~/repos/kol,
  on any tokens.css or Genesis stylesheet, on greenfield work for a Genesis
  surface ("build a landing page for", "new page on genesis.bio", "a deck
  for"), and on "what colour should", "is this accessible", "how many
  series", "which font size", "what is the focus ring". Pair it with
  frontend-design, which builds the interface while this supplies the house
  rules. NOT for non-Genesis web work, and NOT for editing the guidelines
  repository itself (that has its own CLAUDE.md).
categories:
  - design-system
  - web-management
handoff:
  accepts_from:
    - web-presence-manager
    - frontend-design
---

# Genesis design system

**Black, white, grey and photography carry the identity. Hue carries meaning.**

Everything else is downstream of that one rule. An agent that internalises only this
sentence gets most calls right.

The design system lives in **`~/repos/aesthetic_and_web_guidelines`**. This skill routes you
there. It deliberately does not restate the rules — a paraphrase here would be a second
source of truth that drifts from the first, which is the exact failure that repository
exists to prevent. **If you are answering a design question without having read the
document, stop and read it.**

## Preflight

If `~/repos/aesthetic_and_web_guidelines` is missing, say so and stop. Do not reconstruct
the rules from memory.

Several of the documents below — including `docs/rules.md`, `docs/INDEX.md` and
`tools/check_surface.py` — were **untracked** in that repository as of 2026-08-02. A fresh
clone or a `git stash` there will not have them. Check the working tree, not `git show`.

## Read first

Absolute paths are relative to `~/repos/aesthetic_and_web_guidelines`.

| situation | read |
|---|---|
| any Genesis surface | `docs/00-foundations.md`, then the matching guide in `surfaces/` |
| a specific question | `docs/INDEX.md` — routes by question, down to the section anchor |
| adopting the tokens | `docs/09-adoption.md` |
| what a rule is called | `docs/rules.md` |

`docs/rules.md` held **42 rules across 12 prefixes** on 2026-08-02 — COLOR, TYPE, LAYER,
LAYOUT, MOTION, A11Y, FIG, PROV, VOICE, IMG, INST, ADOPT. If your count differs, the
repository is right and this skill is stale: re-read `docs/INDEX.md` before relying on
anything here.

## The tools

Four stdlib-only Python 3 scripts in `tools/`. All take `--json`. Two observe other
repositories and none of them ever write to one. A whole-repo check takes about 125 ms —
run it rather than reasoning about what it would say.

| | |
|---|---|
| `palette_audit.py` | the colour maths. `--selftest` proves it against 26 Sharma reference pairs; `--emit-tokens` re-solves the palette |
| `check_tokens.py` | inward: do the artefacts, the maths and the prose agree? Twelve checks |
| `check_surface.py` | outward: does a surface obey the rules? Read-only, reports by rule ID |
| `audit_surfaces.py` | outward: how far apart is everything today? `--html` for the drift dashboard |

`check_surface.py` **exits 1 when it finds something** — that is a result, not a crash.
Branch on whether stdout parses as JSON. It also needs a **directory**: hand it a single
file and it reports zero findings, which is indistinguishable from clean. `/genesis-check`
wraps both traps.

Every finding carries a `guideline` anchor, so you can go from a finding to the paragraph
arguing for the rule in one read. Do that instead of restating the rule.

## Five rules an agent otherwise breaks confidently

- **Never hand-compute a contrast ratio.** Get it from `palette_audit.py`. Do not estimate,
  do not recall, do not reason about luminance.
- **Never write a number into prose the tools cannot reproduce.** In the guidelines
  repository this is enforced — an unanchored ratio fails the build.
- **Name the ground beside any ratio.** There are three, not two: `#ffffff` light reading,
  `#0a0a0b` dark *reading*, `#000000` *display*. "On black" is ambiguous between the last
  two by enough to flip a pass into a fail.
- **Cite rule IDs.** `COLOR-04`, not a paraphrase. A paraphrase drops the reason first, and
  the reason is what stops the next person overturning the rule.
- **Do not edit a surface repository from the guidelines repository, or the reverse,**
  without being asked in that repository for that change. Adoption is the surface owner's
  decision and the owner is David.

## Enough to answer without opening a file

Verified by `check_tokens.py` as of 2026-08-02. Where this and the repository disagree, the
repository is right.

- **Two token tiers.** `--gx-*` is the raw palette and component CSS never names it.
  Semantic tokens — `--ink`, `--accent`, `--data-1` — are the only ones a stylesheet may use.
- **One declaration per token.** Each scheme-varying token is a single `light-dark(l, d)`.
  Switching scheme is switching `color-scheme` and nothing else. There is no
  `prefers-color-scheme` block to add a token to.
- **Seven reading sizes** — micro, label, fine, body, lead, title, hero — plus two
  instrument sizes only the bioreactor guide may use.
- **No bold.** SN Leif ships one weight (OS/2 class 600) and its italic. Emphasis is size,
  case, tracking, ink and rules.
- **Four data series**, re-derived on every check rather than chosen. Past four:
  direct-label and mute, facet, or aggregate. Diverging scales are blue-to-orange.
- **One accent per surface, one meaning**, declared in a machine-readable comment. On a
  landing page the answer is `none` — the photography carries it.
- **Layer order:** `tokens, fonts, base, layout, components, surface, state, print`. Append
  domain layers after `components`; never reorder; name every layer in the statement.
- **Floors:** body prose 7:1, supporting text and accent 4.5:1, plotted lines and control
  boundaries 3:1, hit targets 44×44 and 60px on an instrument. WCAG 2.2 AA, AAA on body.
- **No AI-generated images, anywhere.** Standing instruction, not a preference.
- **The company is Genesis**, legally Organogenesis Labs, Inc. Never *Spacetime Bio*, *Non
  Fiction* or *Organogenesis Corporation*.

Six questions are still open in `research/open-questions.md`, including whether the margin
column earns its width — which most of the layout rests on and has never been tested with a
reader. Do not present those as settled.

## Delegation

Reading the whole design system into the main thread costs more than the task. For a broad
review — "does this site follow the guidelines" — run `/genesis-check` first and read the
findings; open a guideline only when a finding needs its reasoning. For a narrow question,
`docs/INDEX.md` gives a section anchor directly.

## Not for

- **Non-Genesis web work.** A generic React app is not a Genesis surface.
- **Editing `~/repos/aesthetic_and_web_guidelines` itself** — it has its own CLAUDE.md and a
  stricter contract, including a pre-commit checker you must not bypass.
- **Chart form.** Chart type, encoding and axis treatment stay with `plotting-advisor` for
  Python figures and `dataviz` for everything else. What comes from here is the palette,
  the four-series limit, FIG-01…FIG-04, `genesis.mplstyle` and `genesis_palette.py` — so
  where their palette guidance and ours disagree on a Genesis surface, ours wins.
- **Scientific content.** That is `essay-pipeline`, `latex-document-manager` or
  `papers-library`.
