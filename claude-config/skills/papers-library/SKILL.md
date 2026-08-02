---
name: papers-library
version: 1.1
last_updated: 2026-07-31
description: >
  Use when adding a paper to the library, citing a paper while writing, or scanning a
  manuscript for citations. Triggers on "add this paper to the library", "summarize this
  paper for the library", "cite X here", "scan this grant for citations", a DOI or PDF
  offered for ingestion, or any work inside ~/repos/papers. NOT for reading a paper to
  understand it without recording it (use paper-reader), sourcing parameters with citations
  (use researcher), or multi-stage literature reviews (use lit-pm).
categories:
  - research
  - reading
  - writing
prerequisites:
  - The library at ~/repos/papers
success_criteria:
  - Summaries follow RECIPE.md and pass tools/validate.py with zero errors
  - hmp and flags are computed by tools/hmp.py, never written by hand
  - Cited For entries are appended, never deleted
---

# Papers Library

Announce: "I'm using the papers-library skill."

The library lives at `~/repos/papers`. **Read `~/repos/papers/CLAUDE.md` before acting** —
it holds the full protocol. `~/repos/papers/RECIPE.md` holds the summary format. This skill
exists so those are reachable from any repo; it does not duplicate them.

There is no bare `python` on this machine. Run the tools through the repo venv, which works
from any directory:

```bash
~/repos/papers/.venv/bin/python ~/repos/papers/tools/validate.py
```

## Which flow are you in?

**The user handed you a paper** ("add this", "summarize this for the library", a bare DOI or
PDF path) → ingest. Resolve the DOI, fetch metadata from OpenAlex rather than transcribing
it, work down the fetch ladder in `CLAUDE.md` until you have the best artifact available,
write the sections per the recipe, run `tools/hmp.py`, then `tools/validate.py`. If what you
read was not a PDF held in `pdfs/`, say which artifact it was under `## Full Text Source`.

**The user wants an entry checked** ("does this still hold", "check these against their
PDFs") → verify. Re-run the fetch ladder, read the new artifact, and correct the entry
against it — figure panels first, since a rendering hides the p-values printed on them. Every
factual correction this library has made came out of that step.

**The user is writing and wants a citation** ("cite Igarashi for the OSM point") → this is
the flow that makes the library worth having. Insert the citation into the manuscript in its
native format *and* append the `Cited For` entry to the library file in the same step, with
`cited_by_count` and `cited_by_sources` updated to match — the validator errors when they
disagree. Don't defer it. If the paper isn't in the library, create a stub and say so.

**The user points at a manuscript** ("scan my grant") → `tools/scan.py <path>`. Idempotent;
safe to re-run.

**The user wants an entry changed** ("apply the requests on X", "fix the stats block") → edit
in place, bump `updated:`, leave `## Cited For` untouched, and re-run `hmp.py` if the stats
block moved.

## The rules that matter most

Never hand-write `hmp` or `flags` — run `tools/hmp.py <file>` on the named file. A guessed
harmonic mean p-value looks like knowledge and isn't, and `validate.py` will catch it anyway.

Never delete anything under `## Cited For`. It is append-and-update only.

Never crawl reference lists. The user curates this library. Reading a paper may surface at
most two papers its own authors flag as essential; record them in `suggested:` and create
nothing.

Never claim to have read a full text you couldn't get — that's what `status: abstract-only`
is for.

Finish by running `.venv/bin/python tools/validate.py` in the repo. Zero errors, or the work
isn't done.

## Finding things

```bash
grep -rl "oncostatin" ~/repos/papers/papers/
grep -rn "equivalence_from_nonsignificance" ~/repos/papers/papers/
```

## Adjacent skills

- `paper-reader` — understanding a paper without recording it in the library.
- `researcher` — extracting parameters with citations for a model or calculation.
- `lit-pm` — a comprehensive review across many papers.
