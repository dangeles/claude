# Output Formats

Templates for the two structured outputs technical-pm produces during a coordination session: the work plan (task breakdown) and the progress dashboard. Both are Markdown, both are designed to fit in a single readable scroll.

## Task Breakdown Format

Produce this at the start of a new initiative, after requirements are clear but before agents are dispatched.

```markdown
# Work Plan: [Initiative Name]

**Goal**: [What we're trying to accomplish]
**Timeline**: [Target completion or sprint]
**Status**: [Not Started / In Progress / Blocked / Complete]

## Dependencies
```
[Task A] ──► [Task B] ──► [Task C]
                    └──► [Task D]
```

## Tasks

### 1. [Task Name]
- **Assigned to**: [Agent]
- **Status**: [Pending / In Progress / Complete / Blocked]
- **Depends on**: [Other task IDs, if any]
- **Deliverable**: [What this task produces]
- **Notes**: [Any context]

### 2. ...

## Blockers
| Blocker | Blocking | Resolution Path | Owner |
|---------|----------|-----------------|-------|
| ... | [Task ID] | ... | ... |

## Handoffs Required
- [ ] [Task X] → [Agent Y] when [condition]
- [ ] ...

## Decisions Needed from User
1. [Decision with context]
```

## Tracking Dashboard Format

Update this at every progress checkpoint (per the cadence in the KPI Tracking section of `SKILL.md`).

```markdown
# Progress Dashboard

**As of**: [Date/Time]
**Sprint/Phase**: [Name]

## Status Summary
| Track | Progress | Status | Next Milestone |
|-------|----------|--------|----------------|
| [Track 1] | ███░░ 60% | On Track | [Milestone] |
| [Track 2] | ██░░░ 40% | Blocked | [Blocker] |

## Recent Completions
- [Date]: [What was completed]

## Current Focus
- [Agent]: Working on [Task]

## Upcoming Handoffs
- [Agent A] → [Agent B]: [Document/deliverable]

## Risks
- [Risk and mitigation]
```

## Notes

- Task IDs should be stable across the lifetime of a work plan; do not renumber when reorganizing.
- The dependencies ASCII block is rendered as a code block so the indentation survives copy/paste.
- Progress bars use Unicode block characters (`█`, `░`); they degrade gracefully in plain-text viewers.
- For parallel execution patterns and per-task dispatch templates, see `task-templates.md`.
