# Claude Code config repo

Source of truth for global `~/.claude/` configuration: 54 skills, custom agents, hooks, plugin enable flags. Synced via `sync-config.py`.

## Core invariant

**Modify `claude-config/`, never `~/.claude/` directly.** The repo drives the live system, not the other way around. Manual edits to `~/.claude/` get clobbered on the next `push`.

## Standard workflow (7 steps)

```
1. git status                                          # clean working tree
2. ./sync-config.py plan --title "..."                 # create planning entry
3. edit claude-config/...                              # repo first
4. jq . claude-config/settings.json                    # validate
5. ./sync-config.py push --dry-run                     # preview
6. ./sync-config.py push                               # apply to ~/.claude/
7. git add claude-config/ planning/ && git commit      # only after testing
```

Full workflow with rollback procedures: see `CONFIG_MANAGEMENT.md`.

## Push-time gates (automatic)

`sync-config.py push` blocks on:
- Invalid frontmatter in any `skills/*/SKILL.md` (missing `name:` or `description:`)
- JSON syntax errors in `settings.json`

If push refuses, fix the source, don't bypass.

## Where things live

| Path | Purpose |
|---|---|
| `claude-config/skills/` | All skills (synced bidirectionally with `~/.claude/skills/`) |
| `claude-config/agents/` | Custom agents (synced with `~/.claude/agents/`) |
| `claude-config/settings.json` | Plugin enable/disable flags |
| `claude-config/plugins/installed_plugins.json` | Plugin version pinning |
| `claude-config/skills/CHANGELOG.md` | Skill-level change history |
| `planning/$(hostname -s)/YYYY-MM-DD-title.md` | Per-change planning entry, linked to commit SHA |
| `claude-config/project-configs/` | Per-project overrides (currently empty) |
| `.claude/settings.local.json` | Project-local permissions; **gitignored**, not synced |
| `docs/` | Long-form guides (CLAUDE_CONFIG_GUIDE, SYNC_WORKFLOW) |

## Planning entries

- Use `hostname -s` (short form) for the directory, not `hostname` (FQDN). Canonical here: `planning/mac/`.
- Template: `planning/.template.md`.
- After commit, append the SHA to the planning entry's "Related Commits" section.

## What is NOT in scope here

- `~/.claude/settings.local.json` — user-wide local permissions, machine-specific, never synced.
- `~/.claude/plans/`, `history.jsonl`, `cache/`, `*-workspace/` — ephemeral, never synced. See `sync.config.yaml` for the full exclusion list.

## Skills with strong opinions

A few skills explicitly override defaults — respect their stated precedence:
- `skill-editor` says "always prefer this over `skill-creator` / `plugin-dev:skill-development`" for any skill in `claude-config/`.
- `superpowers:using-superpowers` requires invoking applicable skills before responding.

## Anti-patterns

- Editing `~/.claude/` and expecting it to survive (it won't — next push overwrites it).
- Skipping the planning entry. Untracked config drift is the single biggest pain point in this repo's history.
- Using `--delete` on push without a `--dry-run` first.
- Committing before `sync push && test in a real session`.
