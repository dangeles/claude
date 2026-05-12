# add-frontmatter-colon-space-hint-to-validator

**Date**: 2026-05-12
**Machine**: mac
**Status**: Success

## Objective

Improve the diagnostic from `sync-config.py validate_frontmatter` when an unquoted `: ` (colon followed by space) appears inside a YAML scalar value. The existing validator already catches the error via `yaml.YAMLError` ("mapping values are not allowed here"), but the message is cryptic and doesn't point at the fix. We hit this in Tier 3 when rewriting two descriptions: `research-pipeline` had `... fixed 5 stages: researcher -> ...` and `calculator` had `... feasibility checks are needed: order-of-magnitude ...`. PyYAML interpreted `stages:` and `needed:` as nested mapping keys.

This is exactly the "lesson learned" entry from the Tier 3 planning doc.

## Changes Planned

- [ ] Add `_detect_unquoted_colon_in_value(yaml_block)` helper to `ConfigSync` that scans frontmatter lines for the pattern: unquoted top-level scalar containing `: ` after the key/value separator.
- [ ] Modify `_parse_frontmatter` so when `yaml.YAMLError` is raised, the helper runs and any hint is appended to the error message.
- [ ] Smoke test: create a temp file with `: ` inside an unquoted description, run `./sync-config.py validate`, verify the message includes the hint. Then revert.
- [ ] Plan entry → modify → smoke test → commit (no sync push needed; sync-config.py is at repo root, not synced).

## Expected Outcome

Before:
```
* claude-config/skills/foo/SKILL.md: invalid YAML in frontmatter: mapping values are not allowed here
  in "<unicode string>", line 3, column 89:
     ... description here: with embedded colon ...
                                          ^
```

After:
```
* claude-config/skills/foo/SKILL.md: invalid YAML in frontmatter: mapping values are not allowed here ... (hint: line 3 'description' value contains unquoted ': ' — wrap in double quotes)
```

The cryptic pyyaml message stays (it's authoritative), but the appended hint tells the user exactly what to fix.

## Actual Outcome

Helper `_hint_unquoted_colon_in_value` added as static method on `ConfigSync`. `_parse_frontmatter` calls it from the `yaml.YAMLError` branch and appends `(hint: ...)` to the error message when the pattern matches.

Smoke-tested with 8 edge cases — all behave correctly:
- detects: unquoted scalar with `: ` mid-value
- ignores: double-quoted, single-quoted, block scalars (`|`, `>`), comment lines, list items, lines without `: `, and bare colons (URLs like `https://example.com`)

Regression check: `./sync-config.py validate` passes against current 56 skills + agents. No false positives introduced.

## Assessment

**Result**: Success

## Related Commits

- [pending]: add frontmatter colon-space hint

## Next Steps

If misdiagnoses surface for other YAML pitfalls, extend the helper. For now, this single pattern is the only documented one in our history.
