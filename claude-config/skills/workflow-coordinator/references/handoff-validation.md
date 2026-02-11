# Handoff Validation Reference

**Version**: 1.0
**Schema**: Universal Handoff Schema v3.0
**Date**: 2026-02-07

---

## Overview

The handoff validation system ensures handoff integrity through a three-stage pipeline:

1. **Schema validation** -- Enforce JSON Schema compliance
2. **Token counting** -- Enforce 2000 token limit (research-recommended)
3. **Semantic validation** -- Check references, consistency, circular handoffs

**Key Principle**: Fail fast with clear error messages, provide repair suggestions.

---

## Pilot Scope Limitations

This is a documentation-only pilot. The following are **guidance only** (not enforced automatically):

- **Token counting**: Uses character/4 heuristic estimation, not tiktoken
- **Schema validation**: Requires manual invocation (python3 one-liner)
- **Circular handoff detection**: Advisory, not automated
- **File reference validation**: Not checked at runtime
- **Repair prompts**: Manual application (no CLI tool)

Automated enforcement planned for v1.1 (Python library + CLI tools).

---

## Validation Pipeline

```
Handoff Created
    |
    v
+---------------------+
| Schema Validation   |  <-- JSON Schema validator
+---------------------+
    | FAIL --> Repair Prompt --> Retry (max 3)
    | PASS
    v
+---------------------+
| Token Counting      |  <-- Count tokens in payload
+---------------------+
    | FAIL (>2000) --> Compression --> Retry
    | PASS
    v
+---------------------+
| Semantic Validation |  <-- Check references, consistency
+---------------------+
    | FAIL --> Warnings (non-blocking)
    | PASS
    v
VALID HANDOFF
```

---

## 1. Schema Validation

### Schema Validation (Pilot Approach)

Validate handoff examples against the v3.0 schema:

```bash
python3 -c "
import json, jsonschema
schema = json.load(open('references/universal-handoff-schema-v3.0.json'))
example = json.load(open('examples/skill-editor-to-programming-pm.json'))
jsonschema.validate(example, schema)
print('PASS')
"
```

### Validation Checklist

- [ ] JSON parses without errors
- [ ] Passes JSON Schema validation against v3.0
- [ ] `meta.token_count` within 30% of character/4 estimate
- [ ] `trace_id` is valid UUID v4
- [ ] `timestamp` is valid ISO8601 UTC (ends with Z)
- [ ] `handoff_chain` does not contain target skill name

### Common Validation Errors

| Error | Cause | Repair |
|-------|-------|--------|
| Missing required field | Field not present | Add field with default value |
| Invalid type | Wrong data type | Convert to correct type |
| Exceeds maxLength | String too long (summary >500 chars) | Truncate and add "..." |
| Invalid format | Wrong format (e.g., date not ISO8601) | Convert to ISO8601 UTC with Z suffix |
| Invalid enum value | Value not in allowed list | Use closest match or default |
| `unevaluatedProperties` error | Extra property at root or handoff level | Remove unknown property or check spelling |

---

## 2. Token Budget

### Token Budget Breakdown

```
Total budget: 2000 tokens
+-- Metadata (version, timestamps, trace_id):   ~100 tokens
+-- Source/Target info:                          ~50 tokens
+-- Context (summary, problem_type, etc.):       ~200 tokens
+-- Payload.working:                             ~500 tokens
+-- Payload.session:                             ~800 tokens
+-- Payload.references:                          ~200 tokens
+-- Meta & tracing:                              ~150 tokens
```

### Token Estimation (Pilot Approach)

Without tiktoken, estimate tokens using character count:

**Formula**: `token_count = len(json_string) / 4`

This overestimates by ~20% (conservative, safe for budget enforcement).

Quick estimation command:

```bash
wc -c < handoff.json | awk '{printf "Estimated tokens: %d\n", $1/4}'
```

When populating `meta.token_count` in handoff payloads, use this heuristic and note the method. Actual tiktoken counting replaces this in v1.1.

### Token Limit Enforcement

1. After creating handoff JSON, compute `len(json_string) / 4`
2. If estimate exceeds 2000, apply compression strategies (see below)
3. Update `meta.token_count` with the estimate
4. If still over limit after all compression strategies, escalate to user

---

## 3. Compression Strategies

When a handoff payload exceeds the 2000 token budget, apply these strategies in order:

### Strategy 1: Remove Optional Fields

Remove optional fields to reduce token count:

- `handoff.expires_at`
- `handoff.target.invocation`
- `handoff.target.expected_phase`
- `handoff.context.original_prompt`
- `handoff.context.success_criteria`
- `handoff.payload.session` (keep only working + references)
- `handoff.tracing.tags`

### Strategy 2: Move Session to References

Move large session data to file references:

1. Write `payload.session` to a file at `{session_path}/handoffs/session-context.json`
2. Replace with a file reference in `payload.references.files`
3. Recalculate token count

### Strategy 3: Summarize Long Text

Summarize long text fields:

- `context.summary`: Condense to under 100 tokens
- `context.original_prompt`: Condense to key phrases
- `payload.working.*`: Summarize any string field over 500 characters

---

## 4. Repair Prompts

### Auto-Repair Defaults

When a required field is missing, use these defaults:

| Field | Default Value |
|-------|---------------|
| `version` | `"3.0"` |
| `schema_type` | `"universal"` |
| `trace_id` | New UUID v4 |
| `timestamp` | Current UTC time (ISO8601) |
| `confidence` | `"medium"` |
| `token_count` | `0` (recalculate after repair) |

### LLM-Assisted Repair Prompt

If auto-repair fails, prompt the LLM with:

```
The following handoff failed validation:

{handoff JSON}

Validation errors:
- {error 1}
- {error 2}

Please repair the handoff to comply with universal handoff schema v3.0.
Return only the repaired JSON (no explanation).
```

### Repair Retry Policy

**Maximum retries: 3**

1. **First repair attempt**: Auto-fix using default values
2. **Second repair attempt**: Apply compression strategy
3. **Third repair attempt**: Simplified payload (working context only, remove session and references)
4. **After 3 failures**: Escalate to user with detailed error report

---

## 5. Semantic Validation

### Circular Handoff Prevention

Rules:

1. A workflow MUST NOT hand off to itself (self-loop)
2. If the target workflow appears in `handoff_chain`, the handoff MUST include `meta.handoff_reason` explaining why re-entry is intentional
3. If `handoff_chain` length exceeds 5, require `meta.user_approved: true`
4. Maximum recommended chain length: 10 (truncate middle entries if exceeded)

### Reference Resolution

Check that all file references in `payload.references.files` point to existing files.

- If a file does not exist, log a warning (non-blocking)
- If a file is in `/tmp/`, note that it may be ephemeral

### Consistency Checks

- `meta.token_count` within 30% of actual character/4 estimate
- `trace_id` matches UUID v4 format pattern
- `timestamp` precedes `expires_at` (if both present)
- `source.skill` matches first entry in `meta.handoff_chain`

---

## 6. Graceful Degradation

| Failure | Detection | Fallback | Warning |
|---------|-----------|----------|---------|
| Schema file missing | FileNotFoundError | Basic structural validation | "Schema not found. Basic validation only." |
| Token counting unavailable | ImportError | Character/4 estimation | "tiktoken unavailable. Using estimation." |
| Frontmatter parse failure | YAML error | Skip workflow in discovery | "Could not parse {skill}. Skipping." |
| Handoff file write failure | IOError | Write to /tmp/ fallback | "Using fallback location." |
| Missing trace_id (legacy) | Field absent | Auto-generate UUID v4 | "Generated trace_id for legacy handoff." |

---

## 7. Failure Modes

| Failure | Component | Detection | Recovery |
|---------|-----------|-----------|----------|
| Handoff JSON corrupted | Schema | Parse error | Re-create handoff |
| Target workflow not found | Registry | Empty discovery | Show available workflows, ask user |
| Token limit exceeded after compression | Validation | Count > 2000 | Manual payload reduction |
| Trace log write failure | Tracing | IOException | Continue without tracing |
| Circular handoff detected | Validation | Chain contains target | Reject with chain visualization |
| Schema version mismatch | Schema | Version != 3.0 | Reject; document adapter path |

---

## 8. Default Values for Optional Fields

| Field | Default | Rationale |
|-------|---------|-----------|
| `expires_at` | `null` | Most handoffs don't expire |
| `target.invocation` | `null` | Discovered at runtime |
| `target.expected_phase` | `null` | Target decides entry point |
| `context.focus_areas` | `[]` | No specific focus |
| `context.known_gaps` | `[]` | No known gaps |
| `context.success_criteria` | `[]` | Determined by target |
| `payload.session` | `null` | Not always needed |
| `payload.references` | `null` | Not always needed |
| `meta.handoff_chain` | `[source.skill]` | Minimum chain |
| `meta.user_approved` | `false` | Not approved unless explicit |

---

## 9. Handoff File Location Convention

Standard path:

```
{session_path}/handoffs/{source}-to-{target}-{timestamp}.json
```

Example:

```
/tmp/skill-editor-session/session-123/handoffs/skill-editor-to-programming-pm-20260207T183000Z.json
```

---

## 10. Payload Hash

`payload_hash` is SHA256 of `JSON.stringify(handoff.payload)` in compact form (no whitespace).

Only `sha256` algorithm is supported. Format: `sha256:{64 hex characters}`.

---

## References

- **Universal Handoff Schema v3.0**: `./universal-handoff-schema-v3.0.json`
- **JSON Schema**: https://json-schema.org/
- **Research findings**: Context bloat anti-pattern (2000 token limit)
