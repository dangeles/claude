# Distributed Tracing for Workflow Handoffs

**Version**: 1.0
**Schema**: Universal Handoff Schema v3.0
**Date**: 2026-02-07

---

## Overview

Distributed tracing enables:

1. **Cross-workflow visibility** -- Track handoff chains across multiple workflows
2. **Debugging** -- Identify where handoffs fail or get stuck
3. **Performance analysis** -- Measure handoff latency and bottlenecks
4. **Audit trail** -- Full history of workflow transitions

**Key Principle**: Lightweight tracing that does not add significant overhead (<10ms per handoff).

---

## Standards Alignment

This tracing system is designed for compatibility with:

- **W3C Trace Context** -- trace_id format: UUID v4 = 128-bit identifier
- **OpenTelemetry** -- span hierarchy, event model
- **Google Dapper** -- trace propagation through handoff chain

For future OTel integration, trace_id maps to W3C traceparent:

```
traceparent: 00-{trace_id_no_hyphens}-{span_id_no_hyphens}-01
```

Note: UUID v4 format is used for the pilot. Migration to 32-hex-char W3C format may occur in v1.1 for direct OTel tool integration.

---

## Trace ID Format

### UUID v4 (Standard)

Generate at workflow start:

```bash
python3 -c "import uuid; print(uuid.uuid4())"
# Example: 550e8400-e29b-41d4-a716-446655440000
```

**Properties**:

- 128-bit unique identifier
- Globally unique (no coordination needed)
- Standardized format (RFC 4122)
- Human-readable (hyphenated hex)
- Pattern: `^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$`

---

## Trace ID Propagation

```
User initiates workflow
    |
    v
Generate trace_id (UUID v4)
    |
    v
skill-editor (trace_id: 550e8400-...)
    | handoff
    v
programming-pm (trace_id: 550e8400-...)  <-- same trace_id
    | handoff
    v
senior-developer (trace_id: 550e8400-...) <-- same trace_id
    |
    v
Complete
```

### Trace ID Lifecycle

1. Generated at workflow initialization (source workflow)
2. Stored in session state (e.g., `session-state.json` or `trace-id.txt`)
3. Included in every handoff payload as `handoff.trace_id`
4. Logged with every trace event
5. Persisted in trace log for post-mortem analysis

### Propagation Rules

- **Same trace, same ID**: All workflows in a handoff chain share the same trace_id
- **New trace for new user request**: Each distinct user request gets a new trace_id
- **Never reuse**: Do not reuse a trace_id from a previous workflow run
- **Immutable**: Once generated, trace_id does not change within a trace

---

## Event Model

### Event Types

| Event Type | Description | Workflow | When Logged |
|-----------|-------------|----------|-------------|
| `workflow_started` | Workflow initialization | Source | At workflow start |
| `phase_started` | Phase/stage begins | Source | At phase boundary |
| `phase_completed` | Phase/stage ends | Source | At phase completion |
| `handoff_created` | Handoff payload created | Source | Before transmission |
| `handoff_validated` | Handoff validated against schema | Source | After validation |
| `handoff_transmitted` | Handoff sent to target | Source | After transmission |
| `handoff_received` | Handoff received | Target | On receipt |
| `handoff_processed` | Handoff processing complete | Target | After processing |
| `error` | Error occurred | Any | On error |
| `workflow_completed` | Workflow finished | Any | At workflow end |

### Event Structure

```json
{
  "trace_id": "550e8400-e29b-41d4-a716-446655440000",
  "span_id": "7e4a6b8c-1234-5678-90ab-cdef12345678",
  "parent_span_id": null,
  "timestamp": "2026-02-07T18:30:00.123Z",
  "event_type": "handoff_transmitted",
  "workflow": "skill-editor",
  "phase": "Phase 3: Decision & Review",
  "details": {
    "target": "programming-pm",
    "schema_version": "3.0",
    "token_count": 1847,
    "validation": "passed"
  },
  "duration_ms": 342
}
```

---

## Span Model (OpenTelemetry-Style)

### Span Hierarchy

```
Trace: 550e8400-e29b-41d4-a716-446655440000
|
+-- Span: skill-editor-session (7e4a6b8c-...)
|   +-- Span: phase-1-refinement (a1b2c3d4-...)
|   +-- Span: phase-2-analysis (e5f6g7h8-...)
|   +-- Span: phase-3-decision (i9j0k1l2-...)
|   +-- Span: handoff-to-programming-pm (m3n4o5p6-...)
|
+-- Span: programming-pm-session (q7r8s9t0-...)
    +-- Span: phase-0-archival (u1v2w3x4-...)
    +-- Span: phase-1-requirements (y5z6a7b8-...)
    +-- ...
```

### Span Properties

- `span_id`: Unique identifier for this span (UUID v4)
- `parent_span_id`: Parent span (`null` if root)
- `trace_id`: Trace this span belongs to
- `start_time`: When span started (ISO8601 UTC)
- `end_time`: When span ended (ISO8601 UTC)
- `duration_ms`: Span duration in milliseconds
- `tags`: Metadata (key-value string pairs)

---

## Logging Infrastructure

### Structured Logging (JSON Lines)

**Format**: JSON Lines (one JSON object per line, `.jsonl` extension)

```
{"trace_id":"550e8400-...","timestamp":"2026-02-07T15:00:00.123Z","level":"info","workflow":"skill-editor","event":"workflow_started","details":{}}
{"trace_id":"550e8400-...","timestamp":"2026-02-07T15:05:23.456Z","level":"info","workflow":"skill-editor","event":"phase_started","details":{"phase":"Phase 1"}}
{"trace_id":"550e8400-...","timestamp":"2026-02-07T15:45:12.789Z","level":"info","workflow":"skill-editor","event":"phase_completed","details":{"phase":"Phase 1","duration_ms":2389333}}
```

### Log File Locations

```
Per-session logs:
  /tmp/{workflow}-session/{session-id}/trace.jsonl

Consolidated logs (optional):
  ~/.claude/logs/trace-{trace-id}.jsonl
```

### Log Retention Policy

- **Per-session logs**: Deleted on workflow completion (unless error)
- **Per-session error logs**: Preserved for 7 days
- **Consolidated logs**: Rotated daily, kept for 7 days
- **Error logs**: Kept for 30 days

---

## Handling Legacy Handoffs Without trace_id

When receiving a handoff without `trace_id` (from v2.0 or v1.2 workflows):

1. Generate a new UUID v4 as `trace_id`
2. Log warning: `"Generated trace_id for legacy handoff from {source}"`
3. Proceed with normal processing
4. The generated `trace_id` becomes the root of a new trace chain

Do NOT reject handoffs solely because `trace_id` is missing.

---

## Assumptions and Limitations

1. **Single-machine assumption**: All workflows run on the same machine. Timestamps are comparable without clock synchronization.
2. **Sequential handoffs**: One handoff at a time from a given source.
3. **Local file system**: Trace logs stored locally. No distributed aggregation.
4. **All timestamps MUST use UTC with Z suffix** for consistency.
5. **100% sampling** (all events traced). If volume exceeds 100 traces/day, consider head-based sampling (trace 10% of workflows).

---

## Correlation and Business Context

`trace_id` serves both technical tracing and business correlation purposes. If a separate business correlation ID is needed in the future, add a new field (e.g., `correlation_id`) without breaking the `trace_id` contract.

---

## Trace Visualization

### Trace Timeline

```
Timeline for trace 550e8400-e29b-41d4-a716-446655440000
===========================================================

15:00:00  +----------------------------------------------+
          | skill-editor (3.5 hours)                     |
15:30:00  |   Phase 1: Refinement (45 min)               |
16:15:00  |   Phase 2: Parallel Analysis (1.5 hours)     |
17:45:00  |   Phase 2.5: Strategic Review (45 min)       |
18:30:00  |   Phase 3: Decision & Review (45 min)        |
          +----------------------------------------------+
18:30:00                      | handoff
          +----------------------------------------------+
18:30:00  | programming-pm (estimated 30-44 hours)       |
          |   Phase 0: Archival Setup (15 min)           |
          |   Phase 1: Requirements (1 hour)             |
          |   ...                                        |
          +----------------------------------------------+
```

### Handoff Chain Visualization

```
Handoff Chain for trace 550e8400-e29b-41d4-a716-446655440000
===========================================================

skill-editor
    | schema: 3.0
    | tokens: 1847
    | validation: PASSED
    | latency: 342 ms
    v
programming-pm
    | schema: 1.2 (internal)
    | adapter: 3.0 -> 1.2
    | validation: PASSED
    | latency: 156 ms
    v
senior-developer
    | schema: 1.2 (internal)
    | validation: PASSED
    | latency: 89 ms
    v
COMPLETE
```

---

## Integration with Universal Handoff Schema

### Trace Fields in Handoff

The `tracing` section in the handoff schema carries trace data:

```json
{
  "handoff": {
    "trace_id": "550e8400-e29b-41d4-a716-446655440000",
    "tracing": {
      "events": [
        {
          "timestamp": "2026-02-07T18:30:00.123Z",
          "event_type": "transmitted",
          "workflow": "skill-editor",
          "details": {"target": "programming-pm"}
        }
      ],
      "parent_span_id": "7e4a6b8c-1234-5678-90ab-cdef12345678",
      "tags": {
        "priority": "high",
        "estimated_duration": "30-44 hours"
      }
    }
  }
}
```

### Workflow Integration (Bash)

```bash
# Phase 0: Initialize trace
TRACE_ID=$(python3 -c "import uuid; print(uuid.uuid4())")
echo "$TRACE_ID" > ${SESSION_DIR}/trace-id.txt

# Log workflow start
echo "{\"trace_id\":\"$TRACE_ID\",\"timestamp\":\"$(date -u +%Y-%m-%dT%H:%M:%SZ)\",\"event_type\":\"workflow_started\",\"workflow\":\"skill-editor\"}" >> ${SESSION_DIR}/trace.jsonl

# Phase N: Include trace_id in handoff
TRACE_ID=$(cat ${SESSION_DIR}/trace-id.txt)
# ... include in handoff JSON as handoff.trace_id
```

---

## Performance Impact

| Operation | Latency | Notes |
|-----------|---------|-------|
| Generate trace_id | <1 ms | UUID generation |
| Log event (JSON) | 2-5 ms | Append to file |
| Propagate trace_id | 0 ms | Pass string |
| Query trace (100 events) | 10-20 ms | Read and parse log file |
| Analyze trace | 50-100 ms | Full trace analysis |

**Total overhead per handoff**: <10 ms (negligible compared to handoff latency of 100-500ms).

---

## References

- **OpenTelemetry**: https://opentelemetry.io/docs/concepts/observability-primer/
- **Google Dapper**: https://research.google/pubs/pub36356/
- **W3C Trace Context**: https://www.w3.org/TR/trace-context/
- **Universal Handoff Schema v3.0**: `./universal-handoff-schema-v3.0.json`
