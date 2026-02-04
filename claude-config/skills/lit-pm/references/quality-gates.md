# Quality Gates

Per-stage quality validation criteria for the lit-pm pipeline.

---

## Gate Types

### Automated Gates (Programmatic Validation)

These checks can be validated automatically without agent judgment:

| Check Type | Description | Implementation |
|------------|-------------|----------------|
| Count threshold | Paper count, word count, section count | Integer comparison |
| Presence check | Recency survey exists, all sections present | Boolean |
| Pattern match | No placeholders, no "TODO", no "[CITE]" | Regex scan |
| Range check | Word count in range, section balance | Min/max comparison |

### Human Judgment Gates (Agent Assessment)

These require agent evaluation:

| Check Type | Description | Evaluator |
|------------|-------------|-----------|
| Thesis specificity | Is the claim testable/addressable? | fact-checker |
| Narrative flow | Do sections build logically? | lit-synthesizer |
| Cross-cutting themes | Are connections identified? | lit-synthesizer |
| Claim accuracy | Does evidence support claim? | fact-checker |

---

## Per-Stage Quality Gates

### Stage 1: Scope Refinement

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Research question present | Presence | Not empty | Yes | Yes |
| Research question specific | Semantic | Not "what is X?" | No | Yes |
| Success criteria defined | Count | >= 1 | Yes | Yes |
| In-scope defined | Count | >= 1 | Yes | Yes |
| Out-of-scope defined | Count | >= 0 | Yes | No |
| User approval | User action | Approved | Yes | Yes |

### Stage 2: Review Discovery

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Total reviews collected | Count | >= 6 | Yes | No |
| Minimum reviews | Count | >= 4 | Yes | Yes |
| Convergence achieved | Count | >= 2 | Yes | Yes |
| Reviews with annotations | Count | = total | Yes | Yes |
| Coverage of themes | Semantic | All major themes | No | Yes |

### Stage 3: Outline Synthesis

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Section count | Range | 3-5 | Yes | Yes |
| Minimum sections | Count | >= 2 | Yes | Yes |
| Each section has thesis | Presence | All sections | Yes | Yes |
| Theses are specific | Semantic | Testable claims | No | Yes |
| Section balance | Range | 15-40% each | Yes | No |
| Max section size | Range | < 50% total | Yes | Yes |
| User approval (if checkpoint) | User action | Approved | Yes | Yes |

### Stage 4: Introduction

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Introduction present | Presence | Not empty | Yes | Yes |
| Word count | Range | 300-800 | Yes | No |
| Structure preview present | Semantic | Sections mentioned | No | Yes |
| Matches outline | Semantic | Consistent | No | Yes |
| Editor polish applied | Presence | True | Yes | Yes |

### Stage 5: Section Writing

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Primary papers cited | Count | >= 15 | Yes | Yes |
| Maximum papers | Count | <= 30 | Yes | No |
| Recency survey present | Presence | Subsection exists | Yes | Yes |
| Recency papers (6-12 mo) | Count | >= 3 | Yes | Yes |
| Section addresses thesis | Semantic | Agent judgment | No | Yes |
| No contradictions with intro | Semantic | Agent judgment | No | Yes |
| No placeholder text | Pattern | 0 matches | Yes | Yes |
| Word count | Range | 2000-3000 | Yes | No |

**Placeholder Patterns to Detect**:
```regex
TODO|FIXME|\[CITE\]|\[INSERT\]|\[TBD\]|\[PLACEHOLDER\]|XXX
```

### Stage 6a: Quick Validation (Per-Section)

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Paper count | Count | >= 15 | Yes | Yes |
| Recency survey | Presence | True | Yes | Yes |
| Recency paper count | Count | >= 3 | Yes | Yes |
| Citations have dates | Pattern | All citations | Yes | No |
| Thesis addressed | Semantic | Agent judgment | No | Yes |
| No placeholders | Pattern | 0 matches | Yes | Yes |
| Word count range | Range | 1500-3500 | Yes | No |

### Stage 6b: Comprehensive Fact-Check

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Cross-section consistency | Semantic | No contradictions | No | P0 |
| Citation accuracy (spot-check) | Manual | 10 random | No | P0 if major |
| Quantitative verification | Manual | Values match sources | No | P0 if wrong |
| Gap analysis | Semantic | No obvious gaps | No | P1 |
| Methodological context | Semantic | Present for key data | No | P2 |

**Output Priority Levels**:
- **P0 (Critical)**: Must fix before delivery
- **P1 (Important)**: Should fix
- **P2 (Nice-to-have)**: Editor handles

### Stage 7: Synthesis

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| Narrative flow logical | Semantic | Sections build | No | Yes |
| Cross-cutting themes identified | Presence | >= 1 | No | Yes |
| Conclusion present | Presence | Not empty | Yes | Yes |
| Conclusion synthesizes findings | Semantic | All sections referenced | No | Yes |
| Major additions flagged | Presence | If >20% new | Yes | Yes |
| User approval (if checkpoint) | User action | Approved | Yes | Yes |

### Stage 8: Editorial Polish

| Check | Type | Threshold | Automated | Blocking |
|-------|------|-----------|-----------|----------|
| P0 revisions incorporated | Count | All P0s | Yes | Yes |
| P1 revisions incorporated | Count | All P1s | Yes | No |
| Voice consistent | Semantic | Agent judgment | No | Yes |
| Formatting consistent | Pattern | Uniform citations | Yes | No |
| Final read complete | Boolean | True | Yes | Yes |

---

## Quality Floor (Override Protection)

These checks CANNOT be skipped, even with `--full-auto`:

```yaml
quality_floor:
  stage_2:
    minimum_reviews: 4
    minimum_convergence: 1  # At least 1 review found by multiple agents
    on_failure: "HALT: Review discovery below minimum threshold. Found {count} reviews, need 4."

  stage_3:
    minimum_sections: 2
    max_section_imbalance: 50%  # No section > 50% of total
    on_failure: "HALT: Outline severely unbalanced. Largest section is {pct}% of total."

  stage_5:
    minimum_papers_per_section: 10
    recency_survey_required: true
    on_failure: "HALT: Section research insufficient. {section} has only {count} papers."

  stage_6a:
    no_placeholders: true
    on_failure: "HALT: Section contains placeholder text. Found: {placeholders}"
```

Quality floor violations ALWAYS escalate to user, regardless of checkpoint plan or automation settings.

---

## Gate Execution Protocol

### Before Each Stage

```yaml
pre_stage_check:
  1_load_gate: "Load quality gate for Stage {N}"
  2_validate_inputs: "Check all required inputs present"
  3_log: "Beginning Stage {N} with {X} automated checks, {Y} judgment checks"
```

### After Each Stage

```yaml
post_stage_check:
  1_run_automated:
    for check in automated_checks:
      result = execute_check(check)
      log_result(check, result)
      if check.blocking and not result.passed:
        return FAIL

  2_run_judgment:
    for check in judgment_checks:
      result = agent_evaluate(check)
      log_result(check, result)
      if check.blocking and not result.passed:
        return FAIL

  3_compile_result:
    return {
      passed: all_blocking_passed,
      warnings: non_blocking_failures,
      details: all_results
    }
```

### On Gate Failure

```yaml
on_gate_failure:
  if blocking:
    1_log: "Gate failed: {check_name}"
    2_determine_action:
      if escape_hatch_available:
        increment_cycle_count()
        if cycle_count < max_cycles:
          return_to_producer(issues)
        else:
          escalate_to_user(options)
      else:
        return_to_producer(issues)

  if not blocking:
    1_log: "Warning: {check_name} failed (non-blocking)"
    2_add_to_warnings()
    3_continue()
```

---

## Gate Results Schema

```yaml
gate_result:
  stage: integer
  timestamp: ISO8601
  overall_status: enum  # PASS | FAIL | PASS_WITH_WARNINGS
  checks:
    - name: string
      type: enum  # automated | judgment
      blocking: boolean
      passed: boolean
      value: any  # actual value
      threshold: any  # expected value
      details: string | null
  warnings: list
  action_taken: string | null
```

---

## Escape Hatch Protocol

For blocking gates with potential infinite loops (Stage 6a):

```yaml
escape_hatch:
  max_cycles: 3

  cycle_tracking:
    section_id: string
    cycle_count: integer
    issues_per_cycle:
      - cycle: 1
        issues: list
      - cycle: 2
        issues: list

  on_max_reached:
    notify_user: |
      Section {name} has failed fact-check {max_cycles} times.

      Issues (Cycle {max_cycles}):
      {issues}

      Options:
      1. Accept section as-is (waive requirement)
      2. Adjust requirements for this section
      3. Assign different researcher to section
      4. Remove section from outline

    record_decision:
      field: workflow_state.quality_overrides[]
      value:
        section: string
        decision: enum  # accept | adjust | reassign | remove
        timestamp: ISO8601
        issues_waived: list
```

---

## Threshold Calibration

Thresholds are based on:
- Literature review best practices (15-30 papers per section is standard)
- Prior technical-pm experience (word count ranges)
- Error rates in testing (placeholder detection patterns)

### Adjustment Protocol

If thresholds prove too strict/loose in practice:
1. Log threshold violations with context
2. After 5 workflows, review violation patterns
3. Propose threshold adjustments with rationale
4. Update this document with new thresholds
