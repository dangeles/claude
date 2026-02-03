# Dependency Detection

## Overview

When technical-pm receives a multi-part goal, it must determine which tasks can run in parallel vs which must run sequentially. This document defines the detection algorithm.

## Dependency Types

### 1. Data Dependency (Sequential Required)

Task B explicitly needs output from Task A.

**Detection keywords**: "based on", "using results from", "after", "then", "with the"

**Examples**:
- "Research hepatocytes, THEN synthesize findings" - Sequential
- "Calculate membrane area BASED ON literature values" - Sequential
- "Edit the document USING the fact-check results" - Sequential

### 2. Implicit Dependency (Ask User)

Tasks share domain terms but no explicit connection.

**Detection**: Shared technical terms without dependency keywords

**Examples**:
- "Review hepatocyte literature AND calculate hepatocyte oxygen needs" - Ambiguous
- "Analyze drug candidates AND calculate dosing" - Ambiguous

### 3. Independent (Parallel Safe)

Tasks have no connection.

**Detection**: Different domains, no shared terms, no dependency keywords

**Examples**:
- "Research hepatocytes AND calculate manufacturing costs" - Parallel safe
- "Review literature on topic A AND analyze data for topic B" - Parallel safe

## Detection Algorithm

```python
def classify_dependency(task_a, task_b, goal_text):
    """
    Returns: 'sequential', 'parallel', or 'ask_user'
    """

    # Step 1: Check explicit dependency markers in original goal
    explicit_markers = [
        "based on", "using", "after", "then", "with the",
        "from the", "following", "once", "when"
    ]

    goal_lower = goal_text.lower()
    for marker in explicit_markers:
        if marker in goal_lower:
            # Check if marker connects task_a output to task_b
            if references_task_output(goal_text, task_a, task_b, marker):
                return 'sequential'

    # Step 2: Check for shared technical terms
    terms_a = extract_technical_terms(task_a.description)
    terms_b = extract_technical_terms(task_b.description)
    shared_terms = terms_a.intersection(terms_b)

    if len(shared_terms) > 2:
        return 'ask_user'  # Ambiguous - needs clarification

    # Step 3: Check skill-specific dependencies
    # Some skills naturally depend on others
    natural_dependencies = {
        'synthesizer': ['researcher'],      # Synthesizer needs research
        'devils-advocate': ['synthesizer'], # Review needs draft
        'fact-checker': ['researcher', 'synthesizer'],
        'editor': ['fact-checker', 'devils-advocate'],
    }

    if task_b.skill in natural_dependencies:
        if task_a.skill in natural_dependencies[task_b.skill]:
            return 'sequential'

    # Step 4: Check duration - parallel only worthwhile for long tasks
    if task_a.estimated_duration < 5 and task_b.estimated_duration < 5:
        return 'sequential'  # Overhead exceeds benefit

    # Step 5: Default to parallel for truly independent work
    return 'parallel'
```

## User Confirmation (When Ambiguous)

When `ask_user` is returned:

```
I detect potential dependency between these tasks:

Task A: Review hepatocyte literature
Task B: Calculate oxygen delivery requirements

Shared terms: hepatocyte, oxygen, cells

These tasks MIGHT benefit from sequential execution if the calculator
needs literature values. Or they might be truly independent.

Options:
1. [Sequential] - Calculator waits for researcher (safer, slower)
2. [Parallel] - Run simultaneously (faster, may miss integration)
3. [Help me decide] - I'll explain the trade-offs in detail

Your choice:
```

## Conservative Defaults

When uncertain, default to:
- **Sequential if quality matters** - Research workflows, final documents
- **Parallel if speed matters** - Exploratory work, drafts

User can override with explicit flags:
```
technical-pm goal="..." --sequential  # Force all sequential
technical-pm goal="..." --parallel    # Force all parallel where safe
```

## Estimation Heuristics

| Skill | Estimated Duration | Good for Parallel? |
|-------|-------------------|-------------------|
| researcher | 30-60 min | Yes - long, independent |
| calculator | 5-30 min | Yes - long, independent |
| synthesizer | 15-30 min | Maybe - often needs researcher |
| devils-advocate | 5-15 min | No - needs draft to review |
| fact-checker | 5-15 min | No - needs content to check |
| editor | 5-15 min | No - needs near-final content |
| archivist | 1-5 min | No - quick, sequential fine |

## Natural Skill Dependencies

Some skills have inherent dependencies based on their function:

```
researcher ─────────┬──────────────────────────> fact-checker
                    │
                    v
               synthesizer ────> devils-advocate ────> editor
                    ^
                    │
calculator ─────────┘
```

**Always sequential**:
- synthesizer after researcher (needs research to synthesize)
- devils-advocate after synthesizer (needs draft to critique)
- editor after devils-advocate (needs reviewed draft)
- fact-checker after researcher or synthesizer (needs claims to verify)

**Potentially parallel**:
- researcher AND calculator (if working on different aspects)
- Multiple researchers on different topics
- fact-checker AND devils-advocate (both review synthesizer output)

## False Positive Mitigation

Track decision outcomes in session:
- If user frequently chooses "parallel" for similar patterns, learn preference
- Log dependency decisions for retrospective analysis

## False Negative Detection

After parallel tasks complete, before synthesis:
1. Check if outputs reference same entities
2. If Task A output mentions values that Task B should have used:
   - Flag for user review
   - Ask: "These outputs both reference [X]. Should Task B have waited for Task A?"

## Examples

### Example 1: Clear Sequential

**Goal**: "Review the hollow fiber literature and then synthesize the findings into a summary"

**Analysis**:
- Explicit marker: "and then"
- Task A: Review literature (researcher)
- Task B: Synthesize findings (synthesizer)
- Decision: **SEQUENTIAL** (explicit marker + natural dependency)

### Example 2: Clear Parallel

**Goal**: "Research hollow fiber membranes for oxygenation AND calculate the cost of injection molding"

**Analysis**:
- No explicit markers connecting tasks
- Task A: Research membranes (researcher)
- Task B: Calculate costs (calculator)
- Shared terms: None (different domains)
- Decision: **PARALLEL** (independent tasks)

### Example 3: Ambiguous

**Goal**: "Review hepatocyte oxygen consumption literature AND calculate oxygen delivery requirements for the bioreactor"

**Analysis**:
- No explicit markers
- Task A: Review hepatocyte literature (researcher)
- Task B: Calculate oxygen delivery (calculator)
- Shared terms: "oxygen", "hepatocyte" (>2 shared terms)
- Decision: **ASK USER** (calculation might need literature values)

### Example 4: Duration-Based Sequential

**Goal**: "Archive the completed document AND update the table of contents"

**Analysis**:
- No explicit markers
- Task A: Archive document (archivist, ~2 min)
- Task B: Update ToC (archivist, ~2 min)
- Both < 5 minutes
- Decision: **SEQUENTIAL** (parallel overhead exceeds benefit)

## Override Behavior

When user explicitly requests execution mode:

**--sequential flag**:
- Ignore all parallelization opportunities
- Execute tasks in order listed in goal
- Useful for debugging or when integration is critical

**--parallel flag**:
- Parallelize all tasks without natural dependencies
- Still respect natural skill dependencies
- Useful when speed is priority over potential integration issues
