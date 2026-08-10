# Anthropic Guidelines Summary for Skill Development

This document summarizes key Anthropic guidelines and best practices for developing Claude Code skills.

## Contents

- Core Principles
- Skill Design Patterns
- Integration Patterns
- Error Handling Patterns
- Anti-Patterns to Avoid
- Skill Invocation Best Practices
- Testing and Validation
- Model-Specific Considerations
- Accessibility and Documentation

## Core Principles

### 1. Clear, Specific Instructions

**Guideline:** Skills should provide clear, unambiguous instructions that Claude can follow deterministically.

✅ **Good:**
```markdown
## Workflow

### Step 1: Read File
Use the Read tool to read `claude-config/skills/{skill-name}/SKILL.md`

### Step 2: Validate YAML
Extract YAML frontmatter (content between first two `---` markers)
Parse using YAML validator
```

❌ **Bad:**
```markdown
## Workflow

1. Check the file
2. Make sure it's good
3. Continue if okay
```

### 2. Measurable Success Criteria

**Guideline:** Define clear success criteria so Claude knows when a task is complete.

✅ **Good:**
```markdown
## Success Criteria

- [ ] All tests pass (`pytest` returns 0)
- [ ] YAML validates (`yamllint` returns 0)
- [ ] User approves final result
```

❌ **Bad:**
```markdown
## Success Criteria

- Make it better
- Improve quality
- Finish the work
```

### 3. Appropriate Tool Usage

**Guideline:** Specify which Claude tools to use for each step.

**Available Tools:**
- **Read**: Read files (use instead of `cat`)
- **Write**: Create new files
- **Edit**: Modify existing files (exact string replacement)
- **Bash**: Run shell commands (git, npm, pytest, etc.)
- **Glob**: Find files by pattern
- **Grep**: Search file contents
- **Task**: Launch specialized agents
- **AskUserQuestion**: Interactive user input

✅ **Good:**
```markdown
### Step 2: Search for Pattern

**Tools:** Grep

Search for pattern "function.*export" in `src/**/*.ts` files.
Use output_mode: "content" to show matches.
```

❌ **Bad:**
```markdown
### Step 2: Search for Pattern

Find the pattern somehow.
```

### 4. Context Management

**Guideline:** Keep skills focused. Use agents for complex multi-step tasks.

**Skill Types:**

- **Main context skills** (`context: main`):
  - Execute in main conversation
  - Good for: Quick workflows, user interaction
  - Limit: ~3-5 major steps

- **Agent-based workflows**:
  - Use Task tool to launch specialized agents
  - Good for: Complex workflows, parallel execution
  - Agents run in isolated contexts

✅ **Good:**
```markdown
## Workflow

### Phase 1: Refinement
Launch `skill-editor-request-refiner` agent via Task tool.
Agent refines user request interactively.

### Phase 2: Analysis
For a substantial change, launch the independent analysis tracks in one message so they
run concurrently — e.g. `brainstorming-pm` for multi-perspective analysis and
`skill-editor-edge-case-simulator` for failure modes. For a small change, skip this phase.
```

❌ **Bad (for complex workflow):**
```markdown
## Workflow

Step 1: Analyze codebase (100 substeps)
Step 2: Research external sources (50 substeps)
Step 3: Generate plan (30 substeps)
...
Step 50: Verify everything
```

## Skill Design Patterns

### Pattern 1: Interactive Refinement

Use for tasks requiring user clarification:

```markdown
## Workflow

### Step 1: Gather Requirements
Use AskUserQuestion to ask:
- Question 1: Which approach? [Option A, Option B, Option C]
- Question 2: What scope? [Minimal, Standard, Comprehensive]

### Step 2: Confirm Specification
Present refined specification.
Ask user: "Does this match your intent?"
```

**When to use:** User request is ambiguous or has multiple valid interpretations.

### Pattern 2: Multi-Phase with Quality Gates

Use for complex workflows requiring validation:

```markdown
## Workflow

### Phase 1: Analysis
- Analyze codebase
- Quality Gate 1: Analysis complete

### Phase 2: Implementation
- Make changes
- Quality Gate 2: Tests pass

### Phase 3: Verification
- Verify results
- Quality Gate 3: User approval
```

**When to use:** Complex tasks where errors are costly (system config, deployments).

### Pattern 3: Parallel Execution

Use for independent tasks that can run simultaneously:

```markdown
## Workflow

### Step 3: Analysis
Launch one agent per genuinely independent track, in a single message with multiple Task
calls so they run concurrently. Keep the count to the number of distinct tracks — one
agent per track, not several agents per track.
```

**When to use:** several independent, sizeable analyses are needed. If one agent can cover
the ground, use one.

### Pattern 4: Agent Delegation

Use for specialized subtasks:

```markdown
## Workflow

### Step 1: Orchestrator Determines Approach
Analyze user request.
Identify which specialist agents needed.

### Step 2: Delegate to Specialists
- For code changes: Launch `executor` agent
- For research: Launch `researcher` agent
- For review: Launch `reviewer` agent
```

**When to use:** Workflow has distinct specialist phases.

## Integration Patterns

### Pattern 1: Git Integration

```markdown
## Workflow

### Step 5: Commit Changes
Use Bash tool to commit:

\`\`\`bash
git add claude-config/skills/{skill-name}/
git commit -m "$(cat <<'EOF'
feat(skill-name): Add new feature

Detailed description here.

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
\`\`\`

Do NOT use git commit --amend unless explicitly requested.
```

**Important:**
- Always stage specific files (not `git add -A`)
- Use HEREDOC for multi-line commit messages
- Never use `--no-verify` flag (respects pre-commit hooks)
- Never force push to main/master

### Pattern 2: External Script Integration

```markdown
## Workflow

### Step 3: Sync Configuration
Use Bash to run sync-config.py:

\`\`\`bash
# Preview changes
./sync-config.py push --dry-run

# Apply changes (will prompt user for confirmation)
./sync-config.py push
\`\`\`
```

**Important:**
- Prefer dry-run first
- Let scripts handle user prompts (don't add `yes |`)
- Check exit codes: `./script.sh && echo "Success" || echo "Failed"`

### Pattern 3: Planning Journal Integration

```markdown
## Workflow

### Step 1: Create Planning Entry
Use Bash to create entry:

\`\`\`bash
./sync-config.py plan --title "Brief description of change"
\`\`\`

### Step 7: Update Planning Entry
Document outcome in planning entry:
- What changed
- Commit SHA
- Testing results
```

## Error Handling Patterns

### Pattern 1: Graceful Degradation

```markdown
## Workflow

### Step 3: Attempt Validation
Run validation script.

**If validation fails:**
- Document specific errors
- Attempt automatic fixes (if safe)
- If unfixable: Ask user for guidance via AskUserQuestion

**If validation succeeds:**
- Proceed to next step
```

### Pattern 2: Rollback on Failure

```markdown
## Workflow

### Step 5: Sync Changes
Attempt sync via `./sync-config.py push`.

**If sync fails:**
1. Rollback: `git reset --hard HEAD`
2. Document failure in planning journal
3. Ask user whether to retry or abort
```

### Pattern 3: Pre-flight Checks

```markdown
## Workflow

### Step 1: Pre-flight Checks
Before making changes:
- [ ] Check git status (no uncommitted changes)
- [ ] Verify sync status (no divergence)
- [ ] Validate environment (required tools installed)

If any check fails: Stop and ask user to resolve.
```

## Anti-Patterns to Avoid

### ❌ Anti-Pattern 1: Vague Instructions

```markdown
## Workflow

1. Improve the code
2. Make it better
3. Finish up
```

**Why bad:** Claude doesn't know what "better" means.

**Fix:** Define specific, measurable criteria.

### ❌ Anti-Pattern 2: Assuming Context

```markdown
## Workflow

1. Edit the file (you know which one)
2. Update the config (as usual)
```

**Why bad:** Claude doesn't have implicit context.

**Fix:** Specify exact file paths and changes.

### ❌ Anti-Pattern 3: Unstructured Tool Usage

```markdown
## Workflow

Use whatever tools you need to get it done.
```

**Why bad:** No guidance on appropriate tool selection.

**Fix:** Specify which tools to use in each step.

### ❌ Anti-Pattern 4: No Quality Gates

```markdown
## Workflow

1. Make changes
2. Commit
3. Done
```

**Why bad:** No validation before committing.

**Fix:** Add validation steps between phases.

### ❌ Anti-Pattern 5: Monolithic Skills

```markdown
A single skill that does:
- Analysis
- Research
- Planning
- Implementation
- Testing
- Deployment
- Monitoring

(50+ steps in one skill)
```

**Why bad:** Hard to maintain, debug, reuse.

**Fix:** Break into focused skills or use agent delegation.

## Skill Invocation Best Practices

### When Claude Should Use a Skill

Skills should include clear "When to Use" criteria:

```markdown
## When to Use This Skill

Use this skill when ANY of these conditions apply:

1. **User explicitly requests skill**: "/skill-name" or "use the skill-name skill"
2. **Task matches description**: User asks to [specific task]
3. **Keyword trigger**: User mentions "[specific keyword]"

## When NOT to Use This Skill

Do NOT use this skill when:

1. **Different skill exists**: Use [other-skill] instead for [scenario]
2. **Task too simple**: Manual approach is clearer
3. **Scope mismatch**: Task is outside skill's domain
```

### Skill Composition

Skills can invoke other skills:

```markdown
## Workflow

### Step 3: Sub-skill Invocation
If edge case detected:
- Invoke `edge-case-handler` skill via Skill tool
- Pass parameters: `--scenario={scenario}`

Otherwise:
- Continue with standard workflow
```

## Testing and Validation

Before committing a skill, run the objective checks: the YAML frontmatter parses, the skill
invokes without error, and the sync status is clean.

### Validation Commands

```bash
# Validate YAML
python3 << EOF
import yaml
with open('claude-config/skills/{skill-name}/SKILL.md', 'r') as f:
    parts = f.read().split('---', 2)
    yaml.safe_load(parts[1])
EOF

# Test invocation
claude << EOF
/{skill-name} "test input"
EOF

# Verify sync
./sync-config.py status
```

## Model-Specific Considerations

This repo's configured default is Claude Opus 5 (`claude-opus-5` in
`claude-config/settings.json`).

### Claude Opus 5 — default here

- Strengths: complex agentic coding, long-horizon work, deep reasoning
- Use for: architectural decisions, multi-file changes, critical reviews
- Notes: verifies its own work unprompted and delegates readily, so skills should not
  instruct it to double-check or to delegate more. Writes longer output by default;
  prompt for conciseness rather than lowering effort.

### Claude Sonnet 5

- Strengths: near-Opus quality on coding and agentic work, lower cost
- Use for: high-volume or latency-sensitive workflows

### Claude Haiku 4.5

- Strengths: fastest and cheapest
- Use for: simple, well-scoped, speed-critical tasks

### Claude Fable 5 — not usable for life-sciences work

Anthropic's most capable widely released model, but its safety classifiers target
biology and life-sciences content, including lab methods and molecular mechanisms.
Benign research tasks trip them and fall back to Opus 4.8. Do not set it as the default
for this repo's bioinformatics and research skills.

**Agent Model Selection:**

Every agent definition sets `model:` deliberately. Do **not** default to Opus — that
doctrine produced 27 Opus-pinned agents and 99.5% of this repo's lifetime output tokens on
the top tier. Match the tier to the work:

- `opus` — orchestrators, authoring, and final go/no-go judgment
- `sonnet` — review, analysis, structured generation against clear criteria
- `haiku` — mechanical scans and audits against explicit rules
- `inherit` — only when the agent should track whatever the session is running

Weigh fan-out too: an agent dispatched N-at-a-time multiplies its tier's cost N times.

```yaml
# Mechanical execution of an already-approved plan.
name: skill-editor-executor
model: sonnet

# Final go/no-go gate — judgment is the whole product.
name: skill-editor-adversarial-reviewer
model: opus
```

Use the tier aliases (`opus`/`sonnet`/`haiku`/`inherit`), never pinned version IDs; pinned
IDs across this repo went stale twice.

## Accessibility and Documentation

### Clear Section Headers

Use descriptive headers:

✅ **Good:**
```markdown
## When to Use This Skill
## Workflow
## Quality Gates
## Examples
```

❌ **Bad:**
```markdown
## Usage
## Process
## Checks
## Demos
```

### Progressive Disclosure

Structure content from high-level to detailed:

```markdown
# Skill Name

Brief one-sentence overview.

## When to Use This Skill

High-level criteria (3-5 bullets).

## Workflow

Step-by-step process.

### Step 1: High-Level Name

Detailed description.

**Tools:** Specific tools used
**Input:** What this step needs
**Output:** What this step produces
```

### Examples with Context

Provide realistic examples:

```markdown
## Examples

### Example 1: Add Parallel Execution to Researcher

**Context:** User wants to speed up web research by searching 3 sources simultaneously.

**Input:**
\`\`\`
/skill-editor "Add parallel web search to researcher skill"
\`\`\`

**Output:**
- Refined spec: "Modify researcher/SKILL.md to invoke 3 WebSearch calls in parallel"
- Implementation: Updated Step 2 in researcher workflow
- Validation: Tested with sample research request
- Result: 3x faster research phase
```

## Summary

Key Anthropic guidelines for skills:

1. **Clear instructions**: No ambiguity, specific tools, exact paths
2. **Measurable criteria**: Define "done", validation steps
3. **Appropriate tools**: Match tool to task
4. **Context management**: Main for simple, agents for complex
5. **Quality gates**: Validate at each phase
6. **Error handling**: Graceful degradation, rollback
7. **Testing**: Validate before committing
8. **Documentation**: Progressive disclosure, realistic examples

Following these guidelines ensures skills are:
- **Reliable**: Consistent behavior
- **Maintainable**: Easy to update
- **Composable**: Works with other skills
- **User-friendly**: Clear, helpful
