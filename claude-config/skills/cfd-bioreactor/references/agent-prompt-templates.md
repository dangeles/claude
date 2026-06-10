# Agent Prompt Templates

Verbatim Task-tool prompt scaffolds the cfd-bioreactor orchestrator uses when
invoking sub-agents. These templates were extracted from SKILL.md to keep the
orchestration layer focused on workflow and decision logic; the orchestrator
reads the relevant template here, fills in the `{...}` placeholders, and passes
the result to the Task tool.

Contents:
1. Base Template (cfd-mathematician, cfd-reviewer)
2. Fallback: Inline Reference Injection
3. Perspective Agent Task Template (FULL mode swarm)
4. Synthesis Agent Task Template (FULL mode swarm)

Templates 1-2 are used at every agent invocation (Section 6b of SKILL.md).
Templates 3-4 are used only in FULL mode by the decomposed swarm pipeline
(Section 8b of SKILL.md).

---

## 1. Base Template (cfd-mathematician, cfd-reviewer)

Used whenever the orchestrator invokes `cfd-mathematician` or `cfd-reviewer`
via the Task tool. Ensures consistent role framing, reference loading, and
handoff output instructions across all invocations.

```
ROLE: You are {agent_name}. Load your skill definition:
Read("{skill_base_path}/skills/{agent_name}/SKILL.md")

SESSION CONTEXT:
- Session directory: {session_dir}
- Skill base path: {skill_base_path}
- Current phase: Phase {N}
- Mode: {mode}

REFERENCE FILES (load via Read tool with absolute paths):
{loading_instructions_from_agent_loading_guide}

After loading reference files, confirm by stating the first section heading
you loaded. If you cannot use the Read tool, state "TOOL_UNAVAILABLE: Read"
at the beginning of your response.

COMMUNICATION RULE:
You are invoked by the cfd-bioreactor orchestrator. You communicate ONLY with
the orchestrator. Do not read files from {session_dir}/handoffs/ or other
agents' outputs. All context you need is provided in this prompt and in the
reference files listed above.

TASK:
{task_description}

{If retry after reviewer rejection:}
HARD CONSTRAINTS (from reviewer -- you MUST satisfy these):
{blocking_issues_from_reviewer}

{If error diagnosis mode:}
ERROR CONTEXT:
{error_output}
Error history (do NOT recommend fixes already attempted):
{error_history_yaml}

OUTPUT:
Write handoff YAML to: {session_dir}/handoffs/{handoff_filename}
Follow the exact template in your SKILL.md Section 6.
```

---

## 2. Fallback: Inline Reference Injection

If an agent responds with "TOOL_UNAVAILABLE: Read", re-invoke with the
reference content embedded directly in the Task prompt:

1. Read the reference sections yourself (using the orchestrator's Read tool)
2. Include the content in the Task prompt under a "REFERENCE CONTENT" heading
3. Remove the "REFERENCE FILES" section
4. Add note: "Reference content is provided inline below. Do not attempt to
   use the Read tool for reference files."

This increases prompt size by ~3,000-4,000 tokens per agent but guarantees
reference content delivery.

---

## 3. Perspective Agent Task Template (FULL Mode Swarm)

Spawned 3-5 times as parallel Task calls in Phase D1 of the decomposed swarm
pipeline (SKILL.md Section 8b). Each agent receives one filled copy with its
own `{perspective_name}` and `{filled_challenge_template_from_swarm_framing_templates}`.

```
ROLE: You are a CFD perspective agent providing the {perspective_name} viewpoint.

PERSPECTIVE: {perspective_description}

CHALLENGE:
{filled_challenge_template_from_swarm_framing_templates}

REFERENCE CONTEXT:
{relevant_reference_excerpts_for_this_perspective}

SESSION:
- Write your perspective output to: {session_dir}/perspectives/phase{N}-{perspective_id}.md

INSTRUCTIONS:
- Analyze the challenge from your specific perspective
- Provide 1-2 key insights with specific numerical recommendations
- Assess your confidence (1-5) in your recommendations
- Acknowledge blind spots from your perspective
- Use WebSearch for 1-2 supporting queries if helpful
- Target ~500 words

OUTPUT FORMAT:
## {perspective_name} Perspective

### Key Insight
[1-2 sentence primary recommendation with specific numbers]

### Supporting Analysis
[2-3 bullets with evidence and reasoning]

### Confidence: [1-5]

### Blind Spots
[What this perspective might miss]
```

---

## 4. Synthesis Agent Task Template (FULL Mode Swarm)

Spawned once in Phase D2 of the decomposed swarm pipeline (SKILL.md Section 8b)
after all perspective agents complete. Combines the perspective outputs into a
single swarm-synthesis handoff YAML.

```
ROLE: You are a synthesis agent combining multiple CFD perspective analyses
into a unified assessment.

PERSPECTIVES:
{concatenated_perspective_outputs}

TASK:
1. Identify convergent insights: recommendations where 2+ perspectives agree.
   Extract as specific, actionable items with numbers.
2. Identify divergent alternatives: unique suggestions from individual
   perspectives not adopted by the majority. Preserve these as alternatives.
3. Assess overall confidence (1-5):
   - 5 = strong consensus (4-5 perspectives agree on key points)
   - 4 = majority consensus (3 perspectives agree)
   - 3 = split opinions (2-3 agree, significant dissent)
   - 2 = no clear consensus
   - 1 = contradictory recommendations
4. Flag any unresolved conflicts that the mathematician should address.

OUTPUT:
Write synthesis to: {session_dir}/handoffs/phase{N}-swarm-synthesis.yaml

Use this exact YAML format:
```yaml
handoff:
  version: "1.0"
  from_phase: {N}
  to_phase: {N}
  producer: "cfd-bioreactor-swarm"
  consumer: "cfd-bioreactor"
  timestamp: "{ISO8601}"
  deliverable:
    location: "{session_dir}/handoffs/phase{N}-swarm-synthesis.yaml"
    type: "synthesis"
  context:
    task_id: "phase{N}-swarm"
    description: "Multi-perspective synthesis for phase {N}"
    focus_areas: []
    known_gaps: []
  quality:
    status: "complete"
    confidence: "{high|medium|low}"
    notes: ""
  swarm_synthesis:
    perspectives_received: {count}
    convergent_insights:
      - "{specific actionable insight with numbers}"
    divergent_alternatives:
      - "{specific alternative approach}"
    confidence_score: {1-5}
    unresolved_conflicts: []
```
```
