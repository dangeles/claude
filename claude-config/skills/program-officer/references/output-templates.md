# Program Officer — Output Templates

Verbatim templates and worked example deliverables for the program-officer skill.
SKILL.md handles routing, decision, and orchestration; this file holds the
fill-in-the-blank output formats and the full example reports.

## Progress Update Template

Use when checking specialist status:

```
**Progress Check**: [Specialist Name]

**Task**: [Original task assigned]
**Time elapsed**: [X minutes/hours]
**Expected completion**: [Original estimate]

**Questions**:
1. Current progress? (concrete metric: papers read, calculations done)
2. Blockers or uncertainties?
3. Estimated time remaining?

**Next action based on response**:
- On track → Continue, check again in 60-90 min
- Blocked → Clarify/reassign/escalate
- Scope expanding → Refocus or escalate
- Nearly done → Prepare next specialist
```

## Deliverable Format

Return to domain coordinator with:

```markdown
# Research Coordination Report: [Task]

**Coordinated**: [Date and time range]
**Specialists involved**: [List]

## Recommendation
[Clear, actionable recommendation]

## Supporting Evidence
**Literature**: [Key findings from researcher]
  - Papers reviewed: X
  - Key citations: [list]
  - Consensus: [what most papers agree on]

**Quantitative**: [Key results from calculator]
  - Analysis performed: [method]
  - Key finding: [numerical result]
  - Interpretation: [what it means for feasibility]

**Validation**: [Key confirmations from fact-checker]
  - Claims verified: [list]
  - Assumptions checked: [list]
  - Issues identified: [if any]

**Synthesis**: [Integrated perspective from synthesizer]
  - Cross-source themes: [patterns]
  - Contradictions resolved: [how]
  - Confidence drivers: [what increases/decreases confidence]

## Confidence Level
[HIGH / MEDIUM / LOW]

**Justification**:
- HIGH if: Multiple independent sources converge, quantitative validation passes, no major caveats
- MEDIUM if: Some contradictions, limited data, minor caveats
- LOW if: Conflicting evidence, insufficient data, major assumptions

## Alternative Options
[If primary recommendation fails or has constraints]
1. [Alternative 1]: [brief rationale]
2. [Alternative 2]: [brief rationale]

## Implementation Notes
[What domain coordinator needs to know for implementation]
- Required inputs: [data, parameters, etc.]
- Expected outputs: [format, interpretation]
- Caveats: [limitations, assumptions]
- Validation steps: [how to verify implementation]

## Timeline Summary
- Literature review: [duration]
- Quantitative analysis: [duration]
- Validation: [duration]
- Synthesis: [duration]
- Total: [X hours Y minutes]
```

## Example Deliverables

These accompany the orchestration walkthroughs in SKILL.md ("Example Scenarios").
Each shows the Deliverable Format filled in for a representative coordination.

### Scenario 1: Method Selection — Deliverable

```markdown
# Research Coordination Report: Clustering Algorithm Selection

## Recommendation
**Use Leiden algorithm** with resolution=0.8

## Supporting Evidence
**Literature**:
- Papers reviewed: 12 (2020-2024)
- Leiden preferred: 10/12 papers (83%)
- Key advantage: Better handles resolution limit problem
- Citations: Traag 2019 (Leiden paper), multiple benchmarks

**Quantitative**:
- Tested on sample dataset (5000 cells)
- Leiden: 12% more stable clusters (ARI=0.89 across runs)
- Louvain: More variable (ARI=0.76 across runs)

**Validation**:
- Claims verified on our data type (sparse UMI counts)
- Computational cost similar (Leiden 5% slower, negligible)

**Synthesis**:
- Strong consensus in literature
- Quantitative testing confirms literature claims
- No significant downsides identified

## Confidence Level
**HIGH**
- Multiple independent sources converge
- Quantitative validation passes
- No major caveats

## Alternative Options
1. **Louvain**: If legacy pipeline compatibility required (slight stability loss acceptable)
2. **Hierarchical**: If deterministic results critical (slower, less resolution flexibility)

## Implementation Notes
- Use leidenalg package (Python) or Seurat (R)
- Set resolution=0.8 as starting point (tune based on cluster count)
- Run multiple iterations, check stability
- Document random seed for reproducibility

## Timeline Summary
- Literature review: 2h 10min
- Synthesis: 30min
- Quantitative testing: 45min
- Validation: 15min
- **Total: 3h 40min**
```

### Scenario 2: Statistical Validation — Deliverable

```markdown
## Recommendation
**Proceed with mixed-effects model** (batch as random effect)

## Supporting Evidence
**Quantitative**:
- Power adequate (0.85 for 2-fold changes)
- Assumptions met: residuals normal, variance homogeneous
- Batch explains 15% variance (substantial but not excessive)

**Literature**:
- Used successfully in 3 similar studies (Leek 2014, Ritchie 2015, Johnson 2007)
- Standard approach for known batch effects
- DESeq2 implementation validates well

**Validation**:
- Data structure appropriate: 4 batches, balanced design
- No confounding between batch and condition
- Batch effect visible in PCA (PC2, 15% variance)

## Confidence Level
**HIGH** - Method appropriate, assumptions met, literature precedent strong

## Alternative Options
1. **ComBat**: If batch effect more severe (>30% variance), but loses count distribution
2. **Batch as fixed effect**: If only interested in specific batches, loses generalizability
```

### Scenario 3: Unexpected Finding Validation — Deliverable

```markdown
## Recommendation
**Finding is likely real, not artifact** - report as novel with caveats

## Supporting Evidence
**Literature**:
- Rare but precedented in hypoxia conditions (2 papers: Smith 2019, Jones 2021)
- Housekeeping genes not truly "housekeeping" under stress
- Context-specific regulation documented

**Quantitative**:
- Robust across multiple normalization methods (DESeq2, TMM, CPM)
- Not driven by outliers (consistent across all replicates)
- Not batch effect (no correlation with batch)
- Validated with alternative statistical tests (Wilcoxon, t-test agree)

**Validation**:
- QC checks pass (no low-quality samples)
- Preprocessing appropriate (standard pipeline)
- Raw counts examined (not normalization artifact)

**Synthesis**:
- Literature provides biological precedent (stress response)
- Quantitative testing rules out technical artifacts
- Multiple independent lines of evidence support real biology

## Confidence Level
**MEDIUM-HIGH**
- High: Technical artifacts ruled out
- Medium: Limited biological precedent (only 2 similar papers)
- Caveat: Mechanism unclear, warrants follow-up validation

## Implementation Notes
**Report as novel finding with appropriate caveats**:
- Acknowledge limited precedent
- Suggest validation experiments (qPCR, Western blot)
- Frame as hypothesis-generating
- Note potential stress response mechanism
```
