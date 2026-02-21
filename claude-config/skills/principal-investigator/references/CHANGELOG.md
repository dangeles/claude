# Principal Investigator Skill - Changelog

## 2026-01-28: Team-Directed Workflow Update

### Major Changes

**Previous model**: Direct delegation to bioinformatician or program-officer

**New model**: Team feedback → PI decision → Technical-PM delegation

### New Workflow

1. **Gather team feedback** (ordered by task type)
2. **Synthesize input** and make final decision
3. **Delegate via technical-pm** for implementation
4. **Interpret results** with biological context

### Team Structure

PI now directs:
- **bioinformatician**: Data analysis, statistical methods
- **software-developer**: Code implementation, architecture
- **biologist-commentator**: Biological relevance, interpretation
- **calculator**: Quantitative validation, feasibility
- **technical-pm**: Coordinates implementation team
- **program-officer**: Coordinates research tasks (literature, synthesis)

### Feedback Ordering

**Implementation tasks** (least → most technical):
```
biologist-commentator → bioinformatician → calculator → software-developer
```

**Interpretation tasks** (most → least technical):
```
software-developer → calculator → bioinformatician → biologist-commentator
```

**Research tasks**:
```
Skip team feedback → delegate directly to program-officer
```

### PI Authority

Added explicit guidance that PI has **full authority** to:
- Accept all feedback
- Accept some feedback, reject others
- Override technical recommendations
- Synthesize conflicting input
- Make executive decisions

### Key Additions

1. **Leadership Principles**: Authority, synthesis, scientific judgment, pragmatism
2. **When to Disregard Feedback**: Specific scenarios and examples
3. **Example Workflows**: Three detailed examples showing feedback gathering, synthesis, and delegation
4. **Quick Reference**: Feedback order summary for common task types
5. **Updated Quality Checklist**: Pre-delegation and pre-finalization checks

### Integration Updates

**Technical-PM**: For implementation/execution tasks
**Program-Officer**: For research/literature tasks

Clear decision rules for when to use each coordinator.

### Rationale

This update formalizes the team-based research structure where:
- PI leads strategically
- Specialists provide domain expertise
- PI synthesizes and decides
- Technical-PM manages execution
- Everyone contributes their expertise in the right order

### Files Modified

- `SKILL.md`: Complete rewrite of workflow section
- `CHANGELOG.md`: This file (new)

### Backward Compatibility

**Breaking changes**: None - old workflow patterns still work

**Enhanced**: PI skill now has explicit team coordination workflow

**Migration**: No action needed - new workflow is additive
