# Pipeline Invocation Example

This example demonstrates a complete research-pipeline execution from user request to final output.

## User Request

**Goal**: "Research hepatocyte oxygen consumption and write a comprehensive literature review"

## Pipeline Initialization

```yaml
pipeline:
  workflow_id: "research-hep-oxy-20260203"
  goal: "Research hepatocyte oxygen consumption and write a comprehensive literature review"
  topic: "hepatocyte oxygen consumption"
  scope: comprehensive
  constraints:
    max_pages: null  # No limit
    focus_areas: ["oxygen consumption rates", "measurement methods", "culture conditions"]
    time_limit: null
  current_stage: 1
  started_at: "2026-02-03T10:00:00Z"
```

## Stage 1: Researcher

**Invocation**: `Skill(researcher, topic="hepatocyte oxygen consumption")`

**Duration**: 2 hours 15 minutes

**Output**: `docs/literature/hepatocyte-oxygenation/researcher-draft.md`

**Summary of work**:
- Literature search via PubMed: "hepatocyte[TIAB] AND oxygen consumption[TIAB]"
- Reviews consulted: Jiang 2025 (BAL systems), Hochleitner 2019 (bioreactor design)
- 8 papers read in detail
- Quantitative values extracted with full measurement context

**Sample findings**:
| Parameter | Value | Species | Format | Source |
|-----------|-------|---------|--------|--------|
| OCR | 0.3-0.5 nmol/s/10^6 cells | Human primary | Monolayer | Martinez 2022 |
| OCR | 0.5-0.9 nmol/s/10^6 cells | Rat primary | Spheroids | Chen 2020 |
| OCR | 0.7 nmol/s/10^6 cells | Porcine | Perfused tissue | Sussman 2018 |

**Gaps identified**:
- Limited human data for 3D formats
- No standardized measurement protocol across studies
- Species conversion factors uncertain

**Validation**:
- [x] Document exists: YES
- [x] Contains citations: YES (12 inline superscripts)
- [x] Addresses topic: YES
- [x] Minimum length: YES (3,500 words)

**Handoff created**: `/tmp/handoff-research-hep-oxy-20260203-T1.yaml`

---

## Stage 2: Synthesizer

**Invocation**: `Skill(synthesizer, handoff="/tmp/handoff-research-hep-oxy-20260203-T1.yaml")`

**Duration**: 45 minutes

**Output**: `docs/literature/hepatocyte-oxygenation/synthesis-hepatocyte-oxygenation.md`

**Themes identified**:

1. **Measurement variability**: OCR values vary 3-5x across studies due to species, culture format, and measurement timing
2. **Culture format effects**: 3D formats consistently show higher OCR than 2D monolayers (1.5-2x)
3. **Temporal dynamics**: OCR changes significantly over culture duration (peak at day 2-3, decline by day 7)

**Tensions highlighted**:
- Jiang 2025 reports higher values than Martinez 2022 for similar conditions - likely due to different oxygen probe methods
- Whether rat/porcine values should be extrapolated to human systems is contested

**Implications for design**:
- Bioreactor oxygen delivery must accommodate 0.3-0.9 nmol/s/10^6 cells range
- Design should use conservative (high) estimates for safety margin
- Consider measurement method when comparing literature values

**Validation**:
- [x] Document exists: YES
- [x] Adds synthesis value: YES (themes, tensions, implications)
- [x] Maintains citations: YES
- [x] Addresses goal: YES

**Handoff created**: `/tmp/handoff-research-hep-oxy-20260203-T2.yaml`

---

## Stage 3: Devil's Advocate

**Invocation**: `Skill(devils-advocate, handoff="/tmp/handoff-research-hep-oxy-20260203-T2.yaml")`

**Duration**: 30 minutes (2 exchange rounds)

**Output**: `docs/literature/hepatocyte-oxygenation/review-report.md`

**Thesis identified**: "Bioreactor oxygen delivery requirements can be estimated from literature OCR values, with appropriate safety margins for measurement variability"

**Strategic coherence**: PASS - document addresses thesis throughout

**Challenges raised**:

### Challenge 1: Species extrapolation assumption
**Issue**: Document uses rat/porcine data to inform human bioreactor design
**Writer response**: Added section acknowledging limitation; recommends 1.3x correction factor based on Hochleitner 2019
**Resolution**: ADDRESSED

### Challenge 2: Measurement timing not standardized
**Issue**: OCR values from different time points (day 1 vs day 7) compared without normalization
**Writer response**: Added temporal context to all values; noted peak vs steady-state distinction
**Resolution**: ADDRESSED

### Uncertainty documented:
- Whether in vitro OCR values predict in vivo performance remains uncertain (no direct validation studies found)

**Validation**:
- [x] Review report generated: YES
- [x] Challenges addressed: YES (2/2)
- [x] Approval status: APPROVED with noted uncertainty

**Handoff created**: `/tmp/handoff-research-hep-oxy-20260203-T3.yaml`

---

## Stage 4: Fact-Checker

**Invocation**: `Skill(fact-checker, handoff="/tmp/handoff-research-hep-oxy-20260203-T3.yaml")`

**Duration**: 20 minutes

**Output**: `docs/literature/hepatocyte-oxygenation/fact-check-report.md`

**Citations verified**: 12/12

**Verification summary**:
| Citation | Claim | Source Match | Status |
|----------|-------|--------------|--------|
| [1] Martinez 2022 | OCR 0.3-0.5 nmol/s/10^6 | Exact match Table 2 | VERIFIED |
| [2] Chen 2020 | Spheroid OCR 1.5x monolayer | Verified in Figure 3 | VERIFIED |
| [3] Jiang 2025 | K_oA review | Comprehensive coverage | VERIFIED |
| ... | ... | ... | ... |

**Issues found**: 0

**Citation format check**:
- [x] All superscripts present in text
- [x] References section complete
- [x] DOIs resolve correctly

**Validation**:
- [x] Report generated: YES
- [x] All citations verified: YES (12/12)
- [x] Format compliance: YES

**Handoff created**: `/tmp/handoff-research-hep-oxy-20260203-T4.yaml`

---

## Stage 5: Editor

**Invocation**: `Skill(editor, handoff="/tmp/handoff-research-hep-oxy-20260203-T4.yaml")`

**Duration**: 25 minutes

**Output**: `docs/literature/hepatocyte-oxygenation/review-hepatocyte-oxygenation.md`

**CLAUDE.md checklist**:
- [x] Prose over bullets: Converted 3 bullet lists to flowing prose
- [x] Bridge transitions: Added connecting sentences between all sections
- [x] One analogy per concept: Added oxygen consumption analogy ("metabolic breathing rate")
- [x] "So what?" explicit: Added implications sentence to each major section
- [x] Glossary at end: Moved inline definitions to glossary
- [x] Acronyms defined once: OCR defined in first body paragraph
- [x] Complete sentences: Fixed 2 fragment issues
- [x] First-person plural: Changed passive to "We conclude..."
- [x] Citation format: All superscripts consistent
- [x] LaTeX for equations: Formatted OCR equation properly

**Validation**:
- [x] CLAUDE.md compliant: YES
- [x] Document flows: YES
- [x] Ready for publication: YES

---

## Pipeline Complete

```markdown
# Research Pipeline Complete

**Workflow ID**: research-hep-oxy-20260203
**Original Goal**: Research hepatocyte oxygen consumption and write a comprehensive literature review
**Duration**: 4 hours 15 minutes
**Stages Completed**: 5/5

## Final Output
**Location**: docs/literature/hepatocyte-oxygenation/review-hepatocyte-oxygenation.md

## Pipeline Summary
| Stage | Duration | Status | Notes |
|-------|----------|--------|-------|
| Researcher | 2h 15m | Complete | 8 papers reviewed, 12 citations |
| Synthesizer | 45m | Complete | 3 themes, 2 tensions identified |
| Devil's Advocate | 30m | Complete | 2 challenges addressed, 1 uncertainty documented |
| Fact-Checker | 20m | Complete | 12/12 citations verified |
| Editor | 25m | Complete | CLAUDE.md compliant |

## Quality Indicators
- Citations verified: 12/12
- Challenges addressed: 2/2
- Uncertainties documented: 1
- Style compliance: PASS

## Output Files
- Final document: docs/literature/hepatocyte-oxygenation/review-hepatocyte-oxygenation.md
- Researcher draft: docs/literature/hepatocyte-oxygenation/researcher-draft.md
- Synthesis: docs/literature/hepatocyte-oxygenation/synthesis-hepatocyte-oxygenation.md
- Review report: docs/literature/hepatocyte-oxygenation/review-report.md
- Fact-check report: docs/literature/hepatocyte-oxygenation/fact-check-report.md

## Noted Uncertainties
- Whether in vitro OCR values predict in vivo performance (no validation studies found)
- Species conversion factors based on limited comparative data
```

---

## Key Observations

1. **Handoff format worked**: Each stage received proper context and knew what to focus on
2. **Quality gates prevented issues**: Validation at each stage caught potential problems early
3. **Devil's advocate added value**: Identified species extrapolation assumption that needed explicit acknowledgment
4. **Total time reasonable**: 4+ hours for comprehensive review with full QA pipeline
5. **Intermediate outputs preserved**: User can access any stage's output for reference or modification
