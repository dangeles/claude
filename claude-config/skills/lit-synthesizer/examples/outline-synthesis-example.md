# Outline Synthesis Example: Hepatocyte Oxygenation

This example demonstrates **Mode 1: Outline Synthesis** for creating section structure from 6-9 foundational review papers.

## Scenario

lit-pm Stage 3 calls lit-synthesizer after Stage 2 (Review Discovery) completes. The task is to synthesize 7 review papers on hepatocyte oxygenation into a coherent 6-section outline.

**Input** (from lit-pm Stage 3):
```yaml
mode: outline_synthesis
task_id: outline-20260204-1700
output_dir: /scratchpad/lit-synthesizer/outline-20260204-1700/
reviews:
  - path: /scratchpad/lit-pm/stage2/review-1.md
    title: "Allen & Bhatia 2021 - Hepatocyte Function and Oxygenation in 3D Culture"
    priority: 95
    convergence: 1.0  # Found by all 3 search strategies
  - path: /scratchpad/lit-pm/stage2/review-2.md
    title: "Jiang et al. 2024 - Oxygen Delivery in Liver Bioreactors"
    priority: 90
    convergence: 0.33
  - path: /scratchpad/lit-pm/stage2/review-3.md
    title: "Miki et al. 2018 - Oxygen Consumption Rates of Primary Hepatocytes"
    priority: 83
    convergence: 0.67  # Found by 2 strategies
  - path: /scratchpad/lit-pm/stage2/review-4.md
    title: "Morgan 2019 - Perfusion Bioreactors for Hepatocyte Culture"
    priority: 82
    convergence: 0.33
  - path: /scratchpad/lit-pm/stage2/review-5.md
    title: "Hay et al. 2022 - Microfluidic Liver-on-Chip Platforms"
    priority: 88
    convergence: 0.33
  - path: /scratchpad/lit-pm/stage2/review-6.md
    title: "Sakai 2020 - Oxygen Regulation in Liver Organoids"
    priority: 80
    convergence: 0.33
  - path: /scratchpad/lit-pm/stage2/review-7.md
    title: "Chen et al. 2017 - Mass Transfer Limitations in Hepatocyte Bioreactors"
    priority: 78
    convergence: 0.33
research_question: "What are the key challenges in hepatocyte oxygenation for bioreactor applications?"
```

---

## Process

### Step 1: Read All Reviews

Read all 7 review papers to understand the landscape.

**Key themes extracted**:

From **Allen & Bhatia 2021** (convergence 1.0 - foundational):
- Hepatocyte oxygen consumption rates: 0.3-0.7 nmol/s/10⁶ cells
- Function-oxygen relationship: Albumin synthesis peaks at 5-10% O₂
- Zonation: Periportal vs. pericentral hepatocytes have different oxygen needs
- 3D culture platforms reviewed: spheroids, organoids, bioprinted tissues

From **Jiang et al. 2024** (most recent):
- Membrane oxygenators: PDMS, silicone
- Perfluorocarbon oxygen carriers
- Clinical translation challenges: scalability, cost, regulatory

From **Miki et al. 2018** (quantitative foundation):
- OCR measurements across 15 studies: 0.2-0.9 nmol/s/10⁶ cells (high variability)
- Culture conditions affect OCR: serum-free medium reduces OCR 30%
- Donor variability: 2-3× variation in OCR

From **Morgan 2019** (bioreactor design):
- Perfusion vs. static culture: 3-5× higher viability with perfusion
- Shear stress limits: <1 dyne/cm² for hepatocytes
- Residence time distribution affects oxygenation uniformity

From **Hay et al. 2022** (microfluidics):
- Liver-on-chip platforms enable oxygen gradients (zonation mimicry)
- Throughput limitations: typically <10⁴ cells per device
- Integration with sensors: real-time oxygen monitoring

From **Sakai 2020** (organoids):
- Organoid core hypoxia: >200 μm diameter leads to necrotic core
- HIF-1α stabilization in organoid cores
- Vascularization strategies to overcome size limits

From **Chen et al. 2017** (mass transfer):
- Krogh cylinder model for oxygen diffusion
- Péclet number analysis: convection vs. diffusion
- Critical cell density: >10⁷ cells/mL creates oxygen limitation

---

### Step 2: Identify Cross-Cutting Themes

**Themes spanning multiple reviews**:

1. **Fundamental biology** (Allen, Miki, Chen)
   - OCR measurements, function-oxygen relationship, zonation

2. **Measurement challenges** (Miki, Chen, Allen)
   - OCR variability, donor differences, assay inconsistencies

3. **Mass transfer engineering** (Chen, Morgan, Jiang)
   - Diffusion limits, convection, membrane oxygenators

4. **Platform-specific solutions** (Morgan, Hay, Sakai, Jiang)
   - Perfusion bioreactors, microfluidics, organoids, carriers

5. **Scale and translation** (Jiang, Morgan, Sakai)
   - Clinical scale-up, cost, regulatory approval

6. **Emerging trends** (Hay, Sakai, Jiang - all 2020+)
   - Zonation mimicry, vascularization, real-time monitoring

---

### Step 3: Organize into Narrative Arc

**Logical progression**:
1. Start with fundamentals (what hepatocytes need)
2. Then measurement challenges (how we quantify needs)
3. Then engineering principles (how to meet needs)
4. Then platform comparisons (specific implementations)
5. Then scale-up challenges (clinical translation)
6. Finally future directions (emerging solutions)

**Rationale**: Build from basic biology → measurement → engineering → implementation → translation → innovation

This avoids jumping between abstraction levels and creates a teaching narrative.

---

### Step 4: Create Outline with Thesis Statements

**Output file**: `/scratchpad/lit-synthesizer/outline-20260204-1700/output.md`

```markdown
# Literature Review Outline: Hepatocyte Oxygenation Challenges in Bioreactors

## Research Question
What are the key challenges in hepatocyte oxygenation for bioreactor applications?

## Narrative Arc
Progresses from fundamental biology (what hepatocytes need) → measurement challenges (how we quantify) → engineering principles (how to deliver) → platform implementations (specific solutions) → clinical translation (scale-up barriers) → future innovations (emerging approaches). This structure builds foundational knowledge before addressing implementation challenges.

---

## Section 1: Hepatocyte Oxygen Physiology and Functional Requirements

**Thesis**: Primary hepatocytes consume oxygen at rates of 0.2-0.9 nmol/s/10⁶ cells, with functional outputs (albumin synthesis, CYP450 activity) critically dependent on oxygen tension in the 5-10% O₂ range, but hepatocyte subpopulations exhibit zonation-specific oxygen demands that complicate uniform culture approaches.

**Key Questions**:
- What are the oxygen consumption rates of primary hepatocytes?
- How does oxygen tension affect hepatocyte function?
- What is the biological basis of zonation, and why does it matter for culture?

**Coverage** (Source Mapping):
- OCR measurements: Miki et al. 2018 (systematic review of 15 studies)
- Function-oxygen relationship: Allen & Bhatia 2021 (comprehensive)
- Zonation: Allen & Bhatia 2021, Sakai 2020
- Donor variability: Miki et al. 2018

**Expected Length**: 2000-2500 words

---

## Section 2: Challenges in Measuring and Predicting Oxygen Consumption

**Thesis**: High variability in reported hepatocyte OCR (2-3× range) stems from assay inconsistencies, donor differences, and culture condition effects, limiting the predictive power of oxygen consumption models for bioreactor design.

**Key Questions**:
- Why do OCR measurements vary across studies?
- How do culture conditions (medium, serum, 2D vs. 3D) affect OCR?
- What are the limitations of current measurement techniques?

**Coverage** (Source Mapping):
- OCR variability: Miki et al. 2018 (primary source)
- Assay methods: Miki et al. 2018, Allen & Bhatia 2021
- Donor effects: Miki et al. 2018
- Culture condition effects: Allen & Bhatia 2021, Sakai 2020

**Expected Length**: 1500-2000 words

---

## Section 3: Mass Transfer Principles and Engineering Solutions

**Thesis**: Oxygen delivery to hepatocytes is fundamentally constrained by diffusion limitations (Krogh cylinder model) at cell densities >10⁷ cells/mL, requiring engineering strategies that balance convective transport (perfusion, carriers) against shear-induced damage.

**Key Questions**:
- What are the fundamental mass transfer limits?
- How do diffusion and convection interact in bioreactors?
- What engineering parameters govern oxygen delivery?

**Coverage** (Source Mapping):
- Mass transfer fundamentals: Chen et al. 2017 (primary source)
- Krogh cylinder model: Chen et al. 2017, Allen & Bhatia 2021
- Critical cell density: Chen et al. 2017, Morgan 2019
- Convection-diffusion balance: Chen et al. 2017, Morgan 2019
- Shear stress limits: Morgan 2019

**Expected Length**: 2000-2500 words

---

## Section 4: Platform Comparison - Perfusion Bioreactors, Microfluidics, and Organoids

**Thesis**: Perfusion bioreactors, microfluidic liver-on-chip platforms, and organoid systems represent three distinct oxygenation strategies—each optimized for different scales and applications but facing inherent trade-offs between throughput, physiological fidelity, and translational feasibility.

**Key Questions**:
- How do different platforms address oxygenation?
- What are the trade-offs between platforms?
- Which platform is best for which application?

**Coverage** (Source Mapping):
- Perfusion bioreactors: Morgan 2019 (primary), Jiang et al. 2024
- Microfluidics: Hay et al. 2022 (primary)
- Organoids: Sakai 2020 (primary)
- Comparison framework: Allen & Bhatia 2021, Jiang et al. 2024

**Subsections**:
- Perfusion bioreactors (scale, uniformity, cost)
- Microfluidic platforms (zonation mimicry, throughput limits)
- Organoid approaches (size limits, vascularization)
- Trade-off matrix (table comparing platforms)

**Expected Length**: 3000-3500 words (longest section due to platform diversity)

---

## Section 5: Membrane Oxygenators and Oxygen Carriers

**Thesis**: Membrane oxygenators (PDMS, silicone) and perfluorocarbon oxygen carriers enable decoupling of oxygen supply from medium flow, but clinical translation faces barriers including material biocompatibility, long-term stability, and regulatory approval for carrier technologies.

**Key Questions**:
- How do membrane oxygenators work?
- What are perfluorocarbon carriers and their advantages?
- What limits clinical adoption of these technologies?

**Coverage** (Source Mapping):
- Membrane oxygenators: Jiang et al. 2024 (primary), Morgan 2019
- PFC carriers: Jiang et al. 2024 (primary)
- Biocompatibility: Jiang et al. 2024, Allen & Bhatia 2021
- Regulatory challenges: Jiang et al. 2024

**Expected Length**: 1500-2000 words

---

## Section 6: Clinical Translation Challenges and Future Directions

**Thesis**: Translating hepatocyte bioreactor technologies from bench to bedside requires addressing scale-up economics (cost per patient), regulatory pathways for novel oxygenation approaches, and emerging innovations in vascularization, real-time monitoring, and zonation engineering.

**Key Questions**:
- What are the economic barriers to clinical scale?
- How do regulatory requirements shape technology choices?
- What emerging technologies might overcome current limitations?

**Coverage** (Source Mapping):
- Clinical translation: Jiang et al. 2024 (primary source - 2024 perspective)
- Scale-up challenges: Morgan 2019, Jiang et al. 2024
- Cost analysis: Jiang et al. 2024
- Regulatory: Jiang et al. 2024
- Vascularization: Sakai 2020, Jiang et al. 2024
- Real-time monitoring: Hay et al. 2022
- Zonation engineering: Hay et al. 2022, Sakai 2020

**Expected Length**: 2000-2500 words

---

## Total Review Length
Approximately 12,000-15,000 words (6 sections)

## Timeline for Section Writing (Stage 5)
- Section 1: 2-3 days (foundational biology)
- Section 2: 1-2 days (measurement challenges)
- Section 3: 2-3 days (mass transfer engineering)
- Section 4: 3-4 days (platform comparison - longest)
- Section 5: 1-2 days (membrane/carrier technologies)
- Section 6: 2-3 days (translation and future)

**Total Stage 5 Duration**: 11-17 days (parallelizable with 6 literature-researcher agents)
```

---

### Step 5: Create Metadata File

**Output file**: `/scratchpad/lit-synthesizer/outline-20260204-1700/metadata.yaml`

```yaml
mode: outline_synthesis
status: complete
task_id: outline-20260204-1700
timestamp: 2026-02-04T17:15:00Z

sections_created: 6
narrative_arc: "Progresses from fundamental biology (what hepatocytes need) → measurement challenges (how we quantify) → engineering principles (how to deliver) → platform implementations (specific solutions) → clinical translation (scale-up barriers) → future innovations (emerging approaches)"

reviews_used: 7
convergent_review: "Allen & Bhatia 2021"  # Convergence 1.0

source_mapping_complete: true
thesis_statements_provided: true
expected_total_length: 12000-15000 words
parallelization_ready: true  # All 6 sections can be written in parallel

ready_for_user_approval: true
next_stage: "Stage 5: Section Writing (6 parallel literature-researcher agents)"
```

---

## Handoff to lit-pm

lit-pm reads metadata.yaml and determines next steps:

1. **User Approval**: Present outline to user
   - "Outline created with 6 sections following narrative arc: biology → measurement → engineering → platforms → translation → future"
   - "Ready to proceed with Section Writing (Stage 5)?"

2. **If Approved**: Launch Stage 5
   - Spawn 6 parallel literature-researcher agents (Mode 2: Deep Targeted Research)
   - Each agent writes one section with 15-30 papers
   - Estimated duration: 11-17 days (parallel execution)

3. **If Revisions Requested**: Return to lit-synthesizer
   - User may request: merging sections, splitting sections, reordering
   - lit-synthesizer revises outline, creates new metadata

---

## Key Design Decisions

### 1. Six Sections (Not Five or Eight)

**Rationale**:
- 7 reviews cover diverse topics → need sufficient sections to avoid cramming
- Each section targets 2000-2500 words → 6 sections = 12K-15K total (typical comprehensive review length)
- Avoids both under-granularity (too few sections, each becomes sprawling) and over-granularity (too many sections, each becomes superficial)

### 2. Section 4 is Longest (Platform Comparison)

**Rationale**:
- Three distinct platforms (perfusion, microfluidics, organoids) each deserve deep coverage
- Trade-off analysis requires detailed comparison
- This is the "meat" of the review where readers find actionable guidance

### 3. Measurement Challenges Before Engineering Solutions

**Rationale**:
- Can't solve engineering problems without accurate measurements
- OCR variability (Section 2) directly impacts bioreactor design (Section 3)
- Teaching narrative: understand the problem before solving it

### 4. Separate Section for Membrane/Carrier Technologies

**Rationale**:
- Could have merged into Section 3 (engineering) or Section 4 (platforms)
- But these are specific technologies with clinical translation implications
- Jiang et al. 2024 devotes significant coverage → warrants dedicated section

### 5. Convergent Review (Allen & Bhatia 2021) Appears in Most Sections

**Source Mapping**:
- Section 1: ✅ (function-oxygen relationship)
- Section 2: ✅ (OCR measurements)
- Section 3: ✅ (Krogh cylinder)
- Section 4: ✅ (platform comparison framework)
- Section 5: ✅ (biocompatibility)
- Section 6: - (not clinical focus)

**Rationale**: Convergence 1.0 indicates foundational review → should appear throughout

---

## Lessons Learned

1. **Cross-cutting themes ≠ sections**: Initially identified 6 themes, but themes don't always map 1:1 to sections. Theme 4 (platform-specific solutions) became 2 sections (4 and 5).

2. **Narrative arc justification is critical**: Simply listing sections isn't enough. Explaining *why* Section 2 comes before Section 3 helps Stage 5 agents understand context.

3. **Source mapping prevents coverage gaps**: Explicitly noting which reviews support which sections ensures no review is orphaned and no section lacks sources.

4. **Expected length guides effort allocation**: Flagging Section 4 as longest (3000-3500 words) helps lit-pm allocate more time for that agent in Stage 5.

5. **Convergence score informs structure**: Allen & Bhatia 2021 (convergence 1.0) should be foundational → appears in 5 of 6 sections. Reviews with convergence 0.33 are more specialized → appear in 1-2 sections.
