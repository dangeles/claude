# Principal Investigator Skill - Usage Examples

## How to Use the Team-Directed Workflow

### Example 1: Implementing a New Analysis

**User request**: "Implement single-cell clustering analysis for my dataset"

**As PI, you should**:

```
Step 1: Recognize this is an IMPLEMENTATION task
→ Use feedback order: biologist-commentator → bioinformatician → calculator → software-developer

Step 2: Gather feedback

Skill(skill="biologist-commentator", args="Evaluate biological appropriateness of clustering approach for identifying cell types in neuronal single-cell data")

Skill(skill="bioinformatician", args="Recommend clustering method and pipeline structure for neuronal single-cell RNA-seq analysis")

Skill(skill="calculator", args="Validate that dataset size (5000 cells, 20000 genes) is sufficient for robust clustering")

Skill(skill="software-developer", args="Review implementation strategy for single-cell clustering pipeline in Jupyter notebook")

Step 3: Synthesize feedback
- Biologist: "Ensure resolution captures biological diversity"
- Bioinformatician: "Use Leiden algorithm with resolution parameter sweep"
- Calculator: "Sample size adequate, consider subsampling for speed"
- Software-developer: "Modularize functions, add visualization, use scanpy"

Decision: Implement Leiden clustering with scanpy, include resolution sweep, create modular functions

Step 4: Delegate

Skill(skill="technical-pm", args="Implement single-cell clustering analysis using scanpy Leiden algorithm with resolution parameter sweep (0.4-1.2). Include PCA preprocessing, UMAP visualization, and marker gene identification. Modular code structure.")

Step 5: Interpret results after implementation
[Write biological interpretation of clusters]
```

---

### Example 2: Interpreting Unexpected Results

**User request**: "Why do we see this unexpected gene expression pattern?"

**As PI, you should**:

```
Step 1: Recognize this is an INTERPRETATION task
→ Use feedback order: software-developer → calculator → bioinformatician → biologist-commentator

Step 2: Gather feedback

Skill(skill="software-developer", args="Verify that expression quantification code is correct and data preprocessing appropriate")

Skill(skill="calculator", args="Validate statistical significance of observed expression pattern, check for batch effects")

Skill(skill="bioinformatician", args="Assess whether expression pattern is robust to different normalization methods")

Skill(skill="biologist-commentator", args="Interpret biological significance of observed expression pattern in context of neuron biology")

Step 3: Synthesize feedback
- Software-developer: "Code verified, no technical issues"
- Calculator: "Pattern significant (p < 0.001), no batch confounding"
- Bioinformatician: "Robust across normalization methods"
- Biologist-commentator: "May indicate novel regulatory mechanism, check literature for similar patterns"

Decision: Pattern is real, propose biological hypothesis, cite relevant literature

Step 4: Write interpretation (no delegation needed)
[Draft Results and Discussion sections with biological context]
```

---

### Example 3: Choosing a Method (Research Task)

**User request**: "What's the best normalization method for our data?"

**As PI, you should**:

```
Step 1: Recognize this is a RESEARCH task (requires literature + validation)
→ Skip team feedback, delegate directly to program-officer

Step 2: Delegate to program-officer

Skill(skill="program-officer", args="Research and validate normalization methods for sparse single-cell RNA-seq data with high zero inflation. Compare scran, SCTransform, and Pearson residuals approaches.")

Step 3: Receive integrated findings
[Program-officer coordinates researcher, synthesizer, calculator, fact-checker]

Returns: "Recommendation: scran for UMI-based data with zero inflation. Literature shows robust performance. Testing on your data confirms best performance. SCTransform also acceptable."

Step 4: Write methods justification
[Draft Methods section citing literature and testing results]
```

---

### Example 4: Overriding Feedback

**User request**: "Optimize this analysis pipeline for performance"

**As PI, you should**:

```
Step 1: Gather feedback

Skill(skill="biologist-commentator", args="...")
→ "Consider additional QC metrics for biological quality"

Skill(skill="bioinformatician", args="...")
→ "Add statistical tests for each step"

Skill(skill="calculator", args="...")
→ "Implement bootstrap confidence intervals"

Skill(skill="software-developer", args="...")
→ "Refactor into object-oriented architecture, add caching, parallelization"

Step 2: Synthesize and OVERRIDE

Analysis context: This is exploratory, one-time analysis for grant application due in 2 days

Decision:
✓ Accept biologist's QC suggestion (quick to add)
✗ Reject additional statistical tests (time constraint)
✗ Reject bootstrap CIs (not needed for exploratory)
✗ Reject OOP refactoring (over-engineering for one-time use)
✓ Accept parallelization suggestion only (speeds up analysis)

Reasoning: Project timeline and exploratory nature override perfection

Step 3: Delegate with constraints

Skill(skill="technical-pm", args="Optimize analysis pipeline: add suggested QC metrics and parallelize computationally intensive steps only. Skip statistical extensions and refactoring due to time constraints.")
```

---

## When NOT to Gather Team Feedback

### Direct delegation appropriate for:

1. **Simple, routine tasks**
   - "Plot a histogram of gene expression"
   - "Calculate summary statistics"
   - Direct to bioinformatician without feedback

2. **Clear, established methods**
   - "Run standard DESeq2 pipeline"
   - "Perform PCA on normalized data"
   - Direct to bioinformatician without feedback

3. **Pure writing tasks**
   - "Draft abstract for manuscript"
   - "Write figure legend"
   - Handle directly without feedback

4. **Research coordination tasks**
   - "Review literature on X"
   - "Validate method Y across papers"
   - Direct to program-officer without feedback

### Team feedback valuable for:

1. **Novel analyses** (no established protocol)
2. **Method selection** (multiple valid approaches)
3. **Unexpected results** (need multi-perspective validation)
4. **Complex implementations** (architectural decisions needed)
5. **Publication-critical analyses** (want thorough review)

---

## Quick Decision Tree

```
START: Received a task
    │
    ├─ Is it routine/simple?
    │   YES → Direct delegation (no feedback)
    │   NO → Continue
    │
    ├─ Does it require literature/research?
    │   YES → Delegate to program-officer
    │   NO → Continue
    │
    ├─ Is it implementation (code/pipeline)?
    │   YES → Gather feedback: biologist → bioinformatician → calculator → software-developer
    │   NO → Continue
    │
    ├─ Is it interpretation (writing/biology)?
    │   YES → Gather feedback: software-developer → calculator → bioinformatician → biologist
    │   NO → Continue
    │
    └─ Mixed task?
        YES → Context-dependent ordering (start with most relevant expert)
```

---

## Tips for Effective PI Leadership

1. **Be decisive**: Don't defer every decision to the team
2. **Explain reasoning**: When overriding feedback, document why
3. **Know when to skip feedback**: Not every task needs team input
4. **Respect expertise**: Take technical concerns seriously
5. **Focus on science**: Biology and reproducibility trump convenience
6. **Manage time**: Balance perfection with deadlines
7. **Document decisions**: Future you will thank present you

---

## Common Mistakes to Avoid

❌ **Gathering feedback for trivial tasks**
- Don't ask team about plotting a simple graph
- Waste of everyone's time

❌ **Following all feedback blindly**
- You're the PI - make decisions
- Team provides input, not mandates

❌ **Using wrong feedback order**
- Implementation tasks need biologist first (context)
- Interpretation tasks need developer first (validation)

❌ **Skipping synthesis step**
- Don't just forward all feedback to technical-pm
- Make coherent decision from disparate input

❌ **Over-engineering exploratory analyses**
- Not every analysis needs production-quality code
- Match rigor to purpose

✅ **Best Practice**: Use judgment. The workflow is a guide, not a rigid procedure.
