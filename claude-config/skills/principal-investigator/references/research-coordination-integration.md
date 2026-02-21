# Research Coordination Integration Guide

How principal-investigator skill integrates with research coordination skills (program-officer, researcher, calculator, synthesizer, fact-checker).

## Two-Tier Architecture

### Tier 1: Domain Execution (Bioinformatics Skills)
- **principal-investigator**: Scientific writing and interpretation
- **bioinformatician**: Data analysis implementation
- **copilot**: Code review
- **systems-architect**: Software architecture
- **software-developer**: Production code
- **biologist-commentator**: Biological validation

**Location**: `repos/gpcr_exploration/.claude/skills/`

### Tier 2: Research Coordination (General Skills)
- **program-officer**: Multi-agent coordination
- **researcher**: Literature review
- **calculator**: Quantitative analysis
- **synthesizer**: Cross-source integration
- **fact-checker**: Validation

**Location**: `~/.claude/skills/`

## When to Use Which Tier

### Use Tier 1 Only (Direct PI → Bioinformatician)

**Scenario**: Straightforward analysis with clear methods

**Example**: "Analyze differential expression using DESeq2"

**Flow**:
```
Principal Investigator
    ↓ Writes analysis plan
Bioinformatician
    ↓ Implements
Principal Investigator
    ↓ Interprets results
```

**Characteristics**:
- Established methods (DESeq2, scanpy, standard protocols)
- Clear implementation path
- No need for literature synthesis or validation
- Single-domain expertise sufficient

### Use Both Tiers (PI → Program Officer → Specialists)

**Scenario**: Complex task requiring research coordination

**Example**: "Determine best normalization method for low-expression genes"

**Flow**:
```
Principal Investigator
    ↓ Delegates to program-officer
Program Officer
    ├─ Coordinates researcher (literature review)
    ├─ Coordinates synthesizer (compare methods)
    ├─ Coordinates calculator (quantitative testing)
    └─ Coordinates fact-checker (validate claims)
    ↓ Returns integrated findings
Principal Investigator
    ↓ Interprets and writes methods
Bioinformatician
    ↓ Implements validated approach
Principal Investigator
    ↓ Interprets results
```

**Characteristics**:
- Multiple approaches possible
- Needs literature support
- Requires quantitative validation
- Multi-source integration needed

## Integration Patterns

### Pattern 1: Literature-Informed Analysis

**Goal**: Choose analysis method based on literature

**Full workflow**:

```
Step 1: PI Assessment
────────────────────────
User asks: "What's the best clustering algorithm for single-cell data?"

PI identifies:
- Question requires literature review
- Multiple methods exist (Louvain, Leiden, hierarchical)
- Need quantitative comparison
- Requires synthesis of recommendations

PI action: Invoke program-officer

Step 2: Program Officer Coordination
──────────────────────────────────
Program Officer receives delegation and coordinates:

→ Researcher task: "Review recent papers on single-cell clustering algorithms"
  Returns: List of papers with methods (Louvain, Leiden, hierarchical)

→ Synthesizer task: "Compare clustering algorithms from literature"
  Returns: "Leiden outperforms Louvain for resolution parameter tuning,
           hierarchical good for exploratory but computationally expensive"

→ Calculator task: "Test Leiden vs Louvain on sample dataset"
  Returns: "Leiden gives 12% more stable clusters across resolutions,
           runtime comparable"

→ Fact-Checker task: "Verify performance claims from papers on our data"
  Returns: "Claims verified - Leiden resolution stability confirmed"

Step 3: Program Officer Delivers to PI
────────────────────────────────────
Integrated findings:
- Method recommendation: Leiden algorithm
- Literature support: Traag et al. (2019), preferred in >80% recent papers
- Quantitative validation: Tested on sample data, 12% improvement
- Confidence: High

Step 4: PI Interpretation
─────────────────────────
PI writes methods section:

"We used the Leiden algorithm (Traag et al., 2019) for community detection
based on its superior performance in resolution parameter tuning compared to
the Louvain algorithm (Blondel et al., 2008). We tested both methods on a
subset of our data and confirmed that Leiden produces more stable cluster
assignments across resolutions (12% improvement in stability metric), consistent
with benchmarks in the literature."

Step 5: Implementation
──────────────────────
PI delegates to bioinformatician:
"Implement Leiden clustering with resolution=0.8"

Bioinformatician implements validated approach.
```

**When to use this pattern**:
- Choosing between multiple valid methods
- Need to justify method selection with literature
- Want quantitative validation before committing to approach

### Pattern 2: Quantitative Feasibility Check

**Goal**: Validate statistical approach before large-scale analysis

**Full workflow**:

```
Step 1: PI Assessment
────────────────────────
User proposes: "Use mixed-effects model for batch correction"

PI questions:
- Is this appropriate for our data structure?
- Do we have adequate power?
- Are assumptions met?

PI action: Invoke program-officer

Step 2: Program Officer Coordination
──────────────────────────────────
Program Officer coordinates validation:

→ Calculator task: "Run power analysis for mixed-effects model"
  Returns: "With n=50 samples, 4 batches, adequate power (0.85) to detect
           medium effects"

→ Calculator task: "Check mixed-effects model assumptions on sample data"
  Returns: "Residuals approximately normal, variance homogeneous across batches,
           no severe outliers"

→ Researcher task: "Find papers using mixed-effects for batch correction in similar data"
  Returns: "Used successfully in Smith et al. (2023) for multi-batch RNA-seq,
           Patel et al. (2024) for similar experimental design"

→ Fact-Checker task: "Verify our data meets model requirements"
  Returns: "Requirements met: balanced design, sufficient replication,
           batch effects present but not confounded with treatment"

Step 3: Program Officer Delivers to PI
────────────────────────────────────
Validation report:
- Statistical validity: PASS (assumptions met, adequate power)
- Literature support: Used in similar studies (2 recent papers)
- Recommendation: Proceed with mixed-effects model
- Alternative if failed: limma-voom with batch as covariate
- Confidence: High

Step 4: PI Interpretation
─────────────────────────
PI writes methods section:

"To account for batch effects while preserving biological variation, we used
linear mixed-effects models (lme4 package) with batch as a random effect.
Power analysis indicated adequate statistical power (0.85) for our sample size
(n=50), and model diagnostics confirmed that assumptions were satisfied.
This approach has been used successfully for similar multi-batch experimental
designs (Smith et al., 2023; Patel et al., 2024)."

Step 5: Implementation
──────────────────────
PI delegates to bioinformatician with confidence:
"Implement mixed-effects batch correction as validated"
```

**When to use this pattern**:
- Proposing non-standard statistical approach
- Need to validate assumptions before large-scale analysis
- Want to avoid wasted computational time on invalid methods
- Planning to publish and need statistical rigor

### Pattern 3: Multi-Source Validation

**Goal**: Verify unexpected finding across multiple sources

**Full workflow**:

```
Step 1: PI Assessment
────────────────────────
Bioinformatician reports: "Gene X shows 10-fold upregulation in condition A,
                          but literature says it's constitutively expressed"

PI questions:
- Is this finding real or artifact?
- Technical issue (normalization, batch effect)?
- Biological phenomenon not previously reported?

PI action: Invoke program-officer

Step 2: Program Officer Coordination
──────────────────────────────────
Program Officer coordinates investigation:

→ Researcher task: "Review literature on Gene X expression patterns"
  Returns: "Typically constitutively expressed in most tissues, but 2 papers
           report condition-specific regulation in stress conditions (Jones 2022,
           Kim 2023)"

→ Calculator task: "Test alternative normalizations for Gene X"
  Returns: "Pattern robust across 3 normalization methods (TMM, DESeq2, scran),
           not driven by normalization choice"

→ Fact-Checker task: "Verify Gene X preprocessing - any filtering/imputation issues?"
  Returns: "Gene X passes QC: >100 reads/sample, detected in all samples,
           no zeros or outliers, not flagged by QC pipeline"

→ Calculator task: "Check if batch effect could explain Gene X pattern"
  Returns: "Pattern persists within each batch independently, not confounded
           with batch structure"

→ Synthesizer task: "Integrate evidence - is Gene X finding credible?"
  Returns: "High confidence - robust to methods, QC passed, some literature
           support for condition-specific regulation in stress"

Step 3: Program Officer Delivers to PI
────────────────────────────────────
Validation report:
- Finding validity: REAL (not artifact)
- Technical validation: Robust to normalization, passes QC, not batch effect
- Literature context: Rare but precedented (stress-responsive regulation)
- Biological plausibility: Consistent with stress response pathway
- Confidence: High
- Recommendation: Report as novel finding with appropriate caveats

Step 4: PI Interpretation
─────────────────────────
PI writes results section:

"We observed significant upregulation of Gene X in condition A (log2FC=3.3,
padj<0.001), a finding that was robust across multiple normalization methods
and not attributable to batch effects. While Gene X is typically constitutively
expressed, recent studies have reported condition-specific regulation in stress
responses (Jones et al., 2022; Kim et al., 2023), suggesting that our observation
may reflect a previously under-appreciated stress-responsive function."

PI writes discussion:
"The unexpected regulation of Gene X warrants further investigation, particularly
given its canonical role as a housekeeping gene. Our finding adds to emerging
evidence that 'housekeeping' genes may exhibit context-dependent regulation..."
```

**When to use this pattern**:
- Unexpected results that contradict literature
- Need to rule out technical artifacts
- Want to assess biological plausibility
- Planning to report novel findings and need evidence

## Decision Tree

**Question**: Do I need program-officer or can I proceed directly to bioinformatician?

```
┌─────────────────────────────────────────────────────┐
│ Is task straightforward with established methods?  │
└─────────────────┬───────────────────────────────────┘
                  │
        ┌─────────┴─────────┐
        │                   │
       YES                 NO
        │                   │
        ↓                   ↓
┌───────────────┐   ┌────────────────────────┐
│ Use PI →      │   │ Ask these questions:   │
│ Bioinformatician│ └────────┬───────────────┘
└───────────────┘            │
                             ↓
        ┌────────────────────────────────────────────┐
        │ Need literature review of multiple papers? │
        └────────┬───────────────────────────────────┘
                 │
        ┌────────┴────────┐
       YES               NO
        │                 │
        ↓                 ↓
  Invoke program-officer    Continue...
        │
        ↓
┌──────────────────────────────────────────────┐
│ Need quantitative feasibility check?         │
└────────┬─────────────────────────────────────┘
         │
    ┌────┴────┐
   YES       NO
    │         │
    ↓         ↓
Invoke    Continue...
program-officer
    │
    ↓
┌──────────────────────────────────────────────┐
│ Need validation across multiple sources?    │
└────────┬─────────────────────────────────────┘
         │
    ┌────┴────┐
   YES       NO
    │         │
    ↓         ↓
Invoke    Continue...
program-officer
    │
    ↓
┌──────────────────────────────────────────────┐
│ Multiple specialists with dependencies?      │
└────────┬─────────────────────────────────────┘
         │
    ┌────┴────┐
   YES       NO
    │         │
    ↓         ↓
Invoke    Use PI → Bioinformatician
program-officer
```

**Summary decision rule**:

- **Straightforward task** → PI → Bioinformatician
- **Any of these apply** → PI → Program Officer → Specialists:
  - Literature synthesis needed
  - Quantitative validation needed
  - Multi-source verification needed
  - Complex coordination with dependencies

## Examples by Research Phase

### Planning Phase

**Use program-officer when**:

**Scenario 1**: Designing experimental approach
```
User: "We want to identify cell types in our single-cell RNA-seq data"

PI assessment: Multiple methods exist, need literature-informed choice

PI invokes program-officer:
- Researcher: Review cell type identification methods
- Synthesizer: Compare marker-based vs reference-based vs unsupervised
- Calculator: Estimate computational requirements for each
- Fact-Checker: Verify methods appropriate for our organism/tissue

Result: Validated approach with justified method selection
```

**Scenario 2**: Choosing between methods
```
User: "Should we use bulk or single-cell RNA-seq?"

PI assessment: Requires literature review + cost-benefit analysis

PI invokes program-officer:
- Researcher: Review papers comparing bulk vs single-cell for similar questions
- Calculator: Estimate costs (sequencing, compute, time)
- Synthesizer: Compare resolution, power, artifacts, analysis complexity
- Fact-Checker: Verify budget and timeline constraints

Result: Recommendation with trade-offs clearly articulated
```

**Scenario 3**: Estimating sample size
```
User: "How many samples do we need for differential expression analysis?"

PI assessment: Need power analysis + literature benchmarks

PI invokes program-officer:
- Calculator: Run power analysis for expected effect sizes
- Researcher: Find similar studies and their sample sizes
- Fact-Checker: Verify assumptions (variance, effect size realistic)
- Synthesizer: Integrate power analysis + literature + practical constraints

Result: Sample size recommendation with statistical justification
```

### Analysis Phase

**Use direct PI → Bioinformatician when**:

**Scenario 1**: Implementing established methods
```
User: "Run DESeq2 on our RNA-seq data"

PI writes plan: Standard DESeq2 workflow with QC, filtering, normalization, testing

Bioinformatician implements: Follows established protocol

Result: Standard analysis, no coordination needed
```

**Scenario 2**: Following published protocols
```
User: "Reproduce analysis from Smith et al. paper"

PI writes plan: Follow methods from paper section

Bioinformatician implements: Uses same tools/parameters as paper

Result: Straightforward replication, no coordination needed
```

**Scenario 3**: Running standard QC
```
User: "Check quality of sequencing data"

PI writes plan: Standard QC metrics (read depth, duplication, etc.)

Bioinformatician implements: Runs FastQC/MultiQC

Result: Routine QC, no coordination needed
```

### Interpretation Phase

**Use program-officer when**:

**Scenario 1**: Unexpected results need validation
```
Bioinformatician reports: "Housekeeping gene shows differential expression"

PI assessment: Need to verify not artifact, check literature

PI invokes program-officer:
- Calculator: Test alternative methods
- Fact-Checker: Verify QC for this gene
- Researcher: Check if reported in literature
- Synthesizer: Integrate evidence

Result: Validated finding or identified artifact
```

**Scenario 2**: Cross-study comparisons
```
User: "How do our results compare to published studies?"

PI assessment: Need synthesis across multiple papers

PI invokes program-officer:
- Researcher: Extract results from comparable papers
- Synthesizer: Compare findings (overlap, differences)
- Fact-Checker: Verify comparable methods/conditions
- Calculator: Quantify overlap (Fisher's exact test)

Result: Contextualized findings within literature
```

**Scenario 3**: Quantifying biological significance
```
User: "Are these log2 fold changes biologically meaningful?"

PI assessment: Need literature benchmarks + calculations

PI invokes program-officer:
- Researcher: Find typical effect sizes in similar studies
- Calculator: Convert to biological units (fold change → protein abundance)
- Fact-Checker: Verify calculation assumptions
- Synthesizer: Interpret magnitude in biological context

Result: Effect sizes interpreted with biological context
```

### Writing Phase

**Use program-officer when**:

**Scenario 1**: Literature synthesis for introduction
```
PI needs: Comprehensive review of field for introduction

PI invokes program-officer:
- Researcher: Read key papers and extract themes
- Synthesizer: Identify knowledge gaps and organize narrative
- Fact-Checker: Verify citations and claims
- Calculator: Quantify trends if applicable (meta-analysis)

Result: Literature synthesis ready for PI to write introduction
```

**Scenario 2**: Validating methods description
```
PI drafts methods, needs verification

PI invokes program-officer:
- Fact-Checker: Verify all methods accurately described
- Researcher: Check if methods align with field standards
- Calculator: Verify statistical tests correctly reported

Result: Validated methods section ready for submission
```

**Scenario 3**: Checking statistical reporting
```
PI needs to ensure all stats reported correctly

PI invokes program-officer:
- Calculator: Verify all p-values, effect sizes, CIs reported
- Fact-Checker: Check test assumptions stated
- Researcher: Verify statistical reporting follows journal guidelines

Result: Statistics section ready for submission
```

## Key Principles

### 1. Separation of Concerns

**Principal-Investigator** = Scientific brain
- Frames research questions
- Interprets biological significance
- Writes publication-quality narrative
- Makes scientific judgment calls

**Program Officer** = Research coordinator
- Manages information gathering
- Coordinates multiple specialists
- Handles logistics and dependencies
- Delivers integrated findings

**They work together, not in competition.**

### 2. When to Delegate

**Delegate to program-officer when**:
- Task requires input from 2+ specialists
- Literature synthesis needed
- Quantitative validation needed
- Multi-source verification needed

**Don't delegate when**:
- Straightforward implementation
- Established methods
- Single domain of expertise
- Time-sensitive quick analysis

### 3. What Program Officer Delivers

**Deliverables to PI**:
- Synthesis documents (integrated findings)
- Validation reports (claims verified/refuted)
- Computational results (quantitative analyses)
- Literature notes (papers reviewed and summarized)

**What PI does with deliverables**:
- Interprets biological/scientific significance
- Writes publication-quality narrative
- Frames conclusions in research context
- Identifies next steps

### 4. Integration is Seamless

**From user's perspective**:
- Ask principal-investigator to do complex task
- PI automatically delegates to program-officer when needed
- User receives final interpreted results

**User doesn't need to**:
- Decide when to use program-officer
- Coordinate researchers/calculators themselves
- Integrate findings from multiple specialists

**PI handles coordination decisions.**

## Common Questions

### Q: When should I use program-officer vs researcher directly?

**A**: Never invoke researcher directly if you're doing scientific research. Always go through principal-investigator, who will delegate to program-officer if needed.

**Correct workflow**:
```
User → Principal-Investigator → [Program Officer → Researcher]
```

**Incorrect workflow**:
```
User → Researcher  ❌ (bypasses scientific interpretation)
```

### Q: Can I use program-officer for simple literature reviews?

**A**: Yes, but it's often overkill. If you just need to read 1-2 papers, PI can do it directly. If you need to synthesize 5+ papers, delegate to program-officer.

**Simple (PI handles directly)**:
- "What method did Smith et al. use?"
- "Check if this approach has been tried before"

**Complex (delegate to program-officer)**:
- "Compare normalization methods across literature"
- "Synthesize findings from single-cell clustering papers"

### Q: What if I'm unsure whether to delegate?

**A**: Err on the side of delegating. Program Officer can always decide a task is simple and handle it quickly. Better to delegate and have it be quick than to miss needed coordination.

**Rule of thumb**: If you think "I should probably check the literature / validate this / get a quantitative estimate" → delegate to program-officer.

### Q: Can program-officer write the manuscript?

**A**: No. Program Officer delivers findings, PI writes the narrative.

**Program Officer delivers**:
- "Method A preferred in 8/10 papers"
- "Validation passed: assumptions met"
- "Effect size comparable to Smith et al."

**PI writes**:
- "We used Method A based on its widespread adoption in recent studies..."
- "Statistical assumptions were verified prior to analysis..."
- "Our effect sizes were consistent with previous reports..."

### Q: How do I know if program-officer succeeded?

**A**: Program Officer returns integrated findings with:
- Clear recommendation or summary
- Supporting evidence from multiple sources
- Confidence level (high/medium/low)
- Alternative if recommendation failed

**If you receive this, program-officer succeeded. If you receive fragmented information, it may need refinement.**

### Q: Can I use this pattern for other domain skills?

**A**: Yes! This two-tier architecture works for any domain-specific skill that needs research coordination.

**Examples**:
- **Chemistry PI** → program-officer → researchers/calculators for synthesis planning
- **Physics PI** → program-officer → researchers/calculators for experimental design
- **Clinical PI** → program-officer → researchers/fact-checkers for guideline synthesis

**The pattern is domain-agnostic.**

## Summary

**Two-tier architecture**:
- **Tier 1** (domain): Bioinformatics skills execute scientific work
- **Tier 2** (coordination): Research coordination skills manage information

**PI is the bridge**:
- Assesses task complexity
- Delegates to program-officer when needed
- Receives integrated findings
- Interprets and writes

**Program Officer coordinates**:
- Manages researcher, calculator, synthesizer, fact-checker
- Handles dependencies and integration
- Delivers validated findings

**User benefit**:
- Ask PI to do complex research task
- PI automatically coordinates needed specialists
- Receive interpreted, publication-ready results

**The integration is invisible to users but powerful in execution.**

---

Created: 2026-01-28
Author: David Angeles Albores
Version: 1.0
