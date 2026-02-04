# Review Discovery Example: Hepatocyte Oxygenation

This example demonstrates **Mode 1: Review Discovery** for a literature review on hepatocyte oxygenation in biore actors.

## Scenario

lit-pm orchestrator launches 3 parallel literature-researcher agents with diverse search strategies to find 6-9 foundational review papers on hepatocyte oxygenation.

**Input** (Agent 1 - Keyword-based):
```yaml
mode: review_discovery
task_id: batch-20260204-1640-keyword
output_dir: scratchpad/literature-researcher/batch-20260204-1640-keyword/
search_strategy: keyword_based
target_reviews: 6-9  # Total across all 3 agents
topic: "hepatocyte oxygenation bioreactors"
```

---

## Agent 1: Keyword-Based Search

### Search Queries Executed
1. `"hepatocyte oxygenation review"`
2. `"liver bioreactor oxygen review"`
3. `"hepatocyte culture oxygen transport survey"`

### Candidate Reviews Found

#### Review 1: Jiang et al. 2024
**Title**: "Oxygen Delivery in Liver Bioreactors: From Bench to Clinical Translation"
**Journal**: Tissue Engineering Part B: Reviews
**Year**: 2024
**DOI**: 10.1089/ten.teb.2023.0234
**Citations**: ~45 (recent, rapidly accumulating)
**Abstract Summary**: Comprehensive review of oxygen delivery methods for hepatocyte cultures, covering membrane oxygenators, perfluorocarbon oxygen carriers, and microfluidic approaches. Discusses mass transfer limitations and clinical translation challenges.

**Priority Scoring**:
- Recency: 30 pts (2024, very recent)
- Citations: 20 pts (moderate for recent paper)
- Journal quality: 20 pts (Tier 1 - top journal in tissue engineering)
- Relevance: 20 pts (high - exactly on topic)
- **Total: 90/100**

#### Review 2: Allen & Bhatia 2021
**Title**: "Hepatocyte Function and Oxygenation in 3D Culture Systems"
**Journal**: Advanced Drug Delivery Reviews
**Year**: 2021
**DOI**: 10.1016/j.addr.2021.05.012
**Citations**: ~180
**Abstract Summary**: Reviews relationship between oxygen levels and hepatocyte function (albumin synthesis, CYP450 activity, viability) across various 3D culture platforms. Discusses zonation and oxygen gradients.

**Priority Scoring**:
- Recency: 25 pts (3 years old)
- Citations: 30 pts (highly cited)
- Journal quality: 20 pts (Tier 1)
- Relevance: 20 pts (high)
- **Total: 95/100**

#### Review 3: Miki et al. 2018
**Title**: "Oxygen Consumption Rates of Primary Human Hepatocytes in Culture: A Systematic Review"
**Journal**: Biotechnology Progress
**Year**: 2018
**DOI**: 10.1002/btpr.2701
**Citations**: ~120
**Abstract Summary**: Systematic review of measured oxygen consumption rates (OCR) for primary hepatocytes across different culture conditions. Provides reference values for bioreactor design.

**Priority Scoring**:
- Recency: 20 pts (6 years old)
- Citations: 28 pts (well-cited)
- Journal quality: 15 pts (Tier 2)
- Relevance: 20 pts (high)
- **Total: 83/100**

### Agent 1 Output

**File**: `scratchpad/literature-researcher/batch-20260204-1640-keyword/reviews.md`

```markdown
# Review Discovery Results (Keyword-Based Strategy)

## Search Strategy
Keyword-based queries focusing on direct terms: "hepatocyte oxygenation", "liver bioreactor oxygen"

## Reviews Found (3)

### 1. Jiang et al. 2024 - Priority: 90
"Oxygen Delivery in Liver Bioreactors: From Bench to Clinical Translation"
- **Journal**: Tissue Engineering Part B: Reviews
- **Year**: 2024 (very recent)
- **DOI**: 10.1089/ten.teb.2023.0234
- **Citations**: ~45
- **Key topics**: Membrane oxygenators, perfluorocarbon carriers, microfluidics, mass transfer

### 2. Allen & Bhatia 2021 - Priority: 95
"Hepatocyte Function and Oxygenation in 3D Culture Systems"
- **Journal**: Advanced Drug Delivery Reviews
- **Year**: 2021
- **DOI**: 10.1016/j.addr.2021.05.012
- **Citations**: ~180
- **Key topics**: 3D platforms, function-oxygen relationship, zonation, gradients

### 3. Miki et al. 2018 - Priority: 83
"Oxygen Consumption Rates of Primary Human Hepatocytes in Culture: A Systematic Review"
- **Journal**: Biotechnology Progress
- **Year**: 2018
- **DOI**: 10.1002/btpr.2701
- **Citations**: ~120
- **Key topics**: OCR measurements, reference values, bioreactor design

## Strategy Notes
Keyword strategy successfully identified the most recent major review (Jiang 2024) and high-impact functional review (Allen & Bhatia 2021). Miki et al. 2018 provides quantitative foundation.
```

**Handoff YAML**:
```yaml
mode: review_discovery
task_id: batch-20260204-1640-keyword
status: complete
search_strategy: keyword_based
reviews_found: 3
output_file: scratchpad/literature-researcher/batch-20260204-1640-keyword/reviews.md
reviews:
  - doi: "10.1089/ten.teb.2023.0234"
    title: "Oxygen Delivery in Liver Bioreactors: From Bench to Clinical Translation"
    priority: 90
  - doi: "10.1016/j.addr.2021.05.012"
    title: "Hepatocyte Function and Oxygenation in 3D Culture Systems"
    priority: 95
  - doi: "10.1002/btpr.2701"
    title: "Oxygen Consumption Rates of Primary Human Hepatocytes in Culture: A Systematic Review"
    priority: 83
```

---

## Agent 2: Author-Based Search (Parallel)

**Input**:
```yaml
mode: review_discovery
task_id: batch-20260204-1640-author
search_strategy: author_based
```

### Leading Researchers Identified
- Sangeeta N. Bhatia (MIT) - pioneer in liver tissue engineering
- Jeffrey R. Morgan (Brown) - bioreactor design
- Takeshi Sakai (Osaka) - hepatocyte physiology

### Search Results
1. **Bhatia review 2021** → Allen & Bhatia 2021 (CONVERGENCE with Agent 1!)
2. **Morgan review 2019** → "Perfusion Bioreactors for Hepatocyte Culture" (new)
3. **Sakai review 2020** → "Oxygen Regulation in Liver Organoids" (new)

**Convergent review**: Allen & Bhatia 2021 appeared in both keyword and author strategies → High confidence

---

## Agent 3: Citation-Based Search (Parallel)

**Input**:
```yaml
mode: review_discovery
task_id: batch-20260204-1640-citation
search_strategy: citation_based
```

### Highly-Cited Papers Identified
Top papers on "hepatocyte oxygenation" by citation count, then check if they're reviews or what reviews cite them.

### Search Results
1. **Allen & Bhatia 2021** (CONVERGENCE - found by all 3 strategies!)
2. **Hay et al. 2022** → "Microfluidic Liver-on-Chip Platforms: Oxygen Considerations" (new)
3. **Miki et al. 2018** (CONVERGENCE with Agent 1!)

---

## Convergence Analysis (Performed by lit-pm)

After all 3 agents complete, lit-pm collects and deduplicates:

### Unique Reviews Found (7 total)

| Review | Priority | Convergence | Strategies |
|--------|----------|-------------|------------|
| Allen & Bhatia 2021 | 95 | **1.0 (3/3)** | keyword, author, citation |
| Jiang et al. 2024 | 90 | 0.33 (1/3) | keyword |
| Miki et al. 2018 | 83 | **0.67 (2/3)** | keyword, citation |
| Morgan 2019 | 82 | 0.33 (1/3) | author |
| Hay et al. 2022 | 88 | 0.33 (1/3) | citation |
| Sakai 2020 | 80 | 0.33 (1/3) | author |
| *[Additional from other strategies]* |

### Final Selection (Top 6 for Stage 3)

**Convergent Reviews** (High Confidence):
1. **Allen & Bhatia 2021** - Priority 95, Convergence 1.0 (ALL strategies) ⭐
2. **Miki et al. 2018** - Priority 83, Convergence 0.67 (keyword + citation)

**High-Priority Single-Strategy**:
3. **Jiang et al. 2024** - Priority 90 (most recent, keyword)
4. **Hay et al. 2022** - Priority 88 (microfluidics focus, citation)

**Complementary Perspectives**:
5. **Morgan 2019** - Priority 82 (bioreactor design, author)
6. **Sakai 2020** - Priority 80 (organoid perspective, author)

### Convergence Insights

**Allen & Bhatia 2021** is the foundational review - appeared in all 3 search strategies. This indicates:
- Broadly recognized as important by the field
- Comprehensive coverage (found via keywords, authors, citations)
- High confidence for Stage 3 outline synthesis

**Miki et al. 2018** appeared in 2/3 strategies (keyword + citation). This suggests:
- Important quantitative reference
- Well-cited and keyword-relevant
- Provides data foundation for other reviews

**Single-strategy reviews** (Jiang 2024, Hay 2022, Morgan 2019, Sakai 2020) may represent:
- Specialized perspectives (e.g., microfluidics, organoids)
- Recent developments (Jiang 2024, Hay 2022)
- Author-specific contributions (Morgan, Sakai)

---

## Outcome

**Stage 2 Complete**: 6 reviews selected for Stage 3 (Outline Synthesis)

lit-pm passes these 6 reviews to lit-synthesizer for outline creation, highlighting:
- Allen & Bhatia 2021 as **primary foundation** (convergence 1.0)
- Jiang et al. 2024 as **most recent major review**
- Miki et al. 2018 as **quantitative reference**

**Metrics**:
- Total candidates found: 7+ across 3 strategies
- Convergent reviews: 2 (Allen & Bhatia, Miki)
- Final selection: 6 reviews
- Recency coverage: 2 reviews from last 3 years
- Duration: ~60 minutes (parallel execution)

---

## Lessons Learned

1. **Convergence validates importance**: Allen & Bhatia 2021 appeared in all 3 strategies → clear consensus
2. **Diverse strategies complement**: Keyword found most recent (Jiang 2024), author found specialized (Morgan, Sakai), citation found highly-cited
3. **DOI tracking essential**: Enabled automatic deduplication across parallel agents
4. **Priority scoring balances factors**: Recent papers (Jiang 2024) compete with highly-cited older papers (Allen 2021)
