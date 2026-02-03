# PubMed Search Syntax Quick Reference

**Last Updated**: 2026-01-28
**Purpose**: Essential search operators and field tags for effective PubMed/NCBI E-utilities queries

---

## Boolean Operators

| Operator | Syntax | Example | Notes |
|----------|--------|---------|-------|
| **AND** | `term1 AND term2` | `hepatocyte AND oxygen consumption` | Both terms must appear (narrows results) |
| **OR** | `term1 OR term2` | `hepatocyte OR liver cell` | Either term can appear (broadens results) |
| **NOT** | `term1 NOT term2` | `hepatocyte NOT cancer` | Excludes second term |
| **Phrase** | `"exact phrase"` | `"bioartificial liver"` | Exact phrase matching (use quotes) |
| **Grouping** | `(term1 OR term2) AND term3` | `(porcine OR swine) AND hepatocyte` | Parentheses control precedence |

**Important**: Boolean operators must be UPPERCASE in PubMed. Lowercase "and"/"or" are treated as regular search terms, not operators.

---

## Field Tags

Append `[field]` to restrict search to specific metadata fields:

| Field Tag | Searches | Example | Use When |
|-----------|----------|---------|----------|
| `[Title]` | Article title only | `bioreactor[Title]` | Finding papers specifically about topic (not just mentioning it) |
| `[Title/Abstract]` or `[TIAB]` | Title or abstract | `oxygen transport[TIAB]` | Broader than title, narrower than all fields |
| `[Author]` | Author names | `Yarmush ML[Author]` | Finding work by specific researcher |
| `[Journal]` | Journal name | `Biotechnol Bioeng[Journal]` | Limiting to specific journals |
| `[MeSH Terms]` | Medical Subject Headings | `Hepatocytes[MeSH]` | Using controlled vocabulary (most precise) |
| `[Publication Date]` or `[PDAT]` | Publication date | `2020:2024[PDAT]` | Date range filtering |
| `[DP]` | Date of publication (year only) | `2023[DP]` | Specific year |
| `[All Fields]` | Default search (all metadata) | `bioreactor[All Fields]` | Default if no tag specified |

**Pro tip**: Use `[TIAB]` for most searches—it's broad enough to capture relevant papers but excludes spurious matches in references/supplementary data.

---

## Advanced Search Patterns

### 1. Date Range Searches

```
hepatocyte oxygen[TIAB] AND 2020:2024[PDAT]
```
Finds papers on hepatocyte oxygen published between 2020-2024.

**Syntax**: Use `YYYY:YYYY` format with `[PDAT]` or `[DP]` tag.

### 2. Author + Topic Search

```
Gerlach JC[Author] AND bioartificial liver[TIAB]
```
Finds Gerlach's papers specifically about bioartificial liver systems.

### 3. MeSH Term + Subheading

```
Hepatocytes[MeSH:NoExp] AND Oxygen Consumption[MeSH]
```
Uses controlled vocabulary for precise concept matching. `NoExp` prevents automatic expansion to child terms.

**When to use MeSH**: When you want comprehensive recall of a specific concept, even if authors use different terminology. MeSH terms are manually assigned by indexers.

### 4. Review Articles Only

```
bioartificial liver[TIAB] AND Review[Publication Type]
```
Limits to review articles (excludes primary research papers).

**Useful publication types**:
- `Review[PT]`
- `Systematic Review[PT]`
- `Meta-Analysis[PT]`
- `Clinical Trial[PT]`
- `Randomized Controlled Trial[PT]`

### 5. Exclude Animal-Only Studies

```
hepatocyte[TIAB] NOT (rats[MeSH] OR mice[MeSH]) AND humans[MeSH]
```
Focuses on human-relevant studies, excludes rodent-only papers.

**Caveat**: This is imperfect—some comparative studies with both rodent and human data will be excluded. Use judiciously.

### 6. High-Impact Journals Only

```
hollow fiber bioreactor[TIAB] AND (
  Nature[Journal] OR Science[Journal] OR Cell[Journal] OR
  Biotechnol Bioeng[Journal] OR Hepatology[Journal]
)
```
Limits to specific high-tier journals. See `journal-tiers.md` for full lists.

### 7. Forward Citation Tracking

PubMed doesn't directly support forward citation searches via syntax. Use:
- **OpenAlex database skill** for programmatic forward/backward citation tracking
- **Google Scholar** for manual forward citation checks
- **Cited by** links in PubMed web interface (when available)

---

## Search Strategy Workflow

**Step 1: Start Broad**
```
bioartificial liver[TIAB]
```
Get a sense of the literature volume (~5000 results? ~500? ~50?).

**Step 2: Add Key Concepts with OR**
```
(bioartificial liver[TIAB] OR BAL[TIAB] OR "liver support system"[TIAB])
```
Capture terminology variations.

**Step 3: Narrow by Recent Years**
```
(bioartificial liver[TIAB] OR BAL[TIAB]) AND 2015:2024[PDAT]
```
Focus on recent work.

**Step 4: Add Topic-Specific Terms**
```
(bioartificial liver[TIAB] OR BAL[TIAB]) AND 2015:2024[PDAT] AND (oxygen[TIAB] OR oxygenation[TIAB])
```
Zero in on specific research question.

**Step 5: Filter to Reviews if Overwhelming**
```
(bioartificial liver[TIAB] OR BAL[TIAB]) AND 2015:2024[PDAT] AND Review[PT]
```
Start with review articles to map the landscape.

---

## Common Pitfalls

### Pitfall 1: Forgetting Quotes for Phrases
❌ `bioartificial liver` → searches for "bioartificial" OR "liver" (millions of hits)
✅ `"bioartificial liver"` → searches for exact phrase (~hundreds of hits)

### Pitfall 2: Lowercase Boolean Operators
❌ `hepatocyte and oxygen` → searches for three separate words
✅ `hepatocyte AND oxygen` → Boolean AND operation

### Pitfall 3: Over-Narrowing with AND
❌ `hepatocyte AND oxygen AND consumption AND rate AND measurement` → may miss relevant papers using different phrasing
✅ `hepatocyte[TIAB] AND (oxygen consumption[TIAB] OR OCR[TIAB])` → captures terminology variations

### Pitfall 4: Ignoring MeSH Terms
**Problem**: Authors use inconsistent terminology. Papers about the same concept may use different keywords.

**Solution**: Check MeSH term assignments:
1. Run initial keyword search
2. Look at MeSH terms on relevant papers
3. Re-run search using those MeSH terms for comprehensive recall

**Example**: Searching `"liver failure"[TIAB]` misses papers indexed as `Liver Failure[MeSH]` that use terms like "hepatic insufficiency" or "liver dysfunction" in title/abstract.

### Pitfall 5: Not Using TIAB for Quantitative Terms
❌ `oxygen consumption` (default all fields) → matches papers citing OCR in references, not primary OCR studies
✅ `oxygen consumption[TIAB]` → finds papers actually about OCR

---

## Quick Examples by Research Question

| Research Question | PubMed Query |
|-------------------|--------------|
| Recent reviews on bioartificial liver | `"bioartificial liver"[TIAB] AND 2018:2024[PDAT] AND Review[PT]` |
| Human hepatocyte oxygen consumption | `hepatocyte[TIAB] AND ("oxygen consumption"[TIAB] OR OCR[TIAB]) AND human[MeSH]` |
| Hollow fiber bioreactor design papers | `"hollow fiber"[TIAB] AND (bioreactor[TIAB] OR "membrane bioreactor"[TIAB])` |
| Papers citing hepatocyte isolation methods | `(hepatocyte[TIAB] AND isolation[TIAB]) OR "hepatocyte isolation"[TIAB]` |
| Clinical trials of liver support devices | `("liver support"[TIAB] OR "bioartificial liver"[TIAB]) AND Clinical Trial[PT]` |
| Work by specific research group | `(Gerlach JC[Author] OR Neuhaus P[Author]) AND liver[TIAB]` |

---

## Using PubMed with Scientific Skills

**In Claude Code**, use these skills for programmatic access:

- `Skill: "scientific-skills:pubmed-database"` - Query PubMed via NCBI E-utilities API
- `Skill: "scientific-skills:openalex-database"` - Cross-reference with OpenAlex for citation networks
- `Skill: "scientific-skills:perplexity-search"` - Quick landscape scan before targeted PubMed queries

**Workflow**:
1. Start with `perplexity-search` to understand terminology and key authors
2. Use `pubmed-database` skill with refined queries from this reference
3. Use `openalex-database` skill for forward/backward citation tracking
4. Acquire PDFs via PMC (`https://pmc.ncbi.nlm.nih.gov/articles/PMCxxxxxxx/pdf/`)

---

## Resources

- **PubMed Help**: https://pubmed.ncbi.nlm.nih.gov/help/
- **MeSH Database**: https://www.ncbi.nlm.nih.gov/mesh
- **E-utilities API Documentation**: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **Advanced Search Builder** (web interface): https://pubmed.ncbi.nlm.nih.gov/advanced/

---

**Pro tip**: Save successful query patterns in your research notes. When you find a search that yields high-quality results, document the exact syntax for reuse on related topics.
