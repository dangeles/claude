# Fact-Check Tiers Reference

## Three-Tier Verification System

### Tier 1: Inline Verification (Always Active)

- **When**: Automatically, for every factual claim during Stages 3 and 4
- **Who**: essay-fact-checker sub-agent via Task tool
- **Process**: WebSearch for claim, find authoritative source, confirm via WebFetch
- **Output**: Verified/Unverified/Partial + source URL
- **Standard**: Every factual claim in the final essay MUST have a source URL

### Tier 2: Proactive Research Enrichment (Stage 3, On Request)

- **When**: During argument development, when the orchestrator identifies an opportunity to strengthen an argument
- **Who**: essay-fact-checker sub-agent via Task tool
- **Process**: WebSearch for relevant data, statistics, examples
- **Output**: Suggestions with source URLs, clearly labeled as "enrichment, not required"
- **Standard**: User decides whether to incorporate

### Tier 3: Deep Research Escalation (Rare, User-Initiated)

- **When**: Topic requires literature-level research beyond WebSearch
- **Who**: Escalated to literature-researcher skill (if available)
- **Process**: Orchestrator reports the need; user decides whether to escalate
- **Output**: Research summary with citations (from literature-researcher)
- **Standard**: This is an escape hatch, not a standard pipeline step

## Deferred Verification Queue

When a claim cannot be verified immediately (timeout, service unavailable):

1. **Mark as DEFERRED** in the argument map or paragraph log
2. **Add to session state queue**:
   ```yaml
   deferred_verifications:
     - claim_text: "[exact claim text]"
       section: N
       paragraph: M  # null if during Stage 3
       reason: "[timeout | service_unavailable | no_results]"
       added_at: "{ISO8601}"
   ```
3. **Continue pipeline** -- do not block on deferred claims
4. **Resolve before final approval** -- Quality Gate G5 requires ALL deferred claims to be resolved
5. **Resolution options**: Retry verification, user provides source, revise claim, remove claim

## User Override Protocol

When the fact-checker's findings conflict with the user's intended claim:

### Step 1: Present Evidence
Show the user what the fact-checker found, with sources:
- "The fact-checker found that [fact-checker result]. Your claim states [user's claim]."
- Present the source URL and relevant excerpt

### Step 2: Offer Options
- **Revise claim**: Update the claim to match the verified evidence
- **Rephrase**: Keep the spirit but adjust the specifics (e.g., "approximately" instead of exact number)
- **Keep with override**: User insists on original claim; log the override
- **Provide own source**: User supplies a source the fact-checker did not find

### Step 3: Log All Overrides
Every override is recorded in `fact-check-log.md`:
```markdown
## Override: [claim text]
- **Fact-checker result**: [verified/unverified/partial/contested]
- **Fact-checker source**: [URL]
- **User action**: Override -- keeping original claim
- **User rationale**: "[user's explanation]"
- **Logged at**: [timestamp]
```

### Step 4: Final Review
At Quality Gate G5, ALL overrides are surfaced explicitly:
- "You overrode [N] fact-checker findings. Here they are: [list]. Do you want to revisit any before finalizing?"

## Common Knowledge Exemption Policy

The following types of claims do NOT require citation:

| Category | Example | Action |
|----------|---------|--------|
| Basic scientific facts | "DNA has a double helix structure" | Flag as `common_knowledge` |
| Standard definitions | "CRISPR stands for Clustered Regularly Interspaced Short Palindromic Repeats" | Flag as `common_knowledge` |
| Historical events with broad consensus | "CRISPR was first used for gene editing in human cells in 2013" | Flag as `common_knowledge` |
| Mathematical/logical truths | "Exponential growth means doubling at regular intervals" | Flag as `common_knowledge` |

**When in doubt, verify.** It is better to have an unnecessary citation than a missing one.

## Contradictory Sources Protocol

When authoritative sources disagree on a factual claim:

1. **Report both sides**: Present Source A and Source B with their respective claims
2. **Mark as CONTESTED**: The claim status becomes `contested`
3. **Assess the disagreement**:
   - Is one source more recent?
   - Is one a primary source and the other secondary?
   - Do they measure different things that appear to conflict?
   - Is this a genuine scientific disagreement?
4. **Offer framing options** to the user:
   - Acknowledge the range: "Estimates range from X (Source A) to Y (Source B)"
   - Present the debate: "While Source A found X, this has been disputed by Source B, which found Y"
   - Choose a side with acknowledgment: "According to the more recent Source B, Y -- though earlier estimates suggested X"
5. **User decides**: The orchestrator presents options and the user chooses

## Claim Classification Guide

Before requesting verification, classify each claim:

| Category | Definition | Action | Example |
|----------|-----------|--------|---------|
| **Factual** | Objective claim that can be verified against evidence | Verify (Tier 1) | "Off-target rates dropped 50-fold between 2015-2024" |
| **Interpretive** | Analytical or interpretive claim based on judgment | Skip verification | "This represents a paradigm shift in gene therapy" |
| **Personal** | Personal experience, anecdote, or opinion | Exempt | "In my experience working with these systems..." |
| **Policy** | Value judgment or policy recommendation | Exempt | "We should increase funding for delivery research" |
| **Common knowledge** | Widely established fact | Flag, optional citation | "DNA encodes genetic information" |

**Rule of thumb**: If a claim uses specific numbers, dates, names, or "studies show" language, it is factual and must be verified.

## Source Quality Hierarchy

When multiple sources are available, prefer higher-quality sources:

1. **Primary research**: Peer-reviewed journal articles with DOIs
   - Top tier: Nature, Science, Cell, NEJM, Lancet
   - Strong: PNAS, eLife, PLoS Biology, domain-specific top journals
2. **Systematic reviews and meta-analyses**: Especially Cochrane reviews
3. **Official statistics**: Government agencies (NIH, WHO, CDC, FDA, NSF)
4. **Preprints**: bioRxiv, medRxiv (flag as "not yet peer-reviewed")
5. **Science journalism**: Nature News, Science News, Quanta Magazine (secondary sources)
6. **General news**: Major outlets reporting on science (lower confidence)
7. **Other**: Blog posts, press releases, Wikipedia (flag as "low-confidence source")

**Always prefer the primary source over secondary reporting.** If a news article cites a study, find and cite the study directly.
