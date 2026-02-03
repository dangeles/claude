# Crisis Response Template

Use this template when notifying the user of critical incidents requiring immediate attention.

---

## Template Structure

```
═══════════════════════════════════════════════════════════════
🚨 CRISIS ALERT - [SEVERITY LEVEL]
═══════════════════════════════════════════════════════════════

**Incident Type**: [Data Loss / Security Breach / Critical Failure / System Issue]
**Detected**: [Timestamp]
**Severity**: [P0 - Critical / P1 - Major / P2 - Minor]
**Impact**: [Who/what is affected]

═══════════════════════════════════════════════════════════════
PROBLEM SUMMARY
═══════════════════════════════════════════════════════════════

[Clear, concise description of what went wrong]

**What happened**:
[2-3 sentences describing the incident]

**Current state**:
[What's broken, blocked, or at risk right now]

**Downstream impact**:
[What else is affected? What deadlines are at risk?]

═══════════════════════════════════════════════════════════════
BLAST RADIUS
═══════════════════════════════════════════════════════════════

**Affected deliverables**:
- [Deliverable 1]: [Status/Impact]
- [Deliverable 2]: [Status/Impact]

**Blocked agents/tasks**:
- [Agent/Task 1]: [Why blocked]
- [Agent/Task 2]: [Why blocked]

**Data at risk**:
- [File/Data 1]: [Risk level - Lost / Corrupted / At Risk]
- [File/Data 2]: [Risk level]

**Timeline impact**:
- [Milestone X]: Originally [date], now [delayed/at risk]
- [Milestone Y]: [Impact]

═══════════════════════════════════════════════════════════════
ROOT CAUSE ANALYSIS
═══════════════════════════════════════════════════════════════

**Immediate cause**:
[What directly triggered this incident]

**Contributing factors**:
- [Factor 1]: [How it contributed]
- [Factor 2]: [How it contributed]

**Warning signs missed** (if any):
- [Sign 1]: [When it appeared]
- [Sign 2]: [When it appeared]

═══════════════════════════════════════════════════════════════
RESPONSE OPTIONS
═══════════════════════════════════════════════════════════════

### Option 1: [ROLLBACK] (Recommended for data loss/corruption)

**Action**:
[Specific steps to revert to last known good state]

**Pros**:
- [Benefit 1]
- [Benefit 2]

**Cons**:
- [Trade-off 1]
- [Trade-off 2]

**Timeline**: [How long this takes]
**Risk level**: [LOW / MEDIUM / HIGH]

---

### Option 2: [CONTINUE WITH FIXES]

**Action**:
[Specific steps to fix issue and continue]

**Pros**:
- [Benefit 1]
- [Benefit 2]

**Cons**:
- [Trade-off 1]
- [Trade-off 2]

**Timeline**: [How long this takes]
**Risk level**: [LOW / MEDIUM / HIGH]

---

### Option 3: [ESCALATE/ABORT]

**Action**:
[When to stop work entirely and reassess]

**Pros**:
- [Benefit 1]
- [Benefit 2]

**Cons**:
- [Trade-off 1]
- [Trade-off 2]

**Timeline**: [How long this takes]
**Risk level**: [LOW / MEDIUM / HIGH]

═══════════════════════════════════════════════════════════════
RECOMMENDED ACTION
═══════════════════════════════════════════════════════════════

**My recommendation**: Option [N] - [Name]

**Reasoning**:
[2-3 sentences explaining why this option is best given the circumstances]

**Immediate next steps if approved**:
1. [Step 1]
2. [Step 2]
3. [Step 3]

═══════════════════════════════════════════════════════════════
ARTIFACTS AVAILABLE FOR REVIEW
═══════════════════════════════════════════════════════════════

[List files/logs the user can review to understand the issue]

- `[file-path]` - [Description]
- `[file-path]` - [Description]

═══════════════════════════════════════════════════════════════
Your decision (1/2/3 or provide alternative):
═══════════════════════════════════════════════════════════════
```

---

## Example: Data Loss Incident

```
═══════════════════════════════════════════════════════════════
🚨 CRISIS ALERT - P0 CRITICAL
═══════════════════════════════════════════════════════════════

**Incident Type**: Data Loss
**Detected**: 2026-01-29 15:42:30
**Severity**: P0 - Critical
**Impact**: 3 hours of literature synthesis work lost

═══════════════════════════════════════════════════════════════
PROBLEM SUMMARY
═══════════════════════════════════════════════════════════════

**What happened**:
Synthesizer agent overwrote `docs/literature-review.md` with incomplete draft,
losing all section 3-5 content. Previous version contained 12 pages of synthesized
findings from 8 papers. Current version only has sections 1-2 (4 pages).

**Current state**:
- Complete draft lost (sections 3-5 deleted)
- No backup found in scratchpad directory
- Synthesizer agent continued working on corrupted version for 45 minutes before detection

**Downstream impact**:
- Devil's Advocate review blocked (nothing to review)
- Friday milestone "Complete literature synthesis" now at risk
- User presentation prep depends on this deliverable

═══════════════════════════════════════════════════════════════
BLAST RADIUS
═══════════════════════════════════════════════════════════════

**Affected deliverables**:
- `docs/literature-review.md`: Sections 3-5 lost (8 pages, 3 hours work)
- Devil's Advocate review: Blocked, was scheduled for today
- Final presentation deck: At risk, pulls from literature review

**Blocked agents/tasks**:
- Devil's Advocate: Waiting for complete draft
- Editor: Scheduled for tomorrow, now delayed

**Data at risk**:
- `docs/literature-review.md`: Currently corrupted (sections 3-5 missing)
- WORK-LOG.md: Intact (last backup 2 hours ago)

**Timeline impact**:
- Friday synthesis milestone: Originally Jan 31, now Feb 1 (1-day slip)
- Monday presentation prep: At risk if we can't recover by Saturday

═══════════════════════════════════════════════════════════════
ROOT CAUSE ANALYSIS
═══════════════════════════════════════════════════════════════

**Immediate cause**:
Synthesizer agent used Write tool instead of Edit tool, overwriting entire file
instead of appending section 6.

**Contributing factors**:
- No backup checkpoint created before Synthesizer started work
- WORK-LOG.md not updated after section 3-5 completion (would have captured content)
- Synthesizer instructions didn't specify "use Edit tool" for existing file

**Warning signs missed**:
- 15:20: Synthesizer reported "starting fresh draft"—should have caught this
- 15:35: File size dropped from 14KB to 6KB—went unnoticed

═══════════════════════════════════════════════════════════════
RESPONSE OPTIONS
═══════════════════════════════════════════════════════════════

### Option 1: ROLLBACK TO RESEARCHER OUTPUT (Fastest recovery)

**Action**:
Retrieve sections 3-5 content from Researcher agent's scratchpad notes.
Researcher's extraction files contain all the raw data. Synthesizer redoes
synthesis work (~2 hours) with explicit backup protocol.

**Pros**:
- Can recover today (2-3 hours)
- No information actually lost (source data in scratchpad)
- Still meet Friday milestone with buffer

**Cons**:
- Synthesizer must redo 2 hours of work
- Final prose may differ slightly from lost version

**Timeline**: 2-3 hours to recreate sections 3-5
**Risk level**: LOW

---

### Option 2: RECONSTRUCT FROM MEMORY (Risky)

**Action**:
Ask Synthesizer agent to recreate sections 3-5 from memory (context window).
Agent has full conversation history including discussion of those sections.

**Pros**:
- Faster (1 hour if memory accurate)
- May preserve exact prose from lost version

**Cons**:
- HIGH RISK: Memory may be incomplete or inaccurate
- No verification against source papers
- Could perpetuate any errors from first draft

**Timeline**: 1 hour reconstruction + 30 min verification
**Risk level**: HIGH

---

### Option 3: START FROM SCRATCH (Most thorough)

**Action**:
Synthesizer reads all 8 papers again and creates new synthesis.
Implement checkpoint system: backup after each section.

**Pros**:
- Guaranteed complete and accurate
- Fresh perspective may improve quality
- No risk of propagating old errors

**Cons**:
- Takes 4-5 hours (full resynthesis)
- Misses Friday milestone (slips to Saturday)
- Wastes previous Researcher work

**Timeline**: 4-5 hours
**Risk level**: LOW (quality) but MEDIUM (timeline)

═══════════════════════════════════════════════════════════════
RECOMMENDED ACTION
═══════════════════════════════════════════════════════════════

**My recommendation**: Option 1 - Rollback to Researcher Output

**Reasoning**:
The Researcher's scratchpad contains all extracted data organized by paper.
Sections 3-5 were syntheses of papers 4-8. Those extraction files still exist
in `scratchpad/researcher-liver-support/papers-4-8/`. Synthesizer can rebuild
sections 3-5 in 2-3 hours with proper backup protocol. We still make Friday
milestone with 1-day buffer.

**Immediate next steps if approved**:
1. Create backup: `cp docs/literature-review.md backups/lit-review-v4-partial.md`
2. Archive corrupted state for analysis: `cp -r scratchpad/synthesizer-lit-review/ rollback-archive/2026-01-29-overwrite-incident/`
3. Reassign Synthesizer with instructions:
   - Read researcher extraction files for papers 4-8
   - Recreate sections 3-5 (oxygen transfer, clinical trials, cost analysis)
   - Use Edit tool only (never Write on existing file)
   - Backup after each section: `cp docs/literature-review.md backups/lit-review-v5-section-[N].md`
4. Set 30-minute checkpoints: verify progress every 30 min

═══════════════════════════════════════════════════════════════
ARTIFACTS AVAILABLE FOR REVIEW
═══════════════════════════════════════════════════════════════

- `docs/literature-review.md` - Current corrupted version (sections 1-2 only)
- `scratchpad/researcher-liver-support/papers-4-8/` - Source data for lost sections
- `scratchpad/synthesizer-lit-review/progress.md` - Agent progress log showing overwrite
- `backups/` - Directory (currently empty, need to implement)

═══════════════════════════════════════════════════════════════
Your decision (1/2/3 or provide alternative):
═══════════════════════════════════════════════════════════════
```

---

## Example: Security Breach

```
═══════════════════════════════════════════════════════════════
🚨 CRISIS ALERT - P0 CRITICAL
═══════════════════════════════════════════════════════════════

**Incident Type**: Security Breach
**Detected**: 2026-01-29 14:15:00
**Severity**: P0 - Critical
**Impact**: API key exposed in commit and WORK-LOG.md

═══════════════════════════════════════════════════════════════
PROBLEM SUMMARY
═══════════════════════════════════════════════════════════════

**What happened**:
Calculator agent included OpenAI API key in cost analysis output.
Key was copied into WORK-LOG.md and committed to repository.
Key has been exposed for ~20 minutes.

**Current state**:
- API key visible in git history (commit abc123f)
- Key visible in WORK-LOG.md lines 450-455
- Repository is private, but key should be considered compromised

**Downstream impact**:
- Must regenerate API key immediately
- Must audit API usage for unauthorized calls
- Potential cost exposure if key used maliciously

═══════════════════════════════════════════════════════════════
RESPONSE OPTIONS
═══════════════════════════════════════════════════════════════

### Option 1: IMMEDIATE KEY ROTATION (Recommended)

**Action**:
1. Revoke exposed API key immediately via OpenAI dashboard
2. Generate new API key
3. Remove key from WORK-LOG.md and git history (git filter-branch)
4. Audit API usage logs for unauthorized calls in last 20 minutes

**Pros**:
- Stops exposure immediately
- Audit trail in API logs
- Clean git history

**Cons**:
- Breaks any active processes using old key
- Git history rewrite required

**Timeline**: 15 minutes
**Risk level**: LOW

---

### Option 2: MONITOR THEN ROTATE

**Action**:
Monitor API usage for 1 hour, rotate key if suspicious activity detected.

**Pros**:
- Doesn't break active processes immediately

**Cons**:
- HIGH RISK: Key remains active and exploitable
- May incur costs if compromised

**Timeline**: 60+ minutes
**Risk level**: HIGH - NOT RECOMMENDED

═══════════════════════════════════════════════════════════════
RECOMMENDED ACTION
═══════════════════════════════════════════════════════════════

**My recommendation**: Option 1 - Immediate Key Rotation

**Reasoning**:
Even private repository, any key exposure should trigger immediate rotation.
The cost of rotation (15 min downtime) is far less than potential unauthorized
usage costs or data access risks.

**Immediate next steps if approved**:
1. Revoke key: [Link to OpenAI dashboard]
2. Generate new key, update configuration
3. Clean git history: `git filter-branch --tree-filter 'sed -i "/sk-proj-/d" WORK-LOG.md' HEAD`
4. Check API usage logs for period 13:55-14:15
5. Update Common Pitfalls: "Never include API keys in outputs; use placeholders"

═══════════════════════════════════════════════════════════════
Your decision (1/2):
═══════════════════════════════════════════════════════════════
```

---

## Usage Guidelines

**When to use this template**:
- Any P0 (Critical) incident
- P1 (Major) incidents with significant impact
- User explicitly asks for "detailed incident report"

**When NOT to use** (use standard escalation format instead):
- P2 (Minor) incidents like single agent timeouts
- Routine decision escalations
- Progress updates or status reports

**Customization**:
- Severity levels: Adjust based on project context
- Response options: Provide 2-4 options; always include "abort/escalate"
- Artifacts: Link to specific files user should review
- Recommended action: Always provide one; explain reasoning clearly

**Tone**:
- Factual and direct, not defensive
- Focus on solutions, not blame
- Provide complete context so user can make informed decision
- Use clear formatting (headers, lists, tables) for scanability
