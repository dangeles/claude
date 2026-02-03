# Completion Criteria by Domain

This reference guide provides domain-specific completion criteria to help tailor verification to different types of work. Use this to understand what "complete" means across various task domains.

## Code Implementation Tasks

### Core Criteria

**Functionality Complete**
- All specified features implemented
- Happy path works as intended
- User-facing behavior matches requirements
- Integration points function correctly

**Error Handling**
- Invalid inputs handled gracefully
- Error messages are informative and actionable
- Failure modes don't corrupt data or state
- Recovery paths work where applicable

**Testing**
- Unit tests cover new functionality
- Integration tests verify connected systems
- Edge cases have test coverage
- Existing tests still pass (no regressions)
- Test names clearly describe what they verify

**Code Quality**
- Code follows project style guidelines
- No obvious performance issues introduced
- Security considerations addressed (input validation, auth, etc.)
- Technical debt minimized or documented

**Documentation**
- Inline comments explain non-obvious logic
- API documentation matches implementation
- Configuration options documented
- Breaking changes noted in changelog

### Language-Specific Considerations

**Python**
- Type hints added where project uses them
- Virtual environment requirements updated
- Docstrings follow project convention (Google/NumPy/etc.)
- Linting passes (flake8, pylint, or project standard)

**JavaScript/TypeScript**
- Type definitions complete (if using TypeScript)
- Dependencies added to package.json with appropriate versions
- Build succeeds without warnings
- Browser compatibility considered

**Go**
- Exported functions have comments
- Error handling follows Go conventions (explicit returns)
- go fmt applied
- No race conditions (go test -race passes)

**Rust**
- Ownership and borrowing correct (compiles without unsafe unless justified)
- Error handling uses Result/Option appropriately
- Clippy warnings addressed
- Cargo.toml dependencies properly specified

### Framework-Specific Additions

**Web Frameworks (React, Vue, Django, Rails)**
- Responsive design works on mobile and desktop
- Accessibility basics covered (ARIA labels, keyboard nav)
- Loading states and empty states implemented
- Form validation provides good UX

**Database Changes**
- Migrations written and tested
- Rollback migrations exist and work
- Index performance considered
- Data migration handles existing records

**API Development**
- Endpoints follow project conventions
- Request/response schemas documented
- Rate limiting considered if applicable
- Versioning strategy followed

## Research Tasks

### Core Criteria

**Question Answered**
- Original research question directly addressed
- Answer is clear and unambiguous
- Appropriate depth for the question scope
- Assumptions and limitations stated

**Sources Credible**
- Information from reputable sources
- Sources are current (or historical context explained)
- Multiple sources corroborate key findings
- Conflicting information acknowledged and reconciled

**Evidence Quality**
- Claims supported by evidence, not opinion
- Data and statistics include context
- Anecdotal vs. systematic evidence distinguished
- Source bias considered where relevant

**Actionability**
- Findings presented in usable format
- Recommendations are specific and practical
- Next steps clearly identified
- Decision-makers have what they need

**Completeness**
- All aspects of research question covered
- Related questions addressed or explicitly scoped out
- No obvious gaps in coverage
- Synthesis provided, not just information dump

### Research Type Considerations

**Competitive Analysis**
- All relevant competitors identified
- Consistent evaluation criteria across competitors
- Strengths and weaknesses balanced
- Market positioning clear

**Technical Investigation**
- Proof of concept created if needed
- Performance characteristics quantified
- Compatibility issues identified
- Implementation feasibility assessed

**User Research**
- Sample size and selection described
- Research method appropriate to question
- Raw data and interpretation separated
- Participant quotes contextualized

**Literature Review**
- Search strategy documented
- Inclusion/exclusion criteria clear
- Key themes identified and synthesized
- Gaps in existing literature noted

## Analysis Tasks

### Core Criteria

**Data Quality**
- Data sources identified and vetted
- Data cleaning approach documented
- Missing data handled appropriately
- Outliers investigated, not just removed

**Methodology Sound**
- Analysis approach appropriate to question
- Assumptions stated explicitly
- Statistical methods used correctly
- Limitations acknowledged

**Interpretation Accurate**
- Conclusions follow from data
- Correlation vs. causation distinguished
- Alternative explanations considered
- Confidence levels appropriate

**Visualization Effective**
- Charts and graphs enhance understanding
- Visual design supports message (not decorative)
- Axes and labels clear
- Color choices accessible (colorblind-friendly)

**Recommendations Grounded**
- Suggestions tied to analysis findings
- Prioritization includes rationale
- Implementation considerations addressed
- Success metrics proposed

### Analysis Type Considerations

**Data Analysis**
- Exploratory analysis performed before conclusions
- Key metrics calculated correctly
- Trends distinguished from noise
- Data transformations justified

**Performance Analysis**
- Baseline established for comparison
- Bottlenecks identified with evidence
- Load testing conditions realistic
- Optimization recommendations quantified

**Risk Analysis**
- Risks categorized by likelihood and impact
- Mitigation strategies proposed
- Risk interdependencies identified
- Residual risk after mitigation noted

**Root Cause Analysis**
- Timeline of events established
- Multiple contributing factors considered
- 5 Whys or similar method applied
- Preventive measures address root cause, not symptoms

## Documentation Tasks

### Core Criteria

**Information Accurate**
- Technical details match current implementation
- Examples tested and work as shown
- Commands and code snippets execute successfully
- Version-specific information clearly marked

**Structure Clear**
- Logical organization aids navigation
- Table of contents for longer documents
- Headings follow consistent hierarchy
- Related sections cross-referenced

**Audience Appropriate**
- Terminology matches reader expertise level
- Assumed knowledge stated upfront
- Progressive disclosure (basics to advanced)
- Jargon defined or avoided

**Completeness**
- All necessary information included
- Prerequisites listed
- Common issues addressed (troubleshooting)
- Next steps or related resources provided

**Maintainability**
- Easy to update as system evolves
- Generated elements automated where possible
- Ownership clear (who updates this)
- Review schedule established for living docs

### Documentation Type Considerations

**API Documentation**
- All endpoints documented
- Request/response examples included
- Authentication requirements clear
- Error codes and meanings listed

**User Guides**
- Step-by-step instructions numbered
- Screenshots current and clear
- Success criteria for each step
- What to do if something goes wrong

**README Files**
- Quick start gets user running fast
- Installation steps complete
- Basic usage examples included
- Links to deeper documentation

**Architecture Documentation**
- System components and relationships shown
- Design decisions justified
- Technology choices explained
- Deployment architecture covered

**Runbooks/Playbooks**
- Emergency procedures prioritized
- Decision trees for diagnostics
- Commands ready to copy/paste
- Escalation paths defined

## Configuration/Infrastructure Tasks

### Core Criteria

**Infrastructure as Code**
- Configuration changes in version control
- Infrastructure reproducible from code
- Secrets not hardcoded
- Environment-specific values parameterized

**Testing**
- Configuration validated before applying
- Dry-run or plan reviewed
- Rollback plan exists and tested
- Changes applied to non-prod first

**Security**
- Least privilege principle applied
- Access controls configured correctly
- Audit logging enabled
- Compliance requirements met

**Documentation**
- Configuration decisions explained
- Runbook updated with new procedures
- Monitoring and alerting documented
- Disaster recovery plan current

**Reliability**
- High availability considered
- Failure modes analyzed
- Backup strategy implemented
- Recovery time objectives met

### Infrastructure Type Considerations

**Cloud Resources**
- Cost implications understood and approved
- Resource tagging applied consistently
- Autoscaling configured appropriately
- Region/availability zone strategy followed

**CI/CD Pipelines**
- Pipeline stages logical and efficient
- Test coverage adequate before deploy
- Deployment verification automated
- Pipeline failure notifications configured

**Monitoring/Alerting**
- Key metrics identified and tracked
- Alert thresholds set appropriately
- Alert routing to correct teams
- Runbook linked from alerts

**Database Configuration**
- Backups automated and tested
- Connection pooling configured
- Query performance monitored
- Replication/redundancy established

## Bug Fix Tasks

### Core Criteria

**Bug Reproduced**
- Steps to reproduce documented
- Reproduction confirmed before claiming fix
- Root cause identified, not just symptoms
- Edge case variants considered

**Fix Validated**
- Original bug scenario now works
- Fix doesn't break related functionality
- Regression test added to prevent recurrence
- Manual testing confirms user-visible resolution

**Scope Appropriate**
- Fix addresses reported issue completely
- Related bugs fixed or filed separately
- No over-engineering (fix what's broken)
- Side effects investigated

**Communication**
- Bug report updated with resolution
- Users who reported bug can be notified
- Release notes include fix description
- Known limitations documented if any

### Bug Severity Considerations

**Critical Bugs**
- Fix verified in production-like environment
- Rollback plan prepared
- Stakeholders notified of fix timeline
- Monitoring watches for successful resolution

**Performance Bugs**
- Performance metrics captured before/after
- Fix tested under realistic load
- Resource usage monitored
- Optimization doesn't sacrifice correctness

**UI/UX Bugs**
- Visual regression testing performed
- Multiple browsers/devices tested
- Accessibility not broken by fix
- User workflow end-to-end validated

## Refactoring Tasks

### Core Criteria

**Behavior Preserved**
- All existing tests still pass
- User-visible behavior unchanged
- API contracts maintained
- Performance characteristics similar or better

**Code Improvement**
- Refactoring achieves stated goal (readability, maintainability, etc.)
- Complexity reduced or better organized
- Duplication eliminated
- Abstractions appropriate (not over-engineered)

**Safe Changes**
- Changes made incrementally, not "big bang"
- Each commit leaves code in working state
- Refactoring separated from feature work
- Version control history preserves reasoning

**Testing Coverage**
- Test coverage maintained or improved
- Tests updated to reflect new structure
- Edge cases still tested
- Performance tests still passing

## Migration Tasks

### Core Criteria

**Migration Plan**
- Phased approach with milestones
- Rollback strategy at each phase
- Data integrity verification steps
- Downtime estimated and communicated

**Data Integrity**
- All data migrated completely
- Data transformations correct and reversible
- Validation queries confirm integrity
- Checksums or counts match source and destination

**Backward Compatibility**
- Old and new systems work during transition
- Feature flags control cutover if applicable
- API versions maintained through transition
- Users not disrupted during migration

**Verification**
- Smoke tests pass in new environment
- User acceptance testing completed
- Performance benchmarks comparable
- Monitoring confirms healthy operation

**Documentation**
- Migration runbook followed and updated
- New system documentation complete
- Training materials created for changes
- Decommissioning plan for old system

## Cross-Cutting Concerns

Regardless of domain, consider these aspects:

### Time and Resources
- Work completed within agreed timeframe
- Resource usage reasonable and sustainable
- Dependencies on others resolved
- Blockers escalated or cleared

### Communication
- Stakeholders informed of completion
- Changes communicated to affected teams
- Questions from reviews addressed
- Handoff documentation provided

### Quality Standards
- Professional standards met for the domain
- Peer review feedback incorporated
- Industry best practices followed
- Technical debt either avoided or consciously accepted

### Sustainability
- Solution maintainable by team
- Monitoring in place to detect issues
- Runbook or troubleshooting guide exists
- Knowledge shared, not siloed

## Applying These Criteria

When verifying completion:

1. **Start with core criteria** for the domain
2. **Add domain-specific considerations** based on technology/framework
3. **Check cross-cutting concerns** that apply everywhere
4. **Scale rigor to impact** (critical systems deserve deeper verification)
5. **Use professional judgment** (these are guidelines, not checklists to blindly follow)

Remember: "Complete" means ready for the next stage (review, deployment, user handoff), not perfect. Apply criteria that ensure quality without creating analysis paralysis.
