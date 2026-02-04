# Code Review Checklist

Reference checklist for Quality Gate 4 (Code Review Approval) in the programming-pm workflow.

---

## Automated Checks (Must All Pass)

These checks must pass before human review begins.

| Check | Command | Pass Criteria | Notes |
|-------|---------|---------------|-------|
| Linting | `ruff check .` | 0 errors | Warnings acceptable |
| Type Checking | `mypy --strict src/` | 0 errors | Warnings acceptable |
| Tests | `pytest -v` | All pass | No skipped tests without reason |
| Coverage | `pytest --cov=src --cov-fail-under=80` | >= 80% | For new code |

### Running All Automated Checks

```bash
# Full check before review
ruff check . && \
mypy --strict src/ && \
pytest --cov=src --cov-fail-under=80 -v

# Quick check during development
ruff check --fix . && mypy src/ && pytest -x
```

---

## Manual Review Checklist

### Code Quality

- [ ] **Logic Correctness**: Code does what it's supposed to do
  - Trace through key paths mentally
  - Verify boundary conditions handled
  - Check for off-by-one errors

- [ ] **Design Intent**: Code matches architecture specification
  - Follows component boundaries
  - Uses specified patterns
  - No unexpected dependencies introduced

- [ ] **No Code Duplication**: No repeated blocks >10 lines
  - Extract common functionality
  - Use appropriate abstractions

- [ ] **Error Handling**: Appropriate exceptions with useful messages
  - Specific exception types (not bare `except:`)
  - Error messages include context
  - Recovery where appropriate

- [ ] **Security Considerations**:
  - No hardcoded credentials
  - Input validation present
  - No SQL injection vulnerabilities
  - No path traversal vulnerabilities

### Documentation

- [ ] **Module Docstrings**: Each module has a docstring explaining purpose
  ```python
  """User authentication module.

  Provides functions for user login, logout, and session management.
  """
  ```

- [ ] **Function Docstrings**: All public functions have docstrings (Google style)
  ```python
  def authenticate(username: str, password: str) -> User:
      """Authenticate user credentials.

      Args:
          username: The user's login name.
          password: The user's password (plaintext).

      Returns:
          User object if authentication succeeds.

      Raises:
          AuthenticationError: If credentials are invalid.
      """
  ```

- [ ] **Type Hints**: Present on all public function signatures
  ```python
  def process_items(items: list[Item], config: Config) -> ProcessResult:
  ```

- [ ] **Inline Comments**: Complex logic explained
  - Why, not what (code shows what)
  - Algorithm explanations where needed

- [ ] **README Updated**: If public API changed

### Testing

- [ ] **Unit Tests Cover New Functionality**:
  - Happy path tested
  - Each public function has tests

- [ ] **Edge Cases Tested**:
  - Empty input
  - Boundary values
  - Maximum/minimum values
  - From pre-mortem risks

- [ ] **Error Conditions Tested**:
  - Invalid input handling
  - Exception raising verified

- [ ] **Test Names Descriptive**:
  ```python
  # Good
  def test_authenticate_returns_user_for_valid_credentials(self):

  # Bad
  def test_auth(self):
  ```

- [ ] **Tests Independent**: No shared mutable state between tests

### Architecture Compliance

- [ ] **Follows Component Boundaries**: From architecture specification
- [ ] **No Unexpected Dependencies**: Check imports
- [ ] **Data Flow Matches Design**: Inputs and outputs as specified
- [ ] **Interface Contracts Honored**: Function signatures match specification

---

## Review Outcome Decision

### Approve (LGTM)

Code meets all criteria. Minor suggestions may be included but are not blocking.

**Criteria for approval**:
- All automated checks pass
- All required checklist items satisfied
- No security concerns
- Code is maintainable and readable

### Request Changes

Code has issues that must be addressed before merge.

**Always request changes for**:
- Failing automated checks
- Security vulnerabilities
- Missing tests for new functionality
- Logic errors
- Missing required documentation

**Feedback format**:
```markdown
## Code Review: [PR/Task ID]

### Status: CHANGES_REQUESTED

### Required Changes
1. [File:Line] [Issue description]
   - Current: [what's wrong]
   - Suggested: [how to fix]

2. [File:Line] [Issue description]
   - Current: [what's wrong]
   - Suggested: [how to fix]

### Suggestions (not blocking)
- [Optional improvement]

### Questions
- [Clarification needed]
```

---

## Revision Cycle Protocol

### Maximum Cycles

- **junior-developer code**: 3 revision cycles maximum
- **senior-developer code**: 2 revision cycles maximum

After maximum cycles exceeded:
1. Escalate to programming-pm
2. Document blocking issues
3. Consider task redefinition

### Revision Tracking

| Revision | Changes Made | Reviewer | Status |
|----------|--------------|----------|--------|
| 1 | Initial submission | | CHANGES_REQUESTED |
| 2 | Addressed feedback | | CHANGES_REQUESTED |
| 3 | Final changes | | APPROVED / ESCALATE |

---

## Code Review Best Practices

### For Reviewers

1. **Be Specific**: Point to exact lines, suggest fixes
2. **Be Constructive**: Explain why, not just what
3. **Prioritize**: Distinguish blocking vs. non-blocking
4. **Be Timely**: Review within 30 minutes of submission
5. **Acknowledge Good Work**: Positive feedback matters

### For Authors

1. **Self-Review First**: Use this checklist before submitting
2. **Small PRs**: Easier to review, faster feedback
3. **Context in PR Description**: What, why, how to test
4. **Respond to All Comments**: Even if just acknowledging
5. **Don't Take It Personally**: Review is about code, not you

---

## Override Protocol

When deadline pressure requires merging code that doesn't fully pass review:

### Override Allowed For
- Minor documentation gaps
- Non-critical suggestions
- Style preferences

### Override NOT Allowed For
- Failing tests
- Security vulnerabilities
- Missing error handling for critical paths
- Architecture violations

### Override Process
1. programming-pm approval required
2. Document override in PR description:
   ```
   OVERRIDE: Code Review Gate
   Reason: [Why override needed]
   Items Deferred: [What wasn't addressed]
   Follow-up: [Issue number for tech debt]
   ```
3. Create follow-up issue for deferred items
4. Tag merged code as "TECH_DEBT"

---

## Quality Gate 4 Summary

**Type**: Human judgment (senior-developer review)

**Prerequisites**:
- All automated checks pass
- Self-review completed by author

**Pass Criteria**:
- All required checklist items satisfied
- No blocking issues remaining

**Fail Action**:
- Return to developer with specific feedback
- Track revision count

**Override**:
- programming-pm can approve with TECH_DEBT tag
- Only for deadline-critical situations
- Must create follow-up issue
