# Git Workflow Guide

Reference guide for version control integration in the programming-pm workflow (Phase 6).

---

## Branching Strategy (GitHub Flow)

This workflow uses GitHub Flow - simple and effective for most projects.

### Core Principles

1. **Main branch is always deployable**
2. **Feature branches for all work**
3. **Pull requests for all merges**
4. **Delete branches after merge**

### Branch Lifecycle

```
main ─────────────────────────────────────────────> (always deployable)
       │                              │
       └── feature/auth-module ──────┘ (merge via PR, delete branch)
```

---

## Branch Naming Conventions

### Format

`{type}/{task-id}-{brief-description}`

### Types

| Type | Use Case | Example |
|------|----------|---------|
| `feature/` | New functionality | `feature/TASK-001-user-auth` |
| `bugfix/` | Bug fixes | `bugfix/TASK-042-null-check` |
| `hotfix/` | Critical production fixes | `hotfix/security-patch-001` |
| `refactor/` | Code restructuring | `refactor/extract-validation` |
| `docs/` | Documentation only | `docs/api-reference-update` |
| `test/` | Test additions | `test/integration-coverage` |

### Rules

- Use lowercase
- Use hyphens (not underscores)
- Include task ID when applicable
- Keep description brief (<30 chars)

---

## Commit Message Format

Use Conventional Commits format.

### Structure

```
<type>(<scope>): <subject>

<body>

<footer>
```

### Types

| Type | Description |
|------|-------------|
| `feat` | New feature |
| `fix` | Bug fix |
| `docs` | Documentation only |
| `style` | Formatting, no code change |
| `refactor` | Code restructuring |
| `test` | Adding tests |
| `chore` | Maintenance tasks |

### Examples

```bash
# Simple commit
git commit -m "feat(auth): add JWT token validation"

# Multi-line commit (use HEREDOC)
git commit -m "$(cat <<'EOF'
feat(auth): add JWT token validation

- Implement token parsing with jose library
- Add expiration checking
- Include refresh token support

Closes #123
EOF
)"
```

### Rules

- Subject line: imperative mood, <50 chars, no period
- Body: explain what and why (not how)
- Footer: reference issues, breaking changes

---

## Pull Request Workflow

### Creating a PR

1. **Create feature branch from main**
   ```bash
   git checkout main
   git pull origin main
   git checkout -b feature/TASK-001-description
   ```

2. **Implement and test locally**
   ```bash
   # Make changes
   ruff check . && mypy src/ && pytest
   ```

3. **Commit with conventional format**
   ```bash
   git add src/module.py tests/test_module.py
   git commit -m "feat(module): implement feature X"
   ```

4. **Push and create PR**
   ```bash
   git push -u origin feature/TASK-001-description
   gh pr create --title "feat(module): implement feature X" --body "..."
   ```

### PR Description Template

```markdown
## Summary

[1-2 sentences describing what this PR does]

## Changes

- [Change 1]
- [Change 2]

## Test Plan

- [ ] Unit tests added/updated
- [ ] Integration tests pass
- [ ] Manual testing performed

## Risks

[Any risks identified in pre-mortem relevant to this change]

## Checklist

- [ ] Self-review completed
- [ ] Tests pass locally
- [ ] Documentation updated
```

### PR Review Process

1. **Automated checks run** (CI if configured)
2. **Request review** from senior-developer
3. **Address feedback** (max 3 revision cycles)
4. **Approve and merge** (squash commit)
5. **Delete branch**

---

## Merge Strategy

### Squash and Merge (Default)

- Combines all commits into one
- Clean main history
- Use for feature branches

```bash
gh pr merge --squash
```

### Merge Commit

- Preserves all commits
- Shows branch history
- Use for long-running branches

```bash
gh pr merge --merge
```

### Rebase and Merge

- Linear history
- No merge commits
- Use when commits are clean and atomic

```bash
gh pr merge --rebase
```

---

## Rollback Procedures

### Scenario 1: Quality Gate Failure Rollback

When code review rejects implementation:

```bash
# 1. Document rejection reasons in PR

# 2. Create new branch from pre-implementation state
git checkout main
git checkout -b feature/TASK-001-v2-description

# 3. Archive rejected code (optional)
git stash save "rejected-impl-TASK-001"

# 4. Notify programming-pm

# 5. Implement with clarified requirements
```

### Scenario 2: Tests Fail After Merge

When tests fail on main after merge:

```bash
# 1. Immediate revert
git checkout main
git revert HEAD
git push origin main

# 2. Open issue linking to failed test output
gh issue create --title "Regression: [description]" --body "..."

# 3. Assign root cause analysis

# 4. Re-implement on fresh branch
git checkout -b bugfix/TASK-001-regression
```

### Scenario 3: Integration Breaks Existing Functionality

When merge causes regression in unrelated code:

```bash
# 1. Identify breaking commit
git bisect start
git bisect bad HEAD
git bisect good <last-known-good>
# ... bisect process ...

# 2. Revert breaking commit(s)
git revert <breaking-commit>

# 3. Document regression in test suite
# Add test that catches this regression

# 4. Re-plan implementation
```

### Scenario 4: Emergency Rollback

For critical production issues:

```bash
# 1. Revert immediately (no review required)
git revert HEAD --no-edit
git push origin main

# 2. Create incident issue
gh issue create --title "[P0] Emergency rollback: [description]" --label "priority:p0"

# 3. Post-incident review
# Document in retrospective
```

---

## Git Safety Protocol

### ALWAYS

- [ ] Stage specific files by name (not `git add .`)
- [ ] Review staged changes before commit (`git diff --staged`)
- [ ] Use conventional commit format
- [ ] Create PRs for all changes to main
- [ ] Get code review before merge
- [ ] Delete branches after merge

### NEVER

- [ ] Force push to main (`git push --force origin main`)
- [ ] Commit directly to main (always use PR)
- [ ] Use `--amend` on pushed commits
- [ ] Use `--no-verify` to skip hooks
- [ ] Commit sensitive files (.env, credentials)

### Pre-Commit Checklist

Before every commit:

```bash
# 1. Check what's staged
git status
git diff --staged

# 2. Verify no sensitive files
git diff --staged --name-only | grep -E '\.(env|key|pem|secret)'

# 3. Run quality checks
ruff check . && mypy src/ && pytest

# 4. Commit with good message
git commit -m "type(scope): description"
```

---

## Quality Gate 6: PR Merge

### Prerequisites

- [ ] All previous gates passed (1-5)
- [ ] PR created with proper description
- [ ] Code review approved

### Criteria

- [ ] No merge conflicts
- [ ] CI pipeline green (if configured)
- [ ] PR description includes summary and test plan
- [ ] Branch is up to date with main

### Merge Process

```bash
# 1. Ensure branch is current
git checkout feature/TASK-001-description
git fetch origin main
git rebase origin/main

# 2. Resolve any conflicts
# ... resolve conflicts ...
git add .
git rebase --continue

# 3. Push updated branch
git push --force-with-lease origin feature/TASK-001-description

# 4. Merge via GitHub (squash)
gh pr merge --squash

# 5. Delete branch
git branch -d feature/TASK-001-description
git push origin --delete feature/TASK-001-description
```

### Override

- Repository admin can force merge
- Requires documentation: why override, follow-up actions
- Logged for audit

---

## Common Git Commands

### Daily Workflow

```bash
# Start new work
git checkout main && git pull
git checkout -b feature/TASK-XXX-description

# Save progress
git add src/file.py tests/test_file.py
git commit -m "feat(module): work in progress"

# Update branch with main
git fetch origin main
git rebase origin/main

# Push for review
git push -u origin feature/TASK-XXX-description
```

### Troubleshooting

```bash
# Undo last commit (keep changes)
git reset --soft HEAD~1

# Discard unstaged changes
git checkout -- path/to/file

# Stash changes temporarily
git stash save "description"
git stash pop

# View commit history
git log --oneline -10
```

### Cleanup

```bash
# Delete merged local branches
git branch --merged main | grep -v main | xargs git branch -d

# Prune remote tracking branches
git fetch --prune
```
