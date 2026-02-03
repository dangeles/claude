# Testing Strategies for Systematic Debugging

Effective test patterns for hypothesis validation during debugging.

---

## Binary Search Testing

**Use when**: You know code worked before, broke recently, need to find when

**Strategy**: Divide time/commits in half repeatedly until you find the breaking change

### Example: Finding Breaking Commit

```bash
# Known: Code worked at commit A (100 commits ago), broken at commit B (now)

# Binary search:
git checkout HEAD~50  # Middle point (50 commits back)
python test_case.py

# If it works:
# Bug introduced between commit 50 and now
# Next test: HEAD~25 (midpoint of 0-50)

# If it breaks:
# Bug introduced between commit 100 and 50
# Next test: HEAD~75 (midpoint of 50-100)

# Continue until you find exact commit
```

### Git Bisect (Automated Binary Search)

```bash
# Start bisect:
git bisect start
git bisect bad  # Current commit is bad
git bisect good <commit-hash>  # Last known good commit

# Git will checkout middle commit, you test:
python test_case.py

# If test passes:
git bisect good

# If test fails:
git bisect bad

# Git continues binary search automatically
# At end, shows first bad commit
git bisect reset  # Return to original state
```

**Time savings**: O(log n) instead of O(n)
- Testing 1000 commits: ~10 tests instead of 500 average

---

## Isolation Testing

**Use when**: Multiple components, need to identify which one is failing

**Strategy**: Replace components with known-good versions until system works

### Example: Three-Component System

```python
# System: Data Loader → Processor → Writer
# Bug: Output is corrupted

# Test 1: Replace all with known-good versions
data = load_data_GOOD()  # Use reference implementation
processed = process_GOOD(data)
write_GOOD(processed)
# Result: Works → bug is in one of YOUR components

# Test 2: Use your loader, rest good
data = load_data_YOUR()  # Your implementation
processed = process_GOOD(data)
write_GOOD(processed)
# Result: Works → loader is fine

# Test 3: Use your loader + processor, writer good
data = load_data_YOUR()
processed = process_YOUR(data)  # Your implementation
write_GOOD(processed)
# Result: FAILS → bug is in YOUR processor

# Found it: process_YOUR() is the problem
```

### Isolation Checklist

- [ ] **Component A only**: Replace B, C, D with reference versions
- [ ] **Components A + B**: Replace C, D with reference versions
- [ ] **Components A + B + C**: Replace D with reference version
- [ ] **All components**: Use your implementations

First configuration that fails identifies the buggy component.

---

## Differential Testing

**Use when**: "Works here, fails there" - comparing two environments

**Strategy**: Change one difference at a time until behavior changes

### Example: Works on Mac, Fails on Linux

| Factor | Mac (Works) | Linux (Fails) | Test Result |
|--------|-------------|---------------|-------------|
| OS | macOS | Linux | ? |
| Python version | 3.11.2 | 3.10.8 | ? |
| NumPy version | 1.24.0 | 1.23.5 | ? |
| Data path | /Users/... | /home/... | ? |
| File encoding | UTF-8 | UTF-8 | ? |

**Testing order** (change one at a time on Mac):

1. **Test Python 3.10.8 on Mac**: Still works → Python version not the issue
2. **Test NumPy 1.23.5 on Mac**: Still works → NumPy version not the issue
3. **Test with Linux data path on Mac**: FAILS → Data path is the issue!

**Root cause**: Data path contains spaces on Mac (`/Users/my user/data`), doesn't on Linux (`/home/user/data`). Code fails with spaces in path.

### Differential Testing Template

```markdown
## Differential Analysis

**Working environment** (Environment A):
- Python: 3.11.2
- OS: macOS
- Path: /Users/name/project

**Failing environment** (Environment B):
- Python: 3.10.8
- OS: Linux
- Path: /home/name/project

**Test Plan**: Change A to match B, one variable at a time

| Test | Changed Variable | A → B Value | Result |
|------|------------------|-------------|--------|
| 1 | Python version | 3.11.2 → 3.10.8 | ✓ Still works |
| 2 | OS | macOS → Linux | ? |
| 3 | Path | /Users → /home | ? |
| ... | ... | ... | ... |

First ❌ identifies the critical difference.
```

---

## Boundary Testing

**Use when**: Code fails on some inputs but not others

**Strategy**: Test edge cases and boundaries to define failure zone

### Example: Function Fails on Large Inputs

```python
def process(n):
    result = [i**2 for i in range(n)]
    return sum(result)

# Bug report: "Fails on large n"

# Boundary testing:
test_cases = [
    1,      # Minimum
    10,     # Small
    100,    # Medium
    1000,   # Large
    10000,  # Very large
    100000, # Extreme
]

for n in test_cases:
    try:
        result = process(n)
        print(f"n={n}: ✓ Success (result={result})")
    except Exception as e:
        print(f"n={n}: ✗ FAIL ({e})")

# Output:
# n=1: ✓ Success
# n=10: ✓ Success
# n=100: ✓ Success
# n=1000: ✓ Success
# n=10000: ✓ Success
# n=100000: ✗ FAIL (MemoryError)

# Conclusion: Fails at n >= 100000 due to memory
```

### Binary Search for Exact Boundary

```python
# We know: n=10000 works, n=100000 fails
# Find exact threshold:

low, high = 10000, 100000

while high - low > 1:
    mid = (low + high) // 2
    try:
        process(mid)
        low = mid  # Works at mid, try higher
        print(f"n={mid}: ✓ Works")
    except:
        high = mid  # Fails at mid, try lower
        print(f"n={mid}: ✗ Fails")

print(f"Threshold: works at {low}, fails at {high}")

# Output:
# n=55000: ✓ Works
# n=77500: ✗ Fails
# n=66250: ✓ Works
# n=71875: ✗ Fails
# ...
# Threshold: works at 70000, fails at 70001
```

**Use cases**:
- Find maximum file size before crash
- Find exact string length that triggers error
- Find array size that causes memory error

---

## Stress Testing

**Use when**: Intermittent bugs, race conditions, flaky tests

**Strategy**: Run test many times to establish failure rate

### Example: Flaky Test

```python
def test_clustering():
    """This test fails ~20% of the time"""
    result = run_clustering(data)
    assert len(result.clusters) == 8

# Stress test to measure failure rate:
failures = 0
runs = 100

for i in range(runs):
    try:
        test_clustering()
    except AssertionError:
        failures += 1

print(f"Failure rate: {failures}/{runs} = {failures/runs*100:.1f}%")

# Output: Failure rate: 22/100 = 22.0%
```

### Testing Potential Fix

```python
# Hypothesis: Adding random_state will fix flakiness

def test_clustering_FIXED():
    result = run_clustering(data, random_state=42)  # Added seed
    assert len(result.clusters) == 8

# Stress test again:
failures = 0
for i in range(100):
    try:
        test_clustering_FIXED()
    except AssertionError:
        failures += 1

print(f"Failure rate: {failures}/100")

# Output: Failure rate: 0/100 ✓ Fix confirmed!
```

**Guidelines**:
- Run at least 100 iterations for reliable statistics
- If failure rate > 5%, definitely a bug
- If failure rate drops from 20% → 0%, fix is confirmed

---

## Minimal Reproduction Testing

**Use when**: Complex codebase, hard to debug

**Strategy**: Strip away everything unrelated until bug is isolated

See `examples/minimal-reproduction-example.md` for detailed guide.

**Quick steps**:
1. Remove unrelated imports → still fails?
2. Replace complex data with toy data → still fails?
3. Inline function calls → still fails?
4. Remove error handling → still fails?

**Goal**: 10-20 lines that reproduce the exact same bug.

---

## Regression Testing

**Use when**: After fixing a bug, ensure it stays fixed

**Strategy**: Add test case that would have caught the bug

### Example

```python
# Bug found: Division by zero when processing empty DataFrame

# Regression test (add to test suite):
def test_empty_dataframe():
    """Regression test for issue #42: Division by zero on empty DataFrame"""
    df = pd.DataFrame()  # Empty DataFrame
    result = process_data(df)  # Should not crash
    assert result == 0  # Expected behavior for empty input

# This test will fail if the bug is reintroduced in future
```

**Best practice**: Link to issue number in test name/comment

---

## A/B Comparison Testing

**Use when**: "My changes broke something, what changed?"

**Strategy**: Run identical test on old vs new code, compare outputs

```python
# Save old implementation:
def process_OLD(data):
    # ... original code ...
    pass

# New implementation (possibly buggy):
def process_NEW(data):
    # ... modified code ...
    pass

# Compare on test data:
test_data = load_test_data()
result_old = process_OLD(test_data)
result_new = process_NEW(test_data)

# Check if results match:
if np.allclose(result_old, result_new):
    print("✓ Results match - no regression")
else:
    print("✗ Results differ - regression detected")
    diff = result_new - result_old
    print(f"Max difference: {np.max(np.abs(diff))}")
    print(f"Locations of difference:")
    print(np.where(~np.isclose(result_old, result_new)))
```

**Visual comparison**:
```python
import matplotlib.pyplot as plt

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
ax1.imshow(result_old); ax1.set_title('Old')
ax2.imshow(result_new); ax2.set_title('New')
ax3.imshow(result_new - result_old); ax3.set_title('Difference')
plt.show()
```

---

## Statistical Testing (for Numerical Bugs)

**Use when**: Results are "close but not exact"

### Example: Numerical Precision Issue

```python
# Bug: Results differ slightly between runs

results = []
for i in range(100):
    result = run_experiment()
    results.append(result)

# Statistical analysis:
import numpy as np
print(f"Mean: {np.mean(results)}")
print(f"Std: {np.std(results)}")
print(f"Range: {np.min(results)} - {np.max(results)}")

# If std is small (< 1e-6), likely floating point precision
# If std is large, likely a real bug (nondeterminism, race condition)
```

---

## Summary: Which Strategy When

| Symptom | Strategy | Why |
|---------|----------|-----|
| Worked before, broken now | Binary search | Find breaking change efficiently |
| Multi-component system fails | Isolation | Identify faulty component |
| Works here, fails there | Differential | Find critical environmental difference |
| Fails on some inputs | Boundary testing | Define failure zone |
| Intermittent failure | Stress testing | Measure failure rate, confirm fix |
| Complex codebase | Minimal reproduction | Simplify to essential bug trigger |
| After fixing bug | Regression testing | Prevent recurrence |
| Modified code, uncertain impact | A/B comparison | Detect unintended changes |
| Numerical instability | Statistical testing | Distinguish precision from bug |

---

## Pro Tips

1. **Automate tests**: Write scripts, don't rely on manual testing
2. **Save test cases**: Bug-triggering inputs become regression tests
3. **Document strategy**: Record which tests ruled out which hypotheses
4. **Use version control**: Easy to test old vs new code
5. **Start simple**: Try boundary testing before stress testing
6. **One variable**: Change only one thing between tests
