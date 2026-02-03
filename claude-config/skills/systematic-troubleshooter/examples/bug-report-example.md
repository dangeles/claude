# Bug Report Example: DataFrame Memory Leak in Batch Processing

**Date**: 2026-01-29
**Severity**: Major
**Status**: RESOLVED

## Symptoms

Processing 100 CSV files in a loop caused memory usage to grow from 200MB to 8GB, eventually crashing with `MemoryError`. Each file is only 50MB, so expected peak memory should be ~200MB (4 files buffered).

**Error message**:
```
MemoryError: Unable to allocate 1.2 GB for an array with shape (150000000,) and data type float64
```

**Observed behavior**:
- File 1-20: Memory stays ~500MB
- File 21-50: Memory grows to 2GB
- File 51-80: Memory grows to 5GB
- File 81+: Crashes with MemoryError

## Root Cause

**Issue**: Pandas DataFrames were accumulating in memory despite going out of scope because results list maintained references.

**Detailed explanation**:
```python
results = []
for file in files:
    df = pd.read_csv(file)  # Creates DataFrame
    processed = process_data(df)  # Returns processed DataFrame
    results.append(processed.mean())  # BUG: Keeps reference to 'processed'
    # 'df' and 'processed' go out of scope, but garbage collector
    # can't free them because results[i] holds reference to processed.index
    # which holds reference to processed DataFrame
```

The `.mean()` operation returns a Series that internally references the original DataFrame's index. This prevents garbage collection of the entire DataFrame.

## Investigation Process

**Hypothesis 1: Garbage collector not running**
- **Test**: Added `gc.collect()` after each iteration
- **Result**: REJECTED - Memory still grows, GC runs but doesn't free memory
- **Conclusion**: References exist preventing collection

**Hypothesis 2: Memory leak in pandas C extension**
- **Test**: Used objgraph to find reference cycles
- **Result**: REJECTED - No circular references found
- **Conclusion**: Not a reference cycle issue

**Hypothesis 3: Hidden references in results list**
- **Test**: Used `sys.getsizeof()` on results list and inspected with `gc.get_referents()`
- **Result**: CONFIRMED - Each Series in results maintains reference chain to original DataFrame
- **Observation**: results[0] → Series → Index → DataFrame (50MB each)
- **Conclusion**: Results list unintentionally keeping all DataFrames alive

## Solution

Convert Series to basic Python types to break reference chain:

```diff
results = []
for file in files:
    df = pd.read_csv(file)
    processed = process_data(df)
-   results.append(processed.mean())
+   results.append(processed.mean().to_dict())  # Convert to dict, breaks DataFrame reference
```

**Why this works**:
- `.to_dict()` creates plain Python dict with no pandas object references
- DataFrames can now be garbage collected after each iteration
- Memory usage stays constant at ~200MB

## Verification

**Test 1: Original reproduction case**
- Steps: Process 100 files with fix applied
- Result: ✅ PASS - Memory stays at 180-220MB throughout, no crash

**Test 2: Large file set (1000 files)**
- Steps: Process 1000 files to stress test
- Result: ✅ PASS - Memory stable at 200MB, completes successfully

**Test 3: Performance check**
- Before fix: Crashes at file 83
- After fix: Processes all 100 files in 45 seconds
- Result: ✅ IMPROVEMENT - Can actually complete the job

**Test 4: Regression check**
- Steps: Run existing test suite
- Result: ✅ PASS - All 47 tests pass, no regressions

## Prevention

**How to avoid this in the future**:

1. **Use generators instead of lists when possible**:
   ```python
   # Instead of accumulating in list
   results = (process_file(f).mean().to_dict() for f in files)
   ```

2. **Be aware of pandas reference chains**:
   - Series.mean() → Series → Index → DataFrame
   - DataFrame.head() → DataFrame slice → Original DataFrame
   - Always convert to primitive types if long-lived storage needed

3. **Monitor memory during development**:
   ```python
   import tracemalloc
   tracemalloc.start()
   # ... your code ...
   snapshot = tracemalloc.take_snapshot()
   top_stats = snapshot.statistics('lineno')
   for stat in top_stats[:10]:
       print(stat)
   ```

4. **Use context managers for large data**:
   ```python
   with pd.read_csv(file, chunksize=10000) as reader:
       for chunk in reader:
           process_chunk(chunk)
   # DataFrame chunks automatically freed
   ```

## Related Issues

- Stack Overflow: [Pandas DataFrame not releasing memory](https://stackoverflow.com/questions/...)
- GitHub Issue: pandas#12345 - "Hidden references in Index objects"
- Similar bug in our codebase: 2025-11-15 RNA-seq memory leak (fixed by using `.values` instead of Series)

---

## Lessons Learned

1. **Memory leaks aren't always obvious**: The bug was one line that looked harmless
2. **Python's garbage collection is reference-based**: If any reference exists, object won't be freed
3. **Pandas objects can have hidden reference chains**: Always consider what an object references internally
4. **Minimal reproduction is key**: Simplified test case (3 files instead of 100) made root cause obvious
5. **Tools matter**: `gc.get_referents()` and `sys.getsizeof()` were essential for finding the reference chain
