# Hypothesis Testing Example: Systematic Validation

## Scenario: API Returns 500 Error Intermittently

**Symptoms**:
- API endpoint `/api/v1/process` returns 500 error
- Happens ~20% of the time
- Same request parameters produce different results
- Error message: "Internal Server Error" (not helpful)

**Expected**: 200 OK with JSON response
**Actual**: 500 error ~20% of requests, 200 OK ~80%

---

## Hypothesis 1: Database Connection Timeout

**Reasoning**: Intermittent 500 errors often indicate resource exhaustion. Database connections might timeout under load.

**Prediction**: If this is correct, then errors should correlate with slow database queries.

**Test Design**:
```python
# Add logging to measure database query time
import time
import logging

@app.route('/api/v1/process')
def process():
    start = time.time()
    result = db.query("SELECT * FROM data WHERE id = ?", request.args.get('id'))
    duration = time.time() - start
    logging.info(f"Query duration: {duration}s")
    return jsonify(result)
```

**Test Execution**:
1. Send 100 requests to API
2. Correlate errors with query duration in logs

**Result**:
```
Successful requests: query duration 0.01-0.05s
Failed requests: query duration 0.02-0.04s
```

**Observation**: No correlation between query duration and errors.

**Conclusion**: ❌ REJECTED - Database timeout is not the cause.

---

## Hypothesis 2: Race Condition in Cache

**Reasoning**: Intermittent failures + no pattern in timing suggests race condition. API uses caching, might have concurrent access bug.

**Prediction**: If this is correct, then removing caching should eliminate errors.

**Test Design**:
```python
# Temporarily disable caching
@app.route('/api/v1/process')
def process():
    # cache = get_cache(request.args.get('id'))  # DISABLED
    # if cache:
    #     return jsonify(cache)

    result = compute_result(request.args.get('id'))
    # set_cache(request.args.get('id'), result)  # DISABLED
    return jsonify(result)
```

**Test Execution**:
1. Deploy with caching disabled
2. Send 100 requests
3. Observe error rate

**Result**:
```
With caching: 20% error rate
Without caching: 0% error rate
```

**Observation**: Errors disappear when caching is disabled!

**Conclusion**: ✅ CONFIRMED - Cache is involved in the bug. Need to investigate cache implementation.

---

## Hypothesis 2a: Cache Has Concurrent Write Bug

**Reasoning**: Since cache is involved, likely a thread-safety issue when multiple requests access cache simultaneously.

**Prediction**: If this is correct, then errors should correlate with concurrent requests.

**Test Design**:
```python
# Load test with concurrency
import concurrent.futures

def test_concurrent():
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        # Send 10 simultaneous requests
        futures = [executor.submit(call_api, i) for i in range(10)]
        results = [f.result() for f in futures]
        errors = [r for r in results if r.status_code == 500]
        return len(errors)

# Test with varying concurrency
for concurrency in [1, 2, 5, 10, 20]:
    errors = test_concurrent_load(concurrency)
    print(f"Concurrency {concurrency}: {errors} errors")
```

**Test Execution**: Run test with different concurrency levels.

**Result**:
```
Concurrency 1:  0 errors (0%)
Concurrency 2:  1 errors (10%)
Concurrency 5:  4 errors (40%)
Concurrency 10: 7 errors (70%)
Concurrency 20: 18 errors (90%)
```

**Observation**: Error rate increases with concurrency!

**Conclusion**: ✅ STRONG CONFIRMATION - This is a race condition in cache under concurrent access.

---

## Hypothesis 2b: Cache Uses Non-Thread-Safe Dictionary

**Reasoning**: Python's `dict` is not thread-safe for concurrent writes. If cache uses plain dict, concurrent updates will corrupt state.

**Prediction**: If this is correct, then code review will show cache using plain dict without locking.

**Test Design**: Review cache implementation code.

**Code Review**:
```python
# cache.py
_cache = {}  # ❌ Plain dict, not thread-safe!

def get_cache(key):
    return _cache.get(key)

def set_cache(key, value):
    _cache[key] = value  # ❌ Concurrent writes can corrupt dict
```

**Observation**: Cache uses plain dict with no locking mechanism.

**Conclusion**: ✅ ROOT CAUSE IDENTIFIED - Non-thread-safe cache implementation.

---

## Final Hypothesis: Thread-Safe Cache Will Resolve Issue

**Fix Proposed**:
```python
# cache.py
import threading

_cache = {}
_cache_lock = threading.Lock()  # Add lock

def get_cache(key):
    with _cache_lock:
        return _cache.get(key)

def set_cache(key, value):
    with _cache_lock:
        _cache[key] = value  # Thread-safe now
```

**Prediction**: With thread-safe cache, error rate should drop to 0%.

**Test Design**:
1. Apply fix
2. Run concurrency test with max_workers=20
3. Observe error rate

**Test Execution**: Run 100 requests with concurrency=20.

**Result**:
```
Before fix: 18 errors (90% at concurrency 20)
After fix: 0 errors (0%)
```

**Verification**:
- ✅ Original reproduction case (single request): Still works
- ✅ Concurrent requests (20 workers): No errors
- ✅ Load test (1000 requests): No errors
- ✅ Existing test suite: All pass

**Conclusion**: ✅ FIX CONFIRMED - Thread-safe cache resolves the issue.

---

## Summary of Hypothesis Testing Process

| Hypothesis | Method | Result | Time Spent |
|-----------|--------|--------|-----------|
| 1. Database timeout | Correlation analysis | ❌ Rejected | 15 min |
| 2. Cache involved | Disable caching | ✅ Confirmed | 20 min |
| 2a. Concurrent writes | Concurrency testing | ✅ Confirmed | 30 min |
| 2b. Non-thread-safe dict | Code review | ✅ Root cause | 10 min |
| Final: Thread-safe fix | Apply fix + test | ✅ Verified | 20 min |

**Total debugging time**: ~95 minutes from symptom to verified fix

**Key success factors**:
1. **Started broad, narrowed down**: Hypothesis 1 eliminated database, focused on cache
2. **Each test was specific**: Could clearly confirm or reject each hypothesis
3. **Built on previous results**: Hypothesis 2a built on confirmation of Hypothesis 2
4. **One variable at a time**: Each test changed only one thing
5. **Thorough verification**: Tested fix under multiple scenarios before declaring success

---

## Lessons Learned

1. **Intermittent bugs often indicate race conditions**: When bugs only happen sometimes, think concurrency.

2. **Correlation testing is powerful**: Comparing error rate vs concurrency level immediately pointed to race condition.

3. **Code review after behavioral confirmation**: Confirming "cache is involved" through testing, then reviewing cache code, was more efficient than reading all code upfront.

4. **Thread-safety is easy to miss**: Plain Python dict looks fine, but breaks under concurrency.

5. **Load testing reveals concurrency bugs**: Single-threaded tests would never have caught this.

---

## Hypothesis Testing Best Practices

**Good Hypothesis**:
- ✅ Specific: "Cache uses non-thread-safe dict"
- ✅ Testable: "Code review will show plain dict without locking"
- ✅ Falsifiable: Clear criteria for rejection

**Bad Hypothesis**:
- ❌ Vague: "Something's wrong with the cache"
- ❌ Untestable: "It's probably a timing issue"
- ❌ Not falsifiable: "The code might have a bug somewhere"

**Testing Principles**:
1. **Fastest test first**: Code review (10 min) vs setting up distributed load test (hours)
2. **One variable**: Don't change caching + database + server at once
3. **Reproducible**: Same test should give same result
4. **Quantitative when possible**: "90% error rate at concurrency 20" better than "errors happen often"
5. **Document negative results**: Knowing Hypothesis 1 was rejected saves time later
