# Parallel Coordination Example: Multi-Domain Research Task

## Scenario

A user makes the following request:

> "I need help with three things: 1) Research the latest best practices for API rate limiting in 2026, 2) Analyze the performance bottlenecks in our backend service code, and 3) Find examples of effective error handling patterns in Python microservices. These are for different projects I'm working on."

## Initial Analysis

### Task Identification

The request contains three distinct tasks:

1. **Research Task**: Find current best practices for API rate limiting
2. **Code Analysis Task**: Identify performance bottlenecks in existing code
3. **Pattern Search Task**: Locate error handling examples

### Dependency Check

**Task 1 Dependencies**:
- Inputs: None (web research)
- Outputs: Rate limiting recommendations
- Resources: Web search API

**Task 2 Dependencies**:
- Inputs: Backend service code files
- Outputs: Performance bottleneck report
- Resources: Local codebase files (read-only)

**Task 3 Dependencies**:
- Inputs: None (example search)
- Outputs: Error handling pattern examples
- Resources: Web search API, code repositories

**Independence Verification**:
- ✓ No task requires output from another
- ✓ Task 1 and 3 use same resource (web search) but read-only - SAFE
- ✓ Task 2 uses different resource (local files) - SAFE
- ✓ All results can be aggregated independently

**Decision**: Tasks are fully independent. Proceed with parallelization.

## Parallel Execution Plan

### Task 1: Research Agent - API Rate Limiting

**Objective**: Find and summarize 2026 best practices for API rate limiting.

**Execution Steps**:
1. Perform web search for "API rate limiting best practices 2026"
2. Review top results for current recommendations
3. Identify common patterns (token bucket, leaky bucket, sliding window)
4. Note specific implementation considerations
5. Summarize findings with source citations

**Expected Resources**:
- WebSearch tool
- WebFetch tool for detailed documentation

**Estimated Duration**: 45-60 seconds

### Task 2: Analysis Agent - Performance Bottlenecks

**Objective**: Analyze backend service code to identify performance issues.

**Execution Steps**:
1. Use Glob to find relevant backend service files
2. Read service implementation files
3. Look for common bottlenecks:
   - N+1 query patterns
   - Blocking I/O operations
   - Inefficient loops
   - Missing indexes or caching
4. Prioritize findings by impact
5. Generate report with line numbers and recommendations

**Expected Resources**:
- Glob tool for file discovery
- Read tool for code inspection
- Grep tool for pattern matching

**Estimated Duration**: 30-45 seconds

### Task 3: Pattern Search Agent - Error Handling

**Objective**: Find well-documented examples of error handling in Python microservices.

**Execution Steps**:
1. Search for "Python microservices error handling patterns"
2. Look for code examples on GitHub, documentation sites
3. Identify common approaches:
   - Exception hierarchies
   - Error context propagation
   - Retry patterns
   - Circuit breakers
4. Extract concrete code examples
5. Categorize by use case

**Expected Resources**:
- WebSearch tool
- WebFetch tool for code examples

**Estimated Duration**: 45-60 seconds

## Execution: Launching Parallel Tasks

The coordinator creates all three tasks simultaneously:

```
User Request Received → Dependency Analysis Complete → Launch Parallel Execution

┌─────────────────────────────────────────────────────────────┐
│ Parallel Execution Block (Single Function Call Message)     │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  Task 1: Research API Rate Limiting                         │
│  │                                                           │
│  ├─ WebSearch: "API rate limiting 2026"                     │
│  ├─ WebFetch: Documentation pages                           │
│  └─ Synthesize: Best practices summary                      │
│                                                              │
│  Task 2: Analyze Performance Bottlenecks                    │
│  │                                                           │
│  ├─ Glob: Find backend service files                        │
│  ├─ Read: Load service implementation                       │
│  ├─ Grep: Search for antipatterns                           │
│  └─ Generate: Bottleneck report                             │
│                                                              │
│  Task 3: Find Error Handling Patterns                       │
│  │                                                           │
│  ├─ WebSearch: "Python error handling microservices"        │
│  ├─ WebFetch: Code example repositories                     │
│  └─ Extract: Pattern examples                               │
│                                                              │
└─────────────────────────────────────────────────────────────┘
                          │
                          ▼
            All tasks execute concurrently
                          │
                          ▼
                  Results collected
                          │
                          ▼
                 Integration phase
```

## Results Collection

After all three tasks complete, the coordinator collects outputs:

### Task 1 Results: API Rate Limiting Best Practices

**Summary**: Modern API rate limiting in 2026 emphasizes:
- Token bucket algorithm for burst handling
- Distributed rate limiting using Redis
- Per-user and per-endpoint limits
- Graceful degradation with 429 status codes
- Rate limit headers (X-RateLimit-Remaining, X-RateLimit-Reset)

**Key Sources**:
- IETF draft-ietf-httpapi-ratelimit-headers
- Major cloud provider documentation
- Open source implementations (Kong, Nginx)

### Task 2 Results: Performance Bottleneck Analysis

**Critical Issues Found**:

1. **N+1 Query Pattern** (High Impact)
   - Location: `/services/user_service.py`, lines 145-160
   - Issue: Loading user preferences in loop
   - Recommendation: Use batch loading or eager loading

2. **Blocking Database Calls** (High Impact)
   - Location: `/services/order_service.py`, lines 78-92
   - Issue: Synchronous database queries in request handler
   - Recommendation: Use async/await pattern

3. **Missing Cache Layer** (Medium Impact)
   - Location: `/services/product_service.py`, lines 200-230
   - Issue: Repeated queries for static product data
   - Recommendation: Implement Redis caching

4. **Inefficient JSON Parsing** (Low Impact)
   - Location: `/utils/parsers.py`, lines 45-67
   - Issue: Re-parsing same JSON multiple times
   - Recommendation: Parse once, reuse object

**Estimated Performance Gain**: 60-75% reduction in response time if all addressed.

### Task 3 Results: Error Handling Patterns

**Pattern 1: Custom Exception Hierarchy**

```python
class ServiceException(Exception):
    """Base exception for all service errors"""
    def __init__(self, message, error_code=None, context=None):
        self.message = message
        self.error_code = error_code
        self.context = context or {}
        super().__init__(self.message)

class ValidationError(ServiceException):
    """Raised when input validation fails"""
    pass

class ExternalServiceError(ServiceException):
    """Raised when external dependency fails"""
    pass
```

**Pattern 2: Error Context Propagation**

```python
from contextvars import ContextVar

request_id_var = ContextVar('request_id', default=None)

def handle_request(request):
    request_id_var.set(request.headers.get('X-Request-ID'))
    try:
        return process(request)
    except Exception as e:
        logger.error(
            f"Request {request_id_var.get()} failed",
            exc_info=e,
            extra={'request_id': request_id_var.get()}
        )
        raise
```

**Pattern 3: Retry with Exponential Backoff**

```python
import tenacity

@tenacity.retry(
    stop=tenacity.stop_after_attempt(3),
    wait=tenacity.wait_exponential(multiplier=1, min=2, max=10),
    retry=tenacity.retry_if_exception_type(ExternalServiceError)
)
def call_external_service(data):
    response = requests.post(SERVICE_URL, json=data)
    if response.status_code >= 500:
        raise ExternalServiceError("Service unavailable")
    return response.json()
```

**Pattern 4: Circuit Breaker**

```python
from pybreaker import CircuitBreaker

breaker = CircuitBreaker(
    fail_max=5,
    timeout_duration=60,
    exclude=[ValidationError]
)

@breaker
def call_payment_gateway(payment_data):
    return payment_api.process(payment_data)
```

## Result Integration and Synthesis

### Cross-Task Insights

1. **Rate Limiting + Performance**: The performance bottlenecks identified (especially blocking calls) will affect rate limit behavior. When requests are slow, rate limits are hit more easily. Addressing performance issues should be done before fine-tuning rate limits.

2. **Error Handling + Performance**: The error handling patterns (especially retry logic and circuit breakers) are directly relevant to handling the external service calls identified in the performance analysis.

3. **Unified Recommendation**: Implement circuit breakers (from Task 3) for the external service calls identified in performance analysis (Task 2), and ensure rate limiting (Task 1) accounts for retry attempts.

### Integrated Response Structure

The coordinator synthesizes results into a cohesive response:

---

**Response to User**:

I've completed all three research tasks in parallel. Here's what I found:

## 1. API Rate Limiting Best Practices (2026)

[Task 1 results with sources...]

## 2. Performance Bottleneck Analysis

[Task 2 results with specific file locations and recommendations...]

## 3. Python Microservices Error Handling Patterns

[Task 3 results with code examples...]

## Integration Recommendations

Based on insights across all three areas:

1. **Address performance bottlenecks first** - The blocking I/O patterns identified will impact rate limiting effectiveness
2. **Implement circuit breakers** - Use the patterns from the error handling research for the external service calls in your backend
3. **Design rate limits accounting for retries** - Ensure retry logic doesn't cause users to hit rate limits prematurely

---

## Performance Analysis

### Sequential Execution Time
- Task 1: 60 seconds
- Task 2: 45 seconds
- Task 3: 60 seconds
- **Total: 165 seconds (2 minutes 45 seconds)**

### Parallel Execution Time
- All tasks: ~60 seconds (limited by slowest task)
- **Total: 60 seconds (1 minute)**

### Performance Gain
- **Time Saved: 105 seconds (63% reduction)**
- **Speedup Factor: 2.75x**

## Key Success Factors

1. **True Independence**: No task required output from another
2. **I/O-Bound Operations**: All tasks primarily involved I/O (web requests, file reads)
3. **Appropriate Granularity**: Each task substantial enough to justify parallel overhead
4. **Safe Resource Sharing**: Web search API safely shared by multiple tasks
5. **Effective Synthesis**: Results integrated to provide cross-task insights

## Lessons Learned

1. **Early Dependency Analysis Pays Off**: Thorough upfront analysis ensured no race conditions or blocking
2. **Cross-Task Insights Add Value**: Integration phase revealed connections user might not have considered
3. **Clear Task Boundaries**: Well-defined tasks made parallel execution straightforward
4. **Progress Visibility**: Task status tracking helped monitor overall progress

This example demonstrates how the Parallel Coordinator skill can transform a multi-faceted request into an efficient parallel workflow, delivering comprehensive results in a fraction of sequential execution time.
