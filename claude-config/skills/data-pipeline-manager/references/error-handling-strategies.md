# Error Handling Strategies Reference

This document provides comprehensive error handling patterns and strategies for robust data pipelines, including retry logic, checkpointing, logging, and recovery mechanisms.

## Table of Contents

1. [Error Classification](#error-classification)
2. [Retry Patterns](#retry-patterns)
3. [Checkpointing and Resume](#checkpointing-and-resume)
4. [Logging Strategies](#logging-strategies)
5. [Circuit Breakers](#circuit-breakers)
6. [Graceful Degradation](#graceful-degradation)
7. [Error Recovery](#error-recovery)

## Error Classification

### Error Type Hierarchy

```python
class PipelineError(Exception):
    """Base class for all pipeline errors"""
    pass

class TransientError(PipelineError):
    """Temporary error that may succeed on retry"""
    pass

class PermanentError(PipelineError):
    """Non-recoverable error that should not be retried"""
    pass

class ConfigurationError(PermanentError):
    """Error in pipeline configuration"""
    pass

class DataError(PermanentError):
    """Error in input data quality or format"""
    pass

class ResourceError(TransientError):
    """Temporary resource unavailability"""
    pass

class NetworkError(TransientError):
    """Network-related temporary error"""
    pass
```

### Error Classification Function

```python
import errno
import socket
import requests

def classify_error(exception):
    """
    Classify exception as transient or permanent

    Returns: ('transient', retry_delay) or ('permanent', None)
    """

    # Network errors - retry with backoff
    if isinstance(exception, (socket.timeout, requests.Timeout)):
        return ('transient', 5)

    if isinstance(exception, (ConnectionError, requests.ConnectionError)):
        return ('transient', 10)

    # File system errors
    if isinstance(exception, OSError):
        # Disk full - don't retry immediately
        if exception.errno == errno.ENOSPC:
            return ('permanent', None)

        # Permission denied - permanent
        if exception.errno == errno.EACCES:
            return ('permanent', None)

        # File not found - permanent
        if exception.errno == errno.ENOENT:
            return ('permanent', None)

        # Resource temporarily unavailable - retry
        if exception.errno == errno.EAGAIN:
            return ('transient', 2)

    # Memory errors - don't retry
    if isinstance(exception, MemoryError):
        return ('permanent', None)

    # Value errors (bad data) - don't retry
    if isinstance(exception, (ValueError, TypeError)):
        return ('permanent', None)

    # Default: treat as permanent to avoid infinite retries
    return ('permanent', None)
```

## Retry Patterns

### Exponential Backoff

```python
import time
import random
from functools import wraps

def retry_with_exponential_backoff(
    max_retries=3,
    base_delay=1,
    max_delay=60,
    exponential_base=2,
    jitter=True
):
    """
    Decorator for retrying with exponential backoff

    Args:
        max_retries: Maximum number of retry attempts
        base_delay: Initial delay in seconds
        max_delay: Maximum delay between retries
        exponential_base: Base for exponential calculation (typically 2)
        jitter: Add random jitter to avoid thundering herd
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            last_exception = None

            for attempt in range(max_retries + 1):
                try:
                    return func(*args, **kwargs)

                except Exception as e:
                    last_exception = e

                    # Classify error
                    error_type, suggested_delay = classify_error(e)

                    # Don't retry permanent errors
                    if error_type == 'permanent':
                        raise

                    # Don't retry on last attempt
                    if attempt == max_retries:
                        raise

                    # Calculate delay
                    if suggested_delay is not None:
                        delay = suggested_delay
                    else:
                        delay = min(base_delay * (exponential_base ** attempt), max_delay)

                    # Add jitter
                    if jitter:
                        delay = delay * (0.5 + random.random())

                    print(f"Attempt {attempt + 1}/{max_retries} failed: {e}")
                    print(f"Retrying in {delay:.1f}s...")
                    time.sleep(delay)

            # Should never reach here, but just in case
            raise last_exception

        return wrapper
    return decorator

# Usage example
@retry_with_exponential_backoff(max_retries=3, base_delay=1)
def download_file(url, output_path):
    """Download file with automatic retry"""
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    with open(output_path, 'wb') as f:
        f.write(response.content)
```

### Conditional Retry

```python
def retry_with_conditions(
    max_retries=3,
    retry_on=None,
    dont_retry_on=None,
    retry_condition=None
):
    """
    Decorator for conditional retrying

    Args:
        max_retries: Maximum retry attempts
        retry_on: Exception types to retry (tuple)
        dont_retry_on: Exception types to never retry (tuple)
        retry_condition: Function to determine if should retry
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            for attempt in range(max_retries + 1):
                try:
                    return func(*args, **kwargs)

                except Exception as e:
                    # Check if we should not retry
                    if dont_retry_on and isinstance(e, dont_retry_on):
                        raise

                    # Check if we should retry
                    should_retry = False

                    if retry_on and isinstance(e, retry_on):
                        should_retry = True

                    if retry_condition and retry_condition(e):
                        should_retry = True

                    if not should_retry or attempt == max_retries:
                        raise

                    # Retry
                    delay = 2 ** attempt
                    print(f"Retrying after {delay}s (attempt {attempt + 1}/{max_retries})")
                    time.sleep(delay)

        return wrapper
    return decorator

# Usage example
@retry_with_conditions(
    max_retries=3,
    retry_on=(ConnectionError, TimeoutError),
    dont_retry_on=(ValueError, KeyError)
)
def process_data(data_source):
    """Process data with selective retry"""
    # Implementation
    pass
```

### Async Retry Pattern

```python
import asyncio
from functools import wraps

def async_retry(max_retries=3, base_delay=1):
    """Decorator for async functions with retry"""

    def decorator(func):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            for attempt in range(max_retries + 1):
                try:
                    return await func(*args, **kwargs)

                except Exception as e:
                    error_type, _ = classify_error(e)

                    if error_type == 'permanent' or attempt == max_retries:
                        raise

                    delay = base_delay * (2 ** attempt)
                    print(f"Async retry in {delay}s (attempt {attempt + 1}/{max_retries})")
                    await asyncio.sleep(delay)

        return wrapper
    return decorator

# Usage example
@async_retry(max_retries=3)
async def fetch_data(url):
    """Fetch data asynchronously with retry"""
    async with aiohttp.ClientSession() as session:
        async with session.get(url) as response:
            return await response.text()
```

## Checkpointing and Resume

### File-Based Checkpoint

```python
import json
import os
from typing import Any, Dict, List

class FileCheckpoint:
    """Simple file-based checkpoint manager"""

    def __init__(self, checkpoint_file):
        self.checkpoint_file = checkpoint_file
        self._state = self._load()

    def _load(self) -> Dict:
        """Load checkpoint from file"""
        if os.path.exists(self.checkpoint_file):
            with open(self.checkpoint_file) as f:
                return json.load(f)
        return {}

    def _save(self):
        """Save checkpoint to file"""
        # Write to temp file first (atomic write)
        temp_file = f"{self.checkpoint_file}.tmp"
        with open(temp_file, 'w') as f:
            json.dump(self._state, f, indent=2)

        # Rename to actual file (atomic on POSIX)
        os.rename(temp_file, self.checkpoint_file)

    def mark_completed(self, item_id: str, stage: str):
        """Mark item-stage as completed"""
        if item_id not in self._state:
            self._state[item_id] = {'completed_stages': []}

        if stage not in self._state[item_id]['completed_stages']:
            self._state[item_id]['completed_stages'].append(stage)
            self._save()

    def is_completed(self, item_id: str, stage: str) -> bool:
        """Check if item-stage is completed"""
        if item_id not in self._state:
            return False
        return stage in self._state[item_id].get('completed_stages', [])

    def get_completed_stages(self, item_id: str) -> List[str]:
        """Get list of completed stages for item"""
        if item_id not in self._state:
            return []
        return self._state[item_id].get('completed_stages', [])

    def clear(self):
        """Clear checkpoint (after successful completion)"""
        if os.path.exists(self.checkpoint_file):
            os.remove(self.checkpoint_file)
        self._state = {}

    def save_metadata(self, item_id: str, metadata: Dict[str, Any]):
        """Save arbitrary metadata for item"""
        if item_id not in self._state:
            self._state[item_id] = {}
        self._state[item_id]['metadata'] = metadata
        self._save()

    def get_metadata(self, item_id: str) -> Dict[str, Any]:
        """Get metadata for item"""
        if item_id not in self._state:
            return {}
        return self._state[item_id].get('metadata', {})

# Usage example
def process_pipeline_with_checkpoint(items, stages, checkpoint_file):
    """Process pipeline with checkpointing"""

    checkpoint = FileCheckpoint(checkpoint_file)

    for item in items:
        for stage in stages:
            # Skip if already completed
            if checkpoint.is_completed(item, stage):
                print(f"Skipping {item} - {stage} (already completed)")
                continue

            try:
                # Process stage
                print(f"Processing {item} - {stage}")
                result = process_stage(item, stage)

                # Save metadata
                checkpoint.save_metadata(item, {
                    'stage': stage,
                    'timestamp': time.time(),
                    'result_size': len(result) if hasattr(result, '__len__') else None
                })

                # Mark as completed
                checkpoint.mark_completed(item, stage)

            except Exception as e:
                print(f"Error processing {item} - {stage}: {e}")
                raise

    # Clear checkpoint on success
    checkpoint.clear()
```

### Database-Based Checkpoint

```python
import sqlite3
from datetime import datetime

class DatabaseCheckpoint:
    """SQLite-based checkpoint manager for larger pipelines"""

    def __init__(self, db_path):
        self.db_path = db_path
        self._init_db()

    def _init_db(self):
        """Initialize database schema"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS checkpoint (
                item_id TEXT,
                stage TEXT,
                status TEXT,
                timestamp TEXT,
                metadata TEXT,
                error_message TEXT,
                PRIMARY KEY (item_id, stage)
            )
        """)

        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_status
            ON checkpoint(status)
        """)

        conn.commit()
        conn.close()

    def mark_started(self, item_id: str, stage: str):
        """Mark stage as started"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("""
            INSERT OR REPLACE INTO checkpoint
            (item_id, stage, status, timestamp)
            VALUES (?, ?, 'in_progress', ?)
        """, (item_id, stage, datetime.now().isoformat()))

        conn.commit()
        conn.close()

    def mark_completed(self, item_id: str, stage: str, metadata: Dict = None):
        """Mark stage as completed"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        metadata_json = json.dumps(metadata) if metadata else None

        cursor.execute("""
            INSERT OR REPLACE INTO checkpoint
            (item_id, stage, status, timestamp, metadata)
            VALUES (?, ?, 'completed', ?, ?)
        """, (item_id, stage, datetime.now().isoformat(), metadata_json))

        conn.commit()
        conn.close()

    def mark_failed(self, item_id: str, stage: str, error: str):
        """Mark stage as failed"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("""
            INSERT OR REPLACE INTO checkpoint
            (item_id, stage, status, timestamp, error_message)
            VALUES (?, ?, 'failed', ?, ?)
        """, (item_id, stage, datetime.now().isoformat(), error))

        conn.commit()
        conn.close()

    def is_completed(self, item_id: str, stage: str) -> bool:
        """Check if stage is completed"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("""
            SELECT status FROM checkpoint
            WHERE item_id = ? AND stage = ?
        """, (item_id, stage))

        result = cursor.fetchone()
        conn.close()

        return result and result[0] == 'completed'

    def get_pending_items(self, stage: str) -> List[str]:
        """Get items that haven't completed this stage"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("""
            SELECT DISTINCT item_id FROM checkpoint
            WHERE stage = ? AND status != 'completed'
        """, (stage,))

        items = [row[0] for row in cursor.fetchall()]
        conn.close()

        return items

    def get_statistics(self) -> Dict:
        """Get pipeline statistics"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("""
            SELECT status, COUNT(*) FROM checkpoint
            GROUP BY status
        """)

        stats = dict(cursor.fetchall())
        conn.close()

        return stats
```

## Logging Strategies

### Structured Logging

```python
import logging
import json
from datetime import datetime

class StructuredLogger:
    """Structured logger for pipeline events"""

    def __init__(self, name, log_file=None):
        self.logger = logging.getLogger(name)
        self.logger.setLevel(logging.INFO)

        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'
        ))
        self.logger.addHandler(console_handler)

        # File handler for structured logs
        if log_file:
            file_handler = logging.FileHandler(log_file)
            self.logger.addHandler(file_handler)

    def log_event(self, event_type: str, **kwargs):
        """Log structured event"""
        event = {
            'timestamp': datetime.now().isoformat(),
            'event_type': event_type,
            **kwargs
        }
        self.logger.info(json.dumps(event))

    def log_stage_start(self, item_id: str, stage: str):
        """Log stage start"""
        self.log_event('stage_start', item_id=item_id, stage=stage)

    def log_stage_complete(self, item_id: str, stage: str, duration: float, **metrics):
        """Log stage completion"""
        self.log_event(
            'stage_complete',
            item_id=item_id,
            stage=stage,
            duration=duration,
            metrics=metrics
        )

    def log_stage_error(self, item_id: str, stage: str, error: str, traceback: str = None):
        """Log stage error"""
        self.log_event(
            'stage_error',
            item_id=item_id,
            stage=stage,
            error=error,
            traceback=traceback
        )

    def log_validation_failure(self, item_id: str, validation_type: str, details: Dict):
        """Log validation failure"""
        self.log_event(
            'validation_failure',
            item_id=item_id,
            validation_type=validation_type,
            details=details
        )

# Usage example
logger = StructuredLogger('pipeline', 'pipeline.log')

try:
    logger.log_stage_start('sample_001', 'alignment')
    start_time = time.time()

    # Process
    result = align_reads('sample_001')

    duration = time.time() - start_time
    logger.log_stage_complete(
        'sample_001',
        'alignment',
        duration,
        aligned_reads=result['aligned'],
        alignment_rate=result['rate']
    )

except Exception as e:
    logger.log_stage_error(
        'sample_001',
        'alignment',
        str(e),
        traceback.format_exc()
    )
    raise
```

### Context Manager for Logging

```python
from contextlib import contextmanager
import time

@contextmanager
def log_stage(logger, item_id, stage):
    """Context manager for automatic stage logging"""

    logger.log_stage_start(item_id, stage)
    start_time = time.time()
    error = None

    try:
        yield
    except Exception as e:
        error = e
        logger.log_stage_error(item_id, stage, str(e), traceback.format_exc())
        raise
    finally:
        duration = time.time() - start_time
        if error is None:
            logger.log_stage_complete(item_id, stage, duration)

# Usage
with log_stage(logger, 'sample_001', 'alignment'):
    result = align_reads('sample_001')
```

## Circuit Breakers

### Simple Circuit Breaker

```python
from enum import Enum
import threading

class CircuitState(Enum):
    CLOSED = "closed"  # Normal operation
    OPEN = "open"      # Failing, rejecting requests
    HALF_OPEN = "half_open"  # Testing if recovered

class CircuitBreaker:
    """
    Circuit breaker to prevent cascading failures

    Opens circuit after threshold failures, preventing further attempts.
    Automatically tries to close after timeout.
    """

    def __init__(self, failure_threshold=5, timeout=60, half_open_attempts=1):
        self.failure_threshold = failure_threshold
        self.timeout = timeout
        self.half_open_attempts = half_open_attempts

        self.state = CircuitState.CLOSED
        self.failure_count = 0
        self.last_failure_time = None
        self.lock = threading.Lock()

    def call(self, func, *args, **kwargs):
        """Execute function through circuit breaker"""

        with self.lock:
            # Check if we should try to close circuit
            if self.state == CircuitState.OPEN:
                if self._should_attempt_reset():
                    self.state = CircuitState.HALF_OPEN
                    self.failure_count = 0
                else:
                    raise Exception("Circuit breaker is OPEN")

        # Attempt call
        try:
            result = func(*args, **kwargs)

            with self.lock:
                if self.state == CircuitState.HALF_OPEN:
                    self.state = CircuitState.CLOSED
                self.failure_count = 0

            return result

        except Exception as e:
            with self.lock:
                self.failure_count += 1
                self.last_failure_time = time.time()

                if self.failure_count >= self.failure_threshold:
                    self.state = CircuitState.OPEN
                    print(f"Circuit breaker opened after {self.failure_count} failures")

            raise

    def _should_attempt_reset(self):
        """Check if enough time has passed to attempt reset"""
        if self.last_failure_time is None:
            return True
        return (time.time() - self.last_failure_time) >= self.timeout

    def get_state(self):
        """Get current circuit state"""
        return self.state

# Usage example
breaker = CircuitBreaker(failure_threshold=3, timeout=30)

for i in range(10):
    try:
        result = breaker.call(unreliable_function, arg1, arg2)
        print(f"Success: {result}")
    except Exception as e:
        print(f"Failed: {e}")
        time.sleep(1)
```

## Graceful Degradation

### Fallback Pattern

```python
def with_fallback(primary_func, fallback_func, fallback_on=(Exception,)):
    """
    Execute primary function, fall back to secondary on failure

    Args:
        primary_func: Primary function to try
        fallback_func: Fallback function if primary fails
        fallback_on: Exception types that trigger fallback
    """

    @wraps(primary_func)
    def wrapper(*args, **kwargs):
        try:
            return primary_func(*args, **kwargs)
        except fallback_on as e:
            print(f"Primary function failed ({e}), using fallback")
            return fallback_func(*args, **kwargs)

    return wrapper

# Usage example
def download_from_primary(url):
    """Download from primary server"""
    response = requests.get(f"https://primary.com/{url}")
    response.raise_for_status()
    return response.content

def download_from_mirror(url):
    """Download from mirror server"""
    response = requests.get(f"https://mirror.com/{url}")
    response.raise_for_status()
    return response.content

# Create function with fallback
download = with_fallback(
    download_from_primary,
    download_from_mirror,
    fallback_on=(requests.RequestException,)
)
```

### Partial Success Pattern

```python
def process_with_partial_success(items, process_func, continue_on_error=True):
    """
    Process items, collecting successes and failures separately

    Returns: (successes, failures)
    """

    successes = []
    failures = []

    for item in items:
        try:
            result = process_func(item)
            successes.append({'item': item, 'result': result})

        except Exception as e:
            failures.append({
                'item': item,
                'error': str(e),
                'traceback': traceback.format_exc()
            })

            if not continue_on_error:
                break

    return successes, failures

# Usage example
def process_sample(sample_id):
    # Processing logic
    pass

successes, failures = process_with_partial_success(
    sample_ids,
    process_sample,
    continue_on_error=True
)

print(f"Processed: {len(successes)} succeeded, {len(failures)} failed")

if failures:
    print("\nFailed samples:")
    for failure in failures:
        print(f"  {failure['item']}: {failure['error']}")
```

## Error Recovery

### Automatic Recovery

```python
class RecoverablePipeline:
    """Pipeline with automatic recovery strategies"""

    def __init__(self, checkpoint_file, logger):
        self.checkpoint = FileCheckpoint(checkpoint_file)
        self.logger = logger

    def run(self, items, stages):
        """Run pipeline with recovery"""

        for item in items:
            for stage in stages:
                # Skip if completed
                if self.checkpoint.is_completed(item, stage):
                    continue

                # Attempt with retry and recovery
                self._run_stage_with_recovery(item, stage)

    def _run_stage_with_recovery(self, item, stage, max_retries=3):
        """Run stage with automatic recovery"""

        for attempt in range(max_retries):
            try:
                # Run stage
                result = self._run_stage(item, stage)

                # Validate result
                self._validate_result(item, stage, result)

                # Checkpoint
                self.checkpoint.mark_completed(item, stage)
                return result

            except DataError as e:
                # Data error - try to repair
                self.logger.log_event('attempting_repair', item=item, stage=stage, error=str(e))

                if self._try_repair(item, stage, e):
                    continue  # Retry after repair
                else:
                    raise  # Can't repair, fail

            except ResourceError as e:
                # Resource error - wait and retry
                if attempt < max_retries - 1:
                    delay = 2 ** attempt
                    self.logger.log_event('resource_error_retry', item=item, delay=delay)
                    time.sleep(delay)
                    continue
                else:
                    raise

            except Exception as e:
                # Unknown error - log and fail
                self.logger.log_stage_error(item, stage, str(e))
                raise

    def _try_repair(self, item, stage, error):
        """Attempt to repair data error"""

        # Example repairs
        if "corrupt BAM" in str(error):
            # Try to rebuild index
            self._rebuild_bam_index(item)
            return True

        if "incomplete FASTQ" in str(error):
            # Try to re-download
            self._redownload_fastq(item)
            return True

        return False
```

## Summary

Key error handling principles:

1. **Classify errors appropriately** - Distinguish transient from permanent errors
2. **Retry intelligently** - Use exponential backoff with jitter for transient errors
3. **Checkpoint progress** - Enable resume after failures without repeating work
4. **Log comprehensively** - Capture context needed for debugging
5. **Fail gracefully** - Use circuit breakers and fallbacks to prevent cascading failures
6. **Recover automatically** - Implement repair strategies for common failure modes
7. **Report clearly** - Provide actionable error messages with context

These patterns can be combined and adapted for specific pipeline requirements.
