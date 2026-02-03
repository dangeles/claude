# Common Error Patterns and Their Typical Causes

Quick reference for recognizing frequent bug patterns.

---

## Import Errors

### `ModuleNotFoundError: No module named 'X'`

**Typical causes**:
1. Package not installed: `pip install X`
2. Wrong environment: Check `which python`, activate correct virtualenv
3. Missing `__init__.py`: Package directory needs `__init__.py` file
4. Typo in import: `from sklearn import ...` not `from scikit-learn import ...`
5. Circular import: Module A imports Module B which imports Module A

**Quick diagnosis**:
```bash
# Check if package is installed
pip list | grep package-name

# Check Python path
python -c "import sys; print('\n'.join(sys.path))"

# Check for __init__.py
ls -la package_directory/
```

---

## Index Errors

### `IndexError: list index out of range`

**Typical causes**:
1. Empty list: `items[0]` when `items = []`
2. Off-by-one: `for i in range(len(items)+1): items[i]` (tries items[len])
3. Wrong assumption about size: Assumed 10 items, actually 5

**Quick diagnosis**:
```python
print(f"List length: {len(items)}")
print(f"Trying to access index: {index}")
print(f"List contents: {items}")
```

**Common in bioinformatics**:
- 0-based vs 1-based indexing confusion
- Assuming all samples have same number of features

---

## Type Errors

### `TypeError: unsupported operand type(s) for +: 'int' and 'str'`

**Typical causes**:
1. Forgot to convert: `result = "Value: " + 5` should be `"Value: " + str(5)`
2. Wrong type from input: `input()` returns str, need `int(input())`
3. Mixed types in list: `[1, 2, '3', 4].sum()` fails on string

**Quick diagnosis**:
```python
print(f"Type of x: {type(x)}")
print(f"Value of x: {repr(x)}")  # repr shows quotes for strings
```

---

## Key Errors

### `KeyError: 'gene_name'`

**Typical causes**:
1. Typo: Looking for `gene_name`, key is actually `geneName` (case sensitive)
2. Missing key: Not all dictionaries have this key
3. DataFrame vs dict: `df['column']` works, `row['column']` might not if row is Series

**Safe access**:
```python
# Instead of: value = dict['key']
value = dict.get('key', default_value)  # Returns default if key missing

# Or check first:
if 'key' in dict:
    value = dict['key']
```

---

## Attribute Errors

### `AttributeError: 'NoneType' object has no attribute 'X'`

**Typical causes**:
1. Function returned None: `result = func(); result.method()` but `func()` returned None
2. Failed initialization: `obj = Class(); obj.attribute` but `__init__` failed
3. DataFrame operation returned None: `df.dropna(inplace=True); df.head()` (dropna returns None when inplace=True)

**Quick diagnosis**:
```python
print(f"Object is: {obj}")
print(f"Type is: {type(obj)}")
if obj is None:
    print("ERROR: Object is None!")
```

**Common in pandas**:
```python
# Wrong (returns None):
df = df.dropna(inplace=True)  # inplace returns None!

# Right:
df.dropna(inplace=True)  # Modifies df in place
# OR
df = df.dropna()  # Returns new DataFrame
```

---

## Memory Errors

### `MemoryError: Unable to allocate X GB`

**Typical causes**:
1. Loading entire file into memory: Use chunking
2. Creating huge array by mistake: `np.zeros((1000000, 1000000))` is 8TB
3. Memory leak: Objects not being freed (references kept alive)
4. Inefficient algorithm: O(n²) space when O(n) possible

**Quick diagnosis**:
```python
import sys
print(f"Size of object: {sys.getsizeof(obj) / 1e9:.2f} GB")

# For arrays:
import numpy as np
print(f"Array shape: {arr.shape}")
print(f"Array size: {arr.nbytes / 1e9:.2f} GB")
```

**Solutions**:
```python
# Instead of: df = pd.read_csv('huge.csv')
# Use chunking:
for chunk in pd.read_csv('huge.csv', chunksize=10000):
    process(chunk)

# Instead of: big_list = [process(x) for x in huge_data]
# Use generator:
big_generator = (process(x) for x in huge_data)
```

---

## Value Errors

### `ValueError: could not convert string to float: 'NA'`

**Typical causes**:
1. Missing data coded as string: `float('NA')` fails
2. Unexpected format: `float('1,234')` fails (comma), `float('$100')` fails (dollar sign)
3. Whitespace: `float(' 123 \n')` might fail

**Quick diagnosis**:
```python
print(f"Trying to convert: {repr(value)}")
print(f"Type: {type(value)}")
print(f"Stripped: {repr(value.strip())}")
```

**Solutions**:
```python
# Handle NA values:
value = float(value) if value != 'NA' else np.nan

# Pandas handles this automatically:
df = pd.read_csv('file.csv', na_values=['NA', 'N/A', '', 'null'])
```

---

## FileNotFoundError

### `FileNotFoundError: [Errno 2] No such file or directory: 'data.csv'`

**Typical causes**:
1. Wrong working directory: File is in `/home/user/data.csv`, you're in `/home/user/code/`
2. Relative path: `data.csv` vs `../data.csv` vs `/full/path/data.csv`
3. Typo in filename: `data.csv` vs `Data.csv` (case sensitive on Linux/Mac)
4. File doesn't exist: Need to download or generate it first

**Quick diagnosis**:
```python
import os
print(f"Current directory: {os.getcwd()}")
print(f"File exists: {os.path.exists('data.csv')}")
print(f"Files in directory: {os.listdir('.')}")

# Absolute path always works:
abs_path = os.path.abspath('data.csv')
print(f"Absolute path: {abs_path}")
```

---

## Pandas-Specific Errors

### `ValueError: Length of values does not match length of index`

**Typical causes**:
1. Assigning list of wrong size: `df['new'] = [1, 2, 3]` but df has 100 rows
2. After filtering: `df_filtered = df[df.x > 5]; df_filtered['new'] = original_list` (lengths differ)

**Solution**:
```python
# Check lengths match:
print(f"DataFrame length: {len(df)}")
print(f"List length: {len(new_values)}")

# Use .loc to avoid SettingWithCopyWarning:
df.loc[:, 'new_column'] = new_values
```

### `SettingWithCopyWarning`

**Typical cause**: Chained indexing creates ambiguous reference
```python
# Wrong:
df[df.x > 5]['new'] = 10  # Warning!

# Right:
df.loc[df.x > 5, 'new'] = 10
```

---

## Numerical Errors

### `RuntimeWarning: divide by zero encountered`

**Typical causes**:
1. Division by zero: `result = x / y` when y=0
2. Log of zero: `np.log(x)` when x=0
3. Empty array: `arr.mean()` when `arr = []`

**Solutions**:
```python
# Add epsilon to avoid division by zero:
result = x / (y + 1e-10)

# Use log1p for log(1 + x):
np.log1p(x)  # Handles x=0 gracefully

# Check for zero:
if y != 0:
    result = x / y
else:
    result = 0  # or np.inf, or raise error
```

### `RuntimeWarning: invalid value encountered in` (NaN propagation)

**Typical cause**: Operations with NaN produce NaN
```python
# Example:
x = np.array([1, 2, np.nan, 4])
result = x * 2  # [2, 4, nan, 8]
```

**Solution**:
```python
# Remove NaN before operations:
x_clean = x[~np.isnan(x)]

# Or use nanmean, nansum, etc.:
result = np.nanmean(x)  # Ignores NaN values
```

---

## Bioinformatics-Specific Patterns

### Off-by-One Errors (0-based vs 1-based)

**Common in**:
- BED files (0-based, half-open)
- SAM/BAM files (1-based)
- GTF/GFF files (1-based, closed)

**Example bug**:
```python
# BED: chr1  100  200 (means positions 100-199, 0-based)
# GTF: chr1  101  200 (means positions 101-200, 1-based)
# Mixing these up causes off-by-one errors
```

### Strand Confusion (+/- strand)

**Common bug**: Forgetting to reverse complement - strand sequences

**Check**:
```python
if strand == '-':
    sequence = reverse_complement(sequence)
```

### Missing Chromosome Prefix

**Common bug**: Some files use `chr1`, others use `1`

**Check**:
```python
# Standardize:
chrom = 'chr' + chrom if not chrom.startswith('chr') else chrom
```

---

## Performance Anti-Patterns

### Growing lists in loop (O(n²) behavior)

**Slow**:
```python
result = []
for i in range(1000000):
    result.append(i ** 2)  # Can cause reallocations
```

**Fast**:
```python
# Pre-allocate:
result = np.empty(1000000)
for i in range(1000000):
    result[i] = i ** 2

# Or vectorize:
result = np.arange(1000000) ** 2
```

### Using loop instead of vectorization

**Slow**:
```python
for i in range(len(df)):
    df.loc[i, 'result'] = df.loc[i, 'x'] * 2  # Row-by-row
```

**Fast**:
```python
df['result'] = df['x'] * 2  # Vectorized
```

---

## Quick Diagnostic Commands

```python
# Check types:
print(type(obj), isinstance(obj, expected_type))

# Check values:
print(f"{obj = }")  # Python 3.8+ f-string debug syntax

# Check shapes (arrays/DataFrames):
print(arr.shape, df.shape)

# Check memory:
print(sys.getsizeof(obj))

# Check file paths:
print(os.getcwd(), os.path.exists(path))

# Check imports:
print(module.__file__)  # Shows where module is loaded from

# Check pandas dtypes:
print(df.dtypes)

# Check for NaN/None:
print(pd.isna(value), value is None)
```
