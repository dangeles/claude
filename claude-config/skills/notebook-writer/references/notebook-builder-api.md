# notebook_builder API

Many projects provide `src/utils/notebook_builder.py` with helper functions for programmatic notebook creation. If the project has no such module, write the Jupytext Markdown directly.

## Contents

- [create_notebook_markdown](#create_notebook_markdown)
- [create_parameter_sweep_notebook](#create_parameter_sweep_notebook)
- [create_analysis_report_notebook](#create_analysis_report_notebook)
- [validate_notebook](#validate_notebook)

## create_notebook_markdown

```python
create_notebook_markdown(
    title: str,
    cells: List[Dict[str, str]],
    output_path: Path,
    kernelspec: Optional[Dict] = None
) -> Path
```

**Parameters:**
- `title`: Notebook title (becomes H1 header)
- `cells`: List of dicts with `'type'` (`'code'` or `'markdown'`) and `'content'`
- `output_path`: Where to save `.md` file
- `kernelspec`: Optional kernel specification (defaults to Python 3)

**Example:**
```python
from pathlib import Path
from src.utils.notebook_builder import create_notebook_markdown

cells = [
    {'type': 'markdown', 'content': '## Introduction\n\nThis analysis...'},
    {'type': 'code', 'content': 'import numpy as np'},
    {'type': 'code', 'content': 'x = np.linspace(0, 10)\nprint(x)'}
]

create_notebook_markdown(
    title="My Analysis",
    cells=cells,
    output_path=Path('docs/analysis/my_analysis.md')
)
```

## create_parameter_sweep_notebook

```python
create_parameter_sweep_notebook(
    param_name: str,
    param_range: str,
    calculation_code: str,
    output_path: Path
) -> Path
```

Creates a notebook with imports (numpy, matplotlib, pandas), the parameter range definition, your calculation code, and visualization boilerplate.

**Example:**
```python
from pathlib import Path
from src.utils.notebook_builder import create_parameter_sweep_notebook

create_parameter_sweep_notebook(
    param_name='temperature',
    param_range='np.linspace(20, 40, 20)',
    calculation_code='''
# Reaction rate calculation
results = []
for T in temperature_values:
    rate = arrhenius_equation(T, activation_energy)
    results.append(rate)
''',
    output_path=Path('analysis/temperature_sweep.md')
)
```

## create_analysis_report_notebook

```python
create_analysis_report_notebook(
    analysis_title: str,
    sections: List[Dict[str, str]],
    output_path: Path
) -> Path
```

**Section dict keys:**
- `title`: Section heading (required)
- `description`: Explanatory text (optional)
- `code`: Code to execute (optional)
- `interpretation`: Results interpretation (optional)

**Example:**
```python
from src.utils.notebook_builder import create_analysis_report_notebook

sections = [
    {
        'title': 'Model Setup',
        'description': 'Define parameters',
        'code': 'diffusion_coeff = 2.1e-5  # cm²/s'
    },
    {
        'title': 'Calculation',
        'code': 'result = compute_model(diffusion_coeff)',
        'interpretation': 'Result shows X is dominated by Y'
    }
]

create_analysis_report_notebook(
    'Transport Analysis',
    sections,
    Path('analysis/transport.md')
)
```

## validate_notebook

```python
validate_notebook(notebook_path: Path) -> bool
```

Validates `.ipynb` structure using nbformat. Returns `True` if valid, raises an exception if invalid.

```python
from pathlib import Path
from src.utils.notebook_builder import validate_notebook

validate_notebook(Path('analysis/notebook.ipynb'))
# Returns True or raises ValidationError
```
