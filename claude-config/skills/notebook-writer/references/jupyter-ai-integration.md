# AI assistance inside Jupyter

Modern Jupyter environments (JupyterLab 4.0+, JetBrains IDEs) provide AI-powered assistance inside the notebook itself.

## %%ai magic commands

The `%%ai` cell magic enables code generation and analysis directly in notebooks:

```python
# %%
# %load_ext jupyter_ai_magics

# %%
%%ai chatgpt
Generate a function to calculate the Pearson correlation coefficient between two arrays
```

Typical uses: boilerplate code, data-transformation snippets, questions about DataFrames or arrays, error suggestions, docstrings.

## Context that improves results

Provide, in a markdown cell near the request:

1. **API documentation** for specialized libraries (scanpy, pydeseq2, biopython)
   ```python
   # Example: scanpy.pp.filter_cells(data, min_genes=200)
   ```

2. **Dataset description** — shape, columns, data types
   ```python
   # RNA-seq counts matrix: 20,000 genes × 5,000 cells
   # AnnData object: .X (sparse CSR matrix), .obs (cell metadata), .var (gene metadata)
   ```

3. **Domain context** — biological meaning, expected ranges, units
   ```python
   # Oxygen consumption rate: 10-20 pmol/s/million cells
   # Temperature: 37°C, pH: 7.4
   ```

## Chat UI assistance

JupyterLab's chat interface suits exploratory questions ("What's the best way to normalize this data?"), code review ("Does this analysis handle missing values correctly?"), visualization requests ("Create a heatmap of the top 50 variable genes"), and explanations of an existing cell.

## Where in-notebook AI assistance fits

Good fit: boilerplate (imports, data-loading templates), exploratory analysis (quick plots, summary statistics), learning new library syntax, generating test data.

Write manually instead for core analysis logic (hypothesis testing, modeling), publication-quality figures needing fine-grained control, performance-critical sections, and complex domain-specific algorithms.

Generated code is worth checking against library syntax (APIs change frequently), statistical appropriateness, and correct handling of biological data (species, units, measurement context).

## JetBrains AI Assistant

For notebooks in PyCharm/DataSpell:

- **Explain cell**: Alt+Enter -> "Explain"
- **Create visualization**: generate plots from data descriptions
- **Edit cell**: Alt+Enter -> "AI Actions"
- **Fix errors**: suggestions for runtime errors

Access: right-click cell -> "AI Assistant", or the AI chat sidebar.
