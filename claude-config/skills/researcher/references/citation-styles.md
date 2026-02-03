# Citation Styles for Bioreactor Project

**Last Updated**: 2026-01-28
**Purpose**: Quick reference for Nature-style inline citations (project standard)

---

## Nature-Style Inline Citations (Project Standard)

**Format**: Numbered superscript references in text, full bibliography at end

**Key features**:
- References numbered in order of first appearance
- Superscript numbers in text (¹, ², ³)
- Bracketed numbers in tables [1], [2] for readability
- Full citations in "References" section at document end

---

## In-Text Citation Patterns

### Single Source
```markdown
Hepatocyte oxygen consumption is 0.8 nmol/s/10⁶ cells¹.
```

### Multiple Sources (Sequential)
```markdown
Multiple studies report similar values²,³,⁴.
```
**Or** use range notation for 3+ sequential sources:
```markdown
Multiple studies report similar values²⁻⁴.
```

### Multiple Sources (Non-Sequential)
```markdown
This finding is controversial¹,⁵,⁸.
```

### Citation at End of Sentence
```markdown
The bioreactor maintained 90% cell viability for 14 days³.
```
**Punctuation rule**: Superscript comes AFTER punctuation (period, comma).

### Citation Mid-Sentence
```markdown
Primary hepatocytes⁴ show higher OCR than immortalized cell lines⁷ in equivalent culture systems.
```

### Multiple Claims, Different Sources in Same Sentence
```markdown
Porcine hepatocytes show 20% lower OCR than human cells⁶, while rat hepatocytes are 40% lower⁸.
```

---

## Table Citation Format

**Use bracketed numbers in tables** for readability:

```markdown
| Parameter | Value | Source |
|-----------|-------|--------|
| Oxygen consumption rate | 0.8 nmol/s/10⁶ cells | [1] |
| Cell viability (Day 7) | 87 ± 5% | [2] |
| Albumin secretion | 12 μg/day/10⁶ cells | [1,3] |
```

**Why brackets in tables?** Superscripts are hard to read in table cells. Brackets provide visual clarity.

---

## Full Reference Format (Nature Style)

**Journal article** (most common):
```
Author(s). Title. Journal volume, pages (year). DOI
```

**Examples**:

**Single author**:
```
1. Gerlach, J.C. Development of a hybrid liver support system. *Ann. N.Y. Acad. Sci.* **875**, 390-396 (1999). https://doi.org/10.1111/j.1749-6632.1999.tb08520.x
```

**Two authors**:
```
2. Allen, J.W. & Bhatia, S.N. Formation of steady-state oxygen gradients in vitro: application to liver zonation. *Biotechnol. Bioeng.* **82**, 253-262 (2003). https://doi.org/10.1002/bit.10569
```

**3-5 authors** (list all):
```
3. Kidambi, S., Yarmush, R.S., Novik, E., Chao, P. & Yarmush, M.L. Oxygen-mediated enhancement of primary hepatocyte metabolism. *Proc. Natl. Acad. Sci. U.S.A.* **106**, 15714-15719 (2009). https://doi.org/10.1073/pnas.0906820106
```

**6+ authors** (list first author + et al.):
```
4. Wang, L. et al. Oxygen transport limitations in high-density hepatocyte hollow fiber bioreactors. *Biotechnol. Bioeng.* **120**, 1456-1468 (2023). https://doi.org/10.1002/bit.28345
```

**Journal name formatting**:
- Italicize journal name
- Abbreviate according to ISO 4 standard (e.g., *Proc. Natl. Acad. Sci. U.S.A.*, not *Proceedings of the National Academy of Sciences*)
- For common journals, use standard abbreviations (see "Common Journal Abbreviations" below)

**Volume and page formatting**:
- Bold volume number
- Regular font for page range
- Include DOI when available

---

## Common Journal Abbreviations

| Full Name | Abbreviation |
|-----------|-------------|
| Nature | *Nature* (no abbreviation) |
| Science | *Science* (no abbreviation) |
| Cell | *Cell* (no abbreviation) |
| Proceedings of the National Academy of Sciences | *Proc. Natl. Acad. Sci. U.S.A.* |
| Biotechnology and Bioengineering | *Biotechnol. Bioeng.* |
| Hepatology | *Hepatology* (no abbreviation) |
| Gastroenterology | *Gastroenterology* (no abbreviation) |
| The Lancet | *Lancet* |
| New England Journal of Medicine | *N. Engl. J. Med.* |
| Journal of Biological Chemistry | *J. Biol. Chem.* |
| Tissue Engineering Part A | *Tissue Eng. Part A* |
| Biomaterials | *Biomaterials* (no abbreviation) |
| Artificial Organs | *Artif. Organs* |
| Transplantation | *Transplantation* (no abbreviation) |
| PLOS One | *PLoS One* (note capitalization) |

**Unsure of abbreviation?** Check the journal's official abbreviation on its website or use the [NLM Catalog](https://www.ncbi.nlm.nih.gov/nlmcatalog).

---

## Special Cases

### Book Chapter
```
5. Berthiaume, F., Moghe, P.V., Toner, M. & Yarmush, M.L. Effect of extracellular matrix topology on cell structure, function, and physiological responsiveness: hepatocytes cultured in a sandwich configuration. *FASEB J.* **10**, 1471-1484 (1996). https://doi.org/10.1096/fasebj.10.13.8940293
```

### Article with no DOI (older papers)
```
6. Berry, M.N. & Friend, D.S. High-yield preparation of isolated rat liver parenchymal cells. *J. Cell Biol.* **43**, 506-520 (1969). PMID: 4901131
```
Use PubMed ID (PMID) when DOI is unavailable.

### Preprint (bioRxiv, medRxiv)
```
7. Smith, J., Jones, A. & Williams, B. Novel hollow fiber membrane bioreactor design. *bioRxiv* (2024). https://doi.org/10.1101/2024.01.15.575382
```
Include "Preprint" note if citing in critical context (preprints not peer-reviewed).

### Patent (AVOID in citations when possible)
```
8. Inventor, A.B. & Co-Inventor, C.D. Apparatus for bioartificial liver. US Patent 7,628,910 (2009).
```
**Note**: Project policy is to **NOT download patent PDFs**. Cite patents for historical context only, not for quantitative values (use peer-reviewed sources instead).

### Conference Abstract (USE SPARINGLY)
```
9. Researcher, X. High-density hepatocyte culture in dual-lumen fibers. *Hepatology* **70** (Suppl. 1), 245A (2019). [Abstract]
```
**Warning**: Abstracts lack peer review and full methodological detail. Prefer full papers when available.

---

## Citation Self-Check Checklist

Before submitting any literature review or technical document, verify:

- [ ] Every quantitative claim has an inline citation (superscript)
- [ ] Tables use bracketed citations [1] instead of superscripts
- [ ] Citations numbered in order of first appearance
- [ ] All inline citations have corresponding entries in References section
- [ ] Reference format matches Nature style (Author(s). Title. *Journal* **vol**, pages (year). DOI)
- [ ] Journal names italicized and properly abbreviated
- [ ] Volume numbers in bold
- [ ] DOIs included when available (preferred over PMID alone)
- [ ] No uncited numbers, percentages, rates, or measurements

---

## Inline Citation Syntax in Markdown

**How to create superscript citations**:

**Method 1**: Unicode superscript characters (⁰¹²³⁴⁵⁶⁷⁸⁹)
```markdown
Hepatocyte OCR is 0.8 nmol/s/10⁶ cells¹.
```
**Pros**: Renders correctly in all markdown viewers, GitHub, and plain text
**Cons**: Slightly harder to type (use character map or shortcuts)

**Method 2**: HTML `<sup>` tags
```markdown
Hepatocyte OCR is 0.8 nmol/s/10⁶ cells<sup>1</sup>.
```
**Pros**: Easy to type
**Cons**: Doesn't render in plain text, less clean in markdown source

**Project standard**: Use **Unicode superscripts** (Method 1) for consistency.

**Quick reference for typing superscripts**:
- macOS: Character Viewer (Ctrl+Cmd+Space) → search "superscript"
- Linux: Compose key sequences or unicode input
- Windows: Character Map or Alt codes
- **Or copy-paste**: ⁰¹²³⁴⁵⁶⁷⁸⁹

---

## Example Reference Section

```markdown
## References

1. Gerlach, J.C., Encke, J., Hole, O., Müller, C., Ryan, C.J. & Neuhaus, P. Bioreactor for a larger scale hepatocyte in vitro perfusion. *Transplantation* **58**, 984-988 (1994). PMID: 7940731

2. Kidambi, S., Yarmush, R.S., Novik, E., Chao, P., Yarmush, M.L. & Nahmias, Y. Oxygen-mediated enhancement of primary hepatocyte metabolism, functional polarization, gene expression, and drug clearance. *Proc. Natl. Acad. Sci. U.S.A.* **106**, 15714-15719 (2009). https://doi.org/10.1073/pnas.0906820106

3. Wang, L., Chen, Y., Zhang, M. & Liu, H. Oxygen transport limitations in high-density hepatocyte hollow fiber bioreactors. *Biotechnol. Bioeng.* **120**, 1456-1468 (2023). https://doi.org/10.1002/bit.28345

4. Allen, J.W. & Bhatia, S.N. Formation of steady-state oxygen gradients in vitro: application to liver zonation. *Biotechnol. Bioeng.* **82**, 253-262 (2003). https://doi.org/10.1002/bit.10569

5. Hay, P.D., Veitch, A.R., Smith, M.D., Cousins, R.B. & Gaylor, J.D. Oxygen transfer in a diffusion-limited hollow fiber bioartificial liver. *Artif. Organs* **24**, 278-288 (2000). https://doi.org/10.1046/j.1525-1594.2000.06442.x
```

---

## Common Mistakes to Avoid

### Mistake 1: Missing Citations for Quantitative Claims
❌ `Hepatocyte oxygen consumption is approximately 0.8 nmol/s/10⁶ cells.`
✅ `Hepatocyte oxygen consumption is approximately 0.8 nmol/s/10⁶ cells⁴.`

### Mistake 2: Inconsistent Numbering
❌ Using same number for different sources, or skipping numbers
✅ Number references sequentially in order of first appearance

### Mistake 3: Citation Before Punctuation
❌ `The bioreactor maintained viability³.`
✅ `The bioreactor maintained viability.³`
**Actually correct**: `The bioreactor maintained viability³.` (superscript AFTER punctuation in Nature style)

### Mistake 4: Incomplete Reference Information
❌ `Wang et al. Oxygen transport in bioreactors. (2023).`
✅ `Wang, L., Chen, Y., Zhang, M. & Liu, H. Oxygen transport limitations in high-density hepatocyte hollow fiber bioreactors. *Biotechnol. Bioeng.* **120**, 1456-1468 (2023). https://doi.org/10.1002/bit.28345`

### Mistake 5: Not Italicizing Journal Names
❌ `Biotechnol. Bioeng. 120, 1456-1468 (2023).`
✅ `*Biotechnol. Bioeng.* **120**, 1456-1468 (2023).`

---

## Tools and Resources

**Citation management**:
- PubMed: Export citations in various formats (use "Citation" button)
- DOI resolver: https://doi.org/ (paste DOI to get to paper)
- CrossRef: https://crossref.org/ (DOI lookup and metadata)

**Journal abbreviations**:
- NLM Catalog: https://www.ncbi.nlm.nih.gov/nlmcatalog
- ISO 4 abbreviation tool: https://www.issn.org/services/online-services/access-to-the-ltwa/

**Superscript shortcuts**:
- Copy-paste from here: ⁰¹²³⁴⁵⁶⁷⁸⁹
- macOS Character Viewer: Ctrl+Cmd+Space
- Or use HTML `<sup>1</sup>` tags if needed

---

**Remember**: Every quantitative claim in a review document or technical report MUST have a citation. This is non-negotiable for scientific credibility and engineering reliability.
