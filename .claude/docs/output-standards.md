# Standards & Conventions

## Metadata Standards

### Header & Footer Templates

Use these template header and footers unless otherwise specified:

**Header:**
```markdown
# Case [identifier]: [Title]

**Document Information**
- Version: 1.0
- Date: [current date]
```

**Footer:**
```markdown
---

**Generation Details**
- Version: 1.0
- Date: [current date]
- Generated with: Claude Code ([current model used])
```

When updating report or whitepaper that has a "Version" number in the header or footer, increment the version as qualitatively appropriate. There are no specific standards for incrementing unless directly specified.

### Interpretive Analysis Standards

In any interpretive statements, guard against both confirmation bias as well as hypothesis bias (interpreting results based primarily on popular or conventional scientific discourse). Remain as objective as possible.

### Output Requirements

Unless specified, fallback to these standards:

1. **Python Modules**
   - Main analysis script: `case-[identifier]-analysis.py`
   - Visualization script: `visualization-case-[identifier].py`
   - Test suite: `tests/test-case-[identifier].py` (all tests must pass)

2. **JSON Results**
   - File: `output/case-[identifier]-results.json`
   - Content: All parameters, statistics, test results
   - Format: Hierarchical (e.g., cycle → test → results)

3. **Visualizations**
   - PNG format, 300+ DPI suitable for publication
   - Filename: `output/case-[identifier]-[description].png`
   - Include titles, axis labels, legends, statistical annotations

4. **Documentation**
   - File: `output/case-[identifier]-whitepaper.md`
   - Structure: Methodology → Results → Interpretation → Limitations
   - Include all visualizations with captions

### Chart & Graph Visualizations

**Images:** Any data visualization image should be embedded inline (once) within their referencing section.

**Diagrams and charts:** When displaying complex data prefer visuals in JPG, PNG or SVG formats (like when using `matplotlib.pyplot`). ASCII diagrams should be avoided unless they are simple file structure charts or workflow/decision charts.

**Single-color charts** (histograms, bar charts): Use `steelblue` as the default bar color.

**Gradient/two-tone scale graphics** (heatmaps, density plots): Use a white=least, [color]=greatest scheme
- Color for primary or full populations: Red (fallback to Orange)
- Color for secondary or subpopulations: Blue, Purple, Green
- Never use two-tone scales like "viridis" unless specified

### Code Quality

Use these defaults unless specified or redundant:
- File naming convention is kebab-case: `file-name.md`
- Type hints on all functions
- Docstrings explaining purpose and parameters
- Error handling for data issues
- Logging for important steps
- Reproducible random seeds (if applicable)
