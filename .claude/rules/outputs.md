# Outputs Rules

# Case Outputs
- All case outputs go in subdirectories under the main topic directory
- All `.json`, `.png`, and `.md` outputs should go in `output/`
- All `.py` analysis files should go in `src/`
- All `.py` test files should go in `tests/`

## Whitepaper documents
- File name: `case-[identifier]-whitepaper.md`

## Python Modules
- Main analysis script file name: `case-[identifier]-analysis.py`
- Visualization script file name: `visualization-case-[identifier].py`
- Test suite file name: `test-case-[identifier].py` (all tests must pass)

## JSON Results
- File: `case-[identifier]-results.json`
- Content: All parameters, statistics, test results
- Format: Hierarchical (e.g., cycle → test → results)

## Images

- File name: `case-[identifier]-[description].png`
