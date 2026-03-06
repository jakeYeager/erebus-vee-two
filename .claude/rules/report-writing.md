# Report Writing Rules

## Document Structure
- Numbered sections; you may expand as needed to include relevant content within subsections (i.e. related cross-topic content) that improves readability and interpretation: 
  - Abstract
  - Data Source (brief sample description)
  - Methodology 
  - Results 
  - Interpretation or Summary
  - Limitations
  - References (as applicable)

Use canonical naming references when citing other works within this project. See `.claude/rules/project-naming-and-organization.md` for convention.

### Internal Citation Format

**Inline:** Use the canonical identifier in parentheses, e.g. `(A2.A4)` or `(see A2.B6)`.

**Cross-Topic Comparison section lead-lines:** Use the format `**[Case Title] ([Case ID]):**`, e.g. `**Tectonic Regime Stratification (A2.B3):**`. The case title should match the canonical title from the reference list. Inline case ID references within the body text use the standard inline format `(A2.B3)`.

**Reference list entry:**
```
Yeager, J. ([year]). [Canonical ID]: [Title]. erebus-vee-two internal report.
```

Example:
```
Yeager, J. (2026). A2.A4: Aftershock Phase-Preference Analysis. erebus-vee-two internal report.
```

The year should match the date in the cited document's header.

## Template Header & Footers

Use these template header and footers unless otherwise specified:

**Header:**
```markdown
# [canonical name i.e. "A1.B1"]: [Title]

**Document Information**
- Author: Jake Yeager
- Version: 1.0
- Date: [current date, format: mmmm d, yyyy]
```

**Footer:**
```markdown
---

**Generation Details**
- Version: 1.0
- Generated with: Claude Code ([current model used, readable format])
```

## Images

Any data visualization image should be embedded inline (once) within their referencing section.

## Tables

All tables in a whitepaper must carry a sequentially numbered caption in the format **Table N. Description.** placed immediately above the table. Numbering is document-wide and sequential (Table 1, Table 2, …). Tables should be referenced by their number in surrounding prose where relevant (e.g., "see Table 3").


## Update Version on Content Change

When updating a report or whitepaper that has a "Version" number in the header or footer, increment the version as qualitatively appropriate. There are no specific standards for incrementing unless directly specified.

## Interpretive Analysis Standards

In any interpretive or summary statements remain as objective as possible. Guard against both confirmation bias as well as hypothesis bias (interpreting results based primarily on popular or conventional scientific discourse). 
