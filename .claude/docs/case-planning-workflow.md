# Case Planning Workflow

## Phase 1 - Review Previous Summary

An "initial" planning document will be made:
1. It will be located at `topic-<name>/.claude/docs/planning-initial.md`
2. It will have a generalized hypothesis or purpose statement
3. It will use the previous topic `final-summary.md` to inform the purpose statement
4. It will have a generalized case testing suite to address or satisfy the purpose

> **Lifecycle:** This document is an intermediate artifact. It becomes a historical reference once the Phase 3 refinement document is approved and Phase 4 begins. It is not loaded in active work sessions after that point.

## Phase 2 - Academic Publication References

Research will be conducted on existing academic publications that relate to, support or question what is described in Phase 1. In a concise, descriptive manner these academic references will be documented, for later reference at `topic-<name>/.claude/docs/academic-research.md`

## Phase 3 - Extend Refinement

A second "refinement" planning document will be created:
1. It will be located at `topic-<name>/.claude/docs/planning-refinement.md`
2. Its use is to incorporate the academic research
3. Its use is to make iterative changes to the hypothesis or intent statements
4. Its use is to make iterative changes to the case testing suite
5. It will allow the initial case planning document to remain unchanged during the refinement process

> **Lifecycle:** This document is an intermediate artifact. Once approved and Phase 4 begins, it is superseded by the topic CLAUDE.md and individual case documents. It is not loaded in active work sessions after that point.

## Phase 4 - Build Individual Case Working Documents

Once the Phase 3 refinement document is completed and approved, the following steps must be taken:

### Step 0 — Archive Planning Documents

Add the following status header to the top of both `planning-initial.md` and `planning-refinement.md`:

```
> **Status: Superseded.** See `topic-<name>/.claude/CLAUDE.md` for current case status.
```

This prevents future sessions from treating planning documents as active context.

### Step 1 — Update the topic's CLAUDE.md

The use-case for this CLAUDE.md is topic-level operational context. It should carry general information about the topic, and the essentials information for managing the current state of the cases. Update `topic-<name>/.claude/CLAUDE.md` with:
1. A brief topic intent paragraph: the core hypothesis, purpose, and data being used
2. A case reference list — one line per case, with on-demand file link:
   ```
   Only read this file if you need [case title]: `topic-<name>/.claude/spec/case-<identifier>-spec.md`
   ```
3. A standards section: output header/footer templates and interpretive analysis standards
4. Cross-topic reference links (on-demand, same pattern)

Do NOT copy planning document content into CLAUDE.md. The planning docs are intermediate artifacts; CLAUDE.md is the operational context.

### Step 2 — Create a Case Specification File

A "case specification" file will be made for each case:
1. It will be located at `topic-<name>/.claude/spec/case-<identifier>-spec.md`
2. It serves as the prompt entry point for an agent to build the analysis code — it is a write-once implementation document, not ongoing context
3. It should be structured and have all necessary content such that an agent execute it outside of the main context 
   1. It should eager load any data files, as well as any other applicable completed case files that has meaningful information to its analysis or whitepaper comparative sections
4. The final required step of every spec — after all other deliverables are satisfied — must explicitly:
   - Update `topic-<name>/.claude/docs/topic-summary.md` with key results and final status
   - Update the case status in `topic-<name>/.claude/CLAUDE.md`

> **Skill:** Created alongside the detail file by `/new-case`. Invoke separately only if the spec stub needs to be added or recreated independently.

> **Lifecycle:** Spec files are implementation prompts only. Spec files should not be eagerly loaded after case completion.

## Phase 5 - Final Summary

Once all cases have been analyzed and marked as "complete", "blocked", or "abandoned" create a "topic summary statement":
1. It will be created in the file located at `topic-<name>/.claude/docs/topic-summary.md`
2. This will incorporate the key findings from the other cases documented in that file
3. After the summary is written, update the topic CLAUDE.md overall status to: Complete
4. Add a section to the root `README.md` for this topic listing all completed case whitepapers:
   - Use the topic title and brief descriptor as the section heading
   - Link each entry to its `topic-<name>/output/case-<identifier>-whitepaper.md` file using the whitepaper's H1 title as the link text
5. It will be used to inform the planning of the next topic
