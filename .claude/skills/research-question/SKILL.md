---
name: research-question
description: Spawns a focused academic literature research agent for a given question. Returns a structured report with citations, key findings, a summary table, and a log of paywalled or access-blocked sources.
argument-hint: "[research question]"
---

Use the Task tool to run the `research-topic` agent in the **foreground** (not background). Foreground is required to ensure web search tools are available.

Pass the following as the agent's task prompt:

> Research question: $ARGUMENTS

When the agent completes, present its full report output to the user without summarizing.
