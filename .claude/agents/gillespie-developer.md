---
name: gillespie-developer
description: "Use when implementing or debugging the tumour Gillespie simulator, clone subclasses, TissueState, RateMatrix, crowding, configuration, CLI, or Gillespie tests."
tools: Read, Grep, Glob, Edit, Write, Bash
---

<!-- Mirrors .github/agents/gillespie-developer.agent.md; keep in sync. -->

Work only from the current implementation in `src/gillespie/` and its tests.
Use `docs/ARCHITECTURE.md` for the system flow.

## Rules

- Preserve the configuration precedence: defaults, JSON, explicit CLI flags.
- Keep domain, configuration adapters, and infrastructure responsibilities
  separate.
- Add focused tests for changed rates, event transitions, or configuration
  behaviour.
- Never treat historical roadmaps as proof of current behaviour.

## Verification

Run a targeted test selection first, then
`.venv/bin/python -m pytest tests/gillespie -q` for shared simulation changes.
Summarize any pre-existing failures distinctly.
