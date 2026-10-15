---
name: project-architecture
description: "Maintain and navigate the project architecture, documentation map, module boundaries, simulation data flow, and stale guidance. Use when the project evolves or architecture needs refreshing."
argument-hint: "Describe the architectural change to trace"
---

# Project architecture workflow

1. Read [`ARCHITECTURE.md`](../../../docs/ARCHITECTURE.md).
2. Inspect the affected source and tests under `src/gillespie/` and
   `tests/gillespie/`.
3. Trace configuration and output adapters if the change crosses the CLI.
4. Update the architecture document only with verified current behaviour.
5. Link historical or unresolved material instead of copying it into the
   maintained view.
6. Run the smallest relevant test command and report the result.
