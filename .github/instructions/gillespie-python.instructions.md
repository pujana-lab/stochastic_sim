---
name: Gillespie Python
description: "Use when editing Python under src/gillespie or tests/gillespie, especially clone rates, events, configuration, and simulation state."
applyTo: ["src/gillespie/**/*.py", "tests/gillespie/**/*.py"]
---

# Gillespie Python conventions

- Keep domain code independent of CLI and file-output concerns.
- Pass `TissueState` into effective-rate calculations so interactions use the
  current aggregate population map.
- Use `CloneFactory` for clone creation and `SimulationConfig` for model
  parameters.
- Add regression tests for every changed event or propensity.
- Do not replace per-clone event construction with type aggregation unless the
  architecture document and tests are updated together.
