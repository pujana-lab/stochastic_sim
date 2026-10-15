---
paths:
  - "src/gillespie/**/*.py"
  - "tests/gillespie/**/*.py"
---

<!-- Mirrors .github/instructions/gillespie-python.instructions.md; keep in sync. -->

# Gillespie Python conventions

- Keep domain code independent of CLI and file-output concerns.
- Pass `TissueState` into effective-rate calculations so interactions use the
  current aggregate population map.
- Use `CloneFactory` for clone creation and `SimulationConfig` for model
  parameters.
- Add regression tests for every changed event or propensity.
- Do not replace per-clone event construction with type aggregation unless the
  architecture document and tests are updated together.
