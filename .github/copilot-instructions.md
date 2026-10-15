# Project development guidelines

## Source of truth

- The maintained runtime is `src/gillespie/`; do not introduce new work in
  historical or generated artefacts.
- Read [`docs/ARCHITECTURE.md`](../docs/ARCHITECTURE.md) before changing
  simulation flow, configuration, or output.
- Keep historical diagnosis documents clearly separate from current contracts.

## Python and model changes

- Preserve type hints and the existing domain vocabulary (`CloneId`,
  `TissueState`, `Event`, `RateMatrix`, `SimulationConfig`).
- Reuse the JSON/CLI configuration adapters instead of parsing configuration in
  the simulation domain.
- Make model changes explicit in tests; do not silently change rate formulas,
  event semantics, or scaling.

## Verification

- Run the narrowest relevant `pytest` selection, then the full
  `tests/gillespie` suite when shared simulation code changes.
- Prefer deterministic seeds in tests and keep generated outputs out of commits.
- Update `docs/ARCHITECTURE.md` whenever a boundary or data flow changes.
