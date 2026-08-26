# Developer notes

## Communication style
Be concise and direct. English by default.

## Two codebases
- **`main.py` + `preliminar/`** — earlier "Moran" model (YAML config, XLSX output). Kept for reference; not under active development.
- **`src/gillespie/`** — current Gillespie SSA implementation (event-driven, CSV output). This is where new work should happen.

For Gillespie work, use the `gillespie-developer` skill (`.opencode/skills/gillespie-developer/`).

## Essential commands
```
make test                    # run all tests
make test-cov                # run tests with coverage
make gillespie-homeostasis   # Gillespie scenario: homeostasis
make gillespie-tumour-growth # Gillespie scenario: tumor growth
make gillespie-crowding      # Gillespie scenario: crowding
```

Gillespie CLI: `venv/bin/python -m src.gillespie.infrastructure.cli [args]`

## Known issues
- `use_logistic` must always be `True` (disabling it causes unbounded exponential growth)
- The CLI does not yet pass through custom config arguments; it always builds a default `SimulationConfig()`
- Three test files under `tests/gillespie/` are currently broken
- `tests/gillespie/test_clone_factory.py` is a placeholder (empty)
- Some in-code comments and TODOs are still in Spanish; these reflect real open questions and are worth reading, not just noise
