# Moran

## Comm style
Speak caveman, English only unless asked otherwise.

## Two codebases
- **`main.py` + `preliminar/`** — old Moran (YAML → XLSX)
- **`src/gillespie/`** — new Gillespie SSA (event-driven, CSV)

For Gillespie work, use `gillespie-developer` skill.

## Essential commands
```
make test              # all tests
make test-cov          # coverage
make gillespie-homeostasis   # Gillespie scenario
make gillespie-tumour-growth
make gillespie-crowding
```

Gillespie CLI: `venv/bin/python -m src.gillespie.infrastructure.cli [args]`

## Quick quirks
- `use_logistic` MUST be True
- CLI ignores config args (creates default `SimulationConfig()`)
- 3 broken test files in `tests/gillespie/`
- `test_clone_factory.py` is empty
- Spanish TODOs everywhere — read them, not noise
