# Drafts

Retired code kept for reference. Nothing here is part of the maintained
runtime (`src/gillespie/`), the test suite (`tests/gillespie`) or the model
contracts in `docs/`. It is not kept in sync with model changes.

- `omega_convergence.py`: poster figure. Ensemble mean of a standalone toy
  SSA (C, I, E; not `src/gillespie/`) against the ODE of `main.pdf` eq. 4.
  Writes to `/mnt/user-data/outputs/`.
- `scaling_calculator.py`: physical ↔ dimensionless ↔ Gillespie JSON
  conversion for the three-population model (no fibroblast parameters).
- `test_scaling_calculator.py`: its tests
  (`.venv/bin/python -m pytest drafts -q`).
