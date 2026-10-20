# Project architecture

> **Maintained view.** Update this file when a package boundary, entry point,
> configuration contract, or data flow changes. Historical design discussions
> belong in `docs/archive/` or in commit history, not in this document.

## Purpose and scope

`moran` is a Python research simulator for stochastic tumour clone dynamics.
The maintained implementation is the event-driven Gillespie simulator under
`src/gillespie/`. It models four clone types:

`base` -> `mutated`  
`immune` -> `exhausted`

The repository also contains notebooks, generated research artefacts, and
historical refactoring notes. Those are supporting material, not alternate
runtime implementations.

## System flow

```text
JSON file and/or CLI arguments
              |
              v
infrastructure.config
  load_config_json -> nested cell_parameters
  cli_config      -> explicit CLI overrides
              |
              v
SimulationConfig
  validates/normalizes cell parameters
  scales parameters
  selects crowding strategy
              |
              v
TumorSimulation
  owns time, RNG, history, and the run loop
              |
              +--> CloneFactory --> Clone subclasses
              |                    (rates and clone-local state)
              |
              +--> TissueState ---- (clones + population map)
              |
              +--> RateMatrix ----- Event candidates and weighted choice
              |
              +--> CrowdingStrategy (Simple or Adapted)
              |
              v
infrastructure
  csv_output / display
  Parquet history, clone, rate, and debug artefacts
```

## Responsibility map

| Area | Responsibility | Start here |
|---|---|---|
| Domain model | Clone hierarchy, IDs, events, tissue state, rates | `src/gillespie/clone.py`, `event.py`, `tissue_state.py` |
| Simulation orchestration | Initialise state, build rates, sample and apply events | `src/gillespie/tumor_simulation.py` |
| Configuration | Defaults, scaling, JSON mapping, CLI precedence | `simulation_config.py`, `infrastructure/config/` |
| Infrastructure | CLI, Parquet writers, summary output | `src/gillespie/infrastructure/` |
| Experiments | Batch runs and analysis scripts | `scripts/`, `notebooks/` |
| Verification | Unit and integration tests | `tests/gillespie/` |

## Configuration contract

Configuration is merged in this order:

1. `SimulationConfig` defaults
2. JSON supplied with `--config`
3. explicitly supplied CLI flags

Cell-type settings are nested under `cell_parameters`, with the same field
names as [`CellTypeConfig`](../src/gillespie/simulation_config.py). This is the
only supported JSON configuration schema. The JSON loader removes optional
`_calculator` provenance metadata before constructing `SimulationConfig`; it
does not translate legacy flat fields or the former `cell_types` schema.

The dimensional-to-Gillespie calculator is implemented in
[`scaling_calculator.py`](../src/gillespie/application/scaling_calculator.py).
It emits a valid JSON configuration and records the physical and
dimensionless values in an ignored `_calculator` metadata section. `OMEGA`
is the population/system-size scale: interaction coefficients are divided by
`OMEGA` once, while dimensionless capacities are converted to
`OMEGA * K`. `AdaptedCrowding` then inflates `actual_K` by
`1 / (1 - mu / lambda)` to make the requested baseline homeostatic.

## Event lifecycle

For every simulation step:

1. The step reuses the `RateMatrix` of the current state (built at start-up or
   at the end of the previous step; `create_clone()` invalidates it).
2. A waiting time is sampled from the total propensity.
3. `RateMatrix.choose_event()` selects one event proportionally to its rate.
4. The event mutates the selected clone or creates a child clone.
5. Time and instability are updated, then `TumorSimulation._build_rate_matrix()`
   refreshes `pop_map`, asks each living clone for effective rates, creates the
   non-zero `Event` candidates and stores the effective rates in
   `TissueState.effective_rates`.
6. Optional history is recorded.

The per-clone event construction is an important current behaviour. Each
clone's rate is proportional to its own `N` (interactions read the type totals
in `pop_map`), so the events of one type sum to the type-level propensity
without double counting. Do not describe the implementation as type-aggregated
Gillespie until that design is implemented and tested.

`TissueState.snapshot()` reports `rb`, `rd`, `rm` and `re` as the **effective**
per-clone rates of the snapshot's state (the same values in the `RateMatrix`),
not the per-cell parameters; dead clones report `0.0`. The per-cell parameters
remain available on the `Clone` objects and in `save_clones_parquet`.

## Verification flow

Use the smallest relevant check:

```bash
# Full Gillespie test suite
python -m pytest tests/gillespie -q

# CLI/configuration checks
python -m pytest tests/gillespie/test_cli_config.py tests/gillespie/test_json_config.py -q

# Discover the supported command-line contract
python -m src.gillespie.infrastructure.cli --help
```

The CLI currently pauses for interactive confirmation before running. For
non-interactive experiments, construct `TumorSimulation(SimulationConfig(...))`
from Python or remove that pause as a separate, tested change.

## Change checklist

When changing the model:

- update this file if the flow, ownership, or configuration contract changes;
- add or update a focused test under `tests/gillespie/`;
- update the relevant README or command prompt if a user-facing command changes;
- record unresolved modelling decisions as an issue or a short entry in
  `docs/decisions/`, not as executable-looking code in a roadmap;
- keep generated CSV, Parquet, images, logs, and notebook checkpoints out of
  source-control changes unless they are intentional research fixtures.

## Documentation map

- [`README.md`](../README.md): concise setup and user entry points.
- [`docs/ARCHITECTURE.md`](./ARCHITECTURE.md): this maintained system view.
- [`docs/gillespie_full_diagnosis_roadmap.md`](./gillespie_full_diagnosis_roadmap.md):
  historical modelling diagnosis; verify claims against code.
- [`CLONE_TYPE_REFACTORING.md`](../CLONE_TYPE_REFACTORING.md) and
  [`TISSUE_STATE_REFACTOR.md`](../TISSUE_STATE_REFACTOR.md): historical refactor
  records, not current implementation specifications.
