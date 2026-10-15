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
  load_config_json -> flatten_cell_types
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

Cell-type settings are nested under `cell_parameters`. JSON aliases such as
`lambda_0` are normalised by
[`flatten_cell_types`](../src/gillespie/infrastructure/config/json_config.py).
Do not add a second configuration format without updating that adapter and its
tests.

## Event lifecycle

For every simulation step:

1. `TissueState.update_pop_map()` refreshes aggregate populations.
2. `TumorSimulation._build_rate_matrix()` asks each living clone for effective
   rates and creates non-zero `Event` candidates.
3. A waiting time is sampled from the total propensity.
4. `RateMatrix.choose_event()` selects one event proportionally to its rate.
5. The event mutates the selected clone or creates a child clone.
6. Time, instability, and optional history are updated.

The per-clone event construction is an important current behaviour. Do not
describe the implementation as type-aggregated Gillespie until that design is
implemented and tested.

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
