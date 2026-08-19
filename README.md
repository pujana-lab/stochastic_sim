# Tumor Gillespie Simulator

A stochastic Gillespie-based tumor clone simulator for a simplified tissue model. The implementation in `src/gillespie/` tracks four clone types and computes event rates from clone-specific parameters, crowding, and cell-cell interactions.

## Overview

The simulation models four clone classes:

- `base`: healthy (wild type) cells
- `mutated`: tumor cells
- `immune`: immune cells
- `exhausted`: exhausted immune cells

Each clone type has its own base parameters defined in `simulation_config.py`. Effective rates use the full tissue state kept in `TissueState` and may include interaction terms between clone types. Code architecture allows for additional clone classes (such as intermediate pre-neoplastic states) to be implemented in a straightforward way.

## Installation and setup

1. Create and activate a Python virtual environment:

```bash
python -m venv venv
source venv/bin/activate
```

2. Install the package in editable mode:

```bash
pip install -e .
```

3. Install project dependencies:

```bash
pip install -r requirements.txt
```

4. If you prefer the Makefile helpers:

```bash
make create_venv
source venv/bin/activate
make install-notebook
```

The `requirements.txt` file includes test and notebook dependencies.

## Command line usage

Once the environment is ready, run the Gillespie CLI directly:

```bash
python src/gillespie/infrastructure/cli.py --help
```

The CLI accepts the same simulation parameters used by `SimulationConfig`, such as:

```bash
python src/gillespie/infrastructure/cli.py \
  --N0 50 \
  --lambda0 0.25 --mu0 0.25 --nu0 0.0 \
  --T-max 20 \
  --seed 42 \
  --save-history results.csv \
  --save-debug debug.csv
```

This will run a Gillespie trajectory and write CSV output for the history and debug information.

## Makefile usage

The repository includes Make targets for common workflows:

- `make test`
  - run the test suite against `tests`
- `make test-cov`
  - run tests with coverage reporting
- `make create_venv`
  - create the Python virtual environment `venv`
- `make install-notebook`
  - install the package in editable mode
- `make gillespie-homeostasis`
  - run a homeostasis scenario
- `make gillespie-tumour-growth`
  - run the tumor growth scenario with instability and mutation
- `make gillespie-crowding`
  - run the crowding scenario with logistic limits
- `make gillespie-all`
  - run all three Gillespie scenarios
- `make gillespie-plot`
  - run the CLI and the notebook plotting scripts

After installation, use either the CLI or the Makefile targets depending on whether you want a one-off run or a predefined scenario.

## Model equations

### Event types

Each alive clone can participate in the following events:

- `BIRTH`: one cell divides
- `DEATH`: one cell dies
- `MUTATION`: a new mutant clone is created and the parent clone loses one cell
- `EXHAUSTION`: an immune clone becomes exhausted and one exhausted cell is produced

### Base effective rate formulas

For a clone with population `N`, the base effective rates are:

- birth rate:
  `r_B = λ · N · C(t)`
- death rate:
  `r_D = μ · N`
- mutation rate:
  `r_M = ν · N · (1 + instability)`
- exhaustion rate:
  `r_E = ω · N`

The crowding factor `C(t)` is computed by `CrowdingStrategy`. If `use_logistic` is disabled, `C(t) = 1`.

### Crowding

Crowding is implemented as:

```
C(t) = max(0, 1 - N_crowd / K_t)
```

where:

- `N_crowd` is the competitive population for that clone type
- `K_t` is the effective carrying capacity

The available strategies are:

- `SimpleCrowding`
  - computes `K_t = max(K_min, K)`
- `AdaptedCrowding`
  - if `λ > μ`: `K_t = max(K_min, ceil(K / (1 - μ / λ)))`
  - otherwise: `K_t = ∞`

`AdaptedCrowding` is selected when `use_logistic_adapted=True`.

### Clone-specific interaction terms

Some clone types modify the base formula with interactions.

#### `MutatedClone`

Immune killing adds a tumor-immune interaction to death:

```
r_D(mutated) = μ_mutated · N_mutated + θ_I · N_mutated · N_immune
```

#### `ImmuneClone`

Immune birth includes activation by mutated cells:

```
r_B(immune) = λ_immune · N_immune · C(t) + β · N_immune · N_mutated
```

Immune exhaustion scales with tumor burden:

```
r_E(immune) = ω_immune · N_immune · N_mutated
```

#### `ExhaustedClone`

Exhausted clones do not divide:

```
r_B(exhausted) = 0
```

### Instability dynamics (Work in progress)

Each clone advances instability each step using:

```
instability += (base_instability_buildup + buildup) · Δt
```

The current implementation uses `instability` only to scale mutation rate via `1 + instability`. Instability and buildup parameters are set to 0 by default for the time being until full support is implemented.

### Gillespie algorithm

`TumorSimulation` implements the event-driven Gillespie step:

1. Update `TissueState.pop_map`
2. Build the list of event rates in `RateMatrix`
3. Compute total rate `R = Σ rate`
4. Sample waiting time
   `τ = -log(u) / R`, with `u ∼ Uniform(0, 1)`
5. Choose one event weighted by event rate
6. Apply that event to the selected clone
7. Advance simulation time and record history

The run loop stops when the total rate is zero, there are no more events, or `T_max` is reached.

## Code architecture

### `src/gillespie/simulation_config.py`

Defines global simulation parameters and default clone-type parameters via `CellTypeConfig`.

- scales `beta` and `theta_I` by `OMEGA`
- scales carrying capacity `K`
- selects `AdaptedCrowding` or `SimpleCrowding`

For simplicity `OMEGA` is chosen to be the homeostatic equilibrium value for WT population and user input values for `K`  for each clone tyoe correspond to fractions of this value. 

### `src/gillespie/clone.py`

Defines the `Clone` base class and specialized subclasses:

- `WildTypeClone` (`base`)
- `MutatedClone` (`mutated`)
- `ImmuneClone` (`immune`)
- `ExhaustedClone` (`exhausted`)

Each clone defines:

- `birth_rate_effective(tissue_state)`
- `death_rate_effective(tissue_state)`
- `mutation_rate_effective(tissue_state)`
- `exhaustion_rate_effective(tissue_state)`
- `crowding_numerator(tissue_state)`

### `src/gillespie/crowding_strategy.py`

Implements logistic crowding. The `crowding()` method computes the current crowding multiplier and the concrete strategy decides the effective carrying capacity.

### `src/gillespie/clone_factory.py`

Creates clone instances from configuration. It supports registered clone types and a special `mutated_test` path used for development.

### `src/gillespie/tissue_state.py`

Encapsulates current clone populations and provides:

- `pop_map` of clone type counts
- `total_population()`
- `get_clones_by_type()`
- `snapshot()` for history

### `src/gillespie/rate_matrix.py`

Collects candidate events and their rates, computes the total rate, and chooses a single event proportional to the rates.

### `src/gillespie/tumor_simulation.py`

Executes the simulation loop:

- initializes clones and `TissueState`
- builds per-clone events every step
- samples event times and selects events
- updates state and optional history

## Usage example

```python
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation

config = SimulationConfig(
    T_max=1000,
    seed=42,
    use_logistic=True,
    use_logistic_adapted=True,
)

sim = TumorSimulation(config)
results = sim.run()
```

In order to run a quick test and plot the results one can run:
```bash
make gillespie_plot
```
## Important notes

- `TissueState` stores state and should remain separate from simulation logic.
- `RateMatrix` is currently a simple event list; it does not implement tau-leap.
- `MUTATION` events create a new clone of the configured `next_mutation` type.
- `EXHAUSTION` events kill one source clone and increase the exhausted clone count.

## Testing

```bash
python -m pytest tests/gillespie
```

or

```bash
make test
```