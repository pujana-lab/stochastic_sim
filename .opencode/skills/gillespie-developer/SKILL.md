---
name: gillespie-developer
description: Use when working on src/gillespie/ — tumour clone Gillespie SSA simulator. Clone subclasses, RateMatrix, TissueState, crowding strategies, CLI quirks, broken tests, event rates, and simulation pipeline. Do NOT use for the old preliminar/ Moran process code (main.py).
---

# Gillespie Tumour Simulator — Developer Skill

## Two codebases, one repo

| Codebase | Entrypoint | Config | Output | Tests |
|---|---|---|---|---|
| Old Moran (`preliminar/` + `main.py`) | `venv/bin/python main.py --config configs/*.yaml` | YAML | XLSX | `preliminar/tests/test_population_kata_*.py` |
| **New Gillespie** (`src/gillespie/`) | `venv/bin/python -m src.gillespie.infrastructure.cli` | `SimulationConfig` dataclass | CSV | `tests/gillespie/` |

**Do not mix them.** If `main.py` or `preliminar/` is mentioned, redirect to the old system. The Gillespie code lives entirely under `src/gillespie/`.

## Architecture overview

```
TumorSimulation
  ├── CloneFactory          → creates Clone subclasses from type string
  ├── TissueState           → clones dict + pop_map (computed per type)
  ├── RateMatrix            → list of Event objects, total rate, choose_event()
  ├── CrowdingStrategy      → SimpleCrowding | AdaptedCrowding
  └── SimulationConfig      → frozen dataclass, all model parameters
```

### Data flow per step()

1. `_build_rate_matrix()` — update pop_map, iterate alive clones, ask each clone for its 4 effective rates (`birth`, `death`, `mutation`, `exhaustion`). Rates are cached **per clone type** (first clone of a type computes, others reuse). Each non-zero rate becomes an `Event` in the `RateMatrix`.
2. `_sample_waiting_time(total_rate)` — `tau = -ln(u) / R` where `u ~ U(0,1)`
3. Clamp `tau` so `t + tau ≤ T_max`, advance instability for all clones
4. `rate_matrix.choose_event(u)` — inverse-transform sampling via `np.searchsorted`
5. `_apply_event(event)` — dispatch by `EventType`:
   - `BIRTH` → `clone.divide()` (N += 1)
   - `DEATH` → `clone.kill()` (N = max(0, N-1))
   - `MUTATION` → `_introduce_mutation()`: create child clone with N=2 from `clone.next_child_id()`, kill parent (N-1)
   - `EXHAUSTION` → `_induce_exhaustion()`: kill the immune clone, divide the exhausted clone `(-2,)`
6. Snapshot current tissue state into history

### Clone class hierarchy

```
Clone (base — _registry: Dict[str, Type[Clone]], __init_subclass__(clone_type=...))
  ├── WildTypeClone  (clone_type="base")       — WT cells, next_mutation="mutated"
  ├── MutatedClone   (clone_type="mutated")    — tumour cells, immune-killing death term
  ├── ImmuneClone    (clone_type="immune")     — immune cells, activation boost + exhaustion
  └── ExhaustedClone (clone_type="exhausted")  — exhausted T-cells, no division
```

- Registration is automatic: `class Foo(Clone, clone_type="foo")` → `Clone._registry["foo"] = Foo`
- **No enum.** Types are plain strings everywhere: `pop_map: Dict[str, int]`, `clone.get_type() -> str`
- `CloneId = Tuple[int, ...]` from `cloneId.py` — e.g. `()`, `(1,)`, `(1, 2)`. Initial IDs: `()`=base, `(-3,)`=mutated_root, `(-1,)`=immune, `(-2,)`=exhausted

### Per-subclass effective rate formulas

#### WildTypeClone ("base")
- `birth_rate_effective` = `lambda0 * N_base * crowding`
- `death_rate_effective` = `mu0 * N_base`
- `mutation_rate_effective` = `nu0 * N_base * (1 + instability)`
- `exhaustion_rate_effective` = 0 (exhaustion_rate=0)
- `crowding_numerator` = `N_base + N_mutated`

#### MutatedClone ("mutated")
- `birth_rate_effective` = `(lambda0 * (1 + fitness_gain)) * N_mut * crowding`
- `death_rate_effective` = `mu0 * N_mut + N_mut * N_immune * theta_I` (immune killing!)
- `mutation_rate_effective` = 0 (no further mutation by default — `nu0` not inherited)
- `crowding_numerator` = `N_mutated`

#### ImmuneClone ("immune")
- `birth_rate_effective` = `lambda_Immune * N_immune * crowding + N_immune * N_mutated * beta` (activation boost!)
- `death_rate_effective` = 0 (base death_rate = 0.0)
- `exhaustion_rate_effective` = `mu_Immune * N_immune * N_mutated` (cancer-dependent exhaustion!)
- `crowding_numerator` = `N_immune + N_exhausted`

#### ExhaustedClone ("exhausted")
- `birth_rate_effective` = 0.0 (hardcoded)
- `death_rate_effective` = `mu_Exhausted * N_exhausted`
- `crowding_numerator` = 0 (exhausted cells don't compete for space)

### Crowding strategies

- **SimpleCrowding**: `C(t) = max(0, 1 - N_crowd / max(Kmin, K - decline * t))`
- **AdaptedCrowding**: `C(t) = max(0, 1 - N_crowd / max(Kmin, K / (1 - mu/lambda) - decline * t))`

Both only apply when `config.use_logistic=True` (which **must always be True** — "SIEMPRE TIENE QUE SER O EL CRECIMIENTO EXPONENCIAL EXPLOTA").

If `config.decay=False`, K does not decline with time (no `decline * t` term).

### SimulationConfig key parameters (defaults)

```python
@dataclass(frozen=True)
class SimulationConfig:
    N0: int = 100              # initial WT population
    N_mutant: int = 0          # initial mutant population
    N_immune: int = 50         # initial immune population
    N_exhausted: int = 0       # initial exhausted population

    lambda0: float = 0.005     # WT birth rate
    lambda_Immune: float = 0.005
    mu0: float = 0.002         # WT death rate
    mu_Immune: float = 0.003   # immune → exhausted rate
    mu_Exhausted: float = 0.002
    nu0: float = 0.0002        # WT mutation rate

    T_max: float = 2000
    seed: Optional[int] = None
    OMEGA: int = 100

    use_logistic: bool = True       # MUST be True
    use_logistic_adapted: bool = True
    K0: int = OMEGA                 # = 100
    K_immune: int = ceil(OMEGA/2)   # = 50
    K_mutant: int = OMEGA * 2       # = 200
    decline: float = 0.0
    Kmin: float = 1.0

    theta_I: float = 0.0005    # immune killing rate
    beta: float = 0.0004       # immune activation rate
    fitness_gain: float = 0.2  # mutant birth advantage

    # instability system (for mutational burden over time)
    instability_0: float = 0.0
    buildup_0: float = 0.0
    base_instability_buildup: float = 0.0
    mutation_instability_jump: float = 0.0
    mutation_buildup_gain: float = 0.0
    d1_0: float = 0.0
    d2_0: float = 0.0

    verbose: bool = True
    scale: bool = True
    decay: bool = True
    system_size: int = 10000
```

## Known bugs & quirks

### CLI ignores config args
`cli.py:108-109`: `config = SimulationConfig()` always creates defaults. The parsed `args.*` are only used for output paths and `--top`. Commands like `make gillespie-homeostasis` pass `--N0 50` etc. but they have **no effect**. To test with custom params, create a `TumorSimulation(config=SimulationConfig(...))` in Python directly.

### Broken tests
- `tests/gillespie/test_clone.py`: uses `self.tissue_state` in standalone functions (no class); `__str__` assertion expects `"soy un clone WT"` but actual returns `get_type()` which gives `"base"`
- `tests/gillespie/test_clone_factory.py`: **empty file** (0 lines)
- `tests/gillespie/test_crowding.py`: uses `self.tissue_state` in `TestAdaptedCrowding.test_higher_fitness_yields_higher_effective_birth_rate`

### Memory
History (`self.history`) stores a full snapshot dict at every `step()` — can grow unbounded for long simulations.

### Snapshot encoding
`TissueState.snapshot()` returns `Dict[CloneId, dict]` with keys `Type`, `N`, `rb`, `rd`, `rm`, `re`. The CSV writer (`csv_output.py`) only writes `N, rb, rd` — `rm` and `re` are silently dropped.

### Exhaustion event
`_induce_exhaustion` hardcodes `self.tissue_state.clones[(-2,)]` — it always divides the **original** exhausted clone `(-2,)`, not whichever exhausted clone exists. This works for the single-exhausted-clone setup but won't scale.

### Spanish TODOs
The codebase has Spanish comments and TODOs throughout (e.g., `#TODO: Arreglar este lio.`). Do not treat these as noise — they encode the author's intent.

## Test patterns

### Run commands
```bash
make test                              # all tests
make test-cov                          # with HTML coverage
python -m pytest tests/gillespie -v    # Gillespie tests only
python -m pytest tests/gillespie/test_rate_matrix.py -v   # single file
```

### Known working tests (skeleton)
- `tests/gillespie/test_event.py` — Event dataclass creation, equality, tuple type
- `tests/gillespie/test_rate_matrix.py` — RateMatrix add/clear/total/choose_event
- `tests/gillespie/test_io.py` — CSV output roundtrip (uses `tmp_path`)

### Test conventions
- No pytest fixtures (standalone factory functions like `make_clone()` and `make_clone_factory()` instead)
- No conftest.py
- `SimulationConfig` constructed inline in test helpers
- Clone instances created directly via `Clone(clone_id=..., N=..., config=config)` or via `CloneFactory.create_clone()`

## Gillespie CLI commands (makefile)

```bash
make gillespie-homeostasis      # --N0 50 --lambda0 0.25 --mu0 0.25 --nu0 0.00 --T-max 20 --seed 42
make gillespie-tumour-growth    # --N0 20 --lambda0 0.35 --mu0 0.20 --nu0 0.01 ... instability params
make gillespie-crowding         # same as tumour-growth + --use-logistic --K0 500
make gillespie-all              # runs all three
```

Direct invocation:
```bash
venv/bin/python -m src.gillespie.infrastructure.cli --N0 50 --lambda0 0.25 --mu0 0.25 --T-max 20 --seed 42 --save-history results.csv --save-clones clones.csv --top 15
```

## Dependencies
Install with `pip install -r requirements` (pinned, possibly stale).
Key packages: `numpy`, `pandas`, `PyYAML`, `openpyxl`, `pytest`, `pytest-cov`.

## Files you'll touch most

| File | Purpose |
|---|---|
| `src/gillespie/clone.py` | Clone base + 4 subclasses, rate formulas, registry |
| `src/gillespie/clone_factory.py` | Creates clones by type string |
| `src/gillespie/tissue_state.py` | State container, pop_map, snapshot |
| `src/gillespie/tumor_simulation.py` | Main loop, event dispatch, instability |
| `src/gillespie/simulation_config.py` | All model parameters |
| `src/gillespie/crowding_strategy.py` | logistic growth strategies |
| `src/gillespie/rate_matrix.py` | Event list, total rate, event selection |
| `src/gillespie/event.py` | Event dataclass |
| `src/gillespie/event_type.py` | EventType enum |
| `src/gillespie/infrastructure/cli.py` | CLI entrypoint (has config bug) |
| `src/gillespie/infrastructure/csv_output.py` | History/clone CSV writing |
