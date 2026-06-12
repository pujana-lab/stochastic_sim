## Context

Gillespie simulator has ~20 parameters. No automated fitting exists. You already have a proven ABC-SMC pipeline in `cellSim/scripts/abc_final_v4.py` for a C++ cancer simulator (BRCA1 calibration against Kuchenbaecker 2017 clinical data). That script uses:

- **Sequential Importance Sampling** (Toni et al. 2009) — not simple rejection
- **Weighted SSE** with 1/σ² precision weights from clinical CI
- **Bt ≥ 1 stochastic replicates** per particle
- **ThreadPoolExecutor** with OMP_NUM_THREADS=1 to prevent oversubscription
- **ESS diagnostics** + acceptance rate tracking per generation
- **Prior sensitivity analysis** (width sweep)
- **Reproducibility manifest** (git commit, env, seed, config)
- **Resume** from last saved generation

This design ports that methodology to our Python Gillespie simulator, replacing the C++ binary with `TumorSimulation`.

## Goals / Non-Goals

**Goals:**
- ABC-SMC calibration for Gillespie parameters (mirroring `abc_final_v4.py`)
- Reference data loader (CSV time series or milestone ages)
- Configurable priors, summary statistics, precision weights
- Parallel evaluation via `ThreadPoolExecutor`
- Per-generation diagnostics (ESS, acceptance rate, distance quantiles)
- Resume from last generation
- Prior sensitivity analysis
- Reproducibility manifest

**Non-Goals:**
- No MCMC-ABC (keep SMC with importance weights)
- No GPU acceleration
- No automatic summary statistic selection
- No calibration of the old Moran codebase

## Decisions

### Architecture: Hexagonal (Ports & Adapters)

Follows existing project pattern. Mirrors the `abc_final_v4.py` structure:

```
Domain:
  CalibrationConfig    — n_particles, n_generations, alpha, epsilon_0, epsilon_final
  Prior                — distribution per parameter (uniform, lognorm, truncnorm)
  SummaryStatistic     — computes summary vector from sim/reference time series
  DistanceMetric       — weighted SSE (Σ w_i * (sim_i - ref_i)², w_i = 1/σ²_i)
  AbcResult           — accepted params, weights, distances per generation

Ports:
  SimulatorPort        — run(params_dict) → time_series_dict
  ReferenceDataPort    — load() → reference_time_series + precision_weights
  OutputPort           — save(generation, particles, weights, distances, diagnostics)

Adapters:
  GillespieSimulatorAdapter  — wraps TumorSimulation, injects params into SimulationConfig
  CsvReferenceAdapter        — loads reference CSV + optional weights file
  CsvOutputAdapter           — writes gen_NN.csv, epsilon_schedule.csv, diagnostics.csv, manifest.json

Application:
  AbcSmcService        — orchestrates: gen 0 prior sample → gen 1+ importance sample + perturb → weight → accept vs epsilon → diagnostics
  CalibrationCli       — CLI entrypoint
```

### Decision: ABC-SMC over ABC rejection

`abc_final_v4.py` already uses SMC (Toni et al. 2009). Simple rejection wastes simulations. SMC reuses accepted particles to guide sampling toward high-probability regions.

- Gen 0: uniform prior → simulate → accept if distance ≤ ε₀
- Gen 1+: importance sample from previous population → perturb with Gaussian kernel (σ = 2× weighted std of previous population) → simulate → accept → compute importance weights
- Epsilon schedule: ε_{t+1} = quantile(ε_t, alpha=0.6)

### Decision: Distance = Weighted SSE (not raw Euclidean)

Mirrors `abc_final_v4.py`: `SSE_w = Σ (sim_i - obs_i)² × w_i` where `w_i = 1/σ²_i` (precision weights from reference uncertainty).

**Why:** Clinical/experimental data has different certainty per statistic. Weighted SSE naturally handles this.

### Decision: Per-generation CSV output

Mirrors `abc_final_v4.py`:
- `gen_00.csv`, `gen_01.csv`, ... — one file per generation, columns = param_names + weight + distance + optional summary stats
- `epsilon_schedule.csv` — generation → epsilon
- `generation_diagnostics.csv` — ESS, acceptance rate, prior_reject, sim_fail per gen
- `run_manifest.json` — reproducibility (config, priors, seed, git, platform)

### Decision: Resume support

Load last generation from CSV, continue from next. Epsilon from distance quantile of loaded generation. Perturbation sigma from weighted std of loaded particles.

### Decision: Prior as dict of (lo, hi) tuples + optional distribution type

Same as `abc_final_v4.py` PRIORS dict. Uniform by default. Extensible to lognorm/truncnorm.

### Decision: Bt stochastic replicates (n_reps)

Each particle evaluated n_reps times, distance averaged. Mirrors `abc_final_v4.py --n-reps`. Controls simulation stochasticity.

### Decision: Synthetic reference data from ground-truth params

Since no real experimental data exists, generate synthetic reference from known params:
1. Pick ground-truth params per scenario (mirror Makefile targets)
2. Run N=100 replicates with different seeds
3. Interpolate all trajectories to common time grid (every 5 time units)
4. Compute `mean_N` and `std_N` per (time, type)
5. Precision weight = `1 / max(σ², 1e-4)` for each summary statistic

This gives a known target + natural uncertainty. We validate ABC-SMC by checking posteriors contain the ground truth.

### Decision: Manual prior constraint validation

Like `abc_final_v4.py`: hard prior bounds enforced before simulation. No boundary clipping (reject proposals outside prior — the paper-consistent approach). Additional hard constraints: e.g., `high_delta > low_delta`.

## Package structure

```
src/gillespie/calibration/
├── __init__.py
├── domain/
│   ├── __init__.py
│   ├── calibration_config.py    — CalibrationConfig, Prior, PriorDict
│   ├── summary_statistic.py     — compute_summary_stats()
│   ├── distance_metric.py       — weighted_sse()
│   └── abc_result.py            — AbcResult (particles, weights, dists, diagnostics)
├── ports/
│   ├── __init__.py
│   └── calibration_ports.py     — SimulatorPort, ReferenceDataPort, OutputPort
├── adapters/
│   ├── __init__.py
│   ├── gillespie_adapter.py     — wraps TumorSimulation run()
│   ├── csv_reference.py         — load reference CSV + optional weights
│   └── csv_output.py            — save per-generation CSV + manifest
└── application/
    ├── __init__.py
    ├── abc_smc_service.py       — main SMC loop (run_generation, compute_weights)
    └── cli.py                   — calibration CLI entrypoint
```

## Isolation boundaries (what calibration touches vs not)

Goal: zero collateral damage. Compañero no tiene que saber que existe.

| Componente | Lo toca | Lo NO toca |
|---|---|---|
| `clone.py` | — | ✅ intacto |
| `tumor_simulation.py` | — | ✅ intacto |
| `rate_matrix.py` | — | ✅ intacto |
| `event.py`, `event_type.py` | — | ✅ intacto |
| `tissue_state.py` | — | ✅ intacto |
| `crowding_strategy.py` | — | ✅ intacto |
| `clone_factory.py` | — | ✅ intacto |
| `simulation_config.py` | ✅ via `dataclasses.replace()` (read-only) | ✅ no muta defaults |
| `infrastructure/cli.py` | ✅ añade subcomando `calibrate` | ✅ args existentes siguen igual |
| `infrastructure/csv_output.py` | — | ✅ intacto |
| Tests existentes | — | ✅ intactos |
| Tests nuevos | `tests/gillespie/calibration/` | ✅ sin tocar otros test files |

**Dependencias:** calibrarion **importa** `src/gillespie/` (SimulationConfig, TumorSimulation). `src/gillespie/` **nunca** importa calibration. Grafo acíclico.

**Tests:** suite separada `python -m pytest tests/gillespie/calibration/ -v`. CI no la ejecuta por defecto hasta que esté estable.

**Review:** todo el diff son ficheros nuevos en `calibration/` y 2-3 líneas en `cli.py` para registrar el subcomando. Nada más.

## Risks / Trade-offs

| Risk | Mitigation |
|------|------------|
| High sim time per run (10k+ samples × Gillespie = hours) | Start small; log progress per N attempts; add `--dry-run`; control workers |
| Weight degeneracy (ESS collapse) | Track ESS per generation; warn if ESS/N < 0.1 |
| Zero acceptance rate at high epsilon | `abc_final_v4.py` error handling: fail explicitly, suggest fixes |
| Summary stats mismatch between sim and ref | Validate on load; explicit column mapping |
| Reproducibility with parallel workers | Fixed seed per generation; seed_base from RNG per particle |
| C++ binary vs Python Gillespie | Different target data (cell counts vs cancer incidence). Adapt adapter layer |
