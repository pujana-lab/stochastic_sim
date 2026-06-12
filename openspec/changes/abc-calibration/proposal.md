## Why

Simulation parameters hand-tuned. No automated fitting to experimental data. You already have a working ABC-SMC pipeline in `cellSim/scripts/abc_final_v4.py` (Toni et al. SIS weights, weighted SSE distance, parallel workers, ESS diagnostics, prior sensitivity). Need the same for our Python Gillespie simulator.

## What Changes

- New `abc-calibration` module under `src/gillespie/calibration/`
- **ABC-SMC** (Sequential Monte Carlo, not simple rejection):
  - Gen 0: uniform prior sampling
  - Gen 1+ : importance-weighted proposals + Gaussian perturbation kernel
  - Adaptive epsilon schedule (quantile-based decay)
  - Weighted SSE distance with per-statistic precision weights (1/σ²)
- Parallel evaluation via `ThreadPoolExecutor` (+ worker count control)
- Output: per-generation CSV (particles + weights + distances), epsilon schedule, generation diagnostics (ESS, acceptance rate), posterior summary
- Resume support (load last generation and continue)
- Prior sensitivity analysis
- Reproducibility manifest (git commit, env, config)

## Capabilities

### New Capabilities
- `abc-calibration`: ABC-SMC calibration framework for Gillespie parameters. Mirrors `abc_final_v4.py` methodology: sequential importance sampling, weighted SSE distance, ESS diagnostics, prior sensitivity, resume, reproducibility manifest.

### Modified Capabilities
<!-- No existing specs to modify -->

## Impact

- `src/gillespie/calibration/` — new hexagonal package (domain, ports, adapters, application)
- `SimulationConfig` — may need partial config for calibration targets (which params to vary, which to fix)
- Reference data loading — CSV of clinical/experimental time series (cell counts over time or milestone ages)
- CLI — new `calibrate` subcommand
- Dependencies: `scipy` (stats for ESS, kernels), `matplotlib` optional for posterior plots
- Synthetic reference data generator (`scripts/generate_synthetic_reference.py`) — runs N replicates from known params, aggregates to mean±std, used as calibration target
- No changes to core simulation loop, clone classes, or rate logic
- Existing `abc_final_v4.py` (cellSim C++ calibration) is reference implementation — our Gillespie version mirrors its methodology
