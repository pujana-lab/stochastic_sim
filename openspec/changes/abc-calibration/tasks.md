## 0. Synthetic reference data generator

- [x] 0.1 Create `scripts/generate_synthetic_reference.py` — CLI tool
- [x] 0.2 Define ground-truth params per scenario (homeostasis, tumour-growth, immune-response)
- [x] 0.3 Run N replicates per scenario, collect all trajectories
- [x] 0.4 Interpolate trajectories to common time grid
- [x] 0.5 Compute `mean_N` and `std_N` per (time, type) across replicates
- [x] 0.6 Compute precision weights `w_i = 1 / max(σ², 1e-4)` per statistic
- [x] 0.7 Save `reference_<scenario>.csv` (time, type, mean_N, std_N, n_replicates)
- [x] 0.8 Save `reference_<scenario>_weights.json` (precision weights)
- [x] 0.9 Generate diagnostic plot (time vs N, mean ± 1std bands)

## 1. Scaffold calibration package

- [x] 1.1 Create `src/gillespie/calibration/` package with subpackages (domain, ports, adapters, application)
- [x] 1.2 Add `scipy` to `requirements.txt` (stats module for ESS, norm kernel)
- [x] 1.3 Review `abc_final_v4.py` reference implementation for patterns to mirror

## 2. Domain layer — Core types

- [x] 2.1 Implement `Prior` — dataclass with param_name, lo, hi, distribution type, optional constraint function
- [x] 2.2 Implement `CalibrationConfig` — n_particles, n_generations, alpha, epsilon_0, epsilon_final, max_attempts, n_workers, n_reps, seed, output_dir
- [x] 2.3 Implement `SummaryStatistic` — compute_summary_stats(time_series, config) → np.ndarray
- [x] 2.4 Implement `DistanceMetric` — weighted_sse(sim_stats, ref_stats, precision_weights) → float
- [x] 2.5 Implement `AbcResult` — particles, weights, distances, diagnostics per generation

## 3. Ports — Interfaces

- [x] 3.1 Define `SimulatorPort` protocol — run(params_dict) → time_series dict
- [x] 3.2 Define `ReferenceDataPort` protocol — load() → (stats_vector, precision_weights)
- [x] 3.3 Define `OutputPort` protocol — save(gen, particles, weights, distances, diagnostics), save_manifest(config)

## 4. Adapters — Implementations

- [x] 4.1 Implement `GillespieSimulatorAdapter` — creates TumorSimulation, injects params via config.replace(), returns time series
- [x] 4.2 Implement `CsvReferenceAdapter` — loads reference CSV, validates columns, computes predefined summary stats
- [x] 4.3 Implement `CsvOutputAdapter` — writes gen_NN.csv, epsilon_schedule.csv, generation_diagnostics.csv, run_manifest.json

## 5. Application — ABC-SMC service

- [x] 5.1 Implement `sample_prior()` — draw from uniform priors with constraint validation
- [x] 5.2 Implement `sample_from_previous_gen()` — importance-weighted sample + Gaussian perturbation
- [x] 5.3 Implement `compute_weights()` — importance weight formula with numerical stability (log-space)
- [x] 5.4 Implement `_worker()` — run simulation, compute distance, return result
- [x] 5.5 Implement `run_generation()` — propose → simulate → accept/reject loop with ThreadPoolExecutor
- [x] 5.6 Implement `run_abc_smc()` — full multi-generation loop (gen 0 → gen 1+ → adaptive epsilon)
- [x] 5.7 Implement ESS diagnostic (1/Σw²)
- [x] 5.8 Implement resume support — find last gen, load, continue
- [x] 5.9 Implement prior sensitivity analysis — width factors × distance distribution

## 6. CLI integration

- [x] 6.1 Add `calibrate` subcommand to `infrastructure/cli_config.py`
- [x] 6.2 Implement CLI args (--config, --reference, --output, --dry-run, --resume, --prior-sensitivity)
- [x] 6.3 Wire to `AbcSmcService`

## 7. Tests

- [x] 7.1 Test prior sampling — uniform draws within bounds, constraint validation, prior_reject count
- [x] 7.2 Test summary statistic computation against known time series
- [x] 7.3 Test weighted SSE distance with uniform and precision weights
- [x] 7.4 Test compute_weights — numerical stability, sum to 1
- [x] 7.5 Test run_generation — acceptance loop, worker pool, progress logging
- [x] 7.6 Test CSV adapter roundtrip — write gen_NN.csv → read back
- [x] 7.7 Test resume — write gen_00.csv/gen_01.csv, resume → verify continues at gen 02
- [x] 7.8 Test calibration with Gillespie adapter — synthetic reference, small n_particles, 2 gens
- [x] 7.9 Test prior sensitivity output format

## 8. Smoke test with real scenario

- [x] 8.1 Create calibration JSON for homeostasis scenario (calibrate lambda0, mu0 against known data)
- [x] 8.2 Run 10-particle, 2-generation calibration (reduced for CI speed; 50-particle 3-gen tested manually)
- [x] 8.3 Verify posterior contains expected parameter region
- [x] 8.4 Verify all output files exist and are valid
