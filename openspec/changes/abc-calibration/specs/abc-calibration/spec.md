## ADDED Requirements

### Requirement: System generates synthetic reference data from known params

The system SHALL provide a tool `scripts/generate_synthetic_reference.py` that runs N replicates of the Gillespie simulator with known ground-truth parameters, aggregates trajectories, and produces a reference CSV with `mean_N` and `std_N` per time point per cell type, plus a precision weights JSON.

#### Scenario: Generate homeostasis reference
- **WHEN** user runs with `--scenario homeostasis --n-replicates 50 --seed 42`
- **THEN** output CSV SHALL have columns `time, type, mean_N, std_N, n_replicates`
- **THEN** output weights JSON SHALL keys match summary statistic names

#### Scenario: Ground-truth params are recorded
- **WHEN** synthetic data is generated
- **THEN** a `_ground_truth.json` SHALL be saved with the exact parameters used

#### Scenario: Diagnostic plot
- **WHEN** generating synthetic data
- **THEN** a PNG plot SHALL be saved showing mean ± 1std trajectories per cell type

### Requirement: User can define priors for calibrated parameters

The system SHALL accept a dict of `{param_name: (lo, hi)}` defining uniform priors. Each prior SHALL have hard bounds — proposals outside bounds are rejected without simulation. Additional hard constraints (e.g., `high_delta > low_delta`) SHALL be supported. Distribution types may include `uniform`, `lognorm`, `truncnorm`.

#### Scenario: Define uniform prior
- **WHEN** user configures prior `{"lambda0": (0.001, 0.1), "mu0": (0.001, 0.1)}`
- **THEN** proposals SHALL be drawn uniformly from [lo, hi)

#### Scenario: Hard constraint validation
- **WHEN** a proposal violates a hard constraint or prior bound
- **THEN** it SHALL be rejected without simulation (counted as `prior_reject`)

#### Scenario: Invalid distribution raises error
- **WHEN** user specifies `distribution="invalid_dist"`
- **THEN** ValueError SHALL be raised

### Requirement: User can provide reference data as CSV

The system SHALL load reference data from CSV. Supported formats: (a) time series with `time` + population columns, (b) milestone ages with `age` + `cumulative_risk` columns. Optional weights file for per-statistic precision (1/σ²).

#### Scenario: Load time series reference
- **WHEN** CSV has columns `time,base,mutated,immune`
- **THEN** parse into dict of arrays keyed by cell type

#### Scenario: Missing column raises error
- **WHEN** CSV is missing required columns
- **THEN** ValueError describing the missing column

### Requirement: System computes summary statistics

The system SHALL compute summary statistics from simulation output to compare against reference. Default stats: `mean`, `std`, `final`, `integral` per cell type. Alternatively, milestone-style statistics (e.g., time to reach N cells) if reference is milestone data.

#### Scenario: Default time-series summary
- **WHEN** reference is time series with types [base, mutated]
- **THEN** default stats SHALL be `[mean_base, std_base, final_base, integral_base, mean_mut, std_mut, final_mut, integral_mut]`

#### Scenario: Custom summary statistics
- **WHEN** user specifies custom summary stat config
- **THEN** only those stats SHALL be computed

### Requirement: System computes weighted SSE distance

The system SHALL compute `SSE_w = Σ w_i × (sim_i - ref_i)²` where `w_i = 1/σ²_i` (precision weights from reference uncertainty). If no weights provided, default to uniform.

#### Scenario: Weighted SSE with precision weights
- **WHEN** reference data provides precision weights (1/σ²)
- **THEN** distance SHALL be weighted SSE

#### Scenario: Uniform weights fallback
- **WHEN** no weights provided
- **THEN** all w_i = 1.0 (standard SSE)

### Requirement: System runs ABC-SMC (Sequential Monte Carlo)

The system SHALL implement the Toni et al. 2009 ABC-SMC algorithm:

- **Gen 0**: Sample from prior → simulate → accept if `distance ≤ ε₀`
- **Gen 1+**: Importance-sample from previous population (weighted) → perturb with Gaussian kernel (σ = 2× weighted std) → simulate → accept if `distance ≤ ε_t`
- **Epsilon schedule**: ε_{t+1} = quantile(distances_t, α) (adaptive)
- **Weights**: Gen 0: uniform (1/N). Gen 1+: `w_t = prior(θ) / Σ[w_{t-1} × K_t(θ|θ_{t-1})]`

#### Scenario: Generation 0 prior sampling
- **WHEN** t=0 and n_particles=200
- **THEN** sample 200 particles from prior, run simulation, accept if distance ≤ ε₀

#### Scenario: Generation 1+ importance sampling
- **WHEN** t>0
- **THEN** sample parent from previous population proportional to weights, perturb with Gaussian, run simulation, accept if distance ≤ ε_t

#### Scenario: Weight computation (importance weights)
- **WHEN** t>0 and particle θ accepted
- **THEN** weight = prior(θ) / Σ[w_{t-1} × kernel(θ|θ_{t-1})]

#### Scenario: Adaptive epsilon
- **WHEN** generation t completes
- **THEN** ε_{t+1} = quantile(distances_t, α)

#### Scenario: Generation failure
- **WHEN** n_particles cannot be accepted after max_attempts
- **THEN** RuntimeError with explicit failure message (no padding)

### Requirement: System supports parallel evaluation

The system SHALL use ThreadPoolExecutor for parallel simulation evaluation. Number of workers configurable. Per-generation progress SHALL be logged.

#### Scenario: Parallel simulation
- **WHEN** workers=4 and 200 particles needed
- **THEN** up to 4 simulations SHALL run concurrently

#### Scenario: Progress logging
- **WHEN** generation is running
- **THEN** SHALL log: gen, attempts, prior_reject, sim_fail, accepted/N, acceptance_rate, eps, elapsed

### Requirement: System tracks convergence diagnostics per generation

The system SHALL compute and persist per generation: Effective Sample Size (ESS), ESS/N ratio, acceptance rate, prior reject count, sim failure count, mean/min distance.

#### Scenario: ESS computation
- **WHEN** generation t completes
- **THEN** ESS = 1/Σ(w_i²) SHALL be computed and logged

#### Scenario: Diagnostics persisted
- **WHEN** generation completes
- **THEN** diagnostics SHALL be appended to `generation_diagnostics.csv`

### Requirement: System outputs per-generation results

The system SHALL save one CSV per generation with columns: param_names + weight + distance + optional summary stats.

#### Scenario: Save generation CSV
- **WHEN** generation 0 completes
- **THEN** `gen_00.csv` SHALL be created with all particles, weights, distances

### Requirement: System supports resume from previous run

The system SHALL detect saved generations, load the last one, compute epsilon from distance quantile, compute perturbation sigma from weighted std, and continue from the next generation.

#### Scenario: Resume from gen 3
- **WHEN** `gen_00.csv` to `gen_03.csv` exist and resume=True
- **THEN** SHALL load gen 03, compute epsilon from alpha quantile, continue at gen 04

### Requirement: System generates reproducibility manifest

The system SHALL save a `run_manifest.json` with: timestamp, git commit, python version, platform, binary/script path, priors, fixed parameters, full CLI config.

#### Scenario: Manifest creation
- **WHEN** calibration starts
- **THEN** `run_manifest.json` SHALL be written to output directory

### Requirement: User can run calibration from CLI

The system SHALL provide a `calibrate` CLI subcommand.

#### Scenario: Run calibration
- **WHEN** user runs `python -m src.gillespie.infrastructure.cli calibrate --config calib.json --reference data.csv --output ./calib/`
- **THEN** ABC-SMC SHALL run and results saved to `./calib/`

#### Scenario: Dry run
- **WHEN** user adds `--dry-run`
- **THEN** SHALL print config and exit without simulations

### Requirement: Prior sensitivity analysis

The system SHALL support prior width sensitivity: scale prior bounds by factors (e.g., 0.5, 0.8, 1.2, 1.5), sample N particles per factor, compute distance distribution, report acceptance at reference epsilon.

#### Scenario: Prior sensitivity
- **WHEN** --prior-sensitivity flag set
- **THEN** SHALL sample from scaled priors, compute acceptance rate at ε_ref, save `prior_sensitivity.csv`
