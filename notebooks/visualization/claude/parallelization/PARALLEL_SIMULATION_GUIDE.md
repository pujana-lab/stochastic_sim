# Parallel Simulation Guide

Master multi-core parallelization of your Gillespie tumor simulations.

---

## Quick Start

```python
from parallel_simulator import ParallelSimulator, MemoryMode
from src.gillespie.simulation_config import SimulationConfig
import numpy as np

# 1. Create base configuration
base_config = SimulationConfig(T_max=1000, seed=42)

# 2. Create parallel simulator (auto-detects available cores)
parallel = ParallelSimulator(n_cores=8, memory_mode=MemoryMode.LIGHTWEIGHT)

# 3. Sweep parameter
results = parallel.sweep_parameter(
    param_name='mutation_rate',
    param_values=np.linspace(0.01, 0.1, 50),
    base_config=base_config,
    seeds=[42, 123, 456]  # 3 replicates per param value
)

# 4. Print summary
parallel.print_results_summary(results)

# 5. Analyze
summary = ParallelSimulator.summarize_sweep(results, 'mutation_rate')
for param_val, stats in summary.items():
    print(f"μ={param_val:.3f}: extinction_rate={stats['extinction_rate']:.1%}")
```

**Expected output:**
```
μ=0.010: extinction_rate=95.2%
μ=0.020: extinction_rate=87.3%
μ=0.030: extinction_rate=72.1%
...
```

---

## Understanding Parallelization

### Why Parallelize?

```
Sequential:  Config1  Config2  Config3  ...  Config50
Timeline:    [====][====][====]              [====]  (total: 50× simulation time)

Parallel:    Config1     Config14     Config27    Config40
             Config2     Config15     Config28    Config41
             Config3     Config16     Config29    Config42
             ...
Timeline:    [============================]           (total: ~50/8 × simulation time)
```

**Speedup with 8 cores: ~7-8x faster** (not perfect due to overhead, but close)

### When to Parallelize

| Task | Parallelizable? | Speedup |
|------|---|---|
| Parameter sweep (1D) | ✅ Yes | ~n_cores |
| Parameter sweep (2D) | ✅ Yes | ~n_cores |
| Ensemble runs (different seeds) | ✅ Yes | ~n_cores |
| Single long simulation | ❌ No | 1x (internal parallelization needed) |
| Analyzing results | ⚠️ Sometimes | Depends on task |

---

## Usage Patterns

### Pattern 1: Simple Parameter Sweep

```python
from parallel_simulator import ParallelSimulator, MemoryMode
import numpy as np

# Create simulator with all available cores
parallel = ParallelSimulator(memory_mode=MemoryMode.LIGHTWEIGHT)

# Sweep mutation rate from 0.01 to 0.1
results = parallel.sweep_parameter(
    param_name='mutation_rate',
    param_values=np.linspace(0.01, 0.1, 50),
    base_config=base_config
)

# Analyze
summary = ParallelSimulator.summarize_sweep(results, 'mutation_rate')
```

**Memory:** ~50-100 MB total (1-2 MB per sim × 50)
**Time:** ~2 hours sequential → ~15 minutes parallel (on 8 cores)

---

### Pattern 2: Ensemble with Multiple Replicates

```python
# Run each parameter value multiple times with different seeds
results = parallel.sweep_parameter(
    param_name='mutation_rate',
    param_values=np.linspace(0.01, 0.1, 20),
    base_config=base_config,
    seeds=list(range(42, 52))  # 10 replicates each = 200 total sims
)

# Aggregate to get mean ± std
summary = ParallelSimulator.summarize_sweep(results, 'mutation_rate')
for param_val, stats in summary.items():
    extinction_mean = stats['extinction_rate']
    extinction_std = np.std([r.metadata['extinction'] for r in results 
                             if r.metadata['param_value'] == param_val])
    print(f"μ={param_val:.3f}: {extinction_mean:.1%} ± {extinction_std:.1%}")
```

**Recommended for:** Publication-quality results with error bars

---

### Pattern 3: 2D Parameter Sweep

```python
# Sweep two parameters simultaneously
results, grid_info = parallel.sweep_parameters_2d(
    param1_name='mutation_rate',
    param1_values=np.linspace(0.01, 0.1, 15),
    param2_name='exhaustion_rate',
    param2_values=np.linspace(0.001, 0.01, 15),
    base_config=base_config
)

# Results organized in 15×15 grid (225 simulations)
print(f"Grid shape: {grid_info['shape']}")  # (15, 15)
print(f"Total configurations: {len(results)}")  # 225

# Organize into 2D array for plotting
extinction_grid = np.zeros(grid_info['shape'])
for result in results:
    if result.error is None:
        i = result.metadata['param1_idx']
        j = result.metadata['param2_idx']
        extinction_grid[i, j] = result.metadata['extinction']

# Plot heatmap
import matplotlib.pyplot as plt
plt.imshow(extinction_grid, cmap='RdYlGn_r', aspect='auto')
plt.colorbar(label='Extinction Rate')
plt.xlabel('Exhaustion Rate')
plt.ylabel('Mutation Rate')
plt.show()
```

---

### Pattern 4: Ensemble Runs (Same Config, Different Seeds)

```python
# Run same config with 100 different random seeds
n_replicates = 100
configs = [
    SimulationConfig(**vars(base_config), seed=i) 
    for i in range(42, 42 + n_replicates)
]

results = parallel.run_ensemble(configs)
agg = ParallelSimulator.aggregate_results(results)

print(f"Extinction rate: {agg['extinction_rate']:.1%} ± {np.std(agg['final_populations']):.1f}")
print(f"Mean final population: {agg['final_population_mean']:.0f} ± {agg['final_population_std']:.0f}")
```

**Use for:** Characterizing stochastic behavior

---

### Pattern 5: Adaptive Parallel Strategy

```python
from parallel_simulator import ParallelStrategy

# For many short simulations: use async (better load balancing)
if len(configs) > 100:
    strategy = ParallelStrategy.ASYNC_MAP
elif len(configs) > 10:
    strategy = ParallelStrategy.PROCESS_POOL  # Default, usually best
else:
    strategy = ParallelStrategy.CONCURRENT_FUTURES

parallel = ParallelSimulator(
    n_cores=8,
    strategy=strategy,
    memory_mode=MemoryMode.LIGHTWEIGHT
)

results = parallel.run_ensemble(configs)
```

---

## Parallelization Strategies

### Strategy 1: `PROCESS_POOL` (Recommended)

```python
from parallel_simulator import ParallelStrategy, ParallelSimulator

parallel = ParallelSimulator(strategy=ParallelStrategy.PROCESS_POOL)
```

**When to use:** Most cases (default)
**Pros:**
- Simple, predictable
- Good load balancing
- Works well with tqdm progress bar
- Lowest overhead

**Cons:**
- Results returned in order (may wait for slowest job)

---

### Strategy 2: `CONCURRENT_FUTURES`

```python
parallel = ParallelSimulator(strategy=ParallelStrategy.CONCURRENT_FUTURES)
```

**When to use:** When results complete at very different times
**Pros:**
- Returns results as they complete (no waiting for slowest)
- Excellent for heterogeneous job times
- Built-in timeout support

**Cons:**
- Slightly higher overhead
- Requires Python 3.7+

---

### Strategy 3: `ASYNC_MAP` (Streaming)

```python
parallel = ParallelSimulator(strategy=ParallelStrategy.ASYNC_MAP)
```

**When to use:** Very large ensembles (1000+), when you can't fit all results in memory
**Pros:**
- Minimal memory (processes results immediately)
- Best for streaming to disk
- Excellent cache locality

**Cons:**
- Results returned out of order (must resort)
- Harder to debug

---

## Memory Considerations

### Memory + Cores = Scaling Power

```python
# Configuration 1: Conservative
parallel = ParallelSimulator(
    n_cores=4,
    memory_mode=MemoryMode.LIGHTWEIGHT
)
# Memory: ~4-8 MB total (1-2 MB per simulation × 4 cores)
# Good for: 16 GB RAM, shared systems

# Configuration 2: Balanced
parallel = ParallelSimulator(
    n_cores=8,
    memory_mode=MemoryMode.STANDARD
)
# Memory: ~100-400 MB total (10-50 MB per simulation × 8 cores)
# Good for: 32 GB RAM, dedicated machine

# Configuration 3: Aggressive
parallel = ParallelSimulator(
    n_cores=16,
    memory_mode=MemoryMode.FULL
)
# Memory: ~2-16 GB total (100 MB - 1 GB per simulation × 16 cores)
# Good for: 128 GB RAM, high-memory system
```

### Memory Profile

```python
from parallel_simulator import MemoryMode

# These are per-core memory footprints:
# LIGHTWEIGHT: 1-2 MB/core → 16 cores = 16-32 MB total
# STANDARD:    10-50 MB/core → 8 cores = 80-400 MB total  
# FULL:        100 MB-1 GB/core → 4 cores = 400 MB-4 GB total
```

**Rule of thumb:** Max parallel jobs = RAM / (memory_per_sim × 1.5)

---

## Monitoring & Debugging

### Progress Tracking

```python
# Automatic progress bar (default)
results = parallel.run_ensemble(configs, show_progress=True)

# Disable progress bar
results = parallel.run_ensemble(configs, show_progress=False)

# Custom monitoring
for i, config in enumerate(configs):
    print(f"Processing config {i+1}/{len(configs)}: {config.param_name}={config.param_value}")
```

### Error Handling

```python
results = parallel.run_ensemble(configs)

# Check for errors
failed_results = [r for r in results if r.error is not None]
if failed_results:
    print(f"Failed runs: {len(failed_results)}")
    for result in failed_results:
        print(f"  Config {result.config_id}: {result.error}")

# Filter valid results
valid_results = [r for r in results if r.error is None]
print(f"Successful: {len(valid_results)}/{len(results)}")
```

### Debugging a Failed Config

```python
# Identify which config failed
failed_result = [r for r in results if r.error is not None][0]
config_id = failed_result.config_id
config = configs[config_id]

# Re-run sequentially to debug
from tumor_simulation_optimized import TumorSimulation, MemoryMode

sim = TumorSimulation(config, memory_mode=MemoryMode.FULL)
times, history, final_state, rates = sim.run()  # Will show full traceback
```

---

## Result Analysis

### Quick Summary

```python
parallel.print_results_summary(results)
```

Output:
```
============================================================
Ensemble Results Summary
============================================================
Total runs:         50
Successful:         50
Failed:             0

Final Population:
  Mean:             1234.5
  Std:              456.7
  Range:            [100, 5678]

Extinction Rate:    42.0%
============================================================
```

### Extract Trajectories for Plotting

```python
times_list, trajectories = ParallelSimulator.extract_trajectories(results)

import matplotlib.pyplot as plt

plt.figure(figsize=(12, 6))
for times, trajectory in zip(times_list[:20]):  # Plot first 20
    plt.plot(times, trajectory, alpha=0.3)

plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Population Trajectories')
plt.show()
```

### Organize by Parameter Value

```python
summary = ParallelSimulator.summarize_sweep(results, 'mutation_rate')

# Extract data for plotting
param_values = sorted(summary.keys())
extinction_rates = [summary[p]['extinction_rate'] for p in param_values]
means = [summary[p]['final_population_mean'] for p in param_values]
stds = [summary[p]['final_population_std'] for p in param_values]

plt.errorbar(param_values, means, yerr=stds, fmt='o-', capsize=5)
plt.xlabel('Mutation Rate')
plt.ylabel('Final Population')
plt.show()
```

---

## Performance Benchmarking

### Measure Speedup

```python
import time

base_config = SimulationConfig(T_max=1000)
configs = [SimulationConfig(**vars(base_config), seed=i) for i in range(100)]

# Sequential
print("Sequential mode:")
t0 = time.time()
from tumor_simulation_optimized import TumorSimulation, MemoryMode
for config in configs:
    sim = TumorSimulation(config, memory_mode=MemoryMode.LIGHTWEIGHT)
    sim.run()
t_sequential = time.time() - t0
print(f"  Time: {t_sequential:.1f} seconds")

# Parallel
print("Parallel mode (8 cores):")
parallel = ParallelSimulator(n_cores=8, memory_mode=MemoryMode.LIGHTWEIGHT)
t0 = time.time()
results = parallel.run_ensemble(configs)
t_parallel = time.time() - t0
print(f"  Time: {t_parallel:.1f} seconds")

# Speedup
speedup = t_sequential / t_parallel
efficiency = speedup / 8 * 100
print(f"\nSpeedup: {speedup:.1f}x")
print(f"Efficiency: {efficiency:.1f}%")  # Ideal: 100%, realistic: 70-90%
```

### Memory Profiling

```python
import psutil
import os

process = psutil.Process(os.getpid())

# Baseline
mem_before = process.memory_info().rss / 1e6
print(f"Memory before: {mem_before:.1f} MB")

# Run parallel
parallel = ParallelSimulator(n_cores=8, memory_mode=MemoryMode.STANDARD)
results = parallel.run_ensemble(configs)

# Peak
mem_peak = process.memory_info().rss / 1e6
print(f"Memory peak: {mem_peak:.1f} MB")
print(f"Increase: {mem_peak - mem_before:.1f} MB")

# Estimate per-core
per_core = (mem_peak - mem_before) / 8
print(f"Per-core: {per_core:.1f} MB")
```

---

## Troubleshooting

### Issue: "RuntimeError: An attempt has been made to start a new process..."

**On macOS/Windows, add:**
```python
if __name__ == '__main__':
    # Put your parallel code inside here
    parallel = ParallelSimulator()
    results = parallel.run_ensemble(configs)
```

This is required on non-Linux systems.

---

### Issue: Simulations Are Slow/Not Using All Cores

**Check:**
```python
import multiprocessing
print(f"Available cores: {multiprocessing.cpu_count()}")

# Verify parallelism
from parallel_simulator import ParallelSimulator
parallel = ParallelSimulator()
print(f"Using cores: {parallel.n_cores}")
```

**Common causes:**
1. Python's GIL blocking? No, simulations are CPU-bound and don't hold GIL
2. I/O bottleneck? Unlikely for simulations
3. Job too short? If individual simulations < 1 second, overhead dominates
4. Bad load balancing? Try different strategy

---

### Issue: "PicklingError: Can't pickle X"

**Cause:** Some object in your config isn't serializable

**Fix:**
```python
# Avoid lambdas, nested functions, unpickleable objects
# DON'T do this:
config.custom_fn = lambda x: x**2  # Can't pickle lambdas

# DO this:
def my_custom_fn(x):
    return x**2
config.custom_fn = my_custom_fn  # Pickleable function
```

---

### Issue: Out of Memory

**Solutions (in order):**
```python
# 1. Use LIGHTWEIGHT mode
parallel = ParallelSimulator(memory_mode=MemoryMode.LIGHTWEIGHT)

# 2. Reduce n_cores
parallel = ParallelSimulator(n_cores=4)  # Not 8

# 3. Process results incrementally
results = parallel.run_ensemble(configs, show_progress=True)
for result in results:
    # Process immediately, don't keep in memory
    save_to_disk(result)
```

---

## Advanced Patterns

### Pattern: Save Results Incrementally

```python
import pickle

def save_result(result, output_dir):
    """Save individual result to disk"""
    with open(f"{output_dir}/result_{result.config_id}.pkl", 'wb') as f:
        pickle.dump(result, f)

# Run and save
results = parallel.run_ensemble(configs, show_progress=True)
for result in results:
    save_result(result, '/tmp/results')

print(f"Saved {len(results)} results to /tmp/results")
```

---

### Pattern: Real-Time Monitoring

```python
import time

class ProgressMonitor:
    def __init__(self):
        self.start_time = time.time()
        self.results_so_far = []
    
    def process_result(self, result):
        """Called as each result completes"""
        self.results_so_far.append(result)
        elapsed = time.time() - self.start_time
        rate = len(self.results_so_far) / elapsed if elapsed > 0 else 0
        print(f"Completed {len(self.results_so_far)}: "
              f"{rate:.1f} sims/sec, "
              f"ETA {(len(configs)-len(self.results_so_far))/rate:.0f}s")

# Usage with async strategy
results = parallel.run_ensemble(configs, show_progress=True)
```

---

## Integration with Other Tools

### Jupyter Notebook

```python
# In Jupyter, multiprocessing requires special setup
import multiprocessing as mp

if __name__ == '__main__' or 'ipykernel' in mp.get_context().get_start_method():
    mp.set_start_method('fork', force=True)

from parallel_simulator import ParallelSimulator
parallel = ParallelSimulator()
results = parallel.run_ensemble(configs)
```

### Snakemake Workflow

```snakemake
# Integrate with Snakemake for reproducible pipelines
rule sweep_mutation_rate:
    input:
        config = "config.yaml"
    output:
        results = "results/sweep.pkl"
    shell:
        """
        python scripts/run_sweep.py --output {output.results}
        """
```

---

## Summary

| Scenario | Cores | Memory Mode | Expected Speedup |
|----------|-------|-------------|------------------|
| Quick test | 2 | LIGHTWEIGHT | 1.8x |
| Parameter sweep | 8 | LIGHTWEIGHT | 7.5x |
| Detailed analysis | 4 | STANDARD | 3.8x |
| Large study | 16 | LIGHTWEIGHT | 15x |

**Next:** See `PARALLEL_QUICK_REFERENCE.md` for one-page cheat sheet.
