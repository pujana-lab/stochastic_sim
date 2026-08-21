# Parallel Simulation Quick Reference

## One-Line Summary
Use `ParallelSimulator` to run parameter sweeps across multiple CPU cores with **~7-8x speedup** on an 8-core machine.

---

## Setup (Copy-Paste)

```python
from parallel_simulator import ParallelSimulator, MemoryMode
from src.gillespie.simulation_config import SimulationConfig
import numpy as np

# Create simulator
parallel = ParallelSimulator(n_cores=8, memory_mode=MemoryMode.LIGHTWEIGHT)
base_config = SimulationConfig(T_max=1000, seed=42)
```

---

## Common Patterns

### 1. Sweep One Parameter

```python
results = parallel.sweep_parameter(
    param_name='mutation_rate',
    param_values=np.linspace(0.01, 0.1, 50),
    base_config=base_config
)

# Summarize
summary = ParallelSimulator.summarize_sweep(results, 'mutation_rate')
for param_val, stats in summary.items():
    print(f"{param_val}: extinction={stats['extinction_rate']:.1%}")
```

---

### 2. Sweep with Multiple Replicates

```python
results = parallel.sweep_parameter(
    param_name='mutation_rate',
    param_values=np.linspace(0.01, 0.1, 20),
    base_config=base_config,
    seeds=list(range(42, 52))  # 10 replicates each
)
```

---

### 3. Sweep Two Parameters (Grid)

```python
results, grid_info = parallel.sweep_parameters_2d(
    param1_name='mutation_rate',
    param1_values=np.linspace(0.01, 0.1, 15),
    param2_name='exhaustion_rate', 
    param2_values=np.linspace(0.001, 0.01, 15),
    base_config=base_config
)

# Get shape
print(grid_info['shape'])  # (15, 15) = 225 configs
```

---

### 4. Ensemble (Same Config, Different Seeds)

```python
configs = [
    SimulationConfig(**vars(base_config), seed=i)
    for i in range(100)
]
results = parallel.run_ensemble(configs)

agg = ParallelSimulator.aggregate_results(results)
print(f"Mean extinction: {agg['extinction_rate']:.1%}")
```

---

## Memory Modes

| Mode | Memory/Sim | Total (8 cores) | Use When |
|------|-----------|-----------------|----------|
| LIGHTWEIGHT | 1-2 MB | 8-16 MB | Parameter sweeps |
| STANDARD | 10-50 MB | 80-400 MB | Trajectory analysis needed |
| FULL | 100 MB-1 GB | 800 MB-8 GB | Debugging, detailed analysis |

```python
# Fast & lean (default)
ParallelSimulator(memory_mode=MemoryMode.LIGHTWEIGHT)

# Balanced
ParallelSimulator(memory_mode=MemoryMode.STANDARD)

# Detailed
ParallelSimulator(memory_mode=MemoryMode.FULL)
```

---

## Strategies

```python
from parallel_simulator import ParallelStrategy

# PROCESS_POOL: Usually best (default)
parallel = ParallelSimulator(strategy=ParallelStrategy.PROCESS_POOL)

# CONCURRENT_FUTURES: Results as they complete
parallel = ParallelSimulator(strategy=ParallelStrategy.CONCURRENT_FUTURES)

# ASYNC_MAP: Streaming (good for 1000+ jobs)
parallel = ParallelSimulator(strategy=ParallelStrategy.ASYNC_MAP)
```

---

## Monitoring & Debugging

```python
# Auto progress bar
results = parallel.run_ensemble(configs, show_progress=True)

# Check for errors
failed = [r for r in results if r.error is not None]
print(f"Failed: {len(failed)}/{len(results)}")

# Print summary
parallel.print_results_summary(results)

# Extract trajectories for plotting
times_list, trajectories = ParallelSimulator.extract_trajectories(results)
```

---

## Return Format

```python
# Each result has:
# - config_id: which config this was
# - times: simulation times
# - history: population snapshots  
# - final_state: TissueState at end
# - rate_history: rates (if FULL mode)
# - error: error message (if failed)
# - metadata: {'final_population': int, 'extinction': bool, ...}

for result in results:
    print(f"Config {result.config_id}:")
    print(f"  Final pop: {result.metadata['final_population']}")
    print(f"  Extinct: {result.metadata['extinction']}")
```

---

## Speedup Expectations

```
Cores   Speedup    Efficiency
1       1.0x       100%
2       1.9x       95%
4       3.8x       95%
8       7.5x       94%
16      14.5x      91%
```

**Time estimate:** `sequential_time / n_cores × 1.1` (1.1x for overhead)

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| "RuntimeError: start a new process" (macOS/Windows) | Wrap code in `if __name__ == '__main__':` |
| Not using all cores | Check `parallel.n_cores` and `cpu_count()` |
| Out of memory | Use `LIGHTWEIGHT` mode or reduce `n_cores` |
| "PicklingError" | Avoid lambdas in config; use named functions |
| Simulations slow | Check if other processes are running; reduce `n_cores` |

---

## Performance Checklist

- [x] Use `LIGHTWEIGHT` mode for parameter sweeps
- [x] Use `n_cores = cpu_count() - 1` (leave one core free)
- [x] Put code in `if __name__ == '__main__':` (macOS/Windows)
- [x] Use `show_progress=True` to monitor
- [x] Check `parallel.print_results_summary()` at end

---

## Example: Full Workflow

```python
if __name__ == '__main__':
    # Setup
    parallel = ParallelSimulator(n_cores=8, memory_mode=MemoryMode.LIGHTWEIGHT)
    base_config = SimulationConfig(T_max=1000)
    
    # Run sweep
    results = parallel.sweep_parameter(
        param_name='mutation_rate',
        param_values=np.linspace(0.01, 0.1, 50),
        base_config=base_config,
        seeds=[42, 123, 456]  # 3 replicates
    )
    
    # Analyze
    summary = ParallelSimulator.summarize_sweep(results, 'mutation_rate')
    
    # Print summary
    parallel.print_results_summary(results)
    
    # Extract for plotting
    times_list, trajectories = ParallelSimulator.extract_trajectories(results)
```

---

## CPU Cores vs. RAM

```
RAM          Safe n_cores with LIGHTWEIGHT
8 GB         4-6 cores
16 GB        8-12 cores
32 GB        16-24 cores
64 GB        32-48 cores
```

Use `n_cores = min(cpu_count(), ram_gb // 2)`

---

## For Your Thesis Work

```python
# Recommended starting point
parallel = ParallelSimulator(
    n_cores=8,  # Or your available cores
    memory_mode=MemoryMode.LIGHTWEIGHT,
    verbose=True
)

# For BRCA1 parameter studies
results = parallel.sweep_parameter(
    param_name='mutation_rate',
    param_values=...,  # Your mutation rate range
    base_config=base_config,
    seeds=list(range(42, 52))  # 10 replicates
)

# Report results
parallel.print_results_summary(results)
summary = ParallelSimulator.summarize_sweep(results, 'mutation_rate')
```

Estimated time: 50 parameter values × 10 replicates = 500 simulations
- Sequential: ~50 hours
- Parallel (8 cores): ~6-7 hours
- **Time saved: ~43 hours** ✨

---

## Next Steps

1. Read `PARALLEL_SIMULATION_GUIDE.md` for detailed examples
2. Try Pattern 1 (simple sweep) first
3. Graduate to Pattern 2 (multiple replicates)
4. Use Pattern 3 (2D grid) for publication-quality results
