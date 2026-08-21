"""
Example Scripts for Parallel Simulation

Copy and adapt these for your research!
"""

from parallel_simulator import ParallelSimulator, MemoryMode, ParallelStrategy
from src.gillespie.simulation_config import SimulationConfig
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Any
import time


# ── Example 1: Basic Parameter Sweep ──────────────────────────────────────

def example_1_basic_sweep():
    """
    Run a simple parameter sweep: vary mutation rate, run once per value.
    
    Time: ~30 minutes on 8 cores (for 50 parameter values)
    """
    print("\n" + "="*60)
    print("Example 1: Basic Parameter Sweep")
    print("="*60)
    
    # Setup
    base_config = SimulationConfig(T_max=1000, seed=42)
    parallel = ParallelSimulator(n_cores=8, memory_mode=MemoryMode.LIGHTWEIGHT)
    
    # Define sweep
    mutation_rates = np.linspace(0.01, 0.1, 50)
    
    # Run
    print(f"\nRunning {len(mutation_rates)} simulations...")
    t0 = time.time()
    results = parallel.sweep_parameter(
        param_name='mutation_rate',
        param_values=mutation_rates,
        base_config=base_config
    )
    elapsed = time.time() - t0
    
    # Analyze
    parallel.print_results_summary(results)
    print(f"Completed in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    
    # Summarize by parameter
    summary = ParallelSimulator.summarize_sweep(results, 'mutation_rate')
    
    # Plot
    param_vals = sorted(summary.keys())
    extinction_rates = [summary[p]['extinction_rate'] for p in param_vals]
    
    plt.figure(figsize=(10, 6))
    plt.plot(param_vals, extinction_rates, 'bo-', linewidth=2, markersize=6)
    plt.xlabel('Mutation Rate')
    plt.ylabel('Extinction Rate')
    plt.title('Example 1: Parameter Sweep')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('example1_sweep.png', dpi=150)
    print("\nPlot saved: example1_sweep.png")


# ── Example 2: Ensemble with Replicates ───────────────────────────────────

def example_2_ensemble_replicates():
    """
    Run each parameter value multiple times (different seeds) to get statistics.
    
    Time: ~2 hours on 8 cores (for 20 param values × 10 replicates = 200 sims)
    """
    print("\n" + "="*60)
    print("Example 2: Parameter Sweep with Replicates")
    print("="*60)
    
    base_config = SimulationConfig(T_max=1000)
    parallel = ParallelSimulator(n_cores=8, memory_mode=MemoryMode.LIGHTWEIGHT)
    
    # Setup: 20 parameter values × 10 replicates = 200 simulations
    mutation_rates = np.linspace(0.01, 0.1, 20)
    seeds = list(range(42, 52))  # 10 different seeds
    
    print(f"\nRunning {len(mutation_rates)} × {len(seeds)} = {len(mutation_rates)*len(seeds)} simulations...")
    t0 = time.time()
    results = parallel.sweep_parameter(
        param_name='mutation_rate',
        param_values=mutation_rates,
        base_config=base_config,
        seeds=seeds
    )
    elapsed = time.time() - t0
    
    print(f"Completed in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    parallel.print_results_summary(results)
    
    # Analyze with error bars
    summary = ParallelSimulator.summarize_sweep(results, 'mutation_rate')
    
    param_vals = sorted(summary.keys())
    means = [summary[p]['final_population_mean'] for p in param_vals]
    stds = [summary[p]['final_population_std'] for p in param_vals]
    extinctions = [summary[p]['extinction_rate'] for p in param_vals]
    
    # Plot 1: Population with error bars
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    ax1.errorbar(param_vals, means, yerr=stds, fmt='o-', capsize=5, capthick=2)
    ax1.set_xlabel('Mutation Rate')
    ax1.set_ylabel('Final Population')
    ax1.set_title('Mean ± Std')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Extinction rate
    ax2.plot(param_vals, extinctions, 'ro-', linewidth=2, markersize=6)
    ax2.set_xlabel('Mutation Rate')
    ax2.set_ylabel('Extinction Rate')
    ax2.set_title('Extinction Probability')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 1])
    
    plt.tight_layout()
    plt.savefig('example2_replicates.png', dpi=150)
    print("\nPlot saved: example2_replicates.png")
    
    # Print detailed results
    print("\nDetailed Results (10 replicates per value):")
    print(f"{'Mutation Rate':<15} {'Mean Pop':<15} {'Std':<15} {'Extinction %':<15}")
    print("-" * 60)
    for p in param_vals:
        stats = summary[p]
        print(f"{p:<15.4f} {stats['final_population_mean']:<15.1f} "
              f"{stats['final_population_std']:<15.1f} {stats['extinction_rate']*100:<15.1f}")


# ── Example 3: 2D Parameter Sweep ─────────────────────────────────────────

def example_3_2d_sweep():
    """
    Sweep two parameters simultaneously: mutation rate × exhaustion rate.
    
    Creates heatmap of results.
    
    Time: ~3 hours on 8 cores (for 15×15 = 225 sims)
    """
    print("\n" + "="*60)
    print("Example 3: 2D Parameter Sweep")
    print("="*60)
    
    base_config = SimulationConfig(T_max=1000)
    parallel = ParallelSimulator(n_cores=8, memory_mode=MemoryMode.LIGHTWEIGHT)
    
    # Setup 2D grid
    mutation_rates = np.linspace(0.01, 0.1, 15)
    exhaustion_rates = np.linspace(0.001, 0.01, 15)
    
    print(f"\nRunning {len(mutation_rates)} × {len(exhaustion_rates)} = "
          f"{len(mutation_rates)*len(exhaustion_rates)} simulations...")
    
    t0 = time.time()
    results, grid_info = parallel.sweep_parameters_2d(
        param1_name='mutation_rate',
        param1_values=mutation_rates,
        param2_name='exhaustion_rate',
        param2_values=exhaustion_rates,
        base_config=base_config
    )
    elapsed = time.time() - t0
    
    print(f"Completed in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    parallel.print_results_summary(results)
    
    # Organize into 2D arrays
    extinction_grid = np.zeros(grid_info['shape'])
    population_grid = np.zeros(grid_info['shape'])
    
    for result in results:
        if result.error is None:
            i = result.metadata['param1_idx']
            j = result.metadata['param2_idx']
            extinction_grid[i, j] = result.metadata['extinction']
            population_grid[i, j] = result.metadata['final_population']
    
    # Plot heatmaps
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Heatmap 1: Extinction
    im1 = ax1.imshow(extinction_grid, cmap='RdYlGn_r', aspect='auto', origin='lower')
    ax1.set_xlabel('Exhaustion Rate')
    ax1.set_ylabel('Mutation Rate')
    ax1.set_title('Extinction Probability')
    ax1.set_xticks(range(0, len(exhaustion_rates), 3))
    ax1.set_xticklabels([f'{exhaustion_rates[i]:.3f}' for i in range(0, len(exhaustion_rates), 3)], rotation=45)
    ax1.set_yticks(range(0, len(mutation_rates), 3))
    ax1.set_yticklabels([f'{mutation_rates[i]:.3f}' for i in range(0, len(mutation_rates), 3)])
    plt.colorbar(im1, ax=ax1)
    
    # Heatmap 2: Final population
    im2 = ax2.imshow(population_grid, cmap='viridis', aspect='auto', origin='lower')
    ax2.set_xlabel('Exhaustion Rate')
    ax2.set_ylabel('Mutation Rate')
    ax2.set_title('Final Population')
    ax2.set_xticks(range(0, len(exhaustion_rates), 3))
    ax2.set_xticklabels([f'{exhaustion_rates[i]:.3f}' for i in range(0, len(exhaustion_rates), 3)], rotation=45)
    ax2.set_yticks(range(0, len(mutation_rates), 3))
    ax2.set_yticklabels([f'{mutation_rates[i]:.3f}' for i in range(0, len(mutation_rates), 3)])
    plt.colorbar(im2, ax=ax2)
    
    plt.tight_layout()
    plt.savefig('example3_2d_sweep.png', dpi=150)
    print("\nPlot saved: example3_2d_sweep.png")


# ── Example 4: Ensemble Runs (Stochastic Behavior) ──────────────────────

def example_4_ensemble():
    """
    Run same configuration with 100 different seeds to characterize stochasticity.
    
    Time: ~30 minutes on 8 cores (for 100 sims)
    """
    print("\n" + "="*60)
    print("Example 4: Ensemble with Different Seeds")
    print("="*60)
    
    base_config = SimulationConfig(T_max=1000, mutation_rate=0.05)
    parallel = ParallelSimulator(n_cores=8, memory_mode=MemoryMode.LIGHTWEIGHT)
    
    # Create 100 configs with same parameters, different seeds
    n_replicates = 100
    configs = [
        SimulationConfig(**vars(base_config), seed=i)
        for i in range(42, 42 + n_replicates)
    ]
    
    print(f"\nRunning {n_replicates} replicates of same config...")
    t0 = time.time()
    results = parallel.run_ensemble(configs)
    elapsed = time.time() - t0
    
    print(f"Completed in {elapsed:.1f} seconds")
    
    # Aggregate
    agg = ParallelSimulator.aggregate_results(results)
    
    # Extract final populations
    final_pops = agg['final_populations']
    
    # Plot distribution
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram
    ax1.hist(final_pops, bins=20, edgecolor='black', alpha=0.7)
    ax1.axvline(np.mean(final_pops), color='r', linestyle='--', linewidth=2, label=f'Mean: {np.mean(final_pops):.0f}')
    ax1.axvline(np.median(final_pops), color='g', linestyle='--', linewidth=2, label=f'Median: {np.median(final_pops):.0f}')
    ax1.set_xlabel('Final Population')
    ax1.set_ylabel('Frequency')
    ax1.set_title(f'Distribution of Final Populations (n={n_replicates})')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Box plot
    ax2.boxplot(final_pops, vert=True)
    ax2.set_ylabel('Final Population')
    ax2.set_title('Box Plot')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('example4_ensemble.png', dpi=150)
    print("\nPlot saved: example4_ensemble.png")
    
    # Print statistics
    print("\nStatistics:")
    print(f"  Mean:           {np.mean(final_pops):.1f}")
    print(f"  Std:            {np.std(final_pops):.1f}")
    print(f"  Median:         {np.median(final_pops):.1f}")
    print(f"  Min:            {np.min(final_pops):.0f}")
    print(f"  Max:            {np.max(final_pops):.0f}")
    print(f"  Extinction rate: {agg['extinction_rate']:.1%}")


# ── Example 5: Population Trajectories ────────────────────────────────────

def example_5_trajectories():
    """
    Plot population trajectories from multiple runs.
    """
    print("\n" + "="*60)
    print("Example 5: Population Trajectories")
    print("="*60)
    
    base_config = SimulationConfig(T_max=1000, mutation_rate=0.05)
    parallel = ParallelSimulator(
        n_cores=8,
        memory_mode=MemoryMode.STANDARD  # Need history for trajectories
    )
    
    # Run 30 replicates
    configs = [
        SimulationConfig(**vars(base_config), seed=i)
        for i in range(42, 42 + 30)
    ]
    
    print(f"\nRunning {len(configs)} simulations...")
    results = parallel.run_ensemble(configs)
    
    # Extract trajectories
    times_list, trajectories = ParallelSimulator.extract_trajectories(results)
    
    # Plot all trajectories
    plt.figure(figsize=(12, 7))
    
    for i, (times, trajectory) in enumerate(zip(times_list, trajectories)):
        plt.plot(times, trajectory, alpha=0.4, linewidth=1)
    
    # Add mean trajectory
    if len(trajectories) > 0:
        # Interpolate to common time grid
        common_times = np.linspace(0, max(t[-1] for t in times_list), 100)
        mean_trajectory = np.zeros(len(common_times))
        for traj in trajectories:
            interp = np.interp(common_times, times_list[trajectories.index(traj)], traj)
            mean_trajectory += interp
        mean_trajectory /= len(trajectories)
        plt.plot(common_times, mean_trajectory, 'r-', linewidth=3, label='Mean')
    
    plt.xlabel('Time')
    plt.ylabel('Population')
    plt.title('Population Trajectories (30 replicates)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('example5_trajectories.png', dpi=150)
    print("Plot saved: example5_trajectories.png")


# ── Example 6: Strategy Comparison ────────────────────────────────────────

def example_6_strategy_comparison():
    """
    Compare different parallelization strategies.
    """
    print("\n" + "="*60)
    print("Example 6: Strategy Comparison")
    print("="*60)
    
    base_config = SimulationConfig(T_max=500)
    configs = [
        SimulationConfig(**vars(base_config), seed=i, mutation_rate=0.05 + 0.01*np.random.randn())
        for i in range(100)
    ]
    
    strategies = [
        ParallelStrategy.PROCESS_POOL,
        ParallelStrategy.CONCURRENT_FUTURES,
        ParallelStrategy.ASYNC_MAP
    ]
    
    times = {}
    
    for strategy in strategies:
        parallel = ParallelSimulator(n_cores=8, strategy=strategy, memory_mode=MemoryMode.LIGHTWEIGHT)
        t0 = time.time()
        results = parallel.run_ensemble(configs, show_progress=False)
        elapsed = time.time() - t0
        times[strategy.value] = elapsed
        
        successful = len([r for r in results if r.error is None])
        print(f"{strategy.value:20} {elapsed:8.2f} seconds ({successful}/{len(configs)} successful)")
    
    # Plot
    plt.figure(figsize=(10, 6))
    strategies_names = list(times.keys())
    strategy_times = list(times.values())
    
    plt.bar(strategies_names, strategy_times, color=['blue', 'green', 'red'], alpha=0.7)
    plt.ylabel('Time (seconds)')
    plt.title('Parallelization Strategy Comparison (100 sims, 8 cores)')
    plt.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for i, (strategy, t) in enumerate(zip(strategies_names, strategy_times)):
        plt.text(i, t + 2, f'{t:.1f}s', ha='center', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('example6_strategies.png', dpi=150)
    print("\nPlot saved: example6_strategies.png")


# ── Main: Run All Examples ────────────────────────────────────────────────

if __name__ == '__main__':
    print("\n" + "="*60)
    print("Parallel Simulation Examples")
    print("="*60)
    
    try:
        example_1_basic_sweep()
    except Exception as e:
        print(f"Example 1 failed: {e}")
    
    try:
        example_2_ensemble_replicates()
    except Exception as e:
        print(f"Example 2 failed: {e}")
    
    try:
        example_3_2d_sweep()
    except Exception as e:
        print(f"Example 3 failed: {e}")
    
    try:
        example_4_ensemble()
    except Exception as e:
        print(f"Example 4 failed: {e}")
    
    try:
        example_5_trajectories()
    except Exception as e:
        print(f"Example 5 failed: {e}")
    
    try:
        example_6_strategy_comparison()
    except Exception as e:
        print(f"Example 6 failed: {e}")
    
    print("\n" + "="*60)
    print("All examples completed!")
    print("Check *.png files for plots")
    print("="*60)
