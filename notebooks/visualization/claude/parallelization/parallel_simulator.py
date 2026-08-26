"""
Parallel Gillespie Simulator
Runs multiple simulations in parallel across multiple CPU cores.

Perfect for parameter sweeps, ensemble runs, and large-scale studies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple, Any
from enum import Enum
import multiprocessing as mp
from functools import partial
import numpy as np
from tqdm import tqdm

from src.gillespie.simulation_config import SimulationConfig
from tumor_simulation_optimized import TumorSimulation, MemoryMode


class ParallelStrategy(Enum):
    """Parallelization strategy selection"""
    PROCESS_POOL = "process_pool"      # Classic multiprocessing.Pool
    CONCURRENT_FUTURES = "futures"      # concurrent.futures.ProcessPoolExecutor
    ASYNC_MAP = "async_map"             # Pool.imap_unordered for streaming


@dataclass
class ParallelResult:
    """Container for results from a single parallel simulation"""
    config_id: int
    times: List[float]
    history: List[Dict]
    final_state: Any
    rate_history: Optional[List[List[Dict]]] = None
    error: Optional[str] = None
    metadata: Dict = field(default_factory=dict)


@dataclass
class ParameterSweep:
    """Defines a parameter sweep across one or more parameters"""
    param_name: str
    param_values: List[float]
    base_config: SimulationConfig
    
    def __post_init__(self):
        self.n_configs = len(self.param_values)
    
    def get_config(self, idx: int) -> SimulationConfig:
        """Get a config with the parameter set to param_values[idx]"""
        config = self.base_config.__class__(**vars(self.base_config))
        setattr(config, self.param_name, self.param_values[idx])
        return config


@dataclass
class MultiParameterSweep:
    """Defines a multi-dimensional parameter sweep"""
    param_sweeps: List[ParameterSweep]
    base_config: SimulationConfig
    
    def __post_init__(self):
        self.n_configs = np.prod([s.n_configs for s in self.param_sweeps])
    
    def get_config(self, idx: int) -> SimulationConfig:
        """Get config from linear index (Cartesian product)"""
        config = self.base_config.__class__(**vars(self.base_config))
        
        # Convert linear index to multi-dimensional indices
        indices = np.unravel_index(idx, tuple(s.n_configs for s in self.param_sweeps))
        
        for param_sweep, param_idx in zip(self.param_sweeps, indices):
            setattr(config, param_sweep.param_name, param_sweep.param_values[param_idx])
        
        return config


# ── Serializable Worker Function ──────────────────────────────────────────

def _run_simulation_worker(
    config: SimulationConfig,
    config_id: int,
    memory_mode: MemoryMode = MemoryMode.LIGHTWEIGHT,
    return_rate_history: bool = False
) -> ParallelResult:
    """
    Worker function for multiprocessing.
    
    MUST be defined at module level to be pickleable.
    
    Args:
        config: SimulationConfig for this run
        config_id: Identifier for this configuration
        memory_mode: Which memory mode to use
        return_rate_history: Whether to return rate history
    
    Returns:
        ParallelResult with simulation outputs
    """
    try:
        # Create and run simulation
        sim = TumorSimulation(config, memory_mode=memory_mode)
        times, history, final_state, rate_history = sim.run()
        
        # Package results
        return ParallelResult(
            config_id=config_id,
            times=times,
            history=history,
            final_state=final_state,
            rate_history=rate_history if return_rate_history else None,
            metadata={
                'memory_mode': memory_mode.value,
                'seed': config.seed,
                'final_population': final_state.total_population(),
                'extinction': final_state.total_population() == 0
            }
        )
    except Exception as e:
        # Gracefully handle errors in worker processes
        return ParallelResult(
            config_id=config_id,
            times=[],
            history=[],
            final_state=None,
            error=str(e),
            metadata={'error_type': type(e).__name__}
        )


class ParallelSimulator:
    """
    Orchestrates parallel simulation runs across multiple CPU cores.
    
    Usage:
        parallel = ParallelSimulator(n_cores=8, strategy=ParallelStrategy.PROCESS_POOL)
        results = parallel.run_ensemble(configs)
        parallel.aggregate_results(results)
    """
    
    def __init__(
        self,
        n_cores: Optional[int] = None,
        strategy: ParallelStrategy = ParallelStrategy.PROCESS_POOL,
        memory_mode: MemoryMode = MemoryMode.LIGHTWEIGHT,
        verbose: bool = True
    ):
        """
        Args:
            n_cores: Number of CPU cores to use. None = use all available.
            strategy: Parallelization strategy to use
            memory_mode: Memory mode for simulations (LIGHTWEIGHT recommended for parallel)
            verbose: Show progress bars
        """
        self.n_cores = n_cores or mp.cpu_count()
        self.strategy = strategy
        self.memory_mode = memory_mode
        self.verbose = verbose
        
        if self.verbose:
            print(f"ParallelSimulator initialized:")
            print(f"  CPU cores: {self.n_cores}")
            print(f"  Strategy: {self.strategy.value}")
            print(f"  Memory mode: {self.memory_mode.value}")

    # ── Ensemble Methods ──────────────────────────────────────────────────

    def run_ensemble(
        self,
        configs: List[SimulationConfig],
        return_rate_history: bool = False,
        show_progress: bool = True
    ) -> List[ParallelResult]:
        """
        Run multiple simulations in parallel (different seeds or parameters).
        
        Args:
            configs: List of SimulationConfig objects to run
            return_rate_history: Keep rate history (increases memory, set False for LIGHTWEIGHT)
            show_progress: Show progress bar
        
        Returns:
            List of ParallelResult objects
        
        Example:
            configs = [SimulationConfig(seed=i) for i in range(100)]
            parallel = ParallelSimulator(n_cores=8)
            results = parallel.run_ensemble(configs)
        """
        if self.strategy == ParallelStrategy.PROCESS_POOL:
            return self._run_ensemble_pool(configs, return_rate_history, show_progress)
        elif self.strategy == ParallelStrategy.CONCURRENT_FUTURES:
            return self._run_ensemble_futures(configs, return_rate_history, show_progress)
        else:  # ASYNC_MAP
            return self._run_ensemble_async(configs, return_rate_history, show_progress)

    def _run_ensemble_pool(
        self,
        configs: List[SimulationConfig],
        return_rate_history: bool,
        show_progress: bool
    ) -> List[ParallelResult]:
        """Run using multiprocessing.Pool"""
        worker_fn = partial(
            _run_simulation_worker,
            memory_mode=self.memory_mode,
            return_rate_history=return_rate_history
        )
        
        results = []
        with mp.Pool(processes=self.n_cores) as pool:
            # Create tasks with config_id
            tasks = [(config, i) for i, config in enumerate(configs)]
            
            if show_progress and self.verbose:
                # Use imap for progress bar
                iterator = pool.starmap(worker_fn, tasks, chunksize=max(1, len(configs)//self.n_cores))
                pbar = tqdm(iterator, total=len(configs), desc="Running simulations")
                results = list(pbar)
            else:
                # Faster without progress bar
                results = pool.starmap(worker_fn, tasks)
        
        return results

    def _run_ensemble_futures(
        self,
        configs: List[SimulationConfig],
        return_rate_history: bool,
        show_progress: bool
    ) -> List[ParallelResult]:
        """Run using concurrent.futures.ProcessPoolExecutor"""
        from concurrent.futures import ProcessPoolExecutor
        
        worker_fn = partial(
            _run_simulation_worker,
            memory_mode=self.memory_mode,
            return_rate_history=return_rate_history
        )
        
        results = [None] * len(configs)
        
        with ProcessPoolExecutor(max_workers=self.n_cores) as executor:
            # Submit all tasks
            futures = {
                executor.submit(worker_fn, config, i): i 
                for i, config in enumerate(configs)
            }
            
            # Collect results as they complete
            iterator = futures if not (show_progress and self.verbose) else tqdm(
                futures,
                total=len(configs),
                desc="Running simulations"
            )
            
            for future in iterator:
                result = future.result()
                results[result.config_id] = result
        
        return results

    def _run_ensemble_async(
        self,
        configs: List[SimulationConfig],
        return_rate_history: bool,
        show_progress: bool
    ) -> List[ParallelResult]:
        """Run using Pool.imap_unordered for streaming results"""
        worker_fn = partial(
            _run_simulation_worker,
            memory_mode=self.memory_mode,
            return_rate_history=return_rate_history
        )
        
        with mp.Pool(processes=self.n_cores) as pool:
            tasks = [(config, i) for i, config in enumerate(configs)]
            
            if show_progress and self.verbose:
                iterator = pool.imap_unordered(
                    worker_fn,
                    [(c, i) for c, i in tasks],
                    chunksize=max(1, len(configs)//(self.n_cores*4))
                )
                results = list(tqdm(iterator, total=len(configs), desc="Running simulations"))
            else:
                results = list(pool.imap_unordered(worker_fn, tasks))
            
            # Sort by config_id to restore order
            results.sort(key=lambda r: r.config_id)
        
        return results

    # ── Parameter Sweep Methods ───────────────────────────────────────────

    def sweep_parameter(
        self,
        param_name: str,
        param_values: List[float],
        base_config: SimulationConfig,
        seeds: Optional[List[int]] = None,
        return_rate_history: bool = False
    ) -> List[ParallelResult]:
        """
        Run parameter sweep across a single parameter.
        
        Args:
            param_name: Name of parameter to sweep (e.g., 'mutation_rate')
            param_values: List of values to test
            base_config: Base configuration (will modify param_name for each value)
            seeds: Optional list of seeds for ensemble (if None, uses config.seed)
            return_rate_history: Keep rate history (set False for memory efficiency)
        
        Returns:
            List of ParallelResult objects
        
        Example:
            results = parallel.sweep_parameter(
                param_name='mutation_rate',
                param_values=np.linspace(0.01, 0.1, 50),
                base_config=my_config,
                seeds=[42, 123, 456]  # Run each param value 3 times
            )
        """
        configs = []
        config_ids = []
        
        # Create configuration for each parameter value (and seed if provided)
        for param_idx, param_val in enumerate(param_values):
            if seeds is None:
                # Single run per parameter value
                config = base_config.__class__(**vars(base_config))
                setattr(config, param_name, param_val)
                configs.append(config)
                config_ids.append((param_idx, 0))
            else:
                # Multiple runs per parameter value (different seeds)
                for seed_idx, seed in enumerate(seeds):
                    config = base_config.__class__(**vars(base_config))
                    setattr(config, param_name, param_val)
                    config.seed = seed
                    configs.append(config)
                    config_ids.append((param_idx, seed_idx))
        
        if self.verbose:
            print(f"\nParameter sweep: {param_name}")
            print(f"  Values: {len(param_values)}")
            print(f"  Replicates per value: {len(seeds) if seeds else 1}")
            print(f"  Total configurations: {len(configs)}")
            print(f"  Expected time reduction: ~{self.n_cores}x (ideal)")
        
        results = self.run_ensemble(configs, return_rate_history=return_rate_history)
        
        # Attach param info for later aggregation
        for result, (param_idx, seed_idx) in zip(results, config_ids):
            result.metadata['param_name'] = param_name
            result.metadata['param_value'] = param_values[param_idx]
            result.metadata['param_idx'] = param_idx
            result.metadata['replicate'] = seed_idx
        
        return results

    def sweep_parameters_2d(
        self,
        param1_name: str,
        param1_values: List[float],
        param2_name: str,
        param2_values: List[float],
        base_config: SimulationConfig,
        return_rate_history: bool = False
    ) -> Tuple[List[ParallelResult], Dict]:
        """
        Run 2D parameter sweep.
        
        Example:
            results, grid_info = parallel.sweep_parameters_2d(
                param1_name='mutation_rate',
                param1_values=np.linspace(0.01, 0.1, 10),
                param2_name='exhaustion_rate',
                param2_values=np.linspace(0.001, 0.01, 10),
                base_config=my_config
            )
        """
        configs = []
        
        for p1_val in param1_values:
            for p2_val in param2_values:
                config = base_config.__class__(**vars(base_config))
                setattr(config, param1_name, p1_val)
                setattr(config, param2_name, p2_val)
                configs.append(config)
        
        if self.verbose:
            print(f"\n2D Parameter sweep:")
            print(f"  {param1_name}: {len(param1_values)} values")
            print(f"  {param2_name}: {len(param2_values)} values")
            print(f"  Total configurations: {len(configs)}")
        
        results = self.run_ensemble(configs, return_rate_history=return_rate_history)
        
        # Attach metadata
        idx = 0
        for p1_idx, p1_val in enumerate(param1_values):
            for p2_idx, p2_val in enumerate(param2_values):
                results[idx].metadata.update({
                    'param1_name': param1_name,
                    'param1_value': p1_val,
                    'param1_idx': p1_idx,
                    'param2_name': param2_name,
                    'param2_value': p2_val,
                    'param2_idx': p2_idx,
                })
                idx += 1
        
        grid_info = {
            'param1': {'name': param1_name, 'values': param1_values},
            'param2': {'name': param2_name, 'values': param2_values},
            'shape': (len(param1_values), len(param2_values))
        }
        
        return results, grid_info

    # ── Result Aggregation ────────────────────────────────────────────────

    @staticmethod
    def aggregate_results(
        results: List[ParallelResult],
        statistic: str = 'mean'
    ) -> Dict[str, Any]:
        """
        Aggregate results from ensemble runs.
        
        Args:
            results: List of ParallelResult objects
            statistic: 'mean', 'std', 'all' for aggregation
        
        Returns:
            Dict with aggregated statistics
        
        Example:
            agg = ParallelSimulator.aggregate_results(results, statistic='mean')
            print(agg['final_population_mean'])
            print(agg['extinction_rate'])
        """
        # Filter out failed runs
        valid_results = [r for r in results if r.error is None]
        failed_runs = len(results) - len(valid_results)
        
        if not valid_results:
            return {'error': 'All runs failed'}
        
        final_pops = [r.metadata['final_population'] for r in valid_results]
        extinctions = [r.metadata['extinction'] for r in valid_results]
        
        aggregated = {
            'n_total': len(results),
            'n_valid': len(valid_results),
            'n_failed': failed_runs,
            'final_population_mean': np.mean(final_pops),
            'final_population_std': np.std(final_pops),
            'final_population_min': np.min(final_pops),
            'final_population_max': np.max(final_pops),
            'extinction_rate': np.mean(extinctions),
            'final_populations': final_pops,
        }
        
        return aggregated

    @staticmethod
    def summarize_sweep(
        results: List[ParallelResult],
        param_name: str
    ) -> Dict[str, Any]:
        """
        Summarize parameter sweep results.
        
        Returns results organized by parameter value.
        
        Example:
            summary = ParallelSimulator.summarize_sweep(results, 'mutation_rate')
            for param_val, stats in summary.items():
                print(f"Param={param_val}: extinction_rate={stats['extinction_rate']:.2%}")
        """
        # Group by parameter value
        by_param = {}
        for result in results:
            if result.error is not None:
                continue
            
            param_val = result.metadata.get('param_value')
            if param_val not in by_param:
                by_param[param_val] = []
            
            by_param[param_val].append(result)
        
        # Summarize each parameter value
        summary = {}
        for param_val in sorted(by_param.keys()):
            group_results = by_param[param_val]
            group_agg = ParallelSimulator.aggregate_results(group_results)
            summary[param_val] = group_agg
        
        return summary

    # ── Utility Methods ───────────────────────────────────────────────────

    def print_results_summary(self, results: List[ParallelResult]) -> None:
        """Pretty-print summary of results"""
        agg = self.aggregate_results(results)
        
        print("\n" + "=" * 60)
        print("Ensemble Results Summary")
        print("=" * 60)
        print(f"Total runs:         {agg['n_total']}")
        print(f"Successful:         {agg['n_valid']}")
        print(f"Failed:             {agg['n_failed']}")
        print(f"\nFinal Population:")
        print(f"  Mean:             {agg['final_population_mean']:.1f}")
        print(f"  Std:              {agg['final_population_std']:.1f}")
        print(f"  Range:            [{agg['final_population_min']:.0f}, {agg['final_population_max']:.0f}]")
        print(f"\nExtinction Rate:    {agg['extinction_rate']:.1%}")
        print("=" * 60 + "\n")

    @staticmethod
    def extract_trajectories(results: List[ParallelResult]) -> Tuple[List, List]:
        """
        Extract population trajectories from results.
        
        Returns:
            (times_list, populations_list) where each element is from one run
        """
        trajectories = []
        times_list = []
        
        for result in results:
            if result.error is not None:
                continue
            
            # Reconstruct population trajectory from history
            populations = []
            for snapshot in result.history:
                pop = sum(c.get('N', 0) for c in snapshot.values())
                populations.append(pop)
            
            trajectories.append(populations)
            times_list.append(result.times)
        
        return times_list, trajectories
