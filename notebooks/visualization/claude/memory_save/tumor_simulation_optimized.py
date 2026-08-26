from __future__ import annotations

from enum import Enum
from typing import Dict, List, Optional, Tuple

import numpy as np
from tqdm import tqdm
from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone
from src.gillespie.tissue_state import TissueState
from src.gillespie.event import Event
from src.gillespie.event_type import EventType
from src.gillespie.rate_matrix import RateMatrix
from src.gillespie.crowding_strategy import CrowdingStrategy
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone_factory import CloneFactory


class MemoryMode(Enum):
    """
    LIGHTWEIGHT: Only final state. Minimal memory footprint.
    STANDARD: Final state + periodic snapshots. Balanced.
    FULL: All events, rates, and snapshots. Maximum fidelity.
    """
    LIGHTWEIGHT = "lightweight"  # ~1-2 MB even for large sims
    STANDARD = "standard"         # ~10-100 MB depending on save_interval
    FULL = "full"                 # ~100 MB - 1 GB+ for long simulations


class TumorSimulation:
    def __init__(self, config: SimulationConfig, memory_mode: MemoryMode = MemoryMode.STANDARD) -> None:
        self.config = config
        self.memory_mode = memory_mode
        
        # Configure memory strategy based on mode
        self._configure_memory_mode()
        
        self.step_count = 0
        self.rng = np.random.default_rng(config.seed)
        self.clone_factory = CloneFactory(config)
        self.t = 0.0
        
        print(f"Memory mode: {self.memory_mode.value.upper()}")
        print("Starting with parameters:")
        print(config)
        
        # Initialize clones
        clones_dict: Dict[CloneId, Clone] = {
            (): self.clone_factory.create_clone(clone_type="base"),
            (-3,): self.clone_factory.create_clone(clone_type="mutated"),
            (-1,): self.clone_factory.create_clone(clone_type="immune"),
            (-2,): self.clone_factory.create_clone(clone_type="exhausted")
        }
        
        self.times: List[float] = [0.0]
        self.tissue_state: TissueState = TissueState(t=self.t, clones=clones_dict)
        self.crowding_strategy: CrowdingStrategy = self.config.crowding_strategy

        # Initialize history containers based on memory mode
        if self._should_save_history():
            self.history: List[Dict[CloneId, dict]] = [self.tissue_state.snapshot()]
            self.rate_history: List[List[Dict]] = [] if self.memory_mode == MemoryMode.FULL else None
        else:
            self.history: List[Dict[CloneId, dict]] = []
            self.rate_history: Optional[List[List[Dict]]] = None
        
        # Only allocate events list if saving all steps
        self.events: Optional[List[Event]] = [] if self.memory_mode == MemoryMode.FULL else None
        
        rates0 = self._build_rate_matrix()
        print("STARTING RATES AND STATE")
        self.tissue_state.print_pop_map()
        self.print_event_table(rates0.events)

    # ── Memory Configuration ──────────────────────────────────────────────────

    def _configure_memory_mode(self) -> None:
        """Configure save intervals and history tracking based on memory mode."""
        if self.memory_mode == MemoryMode.LIGHTWEIGHT:
            self.save_interval = 0  # Don't save any snapshots except final
            self.save_all_steps = False
        elif self.memory_mode == MemoryMode.STANDARD:
            # Use config value if provided, otherwise reasonable default
            self.save_interval = self.config.save_interval if self.config.save_interval > 0 else 100
            self.save_all_steps = False
        else:  # FULL
            self.save_interval = self.config.save_interval if self.config.save_interval > 0 else 1
            self.save_all_steps = True

    def _should_save_history(self) -> bool:
        """Check if we should save intermediate history snapshots."""
        return self.save_interval > 0

    def _should_save_rates(self) -> bool:
        """Check if we should save rate history."""
        return self.memory_mode == MemoryMode.FULL

    # ── Internal helpers ──────────────────────────────────────────────────────

    def create_clone(
        self,
        clone_id: CloneId,
        clone_type: str = "base",
        N: int = None,
        parent: Optional[Clone] = None,
    ) -> Clone:
        clone = self.clone_factory.create_clone(
            clone_id=clone_id,
            clone_type=clone_type,
            N=N,
            parent=parent,
        )
        self.tissue_state.clones[clone.clone_id] = clone
        self.tissue_state.update_pop_map()
        return clone

    def _advance_all_instability(self, dt: float) -> None:
        for clone in self.tissue_state.clones.values():
            if clone.is_alive():
                clone.advance_instability(dt, self.config.base_instability_buildup)

    def _build_rate_matrix(self) -> RateMatrix:
        self.tissue_state.update_pop_map()
        rate_matrix = RateMatrix()

        for cid, clone in self.tissue_state.clones.items():
            if not clone.is_alive():
                continue
            
            type_rates: tuple = (
                clone.birth_rate_effective(tissue_state=self.tissue_state),
                clone.death_rate_effective(tissue_state=self.tissue_state),
                clone.mutation_rate_effective(tissue_state=self.tissue_state),
                clone.exhaustion_rate_effective(tissue_state=self.tissue_state),
            )

            rb, rd, rm, re = type_rates

            if rb > 0:
                rate_matrix.add_event(
                    Event(kind=EventType.BIRTH, clone_id=cid, rate=rb, clone_type=clone.get_type())
                )
            if rd > 0:
                rate_matrix.add_event(
                    Event(kind=EventType.DEATH, clone_id=cid, rate=rd, clone_type=clone.get_type())
                )
            if rm > 0:
                rate_matrix.add_event(
                    Event(kind=EventType.MUTATION, clone_id=cid, rate=rm, clone_type=clone.get_type())
                )
            if re > 0:
                rate_matrix.add_event(
                    Event(kind=EventType.EXHAUSTION, clone_id=cid, rate=re, clone_type=clone.get_type())
                )
        return rate_matrix

    def _sample_waiting_time(self, total_rate: float) -> float:
        if total_rate <= 0:
            raise ValueError("Total rate must be positive.")
        return -np.log(self.rng.random()) / total_rate

    def _introduce_mutation(self, clone: Clone) -> None:
        assert clone.is_alive(), "Cannot mutate a dead clone."
        assert clone.next_mutation != "", "Clone is not supposed to be mutating"
        
        self.create_clone(
            clone_id=clone.next_child_id(),
            clone_type=clone.next_mutation,
            N=2,
            parent=clone.clone_id,
        )
        clone.kill()
    
    def _induce_exhaustion(self, clone: Clone) -> None:
        self.tissue_state.clones[clone.clone_id].kill()
        self.tissue_state.clones[(-2,)].divide()
    
    def _apply_event(self, event: Event) -> None:
        clone = self.tissue_state.clones[event.clone_id]
        if event.kind == EventType.BIRTH:
            clone.divide()
        elif event.kind == EventType.DEATH:
            clone.kill()
        elif event.kind == EventType.MUTATION:
            self._introduce_mutation(clone)
        elif event.kind == EventType.EXHAUSTION:
            self._induce_exhaustion(clone)
        else:
            raise ValueError(f"Unknown event kind: {event.kind}")

    def print_event_table(self, events: List[Event]) -> None:
        """Print formatted event rate table."""
        header = f"{'Kind':<12}  | {'Type':<12}| {'Clone ID':<15} | {'Rate':<10} | {'N':<15}"
        print(header)
        print("-" * len(header))
        
        for e in events:
            id_str = str(e.clone_id)
            print(f"{e.kind.name:<12} | {e.clone_type:<12} | {id_str:<15} | {e.rate:<10.4f}| {self.tissue_state.clones[e.clone_id].N}")
        
        print("-" * len(header))

    # ── Public API ────────────────────────────────────────────────────────────

    def total_population(self) -> int:
        return self.tissue_state.total_population()

    def step(self, return_matrix: Optional[bool] = False) -> bool:
        """Advance by one Gillespie step. Returns False when simulation should stop."""
        self.step_count += 1
        rate_matrix = self._build_rate_matrix()
  
        if return_matrix:
            print(self.t)
            self.print_event_table(events=rate_matrix.events)
        
        total_rate = rate_matrix.get_total_rate()
        
        # Only record rates if in FULL mode
        if self._should_save_rates():
            step_rates = [
                {
                    "time": self.t,
                    "kind": e.kind.name,
                    "clone_id": e.clone_id,
                    "clone_type": e.clone_type,
                    "rate": e.rate,
                } for e in rate_matrix.events
            ]
            self.rate_history.append(step_rates)
        
        if total_rate <= 0 or not rate_matrix.events:
            return False

        tau = self._sample_waiting_time(total_rate)
        new_t = self.t + tau

        # Handle end-of-simulation
        if new_t > self.config.T_max:
            tau = self.config.T_max - self.t
            self._advance_all_instability(tau)
            self.t = self.config.T_max
            self.times.append(self.t)
            
            # Always save final state
            self.history.append(self.tissue_state.snapshot())
            print("FINAL RATES AND STATE")
            self.tissue_state.print_pop_map()
            
            return False

        self._advance_all_instability(tau)
        self.t = new_t
        
        event = rate_matrix.choose_event(self.rng.random())
        self._apply_event(event)
        self.tissue_state.update_pop_map()
        
        self.times.append(self.t)
        self.tissue_state.t = self.t
        
        # Save snapshots conditionally
        if self._should_save_history() and self.step_count % self.save_interval == 0:
            self.history.append(self.tissue_state.snapshot())
        
        # Save events only in FULL mode
        if self.memory_mode == MemoryMode.FULL:
            self.events.append(event)
                
        return True
        
    def _stopping_cond(self) -> bool:
        """Check stopping conditions."""
        # Total extinction
        if self.total_population() <= 0:
            print("Total extinction")
            return True
        
        return False

    def run(self) -> Tuple[List[float], List[Dict[CloneId, dict]], TissueState, Optional[List[List[Dict]]]]:
        """
        Run simulation to completion.
        
        Returns:
            Tuple of (times, history, tissue_state, rate_history)
            - rate_history is None in LIGHTWEIGHT/STANDARD modes
            - history contains only final state in LIGHTWEIGHT mode
        """
        with tqdm(desc="Simulating", unit=" steps") as pbar:
            while not self._stopping_cond():
                if not self.step():
                    break
                pbar.update(1)
                pbar.set_postfix({"time": f"{self.t:.2f}/{self.config.T_max:.2f}"})
        
        return self.times, self.history, self.tissue_state, self.rate_history

    # ── Memory inspection ─────────────────────────────────────────────────────

    def estimate_memory_usage(self) -> Dict[str, float]:
        """
        Estimate memory usage in MB.
        
        Returns:
            Dict with breakdown of memory usage by component
        """
        import sys
        
        estimates = {
            "times": sys.getsizeof(self.times) / 1e6,
            "history": sys.getsizeof(self.history) / 1e6,
            "tissue_state": sys.getsizeof(self.tissue_state) / 1e6,
        }
        
        if self.rate_history is not None:
            estimates["rate_history"] = sys.getsizeof(self.rate_history) / 1e6
        
        if self.events is not None:
            estimates["events"] = sys.getsizeof(self.events) / 1e6
        
        estimates["total"] = sum(estimates.values())
        return estimates

    def print_memory_summary(self) -> None:
        """Print memory usage summary."""
        usage = self.estimate_memory_usage()
        print("\n" + "=" * 50)
        print("Memory Usage Summary")
        print("=" * 50)
        for component, mb in usage.items():
            print(f"  {component:<20}: {mb:>8.2f} MB")
        print("=" * 50 + "\n")
