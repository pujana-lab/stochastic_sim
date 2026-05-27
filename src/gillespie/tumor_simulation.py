from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone
from src.gillespie.tissue_state import TissueState
from src.gillespie.event import Event
from src.gillespie.event_type import EventType
from src.gillespie.rate_matrix import RateMatrix
from src.gillespie.crowding_strategy import CrowdingStrategy, SimpleCrowding, AdaptedCrowding
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone_factory import CloneFactory

class TumorSimulation:
    def __init__(self, config: SimulationConfig) -> None:
        self.config = config
        self.rng = np.random.default_rng(config.seed)
        self.clone_factory = CloneFactory(config)
        self.t = 0.0
        
        # Initialize clones
        clones_dict: Dict[CloneId, Clone] = {
            (): self.clone_factory.create_clone(clone_id=(),clone_type="wild_type",N=self.config.N0),
            (0,): self.clone_factory.create_clone(clone_id=(0,),clone_type="mutated",N=self.config.N_mutant),
            (-1,): self.clone_factory.create_clone(clone_id=(-1,),clone_type="immune",N=self.config.N_immune),
            (-2,): self.clone_factory.create_clone(clone_id=(-2,),clone_type="exhausted",N=self.config.N_exhausted)
        }
        
        # Encapsulate tissue state
        self.tissue_state: TissueState = TissueState(clones=clones_dict)

        self.times: List[float] = [0.0]

        self.crowding_strategy: CrowdingStrategy = (
            AdaptedCrowding(config)
            if config.use_logistic_adapted
            else SimpleCrowding(config)
        )
        
        #aqui habria que anyadir lo mismo para elegir strategy pero para el tipo de leap. (Binomial, Poisson, Poisson half etc)

        self.history: List[Dict[CloneId, dict]] = [self.tissue_state.snapshot()]

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _advance_all_instability(self, dt: float) -> None:
        for clone in self.tissue_state.clones.values():
            if clone.is_alive():
                clone.advance_instability(dt, self.config.base_instability_buildup)

    def _build_rate_matrix(self) -> RateMatrix:
        total_N = self.total_population()
        rate_matrix = RateMatrix()

        for cid, clone in self.tissue_state.clones.items():
            if not clone.is_alive():
                continue
                
            crowding_value = self.crowding_strategy.crowding(clone, self.t, total_N)
            
            # Pass tissue_state to clones so they can access population counts
            rb = clone.birth_rate_effective(tissue_state=self.tissue_state, crowding=crowding_value)
            rd = clone.death_rate_effective(tissue_state=self.tissue_state)
            rm = clone.mutation_rate_effective()
            re = clone.exhaustion_rate_effective(tissue_state=self.tissue_state)
            
            if rb > 0:
                rate_matrix.add_event(Event(EventType.BIRTH, cid, rb))
            if rd > 0:
                rate_matrix.add_event(Event(EventType.DEATH, cid, rd))
            if rm > 0:
                rate_matrix.add_event(Event(EventType.MUTATION, cid, rm))
            if re > 0:
                rate_matrix.add_event(Event(EventType.EXHAUSTION,cid,re))
        return rate_matrix

    def _sample_waiting_time(self, total_rate: float) -> float:
        if total_rate <= 0:
            raise ValueError("Total rate must be positive.")
        return -np.log(self.rng.random()) / total_rate

    def _introduce_mutation(self, clone: Clone) -> None:
        assert clone.is_alive(), "Cannot mutate a dead clone."
        clone.kill()
        child = self.clone_factory.create_clone(
            clone_id=clone.next_child_id(),
            clone_type="mutated",
            N=1,
            parent=clone.clone_id,
        )
        self.tissue_state.clones[child.clone_id] = child
    
    def _induce_exhaustion(self, clone: Clone) -> None:
        self.tissue_state.clones[(-1,)].kill()
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

    # ── Public API ────────────────────────────────────────────────────────────

    def total_population(self) -> int:
        return self.tissue_state.total_population()

    def step(self) -> bool:
        """Advance by one Gillespie step. Returns False when simulation should stop."""
        rate_matrix = self._build_rate_matrix()
        total_rate = rate_matrix.get_total_rate()

        if total_rate <= 0 or not rate_matrix.events:
            return False

        tau = self._sample_waiting_time(total_rate)
        new_t = self.t + tau

        if new_t > self.config.T_max:
            tau = self.config.T_max - self.t
            self._advance_all_instability(tau)
            self.t = self.config.T_max
            self.times.append(self.t)
            self.history.append(self.tissue_state.snapshot())
            return False

        self._advance_all_instability(tau)
        self.t = new_t
        event = rate_matrix.choose_event(self.rng.random())
        self._apply_event(event)

        self.times.append(self.t)
        self.history.append(self.tissue_state.snapshot())
        return True

    def run(self) -> Tuple[List[float], List[Dict[CloneId, dict]], TissueState]:
        while self.t < self.config.T_max and self.total_population() > 0:
            if not self.step():
                break
        return self.times, self.history, self.tissue_state
