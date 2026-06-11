from __future__ import annotations

from typing import Dict, List, Optional, Tuple

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
        print("arrancando con parametros:")
        print(config)
        
        #aqui igual merece mas la pena guardarlo como array simplemente o eso o hacerlo por tipos pero en ese caso 
        # Initialize clones
        # clones_dict:Dict[CloneType,Clone]= {}
        # TENGO QUE MOVER EL TIPO DE CLON DE LA FACTORY A LA CLASE CLONE
        clones_dict: Dict[CloneId, Clone] = {
            (): self.clone_factory.create_clone(clone_id=(),clone_type="base",N=self.config.N0),
            (-3,): self.clone_factory.create_clone(clone_id=(-3,),clone_type="mutated",N=self.config.N_mutant),
            # (-4,4,4,): self.clone_factory.create_clone(clone_id=(-4,4,4),clone_type="mutated_test",N=0),
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

    def create_clone(
        self,
        clone_id: CloneId,
        clone_type: str = "base",
        N: int = 1,
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

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _advance_all_instability(self, dt: float) -> None:
        for clone in self.tissue_state.clones.values():
            if clone.is_alive():
                clone.advance_instability(dt, self.config.base_instability_buildup)

    def _build_rate_matrix(self) -> RateMatrix:
        self.tissue_state.update_pop_map()
        rate_matrix = RateMatrix()
        type_rates: Dict[str, tuple[float, float, float, float]] = {}

        for cid, clone in self.tissue_state.clones.items():
            if not clone.is_alive():
                continue

            clone_type = clone.get_type()
            if clone_type not in type_rates:
                crowding_value = self.crowding_strategy.crowding(
                    clone=clone, t=self.t, tissue_state=self.tissue_state
                )
                type_rates[clone_type] = (
                    clone.birth_rate_effective(tissue_state=self.tissue_state, crowding=crowding_value),
                    clone.death_rate_effective(tissue_state=self.tissue_state),
                    clone.mutation_rate_effective(tissue_state=self.tissue_state),
                    clone.exhaustion_rate_effective(tissue_state=self.tissue_state),
                )

            rb, rd, rm, re = type_rates[clone_type]

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

    # TO-DO: Add tests for mutation and exhaustion events, and for the overall step logic.
    def _introduce_mutation(self, clone: Clone) -> None:
        assert clone.is_alive(), "Cannot mutate a dead clone."
        
        # Validate that the next mutation type exists in CloneType
        assert clone.next_mutation is not "", "Clone is not supposed to be mutating"
        
        
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
        self.tissue_state.update_pop_map()

        self.times.append(self.t)
        self.history.append(self.tissue_state.snapshot())
        if self.config.verbose:
            if event.clone_type == "mutant" or event.clone_type == "immune":
                    print("-------------------")
                    print(f"time:{self.t}" )
                    print(f"EVENT: {event.kind.value}")
                    print(f"KIND: {event.clone_type} ")
                    self.tissue_state.print_pop_map()
                    print("-------------------")
                
        return True
    
    def _stopping_cond(self) -> bool:
        """" We will make the system stop if it reaches stability"""

     

        # 2. Extinción total
        if self.total_population() <= 0:
            
            print("total extinction")
            
            return True
                
        
        # if self.tissue_state.pop_map["immune"] <= 0:

        #     print(" Complete tumor escape  ")
        #     return True
        
        # if self.tissue_state.pop_map["mutated"] <=0 and self.tissue_state.pop_map["exhausted"]<=0:
        #     print (" complete tumor control")
        #     return True
            
        
        #
       
        return False

    def run(self) -> Tuple[List[float], List[Dict[CloneId, dict]], TissueState]:
        while not self._stopping_cond():
            if not self.step():
                break
        return self.times, self.history, self.tissue_state
