from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone
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
        
        
        # AHORA MISMO EL SIM SOLO TIENE HARDCODED INIT CLONE COMO WT, deberia incluir posibilidad de elegir distribuciones de initial conditions
        self.clones: Dict[CloneId, Clone] = {
            (): self.clone_factory.create_clone(clone_id=(),clone_type="wild_type",N=self.config.N0),
            (0,): self.clone_factory.create_clone(clone_id=(0,),clone_type="mutated",N=self.config.N_mutant),
            (-1,): self.clone_factory.create_clone(clone_id=(-1,),clone_type="immune",N=self.config.N_immune),
            (-2,): self.clone_factory.create_clone(clone_id=(-2,),clone_type="exhausted",N=self.config.N_exhausted)
        }

        self.times: List[float] = [0.0]

        self.crowding_strategy: CrowdingStrategy = (
            AdaptedCrowding(config)
            if config.use_logistic_adapted
            else SimpleCrowding(config)
        )
        
        #aqui habria que anyadir lo mismo para elegir strategy pero para el tipo de leap. (Binomial, Poisson, Poisson half etc)

        self.history: List[Dict[CloneId, dict]] = [self._snapshot()]

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _snapshot(self) -> Dict[CloneId, dict]:
        total_N = self.total_population()
        return {
            cid: {
                "N": clone.N,
                "rb": clone.birth_rate_effective(
                    crowding=self.crowding_strategy.crowding(clone, self.t, total_N)
                ),
                "rd": clone.death_rate_effective(),
            }
            for cid, clone in self.clones.items()
        }

    def _advance_all_instability(self, dt: float) -> None:
        for clone in self.clones.values():
            if clone.is_alive():
                clone.advance_instability(dt, self.config.base_instability_buildup)

    def _build_rate_matrix(self) -> RateMatrix:
        total_N = self.total_population()
        rate_matrix = RateMatrix()

        for cid, clone in self.clones.items():
            if not clone.is_alive():
                continue
            crowding_value = self.crowding_strategy.crowding(clone, self.t, total_N)
            

            #Necesito una forma de poder pasar a los calculadores de rb rd rm y re las poblaciones. lo mas sencillo seugramente sea pasarle el diccionario de poblaciones N_C N_I N_W y N_E y luego dentro de clone.py pasarle todo el diccionario y que el elija los valores que necesita. la logica del crowding deberia pasarse a dentro del clone type Class en clone.py y cambiar el valor de total_N para que sea la suma de valores que necesita cada tipo.(por ejemplo el numerador del crowding factor para las wildtyp seria N_W+N_C). entonces el crowding_value se calcula dentro de cada Clase directamente. 

            #Para hacer eso necesito poder agrupar los numeros de poblaciones por el tipo de clon y no por el ID (ahora history va por ID)
            #ESTO HAY QUE SACARLO 
            #-----------
            rb = clone.birth_rate_effective(crowding=crowding_value)
            rd = clone.death_rate_effective()
            rm = clone.mutation_rate_effective()
            re = clone.exhaustion_rate_effective()
            #-----------
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

    # TO-DO: Add tests for mutation and exhaustion events, and for the overall step logic.
    def _introduce_mutation(self, clone: Clone) -> None:
        assert clone.is_alive(), "Cannot mutate a dead clone."
        clone.kill()
        child = self.clone_factory.create_clone(
            clone_id=clone.next_child_id(),
            clone_type="mutated",
            N=1,
            parent=clone.clone_id,
        )
        self.clones[child.clone_id] = child
    def _induce_exhaustion(self,clone:Clone) -> None:
        self.clones[(-1,)].kill()
        self.clones[(-2,)].divide()
    def _apply_event(self, event: Event) -> None:
        clone = self.clones[event.clone_id]
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
        return sum(clone.N for clone in self.clones.values())

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
            self.history.append(self._snapshot())
            return False

        self._advance_all_instability(tau)
        self.t = new_t
        event = rate_matrix.choose_event(self.rng.random())
        self._apply_event(event)

        self.times.append(self.t)
        self.history.append(self._snapshot())
        return True

    def run(self) -> Tuple[List[float], List[Dict[CloneId, dict]], Dict[CloneId, Clone]]:
        while self.t < self.config.T_max and self.total_population() > 0:
            if not self.step():
                break
        return self.times, self.history, self.clones
