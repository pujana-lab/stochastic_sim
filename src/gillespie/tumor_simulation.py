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


class TumorSimulation:
    def __init__(self, config: SimulationConfig) -> None:
        self.config = config
        self.rng = np.random.default_rng(config.seed)
        self.t = 0.0

        self.clones: Dict[CloneId, Clone] = {
            (): Clone(
                clone_id=(),
                N=int(config.N0),
                birth_rate=float(config.lambda0),
                death_rate=float(config.mu0),
                mutation_rate=float(config.nu0),
                instability=float(config.instability_0),
                buildup=float(config.buildup_0),
                d1=float(config.d1_0),
                d2=float(config.d2_0),
                parent=None,
            )
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

            rb = clone.birth_rate_effective(crowding=crowding_value)
            rd = clone.death_rate_effective()
            rm = clone.mutation_rate_effective()
            number = 1
            if rb > 0:
                rate_matrix.add_event(Event(EventType.BIRTH, cid, rb, number))
            if rd > 0:
                rate_matrix.add_event(Event(EventType.DEATH, cid, rd, number))
            if rm > 0:
                rate_matrix.add_event(Event(EventType.MUTATION, cid, rm, number))
        return rate_matrix

    def _sample_waiting_time(self, total_rate: float) -> float:
        if total_rate <= 0:
            raise ValueError("Total rate must be positive.")
        return -np.log(self.rng.random()) / total_rate

    def _introduce_mutation(self, clone: Clone) -> None:
        assert clone.is_alive(), "Cannot mutate a dead clone."
        child = clone.mutate(
            fitness_gain=self.config.fitness_gain,
            instability_jump=self.config.mutation_instability_jump,
            buildup_gain=self.config.mutation_buildup_gain,
        )
        self.clones[child.clone_id] = child

    def _apply_event(self, event: Event) -> None:
        clone = self.clones[event.clone_id]
        if event.kind == EventType.BIRTH:
            clone.divide()
        elif event.kind == EventType.DEATH:
            clone.kill()
        elif event.kind == EventType.MUTATION:
            self._introduce_mutation(clone)
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
