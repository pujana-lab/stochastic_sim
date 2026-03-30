from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import numpy as np

CloneId = Tuple[int, ...]


@dataclass
class SimulationConfig:
    N0: int = 100
    lambda0: float = 0.30
    mu0: float = 0.20
    nu0: float = 0.01

    d1_0: float = 0.0
    d2_0: float = 0.0

    instability_0: float = 0.0
    buildup_0: float = 0.0001                    # clone-specific instability buildup
    base_instability_buildup: float = 0.0        # background buildup shared by all clones

    mutation_instability_jump: float = 0.05      # jump in instability at mutation
    mutation_buildup_gain: float = 0.0001        # mutation increases future buildup
    mutation_rate_gain: float = 0.02             # mutation increases nu
    fitness_gain: float = 0.05                   # mutation increases lambda

    T_max: float = 50.0
    seed: Optional[int] = None

    use_logistic: bool = False
    K0: float = 100_000_000
    decline: float = 0.0
    Kmin: float = 1.0


@dataclass
class Clone:
    clone_id: CloneId
    N: int
    birth_rate: float
    death_rate: float
    mutation_rate: float

    # Stored as value at last_update_t
    instability: float = 0.0
    buildup: float = 0.0
    last_update_t: float = 0.0

    d1: float = 0.0
    d2: float = 0.0

    parent: Optional[CloneId] = None
    children_count: int = 0

    def is_alive(self) -> bool:
        return self.N > 0

    def divide(self) -> None:
        self.N += 1

    def die(self) -> None:
        if self.N > 0:
            self.N -= 1

    def next_child_id(self) -> CloneId:
        self.children_count += 1
        return self.clone_id + (self.children_count,)

    # ---------- Lazy time-dependent state ----------

    def instability_at(self, t: float, base_buildup: float) -> float:
        dt = t - self.last_update_t
        if dt < 0:
            raise ValueError("Cannot evaluate instability at a time earlier than last_update_t.")
        return self.instability + (base_buildup + self.buildup) * dt

    def materialize_to(self, t: float, base_buildup: float) -> None:
        self.instability = self.instability_at(t, base_buildup)
        self.last_update_t = t

    # ---------- Derived quantities / hazards ----------

    def mutation_multiplier(self, t: float, base_buildup: float) -> float:
        return 1.0 + self.instability_at(t, base_buildup)

    def birth_hazard(self, t: float, crowding: float, base_buildup: float) -> float:
        _ = t, base_buildup
        return self.birth_rate * self.N * crowding

    def death_hazard(self, t: float, base_buildup: float) -> float:
        _ = t, base_buildup
        return self.death_rate * self.N

    def mutation_hazard(self, t: float, base_buildup: float) -> float:
        return self.mutation_rate * self.N * self.mutation_multiplier(t, base_buildup)

    # ---------- Mutation ----------

    def mutated_child(
        self,
        current_time: float,
        base_buildup: float,
        fitness_gain: float,
        mutation_rate_gain: float,
        instability_jump: float,
        buildup_gain: float,
    ) -> "Clone":
        """
        One cell leaves the parent and becomes a new child clone.

        Parent is first materialized to current_time so the child inherits the
        parent's current instability, then receives an instability jump and a
        larger future buildup rate.
        """
        if self.N <= 0:
            raise ValueError("Cannot mutate a dead clone.")

        self.materialize_to(current_time, base_buildup)

        self.N -= 1
        child_id = self.next_child_id()

        return Clone(
            clone_id=child_id,
            N=1,
            birth_rate=self.birth_rate * (1.0 + fitness_gain),
            death_rate=self.death_rate,
            mutation_rate=self.mutation_rate + mutation_rate_gain,
            instability=self.instability + instability_jump,
            buildup=self.buildup + buildup_gain,
            last_update_t=current_time,
            d1=self.d1,
            d2=self.d2,
            parent=self.clone_id,
        )


@dataclass
class Event:
    kind: str      # "birth", "death", "mutation"
    clone_id: CloneId
    rate: float


class TumorSimulation:
    def __init__(self, config: SimulationConfig):
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
                last_update_t=0.0,
                d1=float(config.d1_0),
                d2=float(config.d2_0),
                parent=None,
            )
        }

        self.times: List[float] = [0.0]
        self.history: List[Dict[CloneId, int]] = [self.snapshot()]

    # ---------- Basic helpers ----------

    def snapshot(self) -> Dict[CloneId, int]:
        return {cid: clone.N for cid, clone in self.clones.items()}

    def total_population(self) -> int:
        return sum(clone.N for clone in self.clones.values())

    def carrying_capacity(self, t: float) -> float:
        cfg = self.config
        return max(cfg.Kmin, cfg.K0 - cfg.decline * t)

    def crowding_factor(self, t: float) -> float:
        if not self.config.use_logistic:
            return 1.0

        Kt = self.carrying_capacity(t)
        if Kt <= 0:
            return 0.0

        return max(0.0, 1.0 - self.total_population() / Kt)

    def materialize_all(self, t: float) -> None:
        base = self.config.base_instability_buildup
        for clone in self.clones.values():
            if clone.is_alive():
                clone.materialize_to(t, base)

    # ---------- Gillespie pieces ----------

    def build_event_table(self) -> Tuple[float, List[Event]]:
        g = self.crowding_factor(self.t)
        base = self.config.base_instability_buildup
        events: List[Event] = []

        for cid, clone in self.clones.items():
            if not clone.is_alive():
                continue

            rb = clone.birth_hazard(self.t, g, base)
            rd = clone.death_hazard(self.t, base)
            rm = clone.mutation_hazard(self.t, base)

            if rb > 0:
                events.append(Event("birth", cid, rb))
            if rd > 0:
                events.append(Event("death", cid, rd))
            if rm > 0:
                events.append(Event("mutation", cid, rm))

        total_rate = sum(event.rate for event in events)
        return total_rate, events

    def sample_waiting_time(self, total_rate: float) -> float:
        if total_rate <= 0:
            raise ValueError("Total rate must be positive.")
        u = self.rng.random()
        return -np.log(u) / total_rate

    def choose_event(self, events: List[Event], total_rate: float) -> Event:
        if total_rate <= 0:
            raise ValueError("Total rate must be positive.")

        rates = np.array([event.rate for event in events], dtype=float)
        cumulative = np.cumsum(rates / total_rate)
        u = self.rng.random()
        idx = np.searchsorted(cumulative, u, side="right")
        if idx >= len(events):
            idx = len(events) - 1
        return events[idx]

    def apply_event(self, event: Event) -> None:
        clone = self.clones[event.clone_id]
        base = self.config.base_instability_buildup

        # Materialize affected clone to the exact event time before changing it
        clone.materialize_to(self.t, base)

        if event.kind == "birth":
            clone.divide()

        elif event.kind == "death":
            clone.die()

        elif event.kind == "mutation":
            if clone.is_alive():
                child = clone.mutated_child(
                    current_time=self.t,
                    base_buildup=base,
                    fitness_gain=self.config.fitness_gain,
                    mutation_rate_gain=self.config.mutation_rate_gain,
                    instability_jump=self.config.mutation_instability_jump,
                    buildup_gain=self.config.mutation_buildup_gain,
                )
                self.clones[child.clone_id] = child

        else:
            raise ValueError(f"Unknown event kind: {event.kind}")

    def step(self) -> bool:
        """
        Advance the simulation by one Gillespie step.
        Returns False if no more events can occur or if the next event would
        occur after T_max.
        """
        total_rate, events = self.build_event_table()
        if total_rate <= 0 or not events:
            return False

        tau = self.sample_waiting_time(total_rate)
        new_t = self.t + tau

        # If the next event would occur after T_max, stop exactly at T_max
        # and materialize continuous-time state there.
        if new_t > self.config.T_max:
            self.t = self.config.T_max
            self.materialize_all(self.t)
            self.times.append(self.t)
            self.history.append(self.snapshot())
            return False

        self.t = new_t
        event = self.choose_event(events, total_rate)
        self.apply_event(event)

        self.times.append(self.t)
        self.history.append(self.snapshot())
        return True

    def run(self) -> Tuple[List[float], List[Dict[CloneId, int]], Dict[CloneId, Clone]]:
        while self.t < self.config.T_max and self.total_population() > 0:
            if not self.step():
                break

        # Make sure returned clone objects reflect the final simulation time
        self.materialize_all(self.t)
        return self.times, self.history, self.clones
    
