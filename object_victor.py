from __future__ import annotations

from dataclasses import dataclass, field
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
    buildup_0: float = 0.0001   # clone-specific buildup rate
    base_instability_buildup: float = 0.00005  # global background buildup

    mutation_instability_jump: float = 0.05    # jump after mutation
    mutation_buildup_gain: float = 0.0001      # mutations make future buildup faster

    T_max: float = 50.0
    seed: Optional[int] = None

    use_logistic: bool = False
    K0: float = 100_000_000
    decline: float = 0.0
    Kmin: float = 1.0

    fitness_gain: float = 0.05

@dataclass
class Clone:
    clone_id: CloneId
    N: int
    birth_rate: float       # lambda
    death_rate: float       # mu
    mutation_rate: float    # nu
    d1: float = 0.0
    d2: float = 0.0
    instability: float =0.0 
    buildup: float =0.0
    parent: Optional[CloneId] = None
    children_count: int = 0

    def is_alive(self) -> bool:
        return self.N > 0
    
    def mutation_multiplier(self) -> float:
        return 1.0 + self.instability

    def birth_prob(self, crowding: float) -> float:
        return self.birth_rate * self.N * crowding

    def death_prob(self) -> float:
        return self.death_rate * self.N

    def mutation_prob(self) -> float:
        return self.mutation_rate * self.N * self.mutation_multiplier()
    def divide(self) -> None:
        self.N += 1

    def die(self) -> None:
        if self.N > 0:
            self.N -= 1

    def next_child_id(self) -> CloneId:
        self.children_count += 1
        return self.clone_id + (self.children_count,)
    def advance_instability(self, dt: float, base_buildup: float) -> None:
        # continuous buildup over time
        self.instability += (base_buildup + self.buildup) * dt
    def mutated_child(
        self,
        fitness_gain: float,
        instability_jump: float,
        buildup_gain: float,
    ) -> "Clone":
        """
        One cell leaves the parent and becomes a new child clone.
        The child is more unstable and has faster future buildup.
        """
        if self.N <= 0:
            raise ValueError("Cannot mutate a dead clone.")

        self.N -= 1
        child_id = self.next_child_id()

        return Clone(
            clone_id=child_id,
            N=1,
            birth_rate=self.birth_rate * (1.0 + fitness_gain),
            death_rate=self.death_rate,
            mutation_rate=self.mutation_rate,
            d1=self.d1,
            d2=self.d2,
            parent=self.clone_id,
            instability= self.instability + instability_jump,
            buildup= self.buildup + buildup_gain,
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
                d1=float(config.d1_0),
                d2=float(config.d2_0),
                parent=None,
            )
        }

        self.times: List[float] = [0.0]
        self.history: List[Dict[CloneId, int]] = [self.snapshot()]
    def advance_all_instability(self, dt: float) -> None:
        for clone in self.clones.values():
            if clone.is_alive():
                clone.advance_instability(dt, self.config.base_instability_buildup)
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

    def build_event_table(self) -> Tuple[float, List[Event]]:
        g = self.crowding_factor(self.t)
        events: List[Event] = []

        for cid, clone in self.clones.items():
            if not clone.is_alive():
                continue
            rb = clone.birth_prob(g)
            rd = clone.death_prob()
            rm = clone.mutation_prob()


            if rb > 0:
                events.append(Event("birth", cid, rb))
            if rd > 0:
                events.append(Event("death", cid, rd))
            if rm > 0:
                events.append(Event("mutation", cid, rm))

        R = sum(event.rate for event in events)
        return R, events
    
    def sample_waiting_time(self, total_rate: float) -> float:
        if total_rate <= 0:
            raise ValueError("Total rate must be positive.")
        u = self.rng.random()
        return -np.log(u) / total_rate

    def choose_event(self, events: List[Event], total_rate: float) -> Event:
        rates = np.array([event.rate for event in events], dtype=float)
        cumulative = np.cumsum(rates / total_rate)
        u = self.rng.random()
        idx = np.searchsorted(cumulative, u, side="right")
        return events[idx]

    def apply_event(self, event: Event) -> None:
        clone = self.clones[event.clone_id]

        if event.kind == "birth":
            clone.divide()

        elif event.kind == "death":
            clone.die()

        elif event.kind == "mutation":
            if clone.is_alive():
                child = clone.mutated_child(
                    fitness_gain=self.config.fitness_gain,
                    instability_jump=self.config.mutation_instability_jump,
                    buildup_gain=self.config.mutation_buildup_gain,
                )
                self.clones[child.clone_id] = child

        else:
            raise ValueError(f"Unknown event kind: {event.kind}")

    def step(self) -> bool:
        """
        Advance the simulation by one step.
        Returns False if no more events can occur.
        """
        total_rate, events = self.build_event_table()
        if total_rate <= 0 or not events:
            return False

        tau = self.sample_waiting_time(total_rate)
        new_t = self.t + tau

        if new_t > self.config.T_max:
            tau = self.config.T_max - self.t
            self.advance_all_instability(tau)
            self.t = self.config.T_max
            self.times.append(self.t)
            self.history.append(self.snapshot())
            return False
        self.advance_all_instability(tau)
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

        return self.times, self.history, self.clones
    
