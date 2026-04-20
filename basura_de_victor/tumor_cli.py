from __future__ import annotations

import argparse
from typing import Dict, List, Tuple

import numpy as np

from src.cloneId import CloneId
from src.event_type import EventType
from src.event import Event
from src.rate_matrix import RateMatrix
from src.clone import Clone
from src.crowding_strategy import CrowdingStrategy, SimpleCrowding, AdaptedCrowding
from src.simulation_config import SimulationConfig
from src.io import save_history_csv, save_clones_csv, print_summary, build_parser


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

        self.crowding_strategy: CrowdingStrategy = (
            AdaptedCrowding(config)
            if config.use_logistic_adapted
            else SimpleCrowding(config)
        )

        self.history: List[Dict[CloneId, int]] = [self.snapshot()]

    def advance_all_instability(self, dt: float) -> None:
        for clone in self.clones.values():
            if clone.is_alive():
                clone.advance_instability(dt, self.config.base_instability_buildup)

    def snapshot(self)-> Dict[CloneId,dict]:
        total_N = self.total_population()
        return {
            cid:{
                 "N": clone.N,
                 "rb": clone.birth_rate_effective(
                    crowding=self.crowding_strategy.crowding(clone, self.t, total_N)
                ),
                 "rd": clone.death_rate_effective()
                 }
            for cid, clone in self.clones.items()
        }
            
 
    def total_population(self) -> int:
        return sum(clone.N for clone in self.clones.values())

    def build_rate_matrix(self) -> RateMatrix:
        total_N = self.total_population()
        rate_matrix=RateMatrix()

        for cid, clone in self.clones.items():
            if not clone.is_alive():
                continue
            crowding_value = self.crowding_strategy.crowding(clone, self.t, total_N)

            rb = clone.birth_rate_effective(crowding = crowding_value)
            rd = clone.death_rate_effective()
            rm = clone.mutation_rate_effective()

            if rb > 0:
                rate_matrix.add_event(Event(EventType.BIRTH, cid, rb))
            if rd > 0:
                rate_matrix.add_event(Event(EventType.DEATH, cid, rd))
            if rm > 0:
                rate_matrix.add_event(Event(EventType.MUTATION, cid, rm))
        return rate_matrix


    def sample_waiting_time(self, total_rate: float) -> float:
        if total_rate <= 0:
            raise ValueError("Total rate must be positive.")
        u = self.rng.random()
        return -np.log(u) / total_rate

    def choose_event(self, rate_matrix:RateMatrix) -> Event:
        u = self.rng.random()
        return rate_matrix.choose_event(u)

    def introduce_mutation(self, clone: Clone) -> None:
        assert clone.is_alive(), "Cannot mutate a dead clone."
        
        child = clone.mutate(
            fitness_gain=self.config.fitness_gain,
            instability_jump=self.config.mutation_instability_jump,
            buildup_gain=self.config.mutation_buildup_gain
        )
        self.clones[child.clone_id] = child

    def apply_event(self, event: Event) -> None:
        clone = self.clones[event.clone_id]

        if event.kind == EventType.BIRTH:
            clone.divide()

        elif event.kind == EventType.DEATH:
            clone.kill()

        elif event.kind == EventType.MUTATION:
            self.introduce_mutation(clone)

        else:
            raise ValueError(f"Unknown event kind: {event.kind}")

    def step(self) -> bool:
        rate_matrix = self.build_rate_matrix()
        total_rate = rate_matrix.get_total_rate()

        if total_rate <= 0 or not rate_matrix.events:
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
        event = self.choose_event(rate_matrix)
        self.apply_event(event)

        self.times.append(self.t)
        self.history.append(self.snapshot())
        return True

    def run(self) -> Tuple[List[float], List[Dict[CloneId, int]], Dict[CloneId, Clone]]:
        while self.t < self.config.T_max and self.total_population() > 0:
            if not self.step():
                break
        return self.times, self.history, self.clones





def config_from_args(args: argparse.Namespace) -> SimulationConfig:
    return SimulationConfig(
        N0=args.N0,
        lambda0=args.lambda0,
        mu0=args.mu0,
        nu0=args.nu0,
        d1_0=args.d1_0,
        d2_0=args.d2_0,
        instability_0=args.instability_0,
        buildup_0=args.buildup_0,
        base_instability_buildup=args.base_instability_buildup,
        mutation_instability_jump=args.mutation_instability_jump,
        mutation_buildup_gain=args.mutation_buildup_gain,
        T_max=args.T_max,
        seed=args.seed,
        use_logistic=args.use_logistic,
        use_logistic_adapted=args.use_logistic_adapted,
        K0=args.K0,
        decline=args.decline,
        Kmin=args.Kmin,
        fitness_gain=args.fitness_gain,
    )




def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    config = config_from_args(args)
    sim = TumorSimulation(config=config)
    times, history, clones = sim.run()

    print_summary(times, clones, args.top)

    if args.save_history is not None:
        save_history_csv(args.save_history, times, history)
        print(f"\nSaved history to {args.save_history}")

    if args.save_clones is not None:
        save_clones_csv(args.save_clones, clones)
        print(f"Saved clones to {args.save_clones}")


if __name__ == "__main__":
    main()