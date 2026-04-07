from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from src.events import CloneId, Event, EventsCollection




@dataclass
class SimulationConfig:
    N0: int = 10
    lambda0: float = 0.50
    mu0: float = 0.20
    nu0: float = 0.00

    d1_0: float = 0.0
    d2_0: float = 0.0

    instability_0: float = 0.0
    buildup_0: float = 0.0000
    base_instability_buildup: float = 0.00000

    mutation_instability_jump: float = 0.05
    mutation_buildup_gain: float = 0.0001

    T_max: float = 5000
    seed: Optional[int] = None

    use_logistic: bool = True
    use_logistic_adapted: bool=True
    K0: float = 100
    decline: float = 0.0
    Kmin: float = 1.0

    fitness_gain: float = 0.05


@dataclass
class Clone:
    clone_id: CloneId
    N: int
    birth_rate: float
    death_rate: float
    mutation_rate: float

    instability: float = 0.0
    buildup: float = 0.0
    d1: float = 0.0
    d2: float = 0.0

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
        self.instability += (base_buildup + self.buildup) * dt

    def mutated_child(
        self,
        fitness_gain: float,
        instability_jump: float,
        buildup_gain: float,
    ) -> "Clone":
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
            instability=self.instability + instability_jump,
            buildup=self.buildup + buildup_gain,
            d1=self.d1,
            d2=self.d2,
            parent=self.clone_id,
        )





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

    # def snapshot(self) -> Dict[CloneId, int]:
    #     return {cid: clone.N for cid, clone in self.clones.items()}
    def snapshot(self)-> Dict[CloneId,dict]:
        return {
            cid:{
                "N": clone.N,
                 "rb": clone.birth_prob(self.crowding_factor_adapted(clone,self.t) if self.config.use_logistic_adapted else self.crowding_factor(self.t)),
                 "rd": clone.death_prob()
                 }
                 for cid, clone in self.clones.items()
                 }
            
 
    def total_population(self) -> int:
        return sum(clone.N for clone in self.clones.values())

    def carrying_capacity(self, t: float) -> float:
        cfg = self.config
        return max(cfg.Kmin, cfg.K0 - cfg.decline * t)
    def carrying_capacity_adapted(self,clone,t:float)->float:
        cfg=self.config
        return max(cfg.Kmin,(cfg.K0/(1-(clone.death_rate/clone.birth_rate))-cfg.decline* t))

    def crowding_factor(self, t: float) -> float:
        if not self.config.use_logistic:
            return 1.0

        Kt = self.carrying_capacity(t)
        if Kt <= 0:
            return 0.0

        return max(0.0, 1.0 - self.total_population() / Kt)
    def crowding_factor_adapted(self,clone,t:float)->float:
        if not self.config.use_logistic:
            return 1.0
        Kt=self.carrying_capacity_adapted(clone,t)
        if Kt <= 0:
            return 0.0
        return max(0.0,1.0-self.total_population()/Kt)

    def build_event_table(self) -> EventsCollection:
        if not self.config.use_logistic_adapted:
          g = self.crowding_factor_adapted(self.t)
        events_collection=EventsCollection()

        for cid, clone in self.clones.items():
            if not clone.is_alive():
                continue
            if self.config.use_logistic_adapted:
                g=self.crowding_factor_adapted(clone,self.t)
            rb = clone.birth_prob(g)
            rd = clone.death_prob()
            rm = clone.mutation_prob()

            if rb > 0:
                events_collection.add_event(Event("birth", cid, rb))
            if rd > 0:
                events_collection.add_event(Event("death", cid, rd))
            if rm > 0:
                events_collection.add_event(Event("mutation", cid, rm))

        return events_collection


    def sample_waiting_time(self, total_rate: float) -> float:
        if total_rate <= 0:
            raise ValueError("Total rate must be positive.")
        u = self.rng.random()
        return -np.log(u) / total_rate

    def choose_event(self, events_collection:EventsCollection) -> Event:
        events = events_collection.events
        
        rates = np.array([event.rate for event in events], dtype=float)
        cumulative = np.cumsum(rates / events_collection.get_total_rate())
        u = self.rng.random()
        idx = np.searchsorted(cumulative, u, side="right")
        if idx >= len(events):
            idx = len(events) - 1
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
        events_collection = self.build_event_table()
        total_rate = events_collection.get_total_rate()
        if total_rate <= 0 or not events_collection.events:
            return False

        tau = self.sample_waiting_time(events_collection.get_total_rate())
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
        event = self.choose_event(events_collection)
        self.apply_event(event)

        self.times.append(self.t)
        self.history.append(self.snapshot())
        return True

    def run(self) -> Tuple[List[float], List[Dict[CloneId, int]], Dict[CloneId, Clone]]:
        while self.t < self.config.T_max and self.total_population() > 0:
            if not self.step():
                break
        return self.times, self.history, self.clones


def clone_id_to_str(clone_id: CloneId) -> str:
    return "root" if len(clone_id) == 0 else ".".join(map(str, clone_id))


# def save_history_csv(path: Path, times: List[float], history: List[Dict[CloneId, int]]) -> None:
#     with path.open("w", newline="") as f:
#         writer = csv.writer(f)
#         writer.writerow(["time", "clone_id", "population"])
#         for t, snap in zip(times, history):
#             for cid, n in snap.items():
#                 writer.writerow([t, clone_id_to_str(cid), n])
def save_history_csv(path: Path, times: List[float], history: List[Dict[CloneId, dict]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time", "clone_id", "N", "rb", "rd"])

        for t, snap in zip(times, history):
            for cid, values in snap.items():
                writer.writerow([
                    t,
                    clone_id_to_str(cid),
                    values["N"],
                    values["rb"],
                    values["rd"],
                ])

def save_clones_csv(path: Path, clones: Dict[CloneId, Clone]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "clone_id",
            "parent",
            "N",
            "birth_rate",
            "death_rate",
            "mutation_rate",
            "instability",
            "buildup",
            "d1",
            "d2",
            "children_count",
        ])
        for cid, clone in sorted(clones.items(), key=lambda x: (len(x[0]), x[0])):
            writer.writerow([
                clone_id_to_str(cid),
                "" if clone.parent is None else clone_id_to_str(clone.parent),
                clone.N,
                clone.birth_rate,
                clone.death_rate,
                clone.mutation_rate,
                clone.instability,
                clone.buildup,
                clone.d1,
                clone.d2,
                clone.children_count,
            ])


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run tumor clone Gillespie simulation.")

    p.add_argument("--N0", type=int, default=10)
    p.add_argument("--lambda0", type=float, default=0.20)
    p.add_argument("--mu0", type=float, default=0.20)
    p.add_argument("--nu0", type=float, default=0.00)

    p.add_argument("--d1-0", dest="d1_0", type=float, default=0.0)
    p.add_argument("--d2-0", dest="d2_0", type=float, default=0.0)

    p.add_argument("--instability-0", dest="instability_0", type=float, default=0.0)
    p.add_argument("--buildup-0", dest="buildup_0", type=float, default=0.0000)
    p.add_argument("--base-instability-buildup", type=float, default=0.00000)

    p.add_argument("--mutation-instability-jump", type=float, default=0.05)
    p.add_argument("--mutation-buildup-gain", type=float, default=0.0001)

    p.add_argument("--T-max", dest="T_max", type=float, default=10.0)
    p.add_argument("--seed", type=int, default=None)

    p.add_argument("--use-logistic", action="store_true")
    p.add_argument("--use-logistic-adapted", action="store_true")
    p.add_argument("--K0", type=float, default=100_000_000)
    p.add_argument("--decline", type=float, default=0.0)
    p.add_argument("--Kmin", type=float, default=1.0)

    p.add_argument("--fitness-gain", type=float, default=0.05)

    p.add_argument("--save-history", type=Path, default=Path("./history.csv"), help="Write long-format history CSV.")
    p.add_argument("--save-clones", type=Path, default=None, help="Write final clone states CSV.")
    p.add_argument("--top", type=int, default=10, help="How many largest final clones to print.")

    return p


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


def print_summary(times: List[float], clones: Dict[CloneId, Clone], top_k: int) -> None:
    living = [c for c in clones.values() if c.N > 0]
    living_sorted = sorted(living, key=lambda c: c.N, reverse=True)

    print(f"final_time: {times[-1]:.6f}")
    print(f"total_clones_created: {len(clones)}")
    print(f"living_clones: {len(living)}")
    print(f"final_population: {sum(c.N for c in clones.values())}")

    if living_sorted:
        biggest = living_sorted[0]
        print(f"largest_clone: {clone_id_to_str(biggest.clone_id)}")
        print(f"largest_clone_size: {biggest.N}")
        print(f"largest_clone_instability: {biggest.instability:.6f}")

    print("\nTop clones:")
    for clone in living_sorted[:top_k]:
        print(
            f"  {clone_id_to_str(clone.clone_id):<12} "
            f"N={clone.N:<8d} "
            f"lambda={clone.birth_rate:.4f} "
            f"mu={clone.death_rate:.4f} "
            f"nu={clone.mutation_rate:.4f} "
            f"instability={clone.instability:.6f} "
            f"buildup={clone.buildup:.6f}"
        )


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    config = config_from_args(args)
    sim = TumorSimulation(config=SimulationConfig)
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