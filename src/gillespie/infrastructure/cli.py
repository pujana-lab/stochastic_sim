from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.infrastructure.csv_output import clone_id_to_str, save_history_csv, save_clones_csv
from src.gillespie.tumor_simulation import TumorSimulation

# TODO: actualizar con nuevos valores del SimulationConfig
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run tumor clone Gillespie simulation.")

    p.add_argument("--N0", type=int, default=10)
    p.add_argument("--lambda0", type=float, default=0.20)
    p.add_argument("--mu0", type=float, default=0.20)
    p.add_argument("--nu0", type=float, default=0.00)

    p.add_argument("--d1-0", dest="d1_0", type=float, default=0.0)
    p.add_argument("--d2-0", dest="d2_0", type=float, default=0.0)

    p.add_argument("--instability-0", dest="instability_0", type=float, default=0.0)
    p.add_argument("--buildup-0", dest="buildup_0", type=float, default=0.0)
    p.add_argument("--base-instability-buildup", type=float, default=0.0)

    p.add_argument("--mutation-instability-jump", type=float, default=0.0)
    p.add_argument("--mutation-buildup-gain", type=float, default=0.0)

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

#TODO: mensaje de arranque y progreso
#TODO: anyadir trazabilidad
#TODO: ver que va guardando en memoria. (explota)
def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    config = SimulationConfig()
    sim = TumorSimulation(config=config)
    times, history, tissue_state = sim.run()

    print_summary(times, tissue_state.clones, args.top)

    if args.save_history is not None:
        save_history_csv(args.save_history, times, history)
        print(f"\nSaved history to {args.save_history}")

    if args.save_clones is not None:
        save_clones_csv(args.save_clones, tissue_state.clones)
        print(f"Saved clones to {args.save_clones}")


if __name__ == "__main__":
    main()
