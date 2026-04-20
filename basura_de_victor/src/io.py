from pathlib import Path
from typing import Dict, List
from src.cloneId import CloneId
from src.clone import Clone
import csv

import argparse


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
