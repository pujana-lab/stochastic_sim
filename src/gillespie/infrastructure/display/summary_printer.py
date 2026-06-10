from __future__ import annotations

import sys
from typing import Dict, List, TextIO

from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone
from src.gillespie.infrastructure.csv_output import clone_id_to_str


def print_summary(
    times: List[float],
    clones: Dict[CloneId, Clone],
    top_k: int,
    file: TextIO = sys.stdout,
) -> None:
    living = [c for c in clones.values() if c.N > 0]
    living_sorted = sorted(living, key=lambda c: c.N, reverse=True)

    print(f"final_time: {times[-1]:.6f}", file=file)
    print(f"total_clones_created: {len(clones)}", file=file)
    print(f"living_clones: {len(living)}", file=file)
    print(f"final_population: {sum(c.N for c in clones.values())}", file=file)

    if living_sorted:
        biggest = living_sorted[0]
        print(f"largest_clone: {clone_id_to_str(biggest.clone_id)}", file=file)
        print(f"largest_clone_size: {biggest.N}", file=file)
        print(f"largest_clone_instability: {biggest.instability:.6f}", file=file)

    print("", file=file)
    for clone in living_sorted[:top_k]:
        print(
            f"  {clone_id_to_str(clone.clone_id):<12} "
            f"N={clone.N:<8d} "
            f"lambda={clone.birth_rate:.4f} "
            f"mu={clone.death_rate:.4f} "
            f"nu={clone.mutation_rate:.4f} "
            f"instability={clone.instability:.6f} "
            f"buildup={clone.buildup:.6f}",
            file=file,
        )
