from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List
from tissue_state import TissueState
from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone

#TODO: HAY QUE REHACER ESTO PARA QUE HAGA MAS USO DE TISSUE STATE Y SU FUNCION SNAPSHOT. Ver que info queremos sacar al csv
def clone_id_to_str(clone_id: CloneId) -> str:
    if len(clone_id) == 0:
        return "root"   
    elif clone_id == (-3,):
        return "mutated_root"
    elif clone_id == (-1,):
        return "immune"
    elif clone_id == (-2,):
        return "exhausted"

    else:
        return ".".join(map(str, clone_id))
    
#TODO: Fix this so that it passes the TissueState Object from which data is pulled
#TODO: Instead of saving the whole history in one go at the end of the simulation we should write it on the go. (makes more sense for big/long simulations i think)
def save_history_csv(path: Path, times: List[float], history: List[Dict[CloneId, dict]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time", "clone_id", "N", "rb", "rd"])
        for t, snap in zip(times, history):
            for cid, values in snap.items():
                writer.writerow([t, clone_id_to_str(cid), values["N"], values["rb"], values["rd"]])


def save_clones_csv(path: Path, clones: Dict[CloneId, Clone]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "clone_id", "parent", "N",
            "birth_rate", "death_rate", "mutation_rate",
            "instability", "buildup", "d1", "d2", "children_count",
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
