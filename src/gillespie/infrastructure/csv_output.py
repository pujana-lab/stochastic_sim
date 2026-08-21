from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone
from src.gillespie.event import Event

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
    
def _normalize_parquet_path(path: Path) -> Path:
    if path.suffix.lower() == ".parquet":
        return path
    return path.with_suffix(".parquet")


#TODO: Fix this so that it passes the TissueState Object from which data is pulled
#TODO: Instead of saving the whole history in one go at the end of the simulation we should write it on the go. (makes more sense for big/long simulations i think)
def save_history_parquet(path: Path, times: List[float], history: List[Dict[CloneId, dict]]) -> Path:
    out_path = _normalize_parquet_path(path)
    rows: list[dict] = []
    for t, snap in zip(times, history):
        for cid, values in snap.items():
            rows.append(
                {
                    "time": t,
                    "type": values["Type"],
                    "clone_id": clone_id_to_str(cid),
                    "N": values["N"],
                    "rb": values["rb"],
                    "rd": values["rd"],
                }
            )

    df = pd.DataFrame(rows, columns=["time", "type", "clone_id", "N", "rb", "rd"])
    df.to_parquet(out_path, index=False)
    return out_path


def save_clones_parquet(path: Path, clones: Dict[CloneId, Clone]) -> Path:
    out_path = _normalize_parquet_path(path)
    rows: list[dict] = []
    for cid, clone in sorted(clones.items(), key=lambda x: (len(x[0]), x[0])):
        rows.append(
            {
                "clone_id": clone_id_to_str(cid),
                "parent": "" if clone.parent is None else clone_id_to_str(clone.parent),
                "N": clone.N,
                "birth_rate": clone.birth_rate,
                "death_rate": clone.death_rate,
                "mutation_rate": clone.mutation_rate,
                "instability": clone.instability,
                "buildup": clone.buildup,
                "d1": clone.d1,
                "d2": clone.d2,
                "children_count": clone.children_count,
            }
        )

    df = pd.DataFrame(
        rows,
        columns=[
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
        ],
    )
    df.to_parquet(out_path, index=False)
    return out_path


def save_debug_history_parquet(
    path: Path,
    times: List[float],
    history: List[Dict[CloneId, dict]],
    events: List[Optional[Event]],
) -> Path:
    out_path = _normalize_parquet_path(path)
    rows: list[dict] = []
    for t, snap, evt in zip(times, history, events):
        event_kind = evt.kind.value if evt is not None else "none"
        event_clone_type = evt.clone_type if evt is not None else "none"
        for cid, values in snap.items():
            rows.append(
                {
                    "time": t,
                    "clone_id": clone_id_to_str(cid),
                    "Type": values.get("Type", ""),
                    "N": values["N"],
                    "rb": values["rb"],
                    "rd": values["rd"],
                    "rm": values.get("rm", ""),
                    "re": values.get("re", ""),
                    "instability": values.get("instability", ""),
                    "buildup": values.get("buildup", ""),
                    "event_kind": event_kind,
                    "event_clone_type": event_clone_type,
                }
            )

    df = pd.DataFrame(
        rows,
        columns=[
            "time",
            "clone_id",
            "Type",
            "N",
            "rb",
            "rd",
            "rm",
            "re",
            "instability",
            "buildup",
            "event_kind",
            "event_clone_type",
        ],
    )
    df.to_parquet(out_path, index=False)
    return out_path


def save_rates_history_parquet(path: Path, rates_history: List[List[Dict]]) -> Path:
    """
    Guarda el histórico de tasas asumiendo que el tiempo 't' ya está en los datos.
    """
    out_path = _normalize_parquet_path(path)
    rows: list[dict] = []
    for step_events in rates_history:
        for event in step_events:
            rows.append(
                {
                    "time": event["time"],
                    "kind": event["kind"],
                    "clone_id": clone_id_to_str(event["clone_id"]),
                    "clone_type": event["clone_type"],
                    "rate": event["rate"],
                }
            )

    df = pd.DataFrame(rows, columns=["time", "kind", "clone_id", "clone_type", "rate"])
    df.to_parquet(out_path, index=False)
    return out_path


# Backward compatibility aliases during transition.
save_history_csv = save_history_parquet
save_clones_csv = save_clones_parquet
save_debug_history_csv = save_debug_history_parquet
save_rates_history_csv = save_rates_history_parquet