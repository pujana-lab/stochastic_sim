from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

DEFAULT_STAT_NAMES = ["mean", "final"]


def _integrate(y: np.ndarray, x: np.ndarray) -> float:
    return float(np.trapezoid(y, x))


def compute_summary_stats(
    pop_by_type: Dict[str, np.ndarray],
    times: np.ndarray,
    stat_names: Optional[List[str]] = None,
    cell_types: Optional[List[str]] = None,
) -> np.ndarray:
    if stat_names is None:
        stat_names = DEFAULT_STAT_NAMES
    if cell_types is None:
        cell_types = sorted(pop_by_type.keys())

    stats = []
    for ct in cell_types:
        arr = pop_by_type[ct]
        for sname in stat_names:
            if sname == "mean":
                stats.append(float(np.mean(arr)))
            elif sname == "std":
                stats.append(float(np.std(arr, ddof=1)))
            elif sname == "final":
                stats.append(float(arr[-1]))
            elif sname == "integral":
                stats.append(_integrate(arr, times))
            else:
                raise ValueError(f"Unknown summary statistic: {sname}")
    return np.array(stats)
