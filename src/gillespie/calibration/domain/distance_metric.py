from __future__ import annotations

import numpy as np


def weighted_sse(
    sim_stats: np.ndarray,
    ref_stats: np.ndarray,
    precision_weights: np.ndarray,
) -> float:
    if sim_stats.shape != ref_stats.shape:
        raise ValueError(
            f"Shape mismatch: sim {sim_stats.shape} vs ref {ref_stats.shape}"
        )
    if precision_weights.shape != sim_stats.shape:
        raise ValueError(
            f"Weights shape {precision_weights.shape} != stats shape {sim_stats.shape}"
        )
    return float(np.sum(precision_weights * (sim_stats - ref_stats) ** 2))
