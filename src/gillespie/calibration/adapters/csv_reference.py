from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from src.gillespie.calibration.domain.summary_statistic import (
    DEFAULT_STAT_NAMES,
    compute_summary_stats,
)
from src.gillespie.calibration.ports.calibration_ports import ReferenceData


class CsvReferenceAdapter:
    def __init__(
        self,
        csv_path: str,
        weights_path: Optional[str] = None,
        stat_names: Optional[List[str]] = None,
    ):
        self.csv_path = Path(csv_path)
        self.weights_path = Path(weights_path) if weights_path else None
        self.stat_names = stat_names or DEFAULT_STAT_NAMES

    def load(self) -> ReferenceData:
        df = pd.read_csv(self.csv_path)

        required = {"time", "type", "mean_N", "std_N"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                f"Reference CSV missing required columns: {missing}. "
                f"Found: {list(df.columns)}"
            )

        cell_types = sorted(df["type"].unique())
        times = df[df["type"] == cell_types[0]]["time"].values

        pop: Dict[str, np.ndarray] = {}
        for ct in cell_types:
            subset = df[df["type"] == ct]
            pop[ct] = subset["mean_N"].values.astype(float)

        ref_stats = compute_summary_stats(pop, times, self.stat_names, cell_types)

        if self.weights_path and self.weights_path.exists():
            with open(self.weights_path) as f:
                raw = json.load(f)
            stat_type_keys = [f"{s}_{ct}" for ct in cell_types for s in self.stat_names]
            precision_weights = np.array([raw.get(k, 1.0) for k in stat_type_keys])
        else:
            precision_weights = np.ones(len(ref_stats))

        return ReferenceData(
            ref_stats=ref_stats,
            precision_weights=precision_weights,
            stat_names=self.stat_names,
            cell_types=cell_types,
        )
