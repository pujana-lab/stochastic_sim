from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation

CELL_TYPES = ["base", "mutated", "immune", "exhausted"]


class GillespieSimulatorAdapter:
    def __init__(
        self,
        base_config: SimulationConfig,
        cell_types: List[str] = CELL_TYPES,
    ):
        self.base_config = base_config
        self.cell_types = cell_types

    def run(
        self, params: Dict[str, float]
    ) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        config_dict = dict(self.base_config.__dict__)
        config_dict.update(params)
        config_dict["seed"] = params.get("seed", self.base_config.seed)
        config = SimulationConfig(**config_dict)

        sim = TumorSimulation(config=config)
        times, history, _ = sim.run()

        t = np.array(times)
        pop: Dict[str, np.ndarray] = {
            ct: np.zeros(len(t), dtype=float) for ct in self.cell_types
        }
        for i, snap in enumerate(history):
            for cid, cd in snap.items():
                typ = cd["Type"]
                if typ in pop:
                    pop[typ][i] += cd["N"]

        return t, pop
