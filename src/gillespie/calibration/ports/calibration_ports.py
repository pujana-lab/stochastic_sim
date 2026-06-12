from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Protocol, Tuple

import numpy as np


@dataclass(frozen=True)
class ReferenceData:
    ref_stats: np.ndarray
    precision_weights: np.ndarray
    stat_names: List[str]
    cell_types: List[str]


class SimulatorPort(Protocol):
    def run(self, params: Dict[str, float]) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        ...


class ReferenceDataPort(Protocol):
    def load(self) -> ReferenceData:
        ...


class OutputPort(Protocol):
    def save_generation(
        self,
        gen: int,
        particles: List[Dict[str, float]],
        weights: np.ndarray,
        distances: np.ndarray,
    ) -> None:
        ...

    def save_epsilon_schedule(self, eps_history: List[float]) -> None:
        ...

    def save_generation_diagnostics(self, diagnostics: List[Dict]) -> None:
        ...

    def save_manifest(
        self,
        calibration_config: dict,
        priors: dict,
    ) -> None:
        ...

    def find_last_generation(self) -> Optional[int]:
        ...

    def load_generation(
        self, gen: int
    ) -> Tuple[List[Dict[str, float]], np.ndarray, np.ndarray]:
        ...
