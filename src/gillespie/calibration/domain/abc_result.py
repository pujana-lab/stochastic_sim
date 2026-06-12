from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np


@dataclass
class AbcResult:
    particles: List[List[Dict[str, float]]] = field(default_factory=list)
    weights: List[np.ndarray] = field(default_factory=list)
    distances: List[np.ndarray] = field(default_factory=list)
    diagnostics: List[Dict] = field(default_factory=list)
    epsilon_history: List[float] = field(default_factory=list)

    def ess(self, gen: int) -> float:
        w = self.weights[gen]
        return 1.0 / float(np.sum(w ** 2))

    def n_generations(self) -> int:
        return len(self.particles)

    def summary(self, gen: int) -> Dict:
        dists = self.distances[gen]
        diag = self.diagnostics[gen]
        w = self.weights[gen]
        return {
            "generation": gen,
            "n_particles": len(self.particles[gen]),
            "epsilon": self.epsilon_history[gen] if gen < len(self.epsilon_history) else None,
            "ess": self.ess(gen),
            "ess_ratio": self.ess(gen) / len(w) if len(w) > 0 else 0.0,
            "mean_distance": float(np.mean(dists)),
            "min_distance": float(np.min(dists)),
            **diag,
        }

    def all_summaries(self) -> List[Dict]:
        return [self.summary(g) for g in range(self.n_generations())]
