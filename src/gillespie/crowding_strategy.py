from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone import Clone

from abc import ABC, abstractmethod

class CrowdingStrategy(ABC):
    @abstractmethod
    def crowding(self, clone: Clone, t: float, total_N: int) -> float:
        ...

class SimpleCrowding(CrowdingStrategy):
    def __init__(self, config: SimulationConfig):
        self.config = config

    def crowding(self, clone: Clone, t: float, total_N: int) -> float:
        if not self.config.use_logistic:
            return 1.0
        Kt = max(clone.K_min, clone.K - self.config.decline * t)
        return max(0.0, 1.0 - total_N / Kt) if Kt > 0 else 0.0

class AdaptedCrowding(CrowdingStrategy):
    def __init__(self, config: SimulationConfig):
        self.config = config

    def crowding(self, clone: Clone, t: float, total_N: int) -> float:
        if not self.config.use_logistic:
            return 1.0
        cfg = self.config
        denom = 1 - clone.death_rate / clone.birth_rate
        if denom <= 0:
            return 0.0
        Kt = max(clone.K_min, clone.K / denom - cfg.decline * t)
        return max(0.0, 1.0 - total_N / Kt) if Kt > 0 else 0.0