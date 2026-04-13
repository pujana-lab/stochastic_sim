from src.simulation_config import SimulationConfig
from src.clone import Clone

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
        Kt = max(self.config.Kmin, self.config.K0 - self.config.decline * t)
        return max(0.0, 1.0 - total_N / Kt) if Kt > 0 else 0.0

class AdaptedCrowding(CrowdingStrategy):
    def __init__(self, config: SimulationConfig):
        self.config = config

    def crowding(self, clone: Clone, t: float, total_N: int) -> float:
        if not self.config.use_logistic:
            return 1.0
        cfg = self.config
        Kt = max(cfg.Kmin, cfg.K0 / (1 - clone.death_rate / clone.birth_rate) - cfg.decline * t)
        return max(0.0, 1.0 - total_N / Kt) if Kt > 0 else 0.0