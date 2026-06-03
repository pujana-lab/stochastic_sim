from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone import Clone
from src.gillespie.tissue_state import TissueState
from abc import ABC, abstractmethod

class CrowdingStrategy(ABC):
    @abstractmethod
    def crowding(self, clone: Clone, t: float, tissue_state: "TissueState") -> float:
        ...

class SimpleCrowding(CrowdingStrategy):
    def __init__(self, config: SimulationConfig):
        self.config = config

    def crowding(self, clone: Clone, t: float, tissue_state: "TissueState") -> float:
        
        if not self.config.use_logistic:
            return 1.0
        numerator_N = clone.crowding_numerator(tissue_state)
        Kt = max(clone.K_min, clone.K - self.config.decline * t)
        return max(0.0, 1.0 - numerator_N / Kt) if Kt > 0 else 0.0

class AdaptedCrowding(CrowdingStrategy):
    def __init__(self, config: SimulationConfig):
        self.config = config

    def crowding(self, clone: Clone, t: float, tissue_state: "TissueState") -> float:

        if not self.config.use_logistic:
            return 1.0
        cfg = self.config
        numerator_N = clone.crowding_numerator(tissue_state)
        if clone.get_type() == "exhausted":
            return 1.0
        denom = 1 - clone.death_rate / clone.birth_rate
        if denom <= 0:
            return 0.0
        Kt = max(clone.K_min, clone.K / denom - cfg.decline * t)
        return max(0.0, 1.0 - numerator_N / Kt) if Kt > 0 else 0.0