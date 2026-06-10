from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone import Clone
from src.gillespie.tissue_state import TissueState
from abc import ABC, abstractmethod

class CrowdingStrategy(ABC):
    
    def crowding(self, clone: Clone, t: float, tissue_state: "TissueState") -> float:
        if not self.config.use_logistic:
            return 1.0
        if clone.get_type() == "exhausted":
            return 1.0
        numerator_N = clone.crowding_numerator(tissue_state)
        if self.config.decay:
            Kt = self.calculate_K(clone,t)
        else:
            Kt = clone.K
        if not Kt:
            return 1.0
        return max(0.0, 1.0 - numerator_N / Kt) if Kt > 0 else 0.0
    @abstractmethod
    def calculate_K(self, clone: Clone, t: float) -> float:
        ...
class SimpleCrowding(CrowdingStrategy):
    def __init__(self, config: SimulationConfig):
        self.config = config

    def calculate_K(self, clone, t):
        return max(clone.K_min, clone.K - self.config.decline * t)
    
    

class AdaptedCrowding(CrowdingStrategy):
    def __init__(self, config: SimulationConfig):
        self.config = config
    
    def calculate_K(self, clone, t):
        if clone.birth_rate <= 0.0:
            return clone.K_min
        denom = 1 - clone.death_rate / clone.birth_rate
        if denom <= 0:
            return clone.K_min
        return max(clone.K_min, clone.K / denom - self.config.decline * t)

       

