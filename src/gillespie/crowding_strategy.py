from abc import ABC, abstractmethod
from typing import TYPE_CHECKING
import math

if TYPE_CHECKING:
    from src.gillespie.simulation_config import SimulationConfig
    from src.gillespie.clone import Clone
    from src.gillespie.tissue_state import TissueState

class CrowdingStrategy(ABC):
    def __init__(self, config: "SimulationConfig"):
        self.config = config
    
    def crowding(self, clone: "Clone", t: float, tissue_state: "TissueState") -> float:
        if not self.config.use_logistic:
            return 1.0
        if clone.get_type() == "exhausted":
            return 1.0
            
        numerator_N = clone.crowding_numerator(tissue_state)
        Kt = self.calculate_K(clone, t)
        
        # Si Kt es infinito, numerator_N / inf = 0.0 -> crowding devuelve 1.0 (sin freno)
        if math.isinf(Kt):
            return 1.0
            
        # Al asegurar numéricamente que Kt >= 1.0 a través de K_min, eliminamos riesgos de ZeroDivisionError
        return max(0.0, 1.0 - numerator_N / Kt)

    @abstractmethod
    def calculate_K(self, clone: "Clone", t: float) -> float:
        """Calcula la capacidad portante microscópica (Kt) para el tiempo dado."""
        ...

class SimpleCrowding(CrowdingStrategy):
    def calculate_K(self, clone: "Clone", t: float) -> float:
        # Aseguramos que clone.K_min sea como mínimo 1.0 en la configuración de tus objetos
        k_floor = max(1.0, getattr(clone, 'K_min', 1.0))
        
        if self.config.decay:
            return max(k_floor, clone.K - self.config.decline * t)
        return max(k_floor, clone.K)

class AdaptedCrowding(CrowdingStrategy):
    def calculate_K(self, clone: "Clone", t: float) -> float:
        # Salvaguarda 1: Si no hay tasa de nacimiento, la población no expande su nicho
        if clone.birth_rate <= 0.0:
            return float('inf')
            
        denom = 1.0 - (clone.death_rate / clone.birth_rate)
        
        # Salvaguarda 2: Población en involución (muerte >= nacimiento). 
        # No satura el volumen por homeostasis, se extingue por su propia tasa lineal.
        if denom <= 0.0:
            return float('inf')
            
        # Flujo estándar: Inflado inverso de la capacidad microscópica
        K_inflated = clone.K / denom
        
        k_floor = max(1.0, getattr(clone, 'K_min', 1.0))
        
        # Salvaguarda 3: Protección frente al decaimiento temporal del tejido
        if self.config.decay:
            return max(k_floor, K_inflated - self.config.decline * t)
        return max(k_floor, K_inflated)
