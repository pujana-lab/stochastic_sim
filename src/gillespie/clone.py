from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING, Dict

from src.gillespie.cloneId import CloneId
from src.gillespie.clone_type import CloneType
from src.gillespie.simulation_config import SimulationConfig

if TYPE_CHECKING:
    from src.gillespie.tissue_state import TissueState


@dataclass
class Clone:
    def __init__(self, clone_id: CloneId, config: SimulationConfig, N: int = 1, parent: Optional[CloneId] = None):
        
        self.clone_id: CloneId = clone_id
        self.config: SimulationConfig = config
        self.N = N
        self.parent= parent
        self.children_count = 0
        self.crowding_value = 0.0
        
        self.birth_rate: float = config.lambda0
        self.death_rate: float = config.mu0
        self.mutation_rate: float = 0.0
        self.exhaustion_rate: float = 0.0

        self.K: int = 0.0
        self.K_min: int = config.Kmin


        self.cell_type: str = ""
        self.instability: float = 0.0
        self.buildup: float = 0.0
        self.d1: float = 0.0
        self.d2: float = 0.0

    def is_alive(self) -> bool:
        return self.N > 0

    def mutation_multiplier(self) -> float:
        return 1.0 + self.instability

    def birth_rate_effective(self) -> float:
        return self.birth_rate * self.N 

    def death_rate_effective(self) -> float:
        return self.death_rate * self.N

    def mutation_rate_effective(self) -> float:
        return self.mutation_rate * self.N * self.mutation_multiplier()

    def exhaustion_rate_effective(self) -> float:
        return self.exhaustion_rate * self.N
        
    def divide(self) -> None:
        self.N += 1

    def kill(self) -> None:
        self.N = max(0, self.N - 1)

    def next_child_id(self) -> CloneId:
        self.children_count += 1
        return self.clone_id + (self.children_count,)

    def advance_instability(self, dt: float, base_buildup: float) -> None:
        self.instability += (base_buildup + self.buildup) * dt       
    def __str__(self) -> str:
        return self.cell_type.value if isinstance(self.cell_type, CloneType) else str(self.cell_type)


# de nuevo aqui hay que pasarle el diccionario de valores N_loquesea y desde ahi pasarle el crowdingstrategy que elijamos y calcule el numero

class WildTypeClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.BASE

    def crowding_numerator(self,popmap: Dict[CloneType,int]) -> int:
        return popmap["mutated"]+popmap["base"]

    
    def birth_rate_effective(self, tissue_state: "TissueState", crowding: float) -> float:
        base_birth = super().birth_rate_effective()
        return base_birth * crowding

class MutatedClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.MUTATED

    def crowding_numerator(self,popmap: Dict[CloneType,int]) -> int:
        return popmap["mutated"]

    def birth_rate_effective(self,tissue_state: "TissueState",crowding:float ) -> float:
        base_birth = super().birth_rate_effective()
        return base_birth * crowding

    
    def death_rate_effective(self, tissue_state: "TissueState") -> float:
        """Death rate increases due to immune cell killing.
        Formula: base_death_rate * N + N * N_immune * theta_I
        (theta_I is a global parameter from config)
        """
        base_death = super().death_rate_effective()
        n_immune = tissue_state.population_by_type("immune")
        immune_killing = self.N * n_immune * self.config.theta_I
        return base_death + immune_killing

    

class ImmuneClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.IMMUNE

    def crowding_numerator(self,popmap: Dict[CloneType,int]) -> int:
        return popmap["immune"]+popmap["exhausted"]
    
    def birth_rate_effective(self, tissue_state: "TissueState", crowding: float) -> float:
        """Birth rate increases through activation by cancer cells.
        Formula: base_birth_rate + N * N_cancer * beta
        (beta is a global parameter from config)
        """
        base_birth = super().birth_rate_effective()
        n_cancer = tissue_state.population_by_type("mutated")
        activation_boost = self.N * n_cancer * self.config.beta
        # crowding_effect= #crowding donde tenemos que el numerador de la logistica es n_exhausted+sel
        return base_birth*crowding + activation_boost
    
    def exhaustion_rate_effective(self, tissue_state: "TissueState") -> float:
        """Exhaustion: exhaustion_rate * N_immune * N_cancer"""
        n_cancer = tissue_state.population_by_type("mutated")
        return self.exhaustion_rate * self.N * n_cancer


class ExhaustedClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.EXHAUSTED 

    def crowding_numerator(self,popmap: Dict[CloneType,int]) -> int:
        return 0
    
    def birth_rate_effective(self, tissue_state: "TissueState", crowding: float) -> float:
        """Exhausted cells cannot divide."""
        return 0.0

        

