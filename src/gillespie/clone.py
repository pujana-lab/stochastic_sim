from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING, Dict, Type

from src.gillespie.cloneId import CloneId
from src.gillespie.simulation_config import SimulationConfig

if TYPE_CHECKING:
    from src.gillespie.tissue_state import TissueState



class Clone:

    _registry: Dict[str, Type["Clone"]] = {}
    _type_name: str = None

    def __init_subclass__(cls, clone_type: str = None, **kwargs):
        super().__init_subclass__(**kwargs)
        if clone_type:
            Clone._registry[clone_type] = cls
            cls._type_name = clone_type

    def __init__(self, clone_id: CloneId, config: SimulationConfig, N: int = 1, parent: Optional[CloneId] = None):
        
        #TODO: todo lo que dependa de la clase y no de la instancia deberia ir guardado dentro de la clase (sacar todos los self. de cosas que no sean dinamicas) y aunque lo sean deberian ir fuera igual y actualizarse para todo el tipo. lo unico que deberia ir dentro son las coasas que dependan del tiempo Y ADEMAS queramos trackear por clon y no por tipo (inestabilidad genomica en un futuro)
        self.clone_id: CloneId = clone_id
        self.config: SimulationConfig = config
        self.N: int = N
        self.parent: CloneId = parent
        self.children_count: int = 0
        self.crowding_value: float = 0.0
        
        self.birth_rate: float = config.lambda0
        self.death_rate: float = config.mu0
        self.mutation_rate: float = 0.0
        self.exhaustion_rate: float = 0.0

        self.K: int = 0.0
        self.K_min: int = config.Kmin

        self.next_mutation: str = None
        self.instability: float = 0.0
        self.buildup: float = 0.0
        self.d1: float = 0.0
        self.d2: float = 0.0

    def get_type(self) -> str:
        return self._type_name
    
    def crowding_numerator(self, tissue_state: "TissueState") -> int:
        return 0

    def is_alive(self) -> bool:
        return self.N > 0

    def mutation_multiplier(self) -> float:
        return 1.0 + self.instability

    def birth_rate_effective(self,tissue_state: "TissueState",crowding:float) -> float:
        return self.birth_rate * tissue_state.pop_map.get(self.get_type(),0) * crowding

    def death_rate_effective(self,tissue_state: "TissueState") -> float:
        return self.death_rate * tissue_state.pop_map.get(self.get_type(),0)

    def mutation_rate_effective(self,tissue_state: "TissueState") -> float:
        return self.mutation_rate * tissue_state.pop_map.get(self.get_type()) * self.mutation_multiplier()

    def exhaustion_rate_effective(self,tissue_state: "TissueState") -> float:
        return self.exhaustion_rate * tissue_state.pop_map.get(self.get_type())
        
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
        return self.get_type()


# de nuevo aqui hay que pasarle el diccionario de valores N_loquesea y desde ahi pasarle el crowdingstrategy que elijamos y calcule el numero

class WildTypeClone(Clone, clone_type = "base"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.next_mutation = "mutated"
        self.mutation_rate = self.config.nu0
        self.K = self.config.K0

    def crowding_numerator(self, tissue_state: "TissueState") -> int:
        pop_map = tissue_state.pop_map
        return pop_map.get("base", 0) + pop_map.get("mutated", 0)

class MutatedClone(Clone, clone_type = "mutated"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.birth_rate = self.config.lambda0 * (1.0 + self.config.fitness_gain)
        self.K = self.config.K_mutant

    def crowding_numerator(self, tissue_state: "TissueState") -> int:
        return tissue_state.pop_map.get("mutated", 0)

    
    def death_rate_effective(self, tissue_state: "TissueState") -> float:
        """Death rate increases due to immune cell killing.
        Formula: base_death_rate * N + N * N_immune * theta_I
        (theta_I is a global parameter from config)
        """
        base_death = super().death_rate_effective(tissue_state)
        self_n = tissue_state.pop_map.get("mutated", 0)
        n_immune = tissue_state.pop_map.get("immune", 0)
        immune_killing = self_n * n_immune * self.config.theta_I
        return base_death + immune_killing

    

class ImmuneClone(Clone, clone_type = "immune"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.birth_rate = self.config.lambda_Immune
        self.death_rate = 0.0
        self.exhaustion_rate = self.config.mu_Immune


    def crowding_numerator(self, tissue_state: "TissueState") -> int:
        pop_map = tissue_state.pop_map
        return pop_map.get("immune", 0) + pop_map.get("exhausted", 0)
    
    def birth_rate_effective(self, tissue_state: "TissueState", crowding: float) -> float:
        """Birth rate increases through activation by cancer cells.
        Formula: base_birth_rate + N * N_cancer * beta
        (beta is a global parameter from config)
        """
        base_birth = super().birth_rate_effective(tissue_state, crowding)
        n_self = tissue_state.pop_map.get("immune", 0)
        n_cancer = tissue_state.pop_map.get("mutated", 0)
        activation_boost = n_self * n_cancer * self.config.beta
        # crowding_effect= #crowding donde tenemos que el numerador de la logistica es n_exhausted+sel
        return base_birth + activation_boost
    
    def exhaustion_rate_effective(self, tissue_state: "TissueState") -> float:
        """Exhaustion: exhaustion_rate * N_immune * N_cancer"""
        base_exhaustion =  super().exhaustion_rate_effective(tissue_state= tissue_state)
        n_cancer = tissue_state.pop_map.get("mutated",0)
        return base_exhaustion * n_cancer


class ExhaustedClone(Clone, clone_type = "exhausted"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.birth_rate = 0.0
        self.death_rate = self.config.mu_Exhausted

    def crowding_numerator(self, tissue_state: "TissueState") -> int:
        return 0
    
    def birth_rate_effective(self, tissue_state: "TissueState", crowding: float) -> float:
        """Exhausted cells cannot divide."""
        return 0.0


        

