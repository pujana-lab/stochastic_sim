from dataclasses import dataclass
from typing import Optional

from src.gillespie.cloneId import CloneId
from src.gillespie.clone_type import CloneType
from src.gillespie.simulation_config import SimulationConfig


@dataclass
class Clone:
    def __init__(self, clone_id: CloneId, config: SimulationConfig, N: int = 1, parent: Optional[CloneId] = None):
        
        self.clone_id: CloneId = clone_id
        self.N = N
        self.parent= parent
        self.children_count= 0
        
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

    def birth_rate_effective(self, crowding: float) -> float:
        return self.birth_rate * self.N * crowding 

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
    def birth_rate_effective(self, crowding):
        return super()

class MutatedClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.MUTATED
    def death_rate_effective(self,N_I: int,killrate:float):
        return super().death_rate_effective() + self.N*N_I*killrate

    

class ImmuneClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.IMMUNE
    def birth_rate_effective(self, crowding,N_C: int,activation:float):
        return super().birth_rate_effective(crowding)+ self.N*N_C*activation
    def exhaustion_rate_effective(self,N_C: int):
        return super().exhaustion_rate_effective()*N_C


class ExhaustedClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.EXHAUSTED  
    def birth_rate_effective(self, crowding):
        return 0.0

        

