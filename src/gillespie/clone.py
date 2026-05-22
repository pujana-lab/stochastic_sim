from dataclasses import dataclass
from typing import Optional

from src.gillespie.cloneId import CloneId
from src.gillespie.clone_type import CloneType
from src.gillespie.clone_factory import CloneFactory
@dataclass
class Clone:
    clone_id: CloneId
    N: int
    # cell_type: CloneType = CloneType.BASE
    birth_rate: float
    death_rate: float
    mutation_rate: float
    exhaustion_rate: float
    cell_type: str = ""
    K: int
    instability: float = 0.0
    buildup: float = 0.0
    d1: float = 0.0
    d2: float = 0.0

    parent: Optional[CloneId] = None
    children_count: int = 0

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

    def mutate(self) -> "Clone":

        if self.N <= 0:
            raise ValueError("Cannot mutate a dead clone.")

        self.kill()
        child = CloneFactory.create_clone(
            clone_id=self.next_child_id(),
            clone_type= "mutated",
            parent= self.clone_id
        )
        
    
        return child        
    def __str__(self) -> str:
        return str(self.cell_type)
class WildTypeClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.BASE

class MutatedClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.MUTATED
    def death_rate_effective(self,N_I: int,activation:float):
        return super().death_rate_effective() + self.N*N_I*activation
        
    

class ImmuneClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.IMMUNE


class ExhaustedClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.EXHAUSTED   

        

