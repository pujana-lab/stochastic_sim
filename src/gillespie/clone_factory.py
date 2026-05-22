from typing import Dict, Optional

from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone import Clone, WildTypeClone,MutatedClone,ImmuneClone,ExhaustedClone
from src.gillespie.cloneId import CloneId
class CloneFactory:
    def __init__(self, config: SimulationConfig):
        self.config = config

    def create_clone(self,clone_id: CloneId, clone_type: str = "wild_type",N: int = 1, parent: Clone | None = None) -> Clone:
        if clone_type == "wild_type":
            x= WildTypeClone(
                clone_id=clone_id,
                N=N,
                mutation_rate=self.config.nu0,
                K= self.config.K0,
            )
        elif clone_type == "mutated":
            x= MutatedClone(
                clone_id=clone_id,
                N=N,
                birth_rate=self.config.lambda0 * (1.0 + self.config.fitness_gain),
                K= self.config.K_mutant,
                parent=parent
            )
        elif clone_type == "immune":
            x= ImmuneClone(
                clone_id=(-1),
                N=N,
                
                birth_rate= self.config.lambda_Immune,
                death_rate= 0.0,
                exhaustion_rate= self.config.exhaustion_rate,
                K= self.config.K_immune,
            )
        elif clone_type == "exhausted":
            x= ExhaustedClone(
                clone_id=(-2),
                N=N,
                birth_rate=0.0,
                death_rate=self.config.mu_Exhausted,
            )
        else:
            raise ValueError(f"Unknown clone type: {clone_type}")
        return x