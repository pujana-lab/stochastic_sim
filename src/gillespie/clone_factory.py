from typing import Optional

from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone import Clone, WildTypeClone, MutatedClone, ImmuneClone, ExhaustedClone
from src.gillespie.cloneId import CloneId


class CloneFactory:
    def __init__(self, config: SimulationConfig):
        self.config = config

    def create_clone(
        self, 
        clone_id: CloneId, 
        clone_type: str = "wild_type",
        N: int = 1, 
        parent: Optional[Clone] = None
    ) -> Clone:
        if clone_type == "wild_type":
            clone = WildTypeClone(
                clone_id=clone_id,
                config=self.config,
                N=N,
                parent=parent
            )
            clone.mutation_rate = self.config.nu0
            clone.K = self.config.K0
            
        elif clone_type == "mutated":
            clone = MutatedClone(
                clone_id=clone_id,
                config=self.config,
                N=N,
                parent=parent
            )
            clone.birth_rate = self.config.lambda0 * (1.0 + self.config.fitness_gain)
            clone.K = self.config.K_mutant
            
        elif clone_type == "immune":
            clone = ImmuneClone(
                clone_id=clone_id,
                config=self.config,
                N=N,
                parent=parent
            )
            clone.birth_rate = self.config.lambda_Immune
            clone.death_rate = 0.0
            clone.exhaustion_rate = self.config.mu_Immune
            
        elif clone_type == "exhausted":
            clone = ExhaustedClone(
                clone_id=clone_id,
                config=self.config,
                N=N,
                parent=parent
            )
            clone.birth_rate = 0.0
            clone.death_rate = self.config.mu_Exhausted
            
        else:
            raise ValueError(f"Unknown clone type: {clone_type}")
            
        return clone