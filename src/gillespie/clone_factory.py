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
                birth_rate=self.config.lambda0,
                death_rate=self.config.mu0,
                mutation_rate=self.config.nu0,
                exhaustion_rate= 0.0,
                K= self.config.K0,
                instability=self.config.instability_0,
                buildup=self.config.buildup_0,
                d1=self.config.d1_0,
                d2=self.config.d2_0,
                parent=parent
            )
        elif clone_type == "mutated":
            x= MutatedClone(
                clone_id=clone_id,
                N=N,
                birth_rate=self.config.lambda0 * (1.0 + self.config.fitness_gain),
                death_rate=self.config.mu0,
                mutation_rate=0.0,
                exhaustion_rate= 0.0,
                K= self.config.K_mutant,
                instability=self.config.instability_0,
                buildup=self.config.buildup_0,
                d1=self.config.d1_0,
                d2=self.config.d2_0,
                parent=parent
            )
        elif clone_type == "immune":
            x= ImmuneClone(
                clone_id=(-1),
                N=N,
                
                birth_rate= self.config.lambda_Immune,
                death_rate= 0.0,
                mutation_rate= 0.0,
                exhaustion_rate= self.config.exhaustion_rate,
                K= self.config.K_immune,
                instability=0.0,
                buildup= 0.0,
                d1= 0.0,
                d2= 0.0,
                parent= None
            )
        elif clone_type == "exhausted":
            x= ExhaustedClone(
                clone_id=(-2),
                N=N,
                birth_rate=0.0,
                death_rate=self.config.mu_Exhausted,
                mutation_rate=0.0,
                exhaustion_rate= 0.0,
                K= 0.0,
                instability=0.0,
                buildup=0.0,
                d1=0.0,
                d2=0.0,
                parent= None
            )
        else:
            raise ValueError(f"Unknown clone type: {clone_type}")
        return x