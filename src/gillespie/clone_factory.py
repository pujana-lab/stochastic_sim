from typing import Dict, Optional

from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone import Clone, WildTypeClone,MutatedClone
from src.gillespie.cloneId import CloneId
class CloneFactory:
    def __init__(self, config: SimulationConfig):
        self.config = config

    def create_clone(self,clone_id: CloneId, clone_type: str = "wild_type", parent: Clone | None = None) -> Clone:
        if clone_type == "wild_type":
            x= WildTypeClone(
                clone_id=clone_id,
                N=self.config.N0,
                birth_rate=self.config.lambda0,
                death_rate=self.config.mu0,
                mutation_rate=self.config.nu0,
                instability=self.config.instability_0,
                buildup=self.config.buildup_0,
                d1=self.config.d1_0,
                d2=self.config.d2_0,
                parent=parent
            )
        elif clone_type == "mutated":
            x= MutatedClone(
                clone_id=clone_id,
                N=self.config.N0,
                birth_rate=self.config.lambda0,
                death_rate=self.config.mu0,
                mutation_rate=self.config.nu0,
                instability=self.config.instability_0,
                buildup=self.config.buildup_0,
                d1=self.config.d1_0,
                d2=self.config.d2_0,
                parent=parent
            )
        else:
            raise ValueError(f"Unknown clone type: {clone_type}")
        return x