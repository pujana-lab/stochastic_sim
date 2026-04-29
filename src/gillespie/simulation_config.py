from dataclasses import dataclass
from typing import Optional

@dataclass(frozen=True)
class SimulationConfig:
    N0: int = 10
    lambda0: float = 0.50
    mu0: float = 0.20
    nu0: float = 0.00

    d1_0: float = 0.0
    d2_0: float = 0.0

    instability_0: float = 0.0
    buildup_0: float = 0.0000
    base_instability_buildup: float = 0.00000

    mutation_instability_jump: float = 0.05
    mutation_buildup_gain: float = 0.0001

    T_max: float = 5000
    seed: Optional[int] = None

    use_logistic: bool = True
    use_logistic_adapted: bool=True
    K0: float = 100
    decline: float = 0.0
    Kmin: float = 1.0

    fitness_gain: float = 0.05

