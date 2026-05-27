from dataclasses import dataclass
from typing import Optional

@dataclass(frozen=True)
class SimulationConfig:
    N0: int = 10
    N_immune: int = 10
    N_exhausted: int = 0
    lambda0: float = 0.50
    lambda_Immune: float = 0.50
    mu0: float = 0.20
    mu_Immune: float = 0.30
    mu_Exhausted: float = 0.2
    exhaustion_rate: float= 0.01
    nu0: float = 0.00
    T_max: float = 5000
    seed: Optional[int] = None

    use_logistic: bool = False
    use_logistic_adapted: bool=True
    K0: int = 100
    K_immune: int = 50
    K_mutant: int = 3000

    decline: float = 0.0
    Kmin: float = 1.0

    fitness_gain: float = 0.05

    # IGNORE FOR NOW ---------------
    d1_0: float = 0.0
    d2_0: float = 0.0

    instability_0: float = 0.0
    buildup_0: float = 0.0000
    base_instability_buildup: float = 0.00000

    mutation_instability_jump: float = 0.0
    mutation_buildup_gain: float = 0.0
    # ----------------------------------


