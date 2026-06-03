from dataclasses import dataclass
from typing import Optional

@dataclass(frozen=True)
class SimulationConfig:
    # TODO: searar parametros por tipo celular. Cambiar pparametros globales aparte. Intentar imlementar escalado por system size.
    N0: int = 10000
    N_immune: int = 100
    N_exhausted: int = 0
    N_mutant: int = 200
    lambda0: float = 0.005
    lambda_Immune: float = 0.005
    mu0: float = 0.002
    mu_Immune: float = 0.003
    mu_Exhausted: float = 0.002
    nu0: float = 0.0002
    T_max: float = 50
    seed: Optional[int] = None

    use_logistic: bool = True # ESTO SIEMPRE TIENE QUE SER TRUE O EL CRECIMIENTO EXPONENCIAL EXPLOTA 
    use_logistic_adapted: bool = True
    K0: int = 100
    K_immune: int = 50
    
    K_mutant: int = 3000

    decline: float = 0.0
    Kmin: float = 1.0

    fitness_gain: float = 0.2
    
    # Interaction parameters
    theta_I: float = 0.0  # Kill rate: immune cells killing mutated cells (N_mutant * N_immune * theta_I)
    beta: float = 0.0     # Activation rate: mutated cells activating immune cells (N_immune * N_mutant * beta)

    # IGNORE FOR NOW ---------------
    d1_0: float = 0.0
    d2_0: float = 0.0

    instability_0: float = 0.0
    buildup_0: float = 0.0000
    base_instability_buildup: float = 0.00000

    mutation_instability_jump: float = 0.0
    mutation_buildup_gain: float = 0.0
    # ----------------------------------

    verbose: bool = True

    scale: bool =True
    system_size: int = 10000
    
    omega = system_size/K0
    #TODO: implementar una forma de escalar por el system size parametros sin tener que hacerlo a mano (relacionado con task de linea 6)


