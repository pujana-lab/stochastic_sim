from dataclasses import dataclass, field
from typing import Optional
import numpy as np
@dataclass(frozen=True)
class CloneParams:
    N0: int = 0
    K: int = 0
    lambda0: float = 0.0
    mu: float = 0.0
    nu: float = 0.0

#TODO: Arreglar este lio. Me gustaria poder meterle los inputs a simulation config de los clones que voy a usar que esten definidos ya fuera y luego al crear las subclases que cada una use el suyo. asi cada subclase tiene sus N,K,birth_rate etc


@dataclass(frozen=True)
class SimulationConfig:

    OMEGA: int = 100 #NUMERO MAXIMO DE CELULAS WT QUE SOPORTA CUANDO NO HAY COMPETICION


    # wt: CloneParams = field(CloneParams(N0=10000, K=10000, birth_rate=0.005, death_rate=0.002, mutation_rate=0.0002))
    # mutant: CloneParams = field( CloneParams(N0=200, K=3000, birth_rate=0.005, death_rate=0.002))
    # immune: CloneParams = field( CloneParams(N0=100, K=50, birth_rate=0.005, death_rate=0.003))
    # exhausted: CloneParams = field(CloneParams(N0=0, death_rate=0.002))
    # TODO: searar parametros por tipo celular. Cambiar pparametros globales aparte. Intentar imlementar escalado por system size.
    N0: int = 100
    N_immune: int = 50
    N_exhausted: int = 0
    N_mutant: int = 0

    # base rates
    lambda0: float = 0.005
    lambda_Immune: float = 0.005
    mu0: float = 0.002
    mu_Immune: float = 0.003
    mu_Exhausted: float = 0.002
    nu0: float = 0.0002

    # simulation control
    T_max: float = 2000
    seed: Optional[int] = None

    # logistic / carrying capacity
    use_logistic: bool = True  # ESTO SIEMPRE TIENE QUE SER TRUE O EL CRECIMIENTO EXPONENCIAL EXPLOTA
    use_logistic_adapted: bool = True
    K0: int = OMEGA
    K_immune: int = np.ceil(OMEGA/2)
    K_mutant: int = OMEGA * 2
    decline: float = 0.0
    Kmin: float = 1.0

    fitness_gain: float = 0.2

    # interaction parameters
    theta_I: float = 0.0005  # Kill rate: immune cells killing mutated cells (N_mutant * N_immune * theta_I)
    beta: float = 0.0004   # Activation rate: mutated cells activating immune cells (N_immune * N_mutant * beta)

    # instability / mutation parameters
    d1_0: float = 0.0
    d2_0: float = 0.0
    instability_0: float = 0.0
    buildup_0: float = 0.0
    base_instability_buildup: float = 0.0
    mutation_instability_jump: float = 0.0
    mutation_buildup_gain: float = 0.0

    # misc
    verbose: bool = True
    scale: bool = True
    decay: bool = True
    
    
    
    #TODO: implementar una forma de escalar por el system size parametros sin tener que hacerlo a mano (relacionado con task de linea 6)


    