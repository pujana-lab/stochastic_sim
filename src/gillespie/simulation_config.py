from dataclasses import dataclass, field
from typing import Optional
import numpy as np
@dataclass(frozen=True)
class CloneParams:
    N: int = 0
    K: int = 0
    lambda0: float = 0.0
    mu: float = 0.0
    nu: float = 0.0

#TODO: Arreglar este lio. Me gustaria poder meterle los inputs a simulation config de los clones que voy a usar que esten definidos ya fuera y luego al crear las subclases que cada una use el suyo. asi cada subclase tiene sus N,K,birth_rate etc


@dataclass(frozen=True)
class SimulationConfig:

    OMEGA: int = 100 #NUMERO MAXIMO DE CELULAS WT QUE SOPORTA CUANDO NO HAY COMPETICION
    
    #TODO: IMPORTANTE: Mover los parametros de las celulas de flat a por tipo. que tengan un string acorde con el tipo y se puedan acceder al crear celulas como my_defaults=self.config.params("cell_type") y luego hacer birth_rate=my_defaults.lambda0
    N0: int = 100
    N_immune: int = 50
    N_exhausted: int = 0
    N_mutant: int = 0

    # base rates
    lambda0: float = 0.005
    lambda_Immune: float = 0.005
    mu0: float = 0.002
    
    mu_Exhausted: float = 0.002
    nu0: float = 0.0002

    # simulation control
    T_max: float = 2000
    seed: Optional[int] = None

    # logistic / carrying capacity
    #TODO: revisar esto. no se si definir Omega como volumen o como numero discreto de celulas. el problema es que si es como numero discreto los K tienen que ser fracciones y si es como volumen los K son enteros. en cualquier caso al multiplicar por Omega luego siempre acaban siendo enteros. pero va variando el tipo de dato lo cual no creo que sea optimo.
    K0: float|int = 1
    K_immune: float|int = np.ceil(1/2)
    K_mutant: float|int = 2
    decline: float = 0.0
    Kmin: float = 1

    fitness_gain: float = 0.2

    # interaction parameters
    theta_I: float = 0.0005  # Kill rate: immune cells killing mutated cells (N_mutant * N_immune * theta_I)
    beta: float = 0.0004   # Activation rate: mutated cells activating immune cells (N_immune * N_mutant * beta)
    mu_Immune: float = 0.003

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
    scale: bool = True # Mantener siempre true a 
    decay: bool = False
    use_logistic: bool = True  # ESTO SIEMPRE TIENE QUE SER TRUE O EL CRECIMIENTO EXPONENCIAL EXPLOTA
    use_logistic_adapted: bool = True
    
    
    

    def __post_init__(self):
        """
        Calcula parámetros derivados basándose en valores base y OMEGA.
        Si se define beta_base, beta se reescala como: beta = beta_base / OMEGA
        Lo mismo para theta_I_base.
        
        Ejemplo JSON:
        {
            "OMEGA": 100,
            "beta_base": 0.04,
            "theta_I_base": 0.05
        }
        
        Resultado:
        - beta = 0.04 / 100 = 0.0004
        - theta_I = 0.05 / 100 = 0.0005
        """
        # Reescalar beta si se proporciona valor base
        if self.scale:
            object.__setattr__(self, 'beta', self.beta / self.OMEGA)      
            object.__setattr__(self, 'theta_I', self.theta_I / self.OMEGA)
            object.__setattr__(self,'mu_Immune', self.mu_Immune / self.OMEGA)
            object.__setattr__(self,'K_0', int(np.ceil(self.K0 * self.OMEGA)))
            object.__setattr__(self,'K_immune',int(np.ceil( self.K_immune * self.OMEGA)))
            object.__setattr__(self,'K_mutant', int(np.ceil(self.K_mutant * self.OMEGA)))
            

    