from dataclasses import dataclass, field
from typing import Optional, TYPE_CHECKING
from src.gillespie.cloneId import CloneId
import numpy as np

if TYPE_CHECKING:
    from src.gillespie.crowding_strategy import CrowdingStrategy
@dataclass(frozen=True)
class CellTypeConfig:
    default_id: CloneId = None
    N: int = 0
    K: float | int  = None 
    lambda0: float = 0.005
    mu: float = 0.002
    nu: float = 0.0
    omega_exhaust: float = 0.0
    next_mutation: str = ""

#TODO: Arreglar este lio. Me gustaria poder meterle los inputs a simulation config de los clones que voy a usar que esten definidos ya fuera y luego al crear las subclases que cada una use el suyo. asi cada subclase tiene sus N,K,birth_rate etc
# Defaults estáticos
_DEFAULT_CELL_PARAMETERS = {
    "base": CellTypeConfig(
        default_id= (),
        N=100,
        K=1,
        nu=0.0002,
        next_mutation= "mutated"
    ),
    "immune": CellTypeConfig(
        default_id= (-1,),
        N=50,
        K=0.5,
        lambda0= 0.005,
        omega_exhaust=0.003,
        mu= 0.0,
    ),
    "mutated": CellTypeConfig(
        default_id= (-3,),
        N=10,
        K=2,
    ),
    "exhausted": CellTypeConfig(
        default_id=(-2,),
        N=0,
        lambda0= 0.0
    ),
}

@dataclass(frozen=True)
class SimulationConfig:

    OMEGA: int = 4000 #NUMERO MAXIMO DE CELULAS WT QUE SOPORTA CUANDO NO HAY COMPETICION
    
    
    #TODO: IMPORTANTE: Mover los parametros de las celulas de flat a por tipo. que tengan un string acorde con el tipo y se puedan acceder al crear celulas como my_defaults=self.config.params("cell_type") y luego hacer birth_rate=my_defaults.lambda0
    cell_parameters: dict = field(default_factory=lambda : dict(_DEFAULT_CELL_PARAMETERS))
   
 
    # simulation control
    T_max: float = 2000
    seed: Optional[int] = None

   
    decline: float = 0.0
    Kmin: float = 1

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
    scale: bool = True # Mantener siempre true a 
    decay: bool = False
    use_logistic: bool = True  # ESTO SIEMPRE TIENE QUE SER TRUE O EL CRECIMIENTO EXPONENCIAL EXPLOTA
    use_logistic_adapted: bool = True
    
    # Crowding strategy (será inicializado en __post_init__)
    crowding_strategy: "CrowdingStrategy" = field(init=False, default=None)
    
    
    def __post_init__(self):
        """Scale interaction parameters by OMEGA and initialize crowding strategy"""
        if self.scale:
            object.__setattr__(self, 'beta', self.beta / self.OMEGA)
            object.__setattr__(self, 'theta_I', self.theta_I / self.OMEGA)
             # Initialize crowding strategy (after K values are calculated)
            from src.gillespie.crowding_strategy import SimpleCrowding, AdaptedCrowding
            strategy_class = AdaptedCrowding if self.use_logistic_adapted else SimpleCrowding
            object.__setattr__(self, 'crowding_strategy', strategy_class(self))
            # Escalar K para cada tipo de célula
            new_cell_params = {}
            for cell_type, config in self.cell_parameters.items():
                if config.K is not None:
                    scaled_K = int(np.ceil(config.K * self.OMEGA))
                    scaled_omega_eshaust = config.omega_exhaust / self.OMEGA if config.omega_exhaust is not 0.0 else config.omega_exhaust
                    # Crear nuevo CellTypeConfig con K escalado
                    new_config = CellTypeConfig(
                        default_id= config.default_id,
                        N=config.N,
                        K=scaled_K,
                        lambda0=config.lambda0,
                        mu=config.mu,
                        nu=config.nu,
                        omega_exhaust=scaled_omega_eshaust,
                        next_mutation=config.next_mutation,
                    )
                    new_cell_params[cell_type] = new_config
                else:
                    new_cell_params[cell_type] = config
            
            object.__setattr__(self, 'cell_parameters', new_cell_params)
        
       
    