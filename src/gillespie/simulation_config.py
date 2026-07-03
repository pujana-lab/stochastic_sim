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
    fitness_gain: float = 0.0
    next_mutation: str = ""

_DEFAULT_CELL_PARAMETERS = {
    "base": CellTypeConfig(
        default_id= (),
        N=100,
        K=1,
        nu=0.00002,
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
        fitness_gain= 0.2, 
        N=2,
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
    T_max: float = 6000
    seed: Optional[int] = None

   
    decline: float = 0.0
    Kmin: float = 1

    

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
        """Initialize configuration: scale parameters and set up strategies."""
        self._initialize_crowding_strategy()
        if self.scale:
            self._scale_interaction_parameters()
            self._scale_cell_parameters()
        
    
    # ─────────────────────────────────────────────────────────────────────────
    # Private initialization methods (called from __post_init__)
    # ─────────────────────────────────────────────────────────────────────────
    
    def _scale_interaction_parameters(self) -> None:
        """Scale beta and theta_I by OMEGA for population-level dynamics."""
        object.__setattr__(self, 'beta', self.beta / self.OMEGA)
        object.__setattr__(self, 'theta_I', self.theta_I / self.OMEGA)
    
    def _scale_cell_parameters(self) -> None:
        """Scale parameters for each cell type.
        
        Processes ALL cell types consistently:
        - Always scales: lambda0 (with fitness gain), omega_exhaust
        - Conditionally scales: K (only if not None)
        
        This ensures scalability: if new parameters need scaling in the future,
        just add the logic to _create_scaled_cell_config().
        """
        scaled_params = {}
        for cell_type, config in self.cell_parameters.items():
            scaled_params[cell_type] = self._create_scaled_cell_config(config)
        
        object.__setattr__(self, 'cell_parameters', scaled_params)
    
    def _create_scaled_cell_config(self, config: CellTypeConfig) -> CellTypeConfig:
        """Create a new CellTypeConfig with scaled parameters.
        
        Scaling logic:
        - lambda0: Always scaled with fitness gain applied
        - omega_exhaust: Always scaled (preserves zero values)
        - K: Only scaled if not None (check happens here)
        
        This separation makes it easy to add new scalable parameters later.
        """
        lambda0 = config.lambda0 * (1 + config.fitness_gain)
        scaled_omega_exhaust = self._scale_omega_exhaust(config.omega_exhaust)
        scaled_K = self._scale_K(config.K) if config.K is not None else None
        
        return CellTypeConfig(
            default_id=config.default_id,
            N=config.N,
            K=scaled_K,
            lambda0=lambda0,
            mu=config.mu,
            nu=config.nu,
            omega_exhaust=scaled_omega_exhaust,
            next_mutation=config.next_mutation,
        )
    
    def _scale_K(self, K: float | int) -> int:
        """Scale carrying capacity by OMEGA.
        
        Called only when K is not None. Separated for clarity and potential
        future modifications to scaling logic.
        """
        return int(np.ceil(K * self.OMEGA))
    
    def _scale_omega_exhaust(self, omega_exhaust: float) -> float:
        """Scale omega_exhaust by OMEGA, preserving zero values.
        
        Always called, but handles the special case where omega_exhaust=0.0
        should remain 0.0 (not divide by OMEGA).
        """
        if omega_exhaust == 0.0:
            return omega_exhaust
        return omega_exhaust / self.OMEGA
    
    def _initialize_crowding_strategy(self) -> None:
        """Initialize crowding strategy based on configuration.
        
        Separated for clarity: strategy initialization is independent of
        parameter scaling and can be tested/modified separately.
        """
        from src.gillespie.crowding_strategy import SimpleCrowding, AdaptedCrowding
        
        strategy_class = AdaptedCrowding if self.use_logistic_adapted else SimpleCrowding
        object.__setattr__(self, 'crowding_strategy', strategy_class(self))
 

