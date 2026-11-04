from dataclasses import dataclass, field, replace
from typing import Optional, TYPE_CHECKING
from src.gillespie.cloneId import CloneId
import numpy as np

if TYPE_CHECKING:
    from src.gillespie.crowding_strategy import CrowdingStrategy
#TODO: poder pasarle el CEllTypeConfig al cli y al json y tal porque ahora mismo esta desactualizado.
#TODO: implementar que si N es inf en el postinit define N como el valor de saturacion homeostatica (corresponde a la scaled_K)
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
    converts_to: str = ""   # target type of the CONVERSION/EXHAUSTION event ("" = none)
    # Crowding pool {T: w_T}: types that compete for this type's niche and their weight,
    # Y_X = Σ_T w_T · N_T. Empty = no crowding (g = 1). Weights are not scaled.
    crowding_pool: dict[str, float] = field(default_factory=dict, hash=False)

# default_id is the id of the root clone of each type: one per type, outside the
# mutation lineage (mutants get parent.clone_id + (k,)). Every other module reads
# the ids from here. The order of this dict is the order of the root clones in
# TumorSimulation, and therefore of the events in the rate matrix: changing it
# changes seeded trajectories.
_DEFAULT_CELL_PARAMETERS: dict[str, CellTypeConfig] = {
    "base": CellTypeConfig(
        default_id= (),
        N=4000,
        K=1,
        nu=0.00002,
        next_mutation= "mutated",
        crowding_pool={"base": 1, "mutated": 1},
    ),
    "mutated": CellTypeConfig(
        default_id= (-3,),
        fitness_gain= 0.2,
        N=0,
        K=2,
        crowding_pool={"mutated": 1},
    ),
    "immune": CellTypeConfig(
        default_id= (-1,),
        N=2000,
        K=0.5,
        lambda0= 0.005,
        omega_exhaust=0.003,
        mu= 0.0,
        converts_to="exhausted",
        crowding_pool={"immune": 1, "exhausted": 1},
    ),
    "exhausted": CellTypeConfig(
        default_id=(-2,),
        N=0,
        lambda0= 0.0
    ),
    # Fibroblasts are OFF by default (N=0): the 4-type model is reproduced exactly.
    "fibroblast": CellTypeConfig(       # F_N, naive / normal fibroblasts
        default_id=(-4,),
        N=0,
        K=1,                            # K_F (fraction of OMEGA)
        lambda0=0.005,                  # λ_F
        mu=0.0,
        converts_to="caf",              # activation F_N -> F_T (rate in FibroblastClone)
        crowding_pool={"fibroblast": 1, "caf": 1, "mutated": 1},
    ),
    "caf": CellTypeConfig(              # F_T, cancer-associated / activated fibroblasts
        default_id=(-5,),
        N=0,
        lambda0=0.0,                    # no proliferation
        mu=0.002,                       # μ_F
    ),
}

@dataclass(frozen=True)
class SimulationConfig:

    OMEGA: int = 4000 #NUMERO MAXIMO DE CELULAS WT QUE SOPORTA CUANDO NO HAY COMPETICION
    
    
    cell_parameters: dict[str, CellTypeConfig] = field(default_factory=lambda : dict(_DEFAULT_CELL_PARAMETERS))
 
    # simulation control
    T_max: float = 1000
    seed: Optional[int] = None

   
    decline: float = 0.0
    Kmin: float = 1

    

    # interaction parameters
    theta_I: float = 0.0005  # Kill rate: immune cells killing mutated cells (N_mutant * N_immune * theta_I)
    beta: float = 0.0004   # Activation rate: mutated cells activating immune cells (N_immune * N_mutant * beta)

    # fibroblast interaction parameters (all 0 -> fibroblast module inert)
    kappa: float = 0.0     # κ: max activation rate F_N -> F_T by mutated cells (Michaelis-Menten, NOT scaled)
    C_half: float = 0.1    # C_½: N_mutated giving half-maximal activation (fraction of OMEGA, x OMEGA)
    omega_F: float = 0.0   # ω_F: exhaustion immune -> exhausted driven by caf (N_immune * N_caf * omega_F)
    alpha_F: float = 0.0   # a_F: caf support of mutated birth, factor (1 + alpha_F * N_caf)
   

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
    save_all_steps: bool = True
    save_interval: int = 0
    
    # Crowding strategy (será inicializado en __post_init__)
    crowding_strategy: "CrowdingStrategy" = field(init=False, default=None)
    def __post_init__(self):
        """Initialize configuration: normalize cell params, scale parameters, and set up strategies."""
        self._normalize_cell_parameters()
        self._initialize_crowding_strategy()
        if self.scale:
            self._scale_interaction_parameters()
            self._scale_cell_parameters()

#TODO: entender esto 
    def _normalize_cell_parameters(self) -> None:
        normalized = dict(_DEFAULT_CELL_PARAMETERS)
        for cell_type, config in self.cell_parameters.items():
            default = _DEFAULT_CELL_PARAMETERS.get(cell_type)
            if isinstance(config, CellTypeConfig):
                cfg = config
                if cfg.default_id is None and default is not None:
                    cfg = replace(cfg, default_id=default.default_id)
            elif isinstance(config, dict):
                payload = dict(config)
                if default is not None:
                    if "crowding_pool" in payload:
                        # merge key by key: only the weights given here change
                        payload["crowding_pool"] = {**default.crowding_pool, **payload["crowding_pool"]}
                    cfg = replace(default, **payload)
                else:
                    cfg = CellTypeConfig(**payload)
            else:
                raise TypeError(...)
            normalized[cell_type] = cfg

        for cell_type, cfg in normalized.items():
            normalized[cell_type] = replace(cfg, crowding_pool=self._clean_crowding_pool(cell_type, cfg, normalized))

        ids = [cfg.default_id for cfg in normalized.values()]
        if None in ids or len(set(ids)) != len(ids):
            raise ValueError(
                "Every cell type needs its own default_id (root clone id): "
                f"{ {t: cfg.default_id for t, cfg in normalized.items()} }"
            )

        object.__setattr__(self, "cell_parameters", normalized)

    @staticmethod
    def _clean_crowding_pool(cell_type: str, cfg: CellTypeConfig, cell_parameters: dict) -> dict[str, float]:
        """Validate a crowding pool and drop zero weights (weight 0 = not in the pool)."""
        unknown = set(cfg.crowding_pool) - set(cell_parameters)
        if unknown:
            raise ValueError(f"crowding_pool of {cell_type} names unknown cell types: {sorted(unknown)}")
        negative = {t: w for t, w in cfg.crowding_pool.items() if w < 0}
        if negative:
            raise ValueError(f"crowding_pool weights of {cell_type} must be >= 0: {negative}")
        return {t: float(w) for t, w in cfg.crowding_pool.items() if w != 0}
    
    # ─────────────────────────────────────────────────────────────────────────
    # Private initialization methods (called from __post_init__)
    # ─────────────────────────────────────────────────────────────────────────
    
    def _scale_interaction_parameters(self) -> None:
        """Scale beta and theta_I by OMEGA for population-level dynamics."""
        object.__setattr__(self, 'beta', self.beta / self.OMEGA)
        object.__setattr__(self, 'theta_I', self.theta_I / self.OMEGA)
        object.__setattr__(self, 'omega_F', self.omega_F / self.OMEGA)
        object.__setattr__(self, 'alpha_F', self.alpha_F / self.OMEGA)
        object.__setattr__(self, 'C_half', self.C_half * self.OMEGA)
    
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
        
        # replace() keeps every unscaled field, so new fields are not lost here
        return replace(
            config,
            K=scaled_K,
            lambda0=lambda0,
            omega_exhaust=scaled_omega_exhaust,
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
 
