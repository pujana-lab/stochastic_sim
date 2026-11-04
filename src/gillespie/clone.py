from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING, Dict, Type

from src.gillespie.cloneId import CloneId
from src.gillespie.event_type import EventType
from src.gillespie.simulation_config import SimulationConfig

if TYPE_CHECKING:
    from src.gillespie.tissue_state import TissueState

# Cells in a new mutant clone: a MUTATION removes 1 parent cell and adds a clone of
# this size. Shared by TumorSimulation._introduce_mutation and Clone.type_delta.
MUTANT_SEED_SIZE: int = 2


class Clone:

    _registry: Dict[str, Type["Clone"]] = {}
    _type_name: str = None

    def __init_subclass__(cls, clone_type: str = None, **kwargs):
        super().__init_subclass__(**kwargs)
        if clone_type:
            Clone._registry[clone_type] = cls
            cls._type_name = clone_type

    def __init__(self, config: SimulationConfig, clone_id: CloneId = None, N: int = None, parent: Optional[CloneId] = None):
        
        #TODO: todo lo que dependa de la clase y no de la instancia deberia ir guardado dentro de la clase (sacar todos los self. de cosas que no sean dinamicas) y aunque lo sean deberian ir fuera igual y actualizarse para todo el tipo. lo unico que deberia ir dentro son las coasas que dependan del tiempo Y ADEMAS queramos trackear por clon y no por tipo (inestabilidad genomica en un futuro)
        
        self.config: SimulationConfig = config
        self.cell_parameters = config.cell_parameters[self.get_type()]
        self.clone_id: CloneId = self.cell_parameters.default_id if clone_id is None else clone_id
        self.N: int = self.cell_parameters.N if N is None else N
        self.parent: CloneId = parent
        self.children_count: int = 0
        self.crowding_value: float = 0.0
        
        self.birth_rate: float = self.cell_parameters.lambda0
        self.death_rate: float = self.cell_parameters.mu
        self.mutation_rate: float = self.cell_parameters.nu
        self.exhaustion_rate: float = self.cell_parameters.omega_exhaust

        self.K: int = self.cell_parameters.K
        self.K_min: int = config.Kmin

        self.next_mutation: str = self.cell_parameters.next_mutation
        self.converts_to: str = self.cell_parameters.converts_to
        self.instability: float = 0.0
        self.buildup: float = 0.0
        self.d1: float = 0.0
        self.d2: float = 0.0

        self.actual_K = config.crowding_strategy.calculate_K(self)
        

    def get_type(self) -> str:
        return self._type_name
    
    def crowding_numerator(self, tissue_state: "TissueState") -> float:
        """Y_X = Σ_T w_T · N_T over this type's crowding_pool (CellTypeConfig)."""
        pop_map = tissue_state.pop_map
        return sum(w * pop_map.get(t, 0) for t, w in self.cell_parameters.crowding_pool.items())

    def is_alive(self) -> bool:
        return self.N > 0

    def mutation_multiplier(self) -> float:
        return 1.0 + self.instability
    def birth_rate_effective(self,tissue_state: "TissueState") -> float:
        return self.birth_rate * self.N * self.config.crowding_strategy.crowding(self,tissue_state=tissue_state) #NO ESTABA MULTIPLICANDO POR EL BIRTH RATE!!!

    def death_rate_effective(self,tissue_state: "TissueState") -> float:
        return self.death_rate * self.N

    def mutation_rate_effective(self,tissue_state: "TissueState") -> float:
        return self.mutation_rate * self.N * self.mutation_multiplier()

    def conversion_rate_effective(self,tissue_state: "TissueState") -> float:
        return self.exhaustion_rate * self.N
        
    def type_delta(self, kind: EventType) -> Dict[str, int]:
        """Change of the type totals N_T caused by one event of this clone (mirrors _apply_event)."""
        t = self.get_type()
        if kind == EventType.BIRTH:
            return {t: 1}
        elif kind == EventType.DEATH:
            return {t: -1}
        elif kind == EventType.MUTATION:
            if not self.next_mutation:
                raise ValueError(f"{t} has no mutation target (next_mutation)")
            target, gain = self.next_mutation, MUTANT_SEED_SIZE
        elif kind == EventType.CONVERSION:
            if not self.converts_to:
                raise ValueError(f"{t} has no conversion target (converts_to)")
            target, gain = self.converts_to, 1
        else:
            raise ValueError(f"Unknown event kind: {kind}")

        delta = {t: -1}
        delta[target] = delta.get(target, 0) + gain   # adds up if target == t
        return {k: v for k, v in delta.items() if v != 0}

    def divide(self) -> None:
        self.N += 1

    def kill(self) -> None:
        self.N = max(0, self.N - 1)

    def next_child_id(self) -> CloneId:
        self.children_count += 1
        return self.clone_id + (self.children_count,)

    def advance_instability(self, dt: float, base_buildup: float) -> None:
        self.instability += (base_buildup + self.buildup) * dt       
    def __str__(self) -> str:
        return self.get_type()


# de nuevo aqui hay que pasarle el diccionario de valores N_loquesea y desde ahi pasarle el crowdingstrategy que elijamos y calcule el numero

class WildTypeClone(Clone, clone_type = "base"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class MutatedClone(Clone, clone_type = "mutated"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def birth_rate_effective(self, tissue_state: "TissueState") -> float:
        """r_B(mutated) = λ_mutated · N_mutated · g(Y_mutated, K_mutated) · (1 + a_F · N_caf)   (a_F already / OMEGA)"""
        base_birth = super().birth_rate_effective(tissue_state)
        n_caf = tissue_state.pop_map.get("caf", 0)
        return base_birth * (1.0 + self.config.alpha_F * n_caf)
    
    def death_rate_effective(self, tissue_state: "TissueState") -> float:
        """Death rate increases due to immune cell killing.
        Formula: base_death_rate * N + N * N_immune * theta_I
        (theta_I is a global parameter from config)
        """
        base_death = super().death_rate_effective(tissue_state)
        self_n = self.N
        n_immune = tissue_state.pop_map.get("immune", 0)
        immune_killing = self_n * n_immune * self.config.theta_I
        return base_death + immune_killing

    

class ImmuneClone(Clone, clone_type = "immune"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def birth_rate_effective(self, tissue_state: "TissueState") -> float:
        """Birth rate increases through activation by cancer cells.
        Formula: base_birth_rate + N * N_cancer * beta
        (beta is a global parameter from config)
        """
        base_birth = super().birth_rate_effective(tissue_state)
        n_self = self.N
        n_cancer = tissue_state.pop_map.get("mutated", 0)
        activation_boost = n_self * n_cancer * self.config.beta
        # crowding_effect= #crowding donde tenemos que el numerador de la logistica es n_exhausted+sel
        return base_birth + activation_boost
    
    def conversion_rate_effective(self, tissue_state: "TissueState") -> float:
        """r_E(immune) = ω · N_immune · N_mutated + ω_F · N_immune · N_caf
(omega_F comes already divided by OMEGA from SimulationConfig when scale=True)"""
        base_exhaustion =  super().conversion_rate_effective(tissue_state= tissue_state)
        n_cancer = tissue_state.pop_map.get("mutated",0)
        n_caf = tissue_state.pop_map.get("caf", 0)
        caf_exhaustion = self.N * n_caf * self.config.omega_F
        return base_exhaustion * n_cancer + caf_exhaustion


class ExhaustedClone(Clone, clone_type = "exhausted"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def birth_rate_effective(self, tissue_state: "TissueState") -> float:
        """Exhausted cells cannot divide."""
        return 0.0


class FibroblastClone(Clone, clone_type = "fibroblast"):
    """F_N: logistic growth with Y_fibroblast = N_fibroblast + N_caf + N_mutated; activated by mutated."""

    def conversion_rate_effective(self, tissue_state: "TissueState") -> float:
        """Activation fibroblast -> caf (Michaelis-Menten in N_mutated):
        r_E(fibroblast) = κ · N_fibroblast · N_mutated / (C_½ + N_mutated)   (C_½ already x OMEGA)
        """
        n_cancer = tissue_state.pop_map.get("mutated", 0)
        denom = self.config.C_half + n_cancer
        if denom <= 0:
            return 0.0
        return self.config.kappa * self.N * n_cancer / denom


class CAFClone(Clone, clone_type = "caf"):
    """F_T: does not divide, r_B(caf) = 0; r_D(caf) = μ_F · N_caf."""

    def birth_rate_effective(self, tissue_state: "TissueState") -> float:
        return 0.0
