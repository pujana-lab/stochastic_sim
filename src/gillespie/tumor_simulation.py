from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
from tqdm import tqdm
from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone, MUTANT_SEED_SIZE
from src.gillespie.tissue_state import TissueState
from src.gillespie.event import Event
from src.gillespie.event_type import EventType
from src.gillespie.rate_matrix import RateMatrix
from src.gillespie.crowding_strategy import CrowdingStrategy
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone_factory import CloneFactory
from src.gillespie.stepper import Stepper, make_stepper

class TumorSimulation:
    def __init__(self, config: SimulationConfig) -> None:
        self.config = config
        self.save_all_steps = config.save_all_steps
        self.save_interval = config.save_interval if config.save_interval > 0 else (1 if config.save_all_steps else 0)
        self.step_count = 0
        self.rng = np.random.default_rng(config.seed)
        self.clone_factory = CloneFactory(config)
        self.t = 0.0
        print("arrancando con parametros:")
        print(config)
        if self.save_interval == 0:
            print(f"Memory mode: minimal (only final state). Use save_all_steps=True to save all steps.")
        
        #aqui igual merece mas la pena guardarlo como array simplemente o eso o hacerlo por tipos pero en ese caso 
        # Initialize clones
        # clones_dict:Dict[CloneType,Clone]= {}
        # TENGO QUE MOVER EL TIPO DE CLON DE LA FACTORY A LA CLASE CLONE
        # One root clone per cell type, keyed by its default_id (see _DEFAULT_CELL_PARAMETERS)
        clones_dict: Dict[CloneId, Clone] = {}
        for clone_type in self.config.cell_parameters:
            clone = self.clone_factory.create_clone(clone_type=clone_type)
            clones_dict[clone.clone_id] = clone
        
        self.times: List[float] = [0.0]
        # Encapsulate tissue state
        self.tissue_state: TissueState = TissueState(t= self.t, clones=clones_dict)

        # Use crowding strategy from config (initialized in SimulationConfig.__post_init__)
        self.crowding_strategy: CrowdingStrategy = self.config.crowding_strategy

        # How the simulation advances (exact SSA, later tau leaping), chosen by name in the config
        self.stepper: Stepper = make_stepper(self.config.stepper)
        # Rates of the current state; rebuilt after every event and reused by the next step
        self._rate_matrix: Optional[RateMatrix] = self._build_rate_matrix()
        rates0 = self._rate_matrix
        self.history: List[Dict[CloneId, dict]] = [self.tissue_state.snapshot()] if self.save_interval > 0 else []
        self.rate_history: List[List[Dict]] = [] if self.save_interval > 0 else []
        self.events: List[Optional[Event]] = [] if self.save_all_steps else []
        print("STARTING RATES AND STATE")
        self.tissue_state.print_pop_map()
        self.print_event_table(rates0.events)


    #TODO: esto huele feo
    def create_clone(
        self,
        clone_id: CloneId,
        clone_type: str = "base",
        N: int = None,
        parent: Optional[Clone] = None,
    ) -> Clone:
        clone = self.clone_factory.create_clone(
            clone_id=clone_id,
            clone_type=clone_type,
            N=N,
            parent=parent,
        )
        self.tissue_state.clones[clone.clone_id] = clone
        self.tissue_state.update_pop_map()
        self._rate_matrix = None
        return clone

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _advance_all_instability(self, dt: float) -> None:
        for clone in self.tissue_state.clones.values():
            if clone.is_alive():
                clone.advance_instability(dt, self.config.base_instability_buildup)

    def _current_rate_matrix(self) -> RateMatrix:
        if self._rate_matrix is None:
            self._rate_matrix = self._build_rate_matrix()
        return self._rate_matrix

    def _build_rate_matrix(self) -> RateMatrix:
        """Build per-clone events and store the effective rates in the TissueState.

        Each clone's rate is proportional to its own N (with type-level
        interactions from pop_map), so summing the events of a type gives the
        type-level rate without double counting.
        """
        self.tissue_state.update_pop_map()
        rate_matrix = RateMatrix()
        effective_rates: Dict[CloneId, Dict[str, float]] = {}

        for cid, clone in self.tissue_state.clones.items():
            if not clone.is_alive():
                continue

            rb = clone.birth_rate_effective(tissue_state=self.tissue_state)
            rd = clone.death_rate_effective(tissue_state=self.tissue_state)
            rm = clone.mutation_rate_effective(tissue_state=self.tissue_state)
            re = clone.conversion_rate_effective(tissue_state=self.tissue_state)
            effective_rates[cid] = {"rb": rb, "rd": rd, "rm": rm, "re": re}

            if rb > 0:
                rate_matrix.add_event(
                    Event(kind=EventType.BIRTH, clone_id=cid, rate=rb, clone_type=clone.get_type())
                )
            if rd > 0:
                rate_matrix.add_event(
                    Event(kind=EventType.DEATH, clone_id=cid, rate=rd, clone_type=clone.get_type())
                )
            if rm > 0:
                rate_matrix.add_event(
                    Event(kind=EventType.MUTATION, clone_id=cid, rate=rm, clone_type=clone.get_type())
                )
            if re > 0:
                rate_matrix.add_event(
                    Event(kind=EventType.CONVERSION, clone_id=cid, rate=re, clone_type=clone.get_type())
                )
        self.tissue_state.effective_rates = effective_rates
        return rate_matrix

    # ── Helpers shared by every stepper (see stepper.py) ──────────────────────
    # Tags: [DYNAMICS] changes the state or the next step; [REPORTING] only records
    # for the analysis (the trajectory is the same without it); [ONLY IF …] condition.

    def _record_rates(self, rate_matrix: RateMatrix) -> None:
        """Rate history: the rates at the start of a step (only when saving history).

        Pure [REPORTING]: read by save_rates_history_parquet, by nothing in the dynamics.
        """
        if self.save_interval > 0:   # [ONLY IF save_interval > 0]
            self.rate_history.append([
                {
                    "time": self.t,
                    "kind": e.kind.name,
                    "clone_id": e.clone_id,
                    "clone_type": e.clone_type,
                    "rate": e.rate,
                } for e in rate_matrix.events
            ])

    def _record(self, event: Optional[Event] = None, final: bool = False) -> None:
        """Record a step once the stepper has set self.t and changed the clones.

        Appends the time, rebuilds the rates of the new state (used by the snapshot
        and by the next step), then saves the snapshot (every save_interval steps,
        always when `final`) and the event (when save_all_steps). Shared by every
        stepper, so the history has the same layout whatever the step method.

        The first block keeps the state consistent and must run; the rest only
        records. They share one function because the snapshot needs the rates of
        the new state, so the rebuild has to come first.
        """
        # [REPORTING] time axis returned by run() (grows every step, even when nothing else is saved)
        self.times.append(self.t)

        # [DYNAMICS] keep the state consistent for the next step. Order matters:
        # t first, because with `decay` the crowding K_t = actual_K − decline·t uses it.
        self.tissue_state.t = self.t
        # [DYNAMICS] pop_map and rates of the new state; the next step reuses this matrix.
        # (It also fills tissue_state.effective_rates, which only the snapshot reads.)
        self._rate_matrix = self._build_rate_matrix()

        # [REPORTING] [ONLY IF final] last state at T_max: always saved, no event to save
        if final:
            self.history.append(self.tissue_state.snapshot())
            print("FINAL RATES AND STATE")
            self.tissue_state.print_pop_map()
            return
        # [REPORTING] [ONLY IF save_interval > 0 and this step is a multiple of it] snapshot
        if self.save_interval > 0 and self.step_count % self.save_interval == 0:
            self.history.append(self.tissue_state.snapshot())
        # [REPORTING] [ONLY IF save_all_steps] event list (debug parquet, --save-debug)
        if self.save_all_steps:
            self.events.append(event)

    # TO-DO: Add tests for mutation and exhaustion events, and for the overall step logic.
    def _introduce_mutation(self, clone: Clone, k: int = 1) -> None:
        """k mutations: k new clones of MUTANT_SEED_SIZE cells each, k parent cells removed."""
        assert clone.is_alive(), "Cannot mutate a dead clone."
        
        # Validate that the next mutation type exists in CloneType
        assert clone.next_mutation != "", "Clone is not supposed to be mutating"
        
        for _ in range(k):
            self.create_clone(
                clone_id=clone.next_child_id(),
                clone_type=clone.next_mutation,
                N=MUTANT_SEED_SIZE,
                parent=clone.clone_id,
            )
        clone.kill(k)
    
    def _induce_conversion(self, clone: Clone, k: int = 1) -> None:
        """Generic conversion: k cells of `clone` become k cells of `clone.converts_to`."""
        if not clone.converts_to:
            raise ValueError(f"{clone.get_type()} has no conversion target (converts_to)")
        target_id = self.config.cell_parameters[clone.converts_to].default_id
        clone.kill(k)
        self.tissue_state.clones[target_id].divide(k)
    
    def _apply_event(self, event: Event, k: int = 1) -> None:
        """Fire `event` k times at once (k = 1 for the SSA; a leap passes its count).

        Refuses to consume more parent cells than the clone has, so a bulk
        conversion or mutation never creates cells from nothing.
        """
        if k < 0:
            raise ValueError(f"k must be >= 0, got {k}")
        if k == 0:
            return
        clone = self.tissue_state.clones[event.clone_id]
        if event.kind.consumes_parent and k > clone.N:
            raise ValueError(
                f"{event.kind.name} x{k} on clone {event.clone_id} with only N = {clone.N} cells"
            )
        if event.kind == EventType.BIRTH:
            clone.divide(k)
        elif event.kind == EventType.DEATH:
            clone.kill(k)
        elif event.kind == EventType.MUTATION:
            self._introduce_mutation(clone, k)
        elif event.kind == EventType.CONVERSION:
            self._induce_conversion(clone, k)
        else:
            raise ValueError(f"Unknown event kind: {event.kind}")
        

    def print_event_table(self,events: List[Event]):
        """
        Imprime una tabla formateada con los eventos de la matriz de tasas.
        """
        # Cabecera de la tabla
        header = f"{'Kind':<12}  | {'Type':<12}| {'Clone ID':<15} | {'Rate':<10} | {'N':<15}"
        print(header)
        print("-" * len(header))
        
        # Filas
        for e in events:
            # e.clone_id suele ser una tupla, la convertimos a string para que quepa bien
            id_str = str(e.clone_id)
            print(f"{e.kind.name:<12} | {e.clone_type:<12} | {id_str:<15} | {e.rate:<10.4f}| {self.tissue_state.clones[e.clone_id].N}")
        
        print("-" * len(header))

    # ── Public API ────────────────────────────────────────────────────────────

    def total_population(self) -> int:
        return self.tissue_state.total_population()

    def step(self, return_matrix: Optional[bool] = False) -> bool:
        """Advance by one step of self.stepper. Returns False when the simulation should stop.

        The step itself (dynamics and recording) is in the stepper, see stepper.py.
        """
        # [REPORTING] [ONLY IF return_matrix] debug print of the rates before the step
        if return_matrix == True:
            print(self.t)
            self.print_event_table(events=self._current_rate_matrix().events)
        return self.stepper.step(self)

    def _stopping_cond(self) -> bool:
        """" We will make the system stop if it reaches stability"""

     

        # 2. Extinción total
        if self.total_population() <= 0:
            
            print("total extinction")
            
            return True
                
        
        # if self.tissue_state.pop_map["immune"] <= 0:

        #     print(" Complete tumor escape  ")
        #     return True
        
        # if self.tissue_state.pop_map["mutated"] <=0 and self.tissue_state.pop_map["exhausted"]<=0:
        #     print (" complete tumor control")
        #     return True
            
        
        #
       
        return False

    def run(self) -> Tuple[List[float], List[Dict[CloneId, dict]], TissueState]:
        with tqdm(desc="Simulating", unit=" steps") as pbar:
            while not self._stopping_cond():
                if not self.step():
                    break
                pbar.update(1)
                pbar.set_postfix({"time": f"{self.t:.2f}/{self.config.T_max:.2f}"})
        return self.times, self.history, self.tissue_state, self.rate_history
