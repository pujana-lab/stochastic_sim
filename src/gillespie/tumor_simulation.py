from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from numba import optional
import numpy as np
from tqdm import tqdm
from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone
from src.gillespie.tissue_state import TissueState
from src.gillespie.event import Event
from src.gillespie.event_type import EventType
from src.gillespie.rate_matrix import RateMatrix
from src.gillespie.crowding_strategy import CrowdingStrategy
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone_factory import CloneFactory

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
        clones_dict: Dict[CloneId, Clone] = {
            (): self.clone_factory.create_clone(clone_type="base"),
            (-3,): self.clone_factory.create_clone(clone_type="mutated"),
            (-1,): self.clone_factory.create_clone(clone_type="immune"),
            (-2,): self.clone_factory.create_clone(clone_type="exhausted")
        }
        
        self.times: List[float] = [0.0]
        # Encapsulate tissue state
        self.tissue_state: TissueState = TissueState(t= self.t, clones=clones_dict)

        # Use crowding strategy from config (initialized in SimulationConfig.__post_init__)
        self.crowding_strategy: CrowdingStrategy = self.config.crowding_strategy

        #aqui habria que anyadir lo mismo para elegir strategy pero para el tipo de leap. (Binomial, Poisson, Poisson half etc)
        self.history: List[Dict[CloneId, dict]] = [self.tissue_state.snapshot()] if self.save_interval > 0 else []
        self.rate_history: List[List[Dict]] = [] if self.save_interval > 0 else []
        self.events: List[Optional[Event]] = [] if self.save_all_steps else []
        rates0=self._build_rate_matrix()
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
        return clone

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _advance_all_instability(self, dt: float) -> None:
        for clone in self.tissue_state.clones.values():
            if clone.is_alive():
                clone.advance_instability(dt, self.config.base_instability_buildup)

    def _build_rate_matrix(self) -> RateMatrix:
        self.tissue_state.update_pop_map()
        rate_matrix = RateMatrix()
        type_rates: Dict[str, tuple[float, float, float, float]] = {}

        for cid, clone in self.tissue_state.clones.items():
            if not clone.is_alive():
                continue
            
            #BUG: aqui estamos mezclando dos logicas. Por una parte estamos asignando por clon y por otra parte por tipo. si hacemos por tipo deberiamos iterar por tipo si no por clon porque si no clones del mismo tipo cuentan varias veces sus rates. 
            # Hay que decidir que hacemos con esto. si hacemos que todos los clones tengan la misma dinamica entonces no sabemos que clon estamos mutando y si diferenciamos tenemos que volver a calcular todo por cada clon. una opcion es hacer un refacotr y en dos partes primero calcular los rates basales y de cada fenotipo y luego ir por cada clon y multiplicarlos por el numero de clones que tenemos. esta logica se podria guardar en tissue stste de forma que los rates basales se calculen y luego el rate matrix pulee de ahi. 
            #TODO: hay que volver a poner los rates por TIPO y simplemente a la hora de aplicar el evento tirar moneda para elegir cual clon de ese tipo prolifera/muere
            
            #TODO: est hay que meterlo a tissue_state para poder pintar en condiciones los rates REALES
            type_rates: tuple[float,float,float,float]= (
                clone.birth_rate_effective(tissue_state=self.tissue_state),
                clone.death_rate_effective(tissue_state=self.tissue_state),
                clone.mutation_rate_effective(tissue_state=self.tissue_state),
                clone.exhaustion_rate_effective(tissue_state=self.tissue_state),
            )

            rb, rd, rm, re = type_rates

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
                    Event(kind=EventType.EXHAUSTION, clone_id=cid, rate=re, clone_type=clone.get_type())
                )
        return rate_matrix

    def _sample_waiting_time(self, total_rate: float) -> float:
        if total_rate <= 0:
            raise ValueError("Total rate must be positive.")
        return -np.log(self.rng.random()) / total_rate

    # TO-DO: Add tests for mutation and exhaustion events, and for the overall step logic.
    def _introduce_mutation(self, clone: Clone) -> None:
        assert clone.is_alive(), "Cannot mutate a dead clone."
        
        # Validate that the next mutation type exists in CloneType
        assert clone.next_mutation != "", "Clone is not supposed to be mutating"
        
        
        self.create_clone(
            clone_id=clone.next_child_id(),
            clone_type=clone.next_mutation,
            N=2,
            parent=clone.clone_id,
        )
        clone.kill()
    
    def _induce_exhaustion(self, clone: Clone) -> None:
        self.tissue_state.clones[clone.clone_id].kill()
        self.tissue_state.clones[(-2,)].divide()
    
    def _apply_event(self, event: Event) -> None:
        clone = self.tissue_state.clones[event.clone_id]
        if event.kind == EventType.BIRTH:
            clone.divide()
        elif event.kind == EventType.DEATH:
            clone.kill()
        elif event.kind == EventType.MUTATION:
            self._introduce_mutation(clone)
        elif event.kind == EventType.EXHAUSTION:
            self._induce_exhaustion(clone)
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
        """Advance by one Gillespie step. Returns False when simulation should stop."""
        self.step_count += 1
        rate_matrix = self._build_rate_matrix()
  
        if return_matrix == True:
            print(self.t)
            self.print_event_table(events = rate_matrix.events)
        total_rate = rate_matrix.get_total_rate()
        
        # Only record rates if saving history
        if self.save_interval > 0:
            step_rates = [
                {
                    "time": self.t,
                    "kind": e.kind.name,
                    "clone_id": e.clone_id,
                    "clone_type": e.clone_type,
                    "rate": e.rate,
                } for e in rate_matrix.events
            ]
            self.rate_history.append(step_rates)
        
        if total_rate <= 0 or not rate_matrix.events:
            return False

        tau = self._sample_waiting_time(total_rate)
        new_t = self.t + tau

        if new_t > self.config.T_max:
            tau = self.config.T_max - self.t
            self._advance_all_instability(tau)
            self.t = self.config.T_max
            self.times.append(self.t)
            # Always save final state
            self.history.append(self.tissue_state.snapshot())
            print("FINAL RATES AND STATE")
            self.tissue_state.print_pop_map()
            # self.print_event_table(events=rate_matrix.events)
            
            return False

        self._advance_all_instability(tau)
        self.t = new_t
        
        event = rate_matrix.choose_event(self.rng.random())
        self._apply_event(event)
        self.tissue_state.update_pop_map()
        
        self.times.append(self.t)
        self.tissue_state.t = self.t
        
        # Save snapshots conditionally
        if self.save_interval > 0 and self.step_count % self.save_interval == 0:
            self.history.append(self.tissue_state.snapshot())
        
        # Save events if recording all steps
        if self.save_all_steps:
            self.events.append(event)
                
        return True
        
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
