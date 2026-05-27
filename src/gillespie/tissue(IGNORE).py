from __future__ import annotations
from typing import Dict, List, Tuple
from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone

class Tissue:
    def __init__(self, initial_clones: Dict[CloneId, Clone] = None) -> None:
        self.clones: Dict[CloneId, Clone] = initial_clones if initial_clones is not None else {}
        
        # Caché interna opcional para evitar recorrer el diccionario en cada evento
        self._populations_by_type: Dict[str, int] = {}
        self._total_n: int = 0
        self._update_cache()

    def _update_cache(self) -> None:
        """Recalcula los agregados de población. Llamar tras mutaciones o cambios estructurales."""
        self._populations_by_type = {"wild_type": 0, "mutated": 0, "immune": 0, "exhausted": 0}
        self._total_n = 0
        for clone in self.clones.values():
            self._populations_by_type[clone.clone_type] = self._populations_by_type.get(clone.clone_type, 0) + clone.N
            self._total_n += clone.N

    def add_clone(self, clone: Clone) -> None:
        self.clones[clone.clone_id] = clone
        self._update_cache()

    def remove_clone(self, clone_id: CloneId) -> None:
        if clone_id in self.clones:
            del self.clones[clone_id]
            self._update_cache()

    def notify_population_change(self) -> None:
        """Llamar cuando clone.divide() o clone.kill() alteren N, sin cambiar el nº de clones."""
        self._update_cache()

    # ── Métricas del Sistema (API para los Clones y Simulación) ──────────────────

    @property
    def total_population(self) -> int:
        return self._total_n

    def get_population_of_type(self, clone_type: str) -> int:
        return self._populations_by_type.get(clone_type, 0)

    def get_tumor_population(self) -> int:
        """Ejemplo: N_W + N_C (WildType + Mutated)"""
        return (self.get_population_of_type("wild_type") + 
                self.get_population_of_type("mutated"))

    def get_immune_population(self) -> int:
        return self.get_population_of_type("immune")

    def get_exhausted_population(self) -> int:
        return self.get_population_of_type("exhausted")