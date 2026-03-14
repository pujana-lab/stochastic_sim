"""
mutation_engine_interface.py - Interfaz abstracta para motores de mutación.
"""
from abc import ABC, abstractmethod

from src.population import Population


class MutationEngineInterface(ABC):
    """Contrato abstracto para motores de mutación."""

    @abstractmethod
    def should_mutate(self, step: int) -> bool:
        """Determina si debe ocurrir una mutación en el paso actual.

        Args:
            step: Número del paso actual de la simulación.

        Returns:
            bool: True si debe ocurrir una mutación.
        """
        ...

    @abstractmethod
    def get_mutation(self, population: Population) -> tuple[str, str, float]:
        """Devuelve los parámetros de la mutación.

        Args:
            population: Estado actual de la población.

        Returns:
            tuple: (victim_group, new_name, new_fitness)
        """
        ...
