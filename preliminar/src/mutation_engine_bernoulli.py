"""
mutation_engine_bernoulli.py - Motor de mutación probabilístico.

La mutación ocurre en cada paso con probabilidad `p` siguiendo una
distribución de Bernoulli.
"""
import random

from preliminar.src.mutation_engine_interface import MutationEngineInterface
from preliminar.src.population import Population


class MutationEngineBernoulli(MutationEngineInterface):
    """Motor de mutación probabilístico para el proceso de Moran.

    En cada paso, la mutación ocurre con probabilidad `p`. En cada mutación:
    - Se elimina un individuo del grupo víctima.
    - Se añade un individuo del grupo mutante (nuevo o existente).

    Args:
        p: Probabilidad de mutación en cada paso. Debe estar en [0, 1].
        victim_group: Nombre del grupo del que se extrae el individuo mutado.
        new_group_name: Nombre del nuevo grupo mutante.
        new_fitness: Fitness del nuevo grupo mutante.
        seed: Semilla para reproducibilidad (opcional).

    Example:
        >>> mutation = MutationEngineBernoulli(
        ...     p=0.1,
        ...     victim_group='weak',
        ...     new_group_name='mutant',
        ...     new_fitness=1.5,
        ...     seed=42,
        ... )
        >>> mutation.get_mutation(population)
        ('weak', 'mutant', 1.5)
    """

    def __init__(
        self,
        p: float,
        victim_group: str,
        new_group_name: str,
        new_fitness: float,
        seed: int | None = None,
    ):
        if not 0.0 <= p <= 1.0:
            raise ValueError(f"La probabilidad p debe estar en [0, 1], recibido: {p}")
        self.p = p
        self.victim_group = victim_group
        self.new_group_name = new_group_name
        self.new_fitness = new_fitness
        self._rng = random.Random(seed)

    def should_mutate(self, step: int) -> bool:
        """Devuelve True con probabilidad p en cada paso.

        Args:
            step: Número del paso actual (no usado, la decisión es independiente).

        Returns:
            bool: True si debe ocurrir una mutación en este paso.
        """
        return self._rng.random() < self.p

    def get_mutation(self, population: Population) -> tuple[str, str, float]:
        """Devuelve los parámetros de la mutación.

        Args:
            population: Estado actual de la población.

        Returns:
            tuple: (victim_group, new_group_name, new_fitness)
        """
        return self.victim_group, self.new_group_name, self.new_fitness
