"""
mutation_engine_deterministic.py - Motor de mutación determinista.

La mutación ocurre cada `every_n_steps` pasos, siempre sobre el mismo
grupo víctima y generando siempre el mismo grupo mutante.
"""
from src.mutation_engine_interface import MutationEngineInterface
from src.population import Population


class MutationEngineDeterministic(MutationEngineInterface):
    """Motor de mutación determinista para el proceso de Moran.

    La mutación ocurre en intervalos fijos de pasos. En cada mutación:
    - Se elimina un individuo del grupo víctima.
    - Se añade un individuo del grupo mutante (nuevo o existente).

    Args:
        every_n_steps: Cada cuántos pasos ocurre una mutación.
        victim_group: Nombre del grupo del que se extrae el individuo mutado.
        new_group_name: Nombre del nuevo grupo mutante.
        new_fitness: Fitness del nuevo grupo mutante.

    Example:
        >>> mutation = MutationEngineDeterministic(
        ...     every_n_steps=10,
        ...     victim_group='weak',
        ...     new_group_name='mutant',
        ...     new_fitness=1.5,
        ... )
        >>> mutation.should_mutate(step=10)
        True
        >>> mutation.should_mutate(step=7)
        False
    """

    def __init__(
        self,
        every_n_steps: int,
        victim_group: str,
        new_group_name: str,
        new_fitness: float,
    ):
        self.every_n_steps = every_n_steps
        self.victim_group = victim_group
        self.new_group_name = new_group_name
        self.new_fitness = new_fitness

    def should_mutate(self, step: int) -> bool:
        """Devuelve True si el paso actual es múltiplo de every_n_steps.

        Args:
            step: Número del paso actual (1-indexed).

        Returns:
            bool: True si debe ocurrir una mutación.
        """
        return step > 0 and step % self.every_n_steps == 0

    def get_mutation(self, population: Population) -> tuple[str, str, float]:
        """Devuelve los parámetros fijos de la mutación.

        Args:
            population: Estado actual de la población (no usado en modo determinista).

        Returns:
            tuple: (victim_group, new_group_name, new_fitness)
        """
        return self.victim_group, self.new_group_name, self.new_fitness
