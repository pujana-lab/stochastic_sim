"""
mutation_engine_deterministic.py - Motor de mutación determinista.

La mutación ocurre cada `every_n_steps` pasos, siempre sobre el mismo
grupo víctima y generando siempre el mismo grupo mutante.
"""
from src.mutation_engine_interface import MutationEngineInterface
from src.population import Population


class MutationEngineDisabled(MutationEngineInterface):
    """Motor de mutación nulo para el proceso de Moran.

        
    Example:
        >>> mutation = MutationEngineDeterministic(
        ...     every_n_steps=10,
        ...     victim_group='weak',
        ...     new_group_name='mutant',
        ...     new_fitness=1.5,
        ... )
        >>> mutation.should_mutate(step=10)
        False
        >>> mutation.should_mutate(step=7)
        False
    """




    def should_mutate(self, step: int) -> bool:
        """Devuelve True si el paso actual es múltiplo de every_n_steps.

        Args:
            step: Número del paso actual (1-indexed).

        Returns:
            bool: False, ya que este motor no muta realmente.
        """
        return False

    def get_mutation(self, population: Population) -> tuple[str, str, float]:
        """Devuelve los parámetros fijos de la mutación.

        Args:
            population: Estado actual de la población (no usado en modo determinista).

        Returns:
            tuple: (victim_group, new_group_name, new_fitness)
        """
        raise Exception("get_mutation no implementado en MutationEngineDeterministic")
