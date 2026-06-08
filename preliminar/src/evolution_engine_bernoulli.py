import random
from preliminar.src.evolution_engine_interface import EvolutionEngineInterface

class EvolutionEngineBernoulli(EvolutionEngineInterface):
    """Motor de evolución no determinista basado en un proceso de Bernoulli para el proceso de Moran."""

    def __init__(self, seed=None):
        """Inicializa el motor de evolución con una semilla opcional para reproducibilidad."""
        self.random = random.Random(seed)
        self.seed = seed

    def get_reproductor_group(self, population) -> str:
        """Selecciona el grupo reproductor de forma no determinista."""
        total_fitness = sum(group.fitness * group.n for group in population.groups.values())
        if total_fitness == 0:
            return None  # Evitar división por cero
        # Seleccionar un grupo basado en la probabilidad proporcional a su contribución al fitness total
        choices = [(group.name, group.fitness * group.n / total_fitness) for group in population.groups.values()]
        groups, probabilities = zip(*choices)
        best_group = self.random.choices(groups, probabilities)[0]
        return best_group
    
    def get_victim_group(self, population) -> str:
        """Selecciona el grupo víctima de forma no determinista."""
        total_fitness = sum(group.fitness * group.n for group in population.groups.values())
        if total_fitness == 0:
            return None  # Evitar división por cero
        # Seleccionar un grupo basado en la probabilidad inversa a su contribución al fitness total
        choices = [(group.name, 1 - (group.fitness * group.n / total_fitness)) for group in population.groups.values()]
        groups, probabilities = zip(*choices)
        worst_group = self.random.choices(groups, probabilities)[0]
        return worst_group