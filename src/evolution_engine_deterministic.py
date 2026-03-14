from src.evolution_engine_interface import EvolutionEngineInterface

class EvolutionEngineDeterministic(EvolutionEngineInterface):
    """Motor de evolución determinista para el proceso de Moran."""

    def get_reproductor_group(self, population):
        """Selecciona el grupo reproductor de forma determinista."""
        total_fitness = sum(group.fitness * group.n for group in population.groups.values())
        if total_fitness == 0:
            return None  # Evitar división por cero
        # Seleccionar el grupo con la mayor contribución al fitness total
        best_group = max(population.groups.values(), key=lambda g: g.fitness * g.n)
        return best_group.name
    
    def get_victim_group(self, population):
        """Selecciona el grupo víctima de forma determinista."""
        # Seleccionar el grupo con la menor contribución al fitness total
        worst_group = min(population.groups.values(), key=lambda g: g.fitness * g.n)
        return worst_group.name