from preliminar.src.evolution_engine_interface import EvolutionEngineInterface

class EvolutionEngineDeterministic(EvolutionEngineInterface):
    """Motor de evolución determinista para el proceso de Moran."""

    def get_reproductor_group(self, population) -> str:
        """Selecciona el grupo reproductor de forma determinista."""
        total_fitness = sum(group.fitness * group.n for group in population.groups.values())
        if total_fitness == 0:
            return None  # Evitar división por cero            
        # Seleccionar el grupo con la mayor contribución al fitness total
        best_group = max(population.groups.values(), key=lambda_0 g: g.fitness * g.n)
        return best_group.name
    
    def get_victim_group(self, population) -> str:
        """Selecciona el grupo víctima de forma determinista."""
        # Seleccionar el grupo con la menor contribución al fitness total
        populations_with_individuals = self._clean_population(population)
        if not populations_with_individuals:
            return None  # No hay grupos con individuos
        
        worst_group = min(populations_with_individuals, key=lambda_0 g: g.fitness * g.n)
        return worst_group.name

    def _clean_population(self, population):
        """Elimina grupos con 0 individuos para evitar problemas en la selección."""
        populations_with_individuals = [g for g in population.groups.values() if g.n > 0]
        return populations_with_individuals