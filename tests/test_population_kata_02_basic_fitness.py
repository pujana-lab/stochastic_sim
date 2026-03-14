"""
Suite de tests para la clase Population.

Tests para la simulación de procesos evolutivos de Moran.
Kata 02: Fitness básico.
"""

from src.population import Population


class TestPopulationFitness:
    """Tests para la inicialización del fitness."""

    def test_population_fitness_by_default(self):
        """Test: la aptitud por defecto es 1.0."""
        groups = {'group1': {'n': 10}}
        pop = Population(groups=groups)
        assert pop.fitness['group1'] == 1.0

    def test_population_fitness_no_default(self):
        """Test: el fitness se puede establecer al crear la población."""
        groups = {'group1': {'n': 10, 'fitness': 2.0}}
        pop = Population(groups=groups)
        assert pop.fitness['group1'] == 2.0