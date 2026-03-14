"""
Suite de tests para la clase Population.

Tests para la simulación de procesos evolutivos de Moran.
Kata 02: Fitness básico.
"""

from src.population import Population
from src.population import SubPopulation

class TestPopulationFitness:
    """Tests para la inicialización del fitness."""

    def test_population_fitness_by_default(self):
        """Test: la aptitud por defecto es 1.0."""
        groups = {'group1': SubPopulation(name='group1', n=10)}
        pop = Population(groups=groups)
        assert pop.fitness['group1'] == 1.0

    def test_population_fitness_no_default(self):
        """Test: el fitness se puede establecer al crear la población."""
        groups = {'group1': SubPopulation(name='group1', n=10, fitness=2.0)}
        pop = Population(groups=groups)
        assert pop.fitness['group1'] == 2.0

    def test_population_fitness_multiple_groups(self):
        """Test: el fitness se establece para múltiples grupos."""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=10, fitness=0.5)
        }
        pop = Population(groups=groups)
        assert pop.fitness['group1'] == 2.0
        assert pop.fitness['group2'] == 0.5

    def test_population_fitness_stored_in_population(self):
        """Test: el fitness se almacena correctamente en la población."""
        groups = {'group1': SubPopulation(name='group1', n=10, fitness=3.0)}
        pop = Population(groups=groups)
        assert 'group1' in pop.fitness
        assert pop.fitness['group1'] == 3.0

    def test_population_find_greatest_fitness(self):
        """Test: encontrar el grupo con mayor fitness."""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=10, fitness=3.0),
            'group3': SubPopulation(name='group3', n=10, fitness=1.0)
        }
        pop = Population(groups=groups)

        greatest_fitness_group = pop.get_greatest_fitness_group()
        assert greatest_fitness_group == 'group2'
        
    def test_population_find_lowest_fitness(self):
        """Test: encontrar el grupo con menor fitness."""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=10, fitness=3.0),
            'group3': SubPopulation(name='group3', n=10, fitness=1.0)
        }
        pop = Population(groups=groups)

        lowest_fitness_group = pop.get_lowest_fitness_group()
        assert lowest_fitness_group == 'group3'

