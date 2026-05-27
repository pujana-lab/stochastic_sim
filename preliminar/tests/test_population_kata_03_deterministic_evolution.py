"""
Suite de tests para la clase Population.

Tests para la simulación de procesos evolutivos de Moran.
Kata 03: Evolución determinista.
"""

from preliminar.src.population import Population
from preliminar.src.population import SubPopulation

class TestPopulationDeterministicEvolution:
    """Tests para la evolución determinista de la población."""

    def test_population_deterministic_evolution(self):
        """Test: evolución determinista con dos grupos."""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=10, fitness=1.0)
        }
        pop = Population(groups=groups)

        pop.append_individual('group1')
        pop.remove_individual('group2')
        assert pop.individuals.count('group1') == 11
        assert pop.individuals.count('group2') == 9
