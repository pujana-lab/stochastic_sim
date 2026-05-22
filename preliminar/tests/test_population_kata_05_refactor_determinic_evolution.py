"""
Suite de tests para la clase Population.

Tests para la simulación de procesos evolutivos de Moran.
Kata 05: Refactorización de la evolución determinista.

"""

from preliminar.src.population import Population
from preliminar.src.population import SubPopulation
from preliminar.src.evolution_engine_deterministic import EvolutionEngineDeterministic

class TestEvolutionEngineDeterministic:
    """Tests para la evolución determinista de la población."""

    def test_population_deterministic_evolution(self):
        """Test: evolución determinista con dos grupos."""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=10, fitness=1.0)
        }
        pop = Population(groups=groups)
        engine = EvolutionEngineDeterministic()

        reproductor_candidate = engine.get_reproductor_group(population=pop)
        victim_candidate = engine.get_victim_group(population=pop)

        assert reproductor_candidate == 'group1'
        assert victim_candidate == 'group2'


    def skip_test_population_determ_when_min_group_is_zero(self):
        """Test: evolución determinista cuando un grupo tiene 0 individuos."""
        groups = {
            'group1': SubPopulation(name='group1', n=0, fitness=2.0),
            'group2': SubPopulation(name='group2', n=10, fitness=1.0)
        }
        pop = Population(groups=groups)
        engine = EvolutionEngineDeterministic()

        reproductor_candidate = engine.get_reproductor_group(population=pop)
        victim_candidate = engine.get_victim_group(population=pop)

        assert reproductor_candidate == 'group2'
        assert victim_candidate == 'group2'



