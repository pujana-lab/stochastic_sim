"""
Suite de tests para la clase Population.

Tests para la simulación de procesos evolutivos de Moran.
Kata 05: Refactorización de la evolución determinista.

"""

from preliminar.src.population import Population
from preliminar.src.population import SubPopulation
from preliminar.src.evolution_engine_bernoulli import EvolutionEngineBernoulli

class TestEvolutionEngineBernoulli:
    """Tests para la evolución no determinista de la población."""

    def test_population_bernoulli_forced_evolution(self):
        """Test: evolución no determinista con dos grupos, donde uno tiene fitness 0... en la práctica es determinista"""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=10, fitness=0)
        }
        pop = Population(groups=groups)
        engine = EvolutionEngineBernoulli(seed=42)

        reproductor_candidate = engine.get_reproductor_group(population=pop)
        victim_candidate = engine.get_victim_group(population=pop)

        assert reproductor_candidate == 'group1'
        assert victim_candidate == 'group2'


    def test_population_bernoulli_random_evolution(self):

        """Test: evolución no determinista con dos grupos, ambos con fitness 1... en la práctica es aleatoria"""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=1),
            'group2': SubPopulation(name='group2', n=10, fitness=1)
        }
        pop = Population(groups=groups)
        

        reproductive_candidates = []
        victim_candidates = []
        for i in range(100):
            engine = EvolutionEngineBernoulli(seed=i)
            reproductor_candidate = engine.get_reproductor_group(population=pop)
            victim_candidate = engine.get_victim_group(population=pop)
            reproductive_candidates.append(reproductor_candidate)
            victim_candidates.append(victim_candidate)

        assert reproductive_candidates.count('group1') > 0
        assert reproductive_candidates.count('group2') > 0
        assert victim_candidates.count('group1') > 0
        assert victim_candidates.count('group2') > 0
        assert reproductive_candidates.count('group1') + reproductive_candidates.count('group2') == 100
        assert victim_candidates.count('group1') + victim_candidates.count('group2') == 100
        assert reproductive_candidates.count('group1') / 100 > 0.3
        assert reproductive_candidates.count('group1') / 100 < 0.7


    def test_population_bernoulli_random_evolution_with_different_fitness(self):
        
        """Test: evolución no determinista con dos grupos, uno con fitness 2 y otro con fitness 1... en la práctica es aleatoria pero con sesgo"""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2),
            'group2': SubPopulation(name='group2', n=10, fitness=1)
        }
        pop = Population(groups=groups)
        

        reproductive_candidates = []
        victim_candidates = []
        for i in range(100):
            engine = EvolutionEngineBernoulli(seed=i)
            reproductor_candidate = engine.get_reproductor_group(population=pop)
            victim_candidate = engine.get_victim_group(population=pop)
            reproductive_candidates.append(reproductor_candidate)
            victim_candidates.append(victim_candidate)

        assert reproductive_candidates.count('group1') > reproductive_candidates.count('group2')
        assert victim_candidates.count('group2') > victim_candidates.count('group1')