"""
Suite de tests para el motor de mutación Bernoulli.

Kata 08: Mutaciones con distribución binomial.
La mutación ocurre en cada paso con probabilidad `p`.
Con seed fijo, el comportamiento es reproducible.
"""
import pytest

from src.evolution_engine_deterministic import EvolutionEngineDeterministic
from src.mutation_engine_bernoulli import MutationEngineBernoulli
from src.population import Population, SubPopulation
from src.simulator import Simulator


class TestMutationEngineBernoulli:
    """Tests para el motor de mutación Bernoulli."""

    @pytest.fixture
    def engine(self):
        """Fixture: motor Bernoulli con p=0.5 y seed fijo."""
        return MutationEngineBernoulli(
            p=0.5,
            victim_group='weak',
            new_group_name='mutant',
            new_fitness=3.0,
            seed=42,
        )

    @pytest.fixture
    def engine_never(self):
        """Fixture: motor Bernoulli con p=0.0, nunca muta."""
        return MutationEngineBernoulli(
            p=0.0,
            victim_group='weak',
            new_group_name='mutant',
            new_fitness=3.0,
            seed=42,
        )

    @pytest.fixture
    def engine_always(self):
        """Fixture: motor Bernoulli con p=1.0, siempre muta."""
        return MutationEngineBernoulli(
            p=1.0,
            victim_group='weak',
            new_group_name='mutant',
            new_fitness=3.0,
            seed=42,
        )

    def test_never_mutates_with_p_zero(self, engine_never):
        """Con p=0.0, should_mutate siempre devuelve False."""
        for step in range(1, 20):
            assert engine_never.should_mutate(step) is False

    def test_always_mutates_with_p_one(self, engine_always):
        """Con p=1.0, should_mutate siempre devuelve True."""
        for step in range(1, 20):
            assert engine_always.should_mutate(step) is True

    def test_get_mutation_returns_correct_values(self, engine):
        """get_mutation devuelve el grupo víctima, nuevo nombre y nuevo fitness."""
        population = Population({
            'dominant': SubPopulation(name='dominant', n=6, fitness=2.0),
            'weak':     SubPopulation(name='weak',     n=4, fitness=1.0),
        })
        victim, new_name, new_fitness = engine.get_mutation(population)
        assert victim == 'weak'
        assert new_name == 'mutant'
        assert new_fitness == 3.0

    def test_same_seed_produces_same_sequence(self):
        """Dos motores con el mismo seed producen la misma secuencia de mutaciones."""
        engine_a = MutationEngineBernoulli(p=0.5, victim_group='weak', new_group_name='mutant', new_fitness=3.0, seed=99)
        engine_b = MutationEngineBernoulli(p=0.5, victim_group='weak', new_group_name='mutant', new_fitness=3.0, seed=99)
        results_a = [engine_a.should_mutate(step) for step in range(1, 20)]
        results_b = [engine_b.should_mutate(step) for step in range(1, 20)]
        assert results_a == results_b

    def test_different_seeds_produce_different_sequences(self):
        """Dos motores con distinto seed producen secuencias distintas."""
        engine_a = MutationEngineBernoulli(p=0.5, victim_group='weak', new_group_name='mutant', new_fitness=3.0, seed=1)
        engine_b = MutationEngineBernoulli(p=0.5, victim_group='weak', new_group_name='mutant', new_fitness=3.0, seed=2)
        results_a = [engine_a.should_mutate(step) for step in range(1, 20)]
        results_b = [engine_b.should_mutate(step) for step in range(1, 20)]
        assert results_a != results_b


class TestMutationBernoulliTracking:
    """Tests de integración: mutación Bernoulli con Simulator y tracking."""

    @pytest.fixture
    def population(self):
        """Fixture: población con dos grupos."""
        return Population({
            'dominant': SubPopulation(name='dominant', n=6, fitness=2.0),
            'weak':     SubPopulation(name='weak',     n=4, fitness=1.0),
        })

    def test_events_column_present_with_bernoulli_mutation(self, population):
        """El tracking incluye la columna 'events' con mutación Bernoulli."""
        mutation = MutationEngineBernoulli(
            p=1.0, victim_group='weak', new_group_name='mutant', new_fitness=3.0, seed=42
        )
        simulator = Simulator(
            population=population,
            evol_engine=EvolutionEngineDeterministic(),
            mutation_engine=mutation,
        )
        simulator.run()
        df = simulator.get_tracking_summary_df()
        assert 'events' in df.columns

    def test_mutation_event_recorded_when_p_one(self, population):
        """Con p=1.0, todos los pasos registran un evento de mutación."""
        mutation = MutationEngineBernoulli(
            p=1.0, victim_group='weak', new_group_name='mutant', new_fitness=3.0, seed=42
        )
        simulator = Simulator(
            population=population,
            evol_engine=EvolutionEngineDeterministic(),
            mutation_engine=mutation,
        )
        for _ in range(3):
            simulator.run()
        df = simulator.get_tracking_summary_df()
        assert (df['events'] != '').all()

    def test_no_mutation_event_recorded_when_p_zero(self, population):
        """Con p=0.0, ningún paso registra evento de mutación."""
        mutation = MutationEngineBernoulli(
            p=0.0, victim_group='weak', new_group_name='mutant', new_fitness=3.0, seed=42
        )
        simulator = Simulator(
            population=population,
            evol_engine=EvolutionEngineDeterministic(),
            mutation_engine=mutation,
        )
        for _ in range(5):
            simulator.run()
        df = simulator.get_tracking_summary_df()
        assert (df['events'] == '').all()
