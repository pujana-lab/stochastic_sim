"""
Suite de tests para el motor de mutación determinista.

Kata 07: Mutación determinista básica.

En esta kata se comprueba:
- should_mutate devuelve True/False según el paso
- get_mutation devuelve los parámetros correctos
- El tracking incluye la columna 'events' (vacía por ahora)
"""
import pytest

from src.population import Population, SubPopulation
from src.evolution_engine_deterministic import EvolutionEngineDeterministic
from src.mutation_engine_deterministic import MutationEngineDeterministic
from src.simulator import Simulator


@pytest.fixture
def population():
    """Fixture: población con dos grupos."""
    return Population({
        'dominant': SubPopulation(name='dominant', n=60, fitness=2.0),
        'weak':     SubPopulation(name='weak',     n=40, fitness=1.0),
    })


@pytest.fixture
def mutation_engine():
    """Fixture: motor de mutación que muta cada 10 pasos."""
    return MutationEngineDeterministic(
        every_n_steps=10,
        victim_group='weak',
        new_group_name='mutant',
        new_fitness=1.5,
    )


class TestMutationEngineDeterministic:
    """Tests para el motor de mutación determinista."""

    def test_should_mutate_on_correct_step(self, mutation_engine):
        """should_mutate devuelve True en múltiplos de every_n_steps."""
        assert mutation_engine.should_mutate(step=10) is True
        assert mutation_engine.should_mutate(step=20) is True
        assert mutation_engine.should_mutate(step=30) is True

    def test_should_not_mutate_on_other_steps(self, mutation_engine):
        """should_mutate devuelve False en pasos que no son múltiplos."""
        assert mutation_engine.should_mutate(step=0)  is False
        assert mutation_engine.should_mutate(step=1)  is False
        assert mutation_engine.should_mutate(step=7)  is False
        assert mutation_engine.should_mutate(step=11) is False

    def test_get_mutation_returns_correct_params(self, mutation_engine, population):
        """get_mutation devuelve siempre los mismos parámetros."""
        victim, new_name, new_fitness = mutation_engine.get_mutation(population)
        assert victim    == 'weak'
        assert new_name  == 'mutant'
        assert new_fitness == 1.5


class TestSimulatorTrackingWithMutation:
    """Tests para el tracking del Simulator con columna 'events'."""

    def test_tracking_summary_has_events_column(self, population, mutation_engine):
        """El summary df incluye la columna 'events'."""
        simulator = Simulator(
            population=population,
            evol_engine=EvolutionEngineDeterministic(),
            mutation_engine=mutation_engine,
        )
        simulator.run()
        df = simulator.get_tracking_summary_df()
        assert 'events' in df.columns

    def test_tracking_events_column_is_empty(self, population, mutation_engine):
        """La columna 'events' está vacía (None/NaN) en todos los pasos."""
        simulator = Simulator(
            population=population,
            evol_engine=EvolutionEngineDeterministic(),
            mutation_engine=mutation_engine,
        )
        for _ in range(5):
            simulator.run()
        df = simulator.get_tracking_summary_df()
        assert df['events'].isna().all() or (df['events'] == '').all()


    def test_tracking_events_column_writes_mutation_events(self, population, mutation_engine):
        """La columna 'events' registra las mutaciones ocurridas."""
        simulator = Simulator(
            population=population,
            evol_engine=EvolutionEngineDeterministic(),
            mutation_engine=mutation_engine,
        )
        for step in range(1, 21):
            simulator.run()
            df = simulator.get_tracking_summary_df()
            if step % 10 == 0:
                # En pasos 10 y 20 debería haber un evento de mutación
                assert df.loc[step-1, 'events'] == 'Mutación: weak -> mutant (fitness=1.5)'
            else:
                # En otros pasos no debería haber eventos
                print(df.loc[step-1, 'events'])
                assert df.loc[step-1, 'events'] is None or df.loc[step-1, 'events'] == ''