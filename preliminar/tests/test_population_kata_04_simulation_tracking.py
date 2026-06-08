"""
Suite de tests para la clase Population.

Tests para la simulación de procesos evolutivos de Moran.
Kata 04: Seguimiento de la simulación.
"""

from preliminar.src.simulator import Simulator
from preliminar.src.population import Population
from preliminar.src.population import SubPopulation
from preliminar.src.evolution_engine_deterministic import EvolutionEngineDeterministic

class TestPopulationSimulationTracking:
    """Tests para el seguimiento de la simulación de la población."""

    def test_population_simulation_tracking_builder(self):
        """Test: seguimiento de la simulación con dos grupos."""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=10, fitness=1.0)
        }
        pop = Population(groups=groups)

        sim = Simulator(population=pop)
        assert len(sim.tracking) == 0
        assert isinstance(sim.evol_engine, EvolutionEngineDeterministic)

    def test_population_simulation_tracking_after_evolution(self):
        """Test: seguimiento de la simulación después de la evolución."""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=10, fitness=1.0)
        }
        pop = Population(groups=groups)

        sim = Simulator(population=pop)
        sim.run()
        sim.run()
        assert len(sim.tracking) == 2
        assert sim.tracking[0].count('group1') == 11
        assert sim.tracking[0].count('group2') == 9
        assert sim.tracking[1].count('group1') == 12
        assert sim.tracking[1].count('group2') == 8


    def test_tracker_to_dataframe(self):
        """Test: el seguimiento se puede convertir a un DataFrame."""
        import pandas as pd

        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=10, fitness=1.0)
        }
        pop = Population(groups=groups)

        sim = Simulator(population=pop)
        sim.run()
        sim.run()

        df = sim.get_tracking_dataframe()
        assert isinstance(df, pd.DataFrame)
        assert df.shape == (2, 20)

        assert df.iloc[0].value_counts()['group1'] == 11
        assert df.iloc[0].value_counts()['group2'] == 9
        assert df.iloc[1].value_counts()['group1'] == 12
        assert df.iloc[1].value_counts()['group2'] == 8

    def test_tracker_dataframe_summary(self):
        """Test: el DataFrame de seguimiento se puede resumir."""

        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=10, fitness=1.0)
        }
        pop = Population(groups=groups)

        sim = Simulator(population=pop)
        sim.run()
        sim.run()

        summary = sim.get_tracking_summary_df()
        assert summary.shape == (2, 4)
        assert summary.columns.tolist() == ['time', 'group1', 'group2', 'events']
        assert summary.loc[0, 'group1'] == 11
        assert summary.loc[0, 'group2'] == 9
        assert summary.loc[1, 'group1'] == 12
        assert summary.loc[1, 'group2'] == 8


    def test_tracker_on_group_lowest_zero(self):
        """Test: seguimiento cuando el grupo con menor fitness tiene 0 individuos."""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=0, fitness=1.0)
        }
        pop = Population(groups=groups)

        sim = Simulator(population=pop)
        sim.run()

        summary = sim.get_tracking_summary_df()

        assert summary.loc[0, 'group1'] == 10

        assert summary.loc[0, 'group2'] == 0
        

    def test_tracker_on_group_lowest_zero_three_groups(self):
        """Test: seguimiento cuando el grupo con menor fitness tiene 0 individuos."""
        groups = {
            'group1': SubPopulation(name='group1', n=10, fitness=2.0),
            'group2': SubPopulation(name='group2', n=6, fitness=1.0),
            'group3': SubPopulation(name='group3', n=0, fitness=1.0)
        }
        pop = Population(groups=groups)

        sim = Simulator(population=pop)
        sim.run()

        summary = sim.get_tracking_summary_df()

        assert summary.loc[0, 'group1'] == 11
        assert summary.loc[0, 'group2'] == 5
        assert summary.loc[0, 'group3'] == 0
        
