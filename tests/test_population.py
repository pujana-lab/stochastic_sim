"""
Suite de tests para la clase Population.

Tests para la simulación de procesos evolutivos de Moran.
"""

import pytest
from src.population import Population


class TestPopulationInitialization:
    """Tests para la inicialización de la población."""

    def test_create_population_with_n_individuals(self):
        """Test: crear una población de n individuos."""
        groups = {'group1': 100}
        pop = Population(groups=groups)
        assert len(pop.individuals) == 100
        assert pop.n == 100
    
    def test_population_size_stored(self):
        """Test: el tamaño se almacena correctamente."""
        groups = {'group1': 50}
        pop = Population(groups=groups)
        assert len(pop.individuals) == 50
        assert pop.n == 50
    
    def test_population_two_groups(self):
        """Test: la población se inicializa con dos grupos."""
        groups = {'group1': 50, 'group2': 50}
        pop = Population(groups=groups)
        assert len(pop.individuals) == 100
        assert pop.n == 100
        assert pop.individuals.count('group1') == 50
        assert pop.individuals.count('group2') == 50
