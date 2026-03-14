import pandas as pd
from src.evolution_engine_interface import EvolutionEngineInterface
from src.evolution_engine_deterministic import EvolutionEngineDeterministic
from src.population import Population


class Simulator:
    """Clase para simular procesos evolutivos de Moran."""

    def __init__(self, population: Population, evol: EvolutionEngineInterface = EvolutionEngineDeterministic()):
        """Inicializa el simulador con una población."""
        self.population = population
        self.evol = evol
        self.tracking = []

    def run(self):
        """Ejecuta un paso del proceso de Moran.

        Si la víctima tiene 0 individuos, el paso no modifica la población
        pero sí registra el estado actual en el tracking.
        """
        group_reproductor = self.evol.get_reproductor_group(self.population)
        group_victim = self.evol.get_victim_group(self.population)

        victim_has_individuals = (
            group_victim is not None
            and self.population.groups[group_victim].n > 0
        )

        if victim_has_individuals:
            self.population.append_individual(group_reproductor)
            self.population.remove_individual(group_victim)

        if self.population.__len__() > self.population.n_init:
            raise ValueError("La población ha crecido más allá del tamaño inicial. Revisa la lógica de reproducción y eliminación.")

        self.tracking.append(self.population.individuals.copy())

    def get_tracking_dataframe(self):
        df = pd.DataFrame(self.tracking)
        return df

    def get_tracking_summary_df(self):
        all_groups = list(self.population.groups.keys())
        df = pd.DataFrame(self.tracking)
        summary_df = df.apply(lambda x: x.value_counts(), axis=1).fillna(0)
        # Garantiza que todos los grupos aparecen como columnas, incluso con n=0
        for group in all_groups:
            if group not in summary_df.columns:
                summary_df[group] = 0
        summary_df = summary_df[all_groups].astype(int)
        summary_df.insert(0, 'time', range(len(summary_df)))
        return summary_df
