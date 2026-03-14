import pandas as pd
from src.evolution_engine_interface import EvolutionEngineInterface
from src.evolution_engine_deterministic import EvolutionEngineDeterministic
from src.population import Population


class Simulator:
    """Clase para simular procesos evolutivos de Moran."""
    pass

    def __init__(self, population: Population , evol: EvolutionEngineInterface = EvolutionEngineDeterministic()):
        """Inicializa el simulador con una población."""
        self.population = population
        self.evol = evol
        self.tracking = []

    def run(self):
        group_reproductor = self.evol.get_reproductor_group(self.population)
        group_victim = self.evol.get_victim_group(self.population)
        self.population.append_individual(group_reproductor)
        self.population.remove_individual(group_victim)

        self.tracking.append(self.population.individuals.copy())

    def get_tracking_dataframe(self):
        df = pd.DataFrame(self.tracking)
        return df

    def get_tracking_summary_df(self):
        df = pd.DataFrame(self.tracking)
        summary_df = df.apply(lambda x: x.value_counts(), axis=1).fillna(0)
        summary_df.insert(0, 'time', range(len(summary_df)))
        return summary_df
