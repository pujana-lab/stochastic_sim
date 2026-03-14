import pandas as pd

class Simulator:
    """Clase para simular procesos evolutivos de Moran."""
    pass

    def __init__(self, population):
        """Inicializa el simulador con una población."""
        self.population = population
        self.tracking = []

    def run(self):
        self.population.evolve_deterministic()
        self.tracking.append(self.population.individuals.copy())

    def get_tracking_dataframe(self):
        df = pd.DataFrame(self.tracking)
        return df

    def get_tracking_summary_df(self):
        df = pd.DataFrame(self.tracking)
        summary_df = df.apply(lambda x: x.value_counts(), axis=1).fillna(0)
        summary_df.insert(0, 'time', range(len(summary_df)))
        return summary_df
