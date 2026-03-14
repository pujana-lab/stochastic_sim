import pandas as pd
from src.evolution_engine_interface import EvolutionEngineInterface
from src.evolution_engine_deterministic import EvolutionEngineDeterministic
from src.population import Population
from src.mutation_engine_disabled import MutationEngineDisabled
from src.mutation_engine_interface import MutationEngineInterface

class Simulator:
    """Clase para simular procesos evolutivos de Moran."""

    def __init__(self, 
                population: Population, 
                evol_engine: EvolutionEngineInterface = EvolutionEngineDeterministic(),
                mutation_engine: MutationEngineInterface = MutationEngineDisabled()):
        """Inicializa el simulador con una población."""
        self.population = population
        self.evol_engine = evol_engine
        self.mutation_engine = mutation_engine
        self.tracking = []
        self.events = []  # Lista para registrar eventos de mutación (vacía por ahora)

    def run(self):
        """Ejecuta un paso del proceso de Moran.

        Si la víctima tiene 0 individuos, el paso no modifica la población
        pero sí registra el estado actual en el tracking.
        """
        group_reproductor = self.evol_engine.get_reproductor_group(self.population)
        group_victim = self.evol_engine.get_victim_group(self.population)

        victim_has_individuals = (
            group_victim is not None
            and self.population.groups[group_victim].n > 0
        )

        if victim_has_individuals:
            self.population.append_individual(group_reproductor)
            self.population.remove_individual(group_victim)

        if self.population.__len__() > self.population.n_init:
            raise ValueError("La población ha crecido más allá del tamaño inicial. Revisa la lógica de reproducción y eliminación.")
        
        if self.mutation_engine.should_mutate(step=len(self.tracking) + 1):
            victim, new_name, new_fitness = self.mutation_engine.get_mutation(self.population)
            if victim in self.population.groups and self.population.groups[victim].n > 0:
                self.population.mutate(
                    victim_group=victim, 
                    new_group_name=new_name, 
                    new_fitness=new_fitness
                )
                self.events.append(f"Mutación: {victim} -> {new_name} (fitness={new_fitness})")
            else:
                self.events.append(f"Mutación fallida: {victim} no tiene individuos para mutar.")

        else:
            self.events.append('')  # Por ahora no hay eventos de mutación, pero se reserva el espacio
        
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
        summary_df.insert(len(summary_df.columns), 'events', self.events)  # Columna para eventos de mutación (vacía por ahora)
        return summary_df
