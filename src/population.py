from dataclasses import dataclass
@dataclass
class SubPopulation:
    def __init__(self, name, n, fitness=1.0):
        """ Inicializa un subgrupo de la población.
        Args:
            name (str): Nombre del subgrupo.
            n (int): Número de individuos en el subgrupo.
            fitness (float): Aptitud del subgrupo (por defecto 1.0).
        """
        self.name = name
        self.n = n
        self.fitness = fitness

class Population:
    def __init__(self, groups):
        """ Inicializa la población con los grupos y sus tamaños.
        Args:
            groups (dict): Diccionario con el nombre del grupo y su tamaño.
        """
        self.groups = groups
        self.individuals = []
        self.fitness = {}
        for group, info in groups.items():
            n = info.n
            self.individuals.extend([group] * n)
            self.fitness[group] = info.fitness
        self.n = len(self.individuals)

    def get_greatest_fitness_group(self):
        """ Devuelve el grupo con la mayor aptitud.
        Returns:
            str: Nombre del grupo con la mayor aptitud.
        """
        greatest_fitness_group = max(self.fitness, key=self.fitness.get)
        return greatest_fitness_group

    def get_lowest_fitness_group(self):
        """ Devuelve el grupo con la menor aptitud.
        Returns:
            str: Nombre del grupo con la menor aptitud.
        """
        lowest_fitness_group = min(self.fitness, key=self.fitness.get)
        return lowest_fitness_group

    def evolve_deterministic(self):
        """ Simula una evolución determinista de la población.
        El grupo con mayor aptitud gana un individuo y el grupo con menor aptitud pierde un individuo.
        """
        greatest_fitness_group = self.get_greatest_fitness_group()
        lowest_fitness_group = self.get_lowest_fitness_group()

        # Aumentar un individuo al grupo con mayor aptitud
        self.individuals.append(greatest_fitness_group)
        self.fitness[greatest_fitness_group] += 1

        # Eliminar un individuo del grupo con menor aptitud
        if lowest_fitness_group in self.individuals:
            self.individuals.remove(lowest_fitness_group)
            self.fitness[lowest_fitness_group] -= 1

        # Actualizar el tamaño total de la población
        self.n = len(self.individuals)