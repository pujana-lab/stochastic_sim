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



    def append_individual(self, group):
        """ Agrega un individuo a un grupo específico.
        Args:
            group (str): Nombre del grupo al que se agregará el individuo.
        """
        self.individuals.append(group)
        self.fitness[group] += 1
        self.n += 1

    def remove_individual(self, group):
        """ Elimina un individuo de un grupo específico.
        Args:
            group (str): Nombre del grupo del que se eliminará el individuo.
        """
        if group in self.individuals:
            self.individuals.remove(group)
            self.fitness[group] -= 1
            self.n -= 1