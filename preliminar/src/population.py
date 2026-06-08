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
        n_total = 0
        for group, info in groups.items():
            n = info.n
            n_total += n
            self.individuals.extend([group] * n)
            self.fitness[group] = info.fitness
        self.n_init = n_total
        self.n = self.__len__()

    def append_individual(self, group):
        """ Agrega un individuo a un grupo específico.
        Args:
            group (str): Nombre del grupo al que se agregará el individuo.
        """
        self.individuals.append(group)
        self.groups[group].n += 1
        self.n += 1

    def remove_individual(self, group):
        """ Elimina un individuo de un grupo específico.
        Args:
            group (str): Nombre del grupo del que se eliminará el individuo.
        """
        if group in self.individuals:
            self.individuals.remove(group)
            self.groups[group].n -= 1
            self.n -= 1

    def mutate(self, victim_group, new_group_name, new_fitness):
        """ Realiza una mutación eliminando un individuo del grupo víctima y agregando uno al nuevo grupo.
        Args:
            victim_group (str): Nombre del grupo del que se eliminará el individuo.
            new_group_name (str): Nombre del nuevo grupo al que se agregará el individuo.
            new_fitness (float): Aptitud del nuevo grupo.
        """
        self.remove_individual(victim_group)
        if new_group_name not in self.groups:
            self.groups[new_group_name] = SubPopulation(name=new_group_name, n=0, fitness=new_fitness)
            self.fitness[new_group_name] = new_fitness
        self.append_individual(new_group_name)

    def __len__(self):
        return len(self.individuals)