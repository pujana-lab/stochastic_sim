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
            n = info['n'] if isinstance(info, dict) else info
            self.individuals.extend([group] * n)
            if isinstance(info, dict) and 'fitness' in info:
                self.fitness[group] = info['fitness']
            else: 
                self.fitness[group] = 1.0  # Fitness por defecto
        self.n = len(self.individuals)