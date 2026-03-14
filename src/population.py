class Population:
    def __init__(self, groups):
        """ Inicializa la población con los grupos y sus tamaños.
        Args:
            groups (dict): Diccionario con el nombre del grupo y su tamaño.
        """
        self.groups = groups
        self.individuals = []
        for group, count in groups.items():
            self.individuals.extend([group] * count)
        self.n = len(self.individuals)