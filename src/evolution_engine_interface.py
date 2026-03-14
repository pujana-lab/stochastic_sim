from abc import ABC
class EvolutionEngineInterface(ABC):
    """Interfaz para el motor de evolución de Moran."""
    pass

    def evolve(self):
        """Evoluciona la población según el proceso de Moran."""
        raise NotImplementedError("El método evolve debe ser implementado por la clase hija.")