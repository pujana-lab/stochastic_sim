from abc import ABC
class EvolutionEngineInterface(ABC):
    """Interfaz para el motor de evolución de Moran."""
    pass

    def get_reproductor_group(self, population) -> str:
        """Selecciona el grupo reproductor de la población."""
        raise NotImplementedError("Este método debe ser implementado por subclases.")

    def get_victim_group(self, population) -> str:
        """Selecciona el grupo víctima de la población."""
        raise NotImplementedError("Este método debe ser implementado por subclases.")