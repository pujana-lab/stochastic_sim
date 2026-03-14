from abc import ABC, abstractmethod
class EvolutionEngineInterface(ABC):
    """Interfaz para el motor de evolución de Moran."""
    pass

    @abstractmethod
    def get_reproductor_group(self, population) -> str:
        """Selecciona el grupo reproductor de la población."""
        ...

    @abstractmethod
    def get_victim_group(self, population) -> str:
        """Selecciona el grupo víctima de la población."""
        