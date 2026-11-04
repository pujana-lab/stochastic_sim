from enum import Enum

class EventType(Enum):
    BIRTH = "birth"
    DEATH = "death"
    MUTATION = "mutation"
    CONVERSION = "conversion"   # one cell becomes one cell of `converts_to` (immune -> exhausted, fibroblast -> caf)

    @property
    def consumes_parent(self) -> bool:
        """True if the event removes one cell from its own clone (all but BIRTH)."""
        return self is not EventType.BIRTH

