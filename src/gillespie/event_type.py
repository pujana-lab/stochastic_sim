from enum import Enum

class EventType(Enum):
    BIRTH = "birth"
    DEATH = "death"
    MUTATION = "mutation"
    EXHAUSTION = "exhaustion"   # generic conversion to `converts_to` (immune -> exhausted, fibroblast -> caf)

    @property
    def consumes_parent(self) -> bool:
        """True if the event removes one cell from its own clone (all but BIRTH)."""
        return self is not EventType.BIRTH

