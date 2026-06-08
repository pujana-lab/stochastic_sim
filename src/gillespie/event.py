from dataclasses import dataclass

from src.gillespie.cloneId import CloneId
from src.gillespie.event_type import EventType


@dataclass(frozen=True)
class Event:
    kind: EventType
    clone_id: CloneId
    rate: float
    clone_type: str = "base"
    reaction_number: int = 0

