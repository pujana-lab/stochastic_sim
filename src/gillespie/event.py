from dataclasses import dataclass

from src.gillespie.cloneId import CloneId
from src.gillespie.clone_type import CloneType
from src.gillespie.event_type import EventType


@dataclass(frozen=True)
class Event:
    kind: EventType
    clone_id: CloneId
    clone_type: CloneType
    rate: float
    reaction_number: int = 0

