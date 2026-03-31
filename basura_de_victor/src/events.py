from dataclasses import dataclass
from typing import Tuple

CloneId = Tuple[int, ...]

@dataclass
class Event:
    kind: str
    clone_id: CloneId
    rate: float



class EventsCollection:
    def __init__(self):
        self.events: list[Event] = []

    def add_event(self, event: Event) -> None:
        self.events.append(event)

    def clear(self) -> None:
        self.events.clear()

    def get_total_rate(self) -> float:
        return sum(event.rate for event in self.events)