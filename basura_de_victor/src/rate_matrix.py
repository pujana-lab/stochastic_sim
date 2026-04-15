import numpy as np
from .event import Event


class RateMatrix:
    def __init__(self):
        self.events: list[Event] = []

    def add_event(self, event: Event) -> None:
        self.events.append(event)

    def clear(self) -> None:
        self.events.clear()

    def get_total_rate(self) -> float:
        return sum(event.rate for event in self.events)

    def choose_event(self, u: float) -> Event:
        events = self.events
        rates = np.array([event.rate for event in events], dtype=float)
        cumulative = np.cumsum(rates / self.get_total_rate())

        idx = np.searchsorted(cumulative, u, side="right")
        if idx >= len(events):
            idx = len(events) - 1
        return events[idx]

