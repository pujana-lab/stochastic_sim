import dataclasses

import numpy as np

from .event import Event


class RateMatrix:
    def __init__(self):
        self.events: list[Event] = []
        self.total_rate: float | None = None

    def add_event(self, event: Event) -> None:
        self.events.append(event)
        self.total_rate = None

    def clear(self) -> None:
        self.events.clear()
        self.total_rate = None

    def get_total_rate(self) -> float:
        if self.total_rate is None:
            self.total_rate = sum(event.rate for event in self.events)
        return self.total_rate
        
    def get_reaction_number_poisson(self, tau: float) -> np.ndarray:
        events = self.events
        rates = np.array([event.rate for event in events], dtype=float)
        K_j = np.random.poisson(lam=rates * tau)
        self.events = [
            dataclasses.replace(event, reaction_number=int(K_j[i]))
            for i, event in enumerate(events)
        ]
        return K_j

    def choose_event(self, u: float) -> Event:
        events = self.events
        rates = np.array([event.rate for event in events], dtype=float)
        total_rate = self.get_total_rate()
        cumulative = np.cumsum(rates / total_rate)

        idx = np.searchsorted(cumulative, u, side="right")
        if idx >= len(events):
            idx = len(events) - 1
        return events[idx]

