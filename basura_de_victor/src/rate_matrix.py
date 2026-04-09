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