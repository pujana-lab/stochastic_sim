from src.gillespie.event import Event
from src.gillespie.event_type import EventType


def test_event_creation():
    e = Event(kind=EventType.BIRTH, clone_id=(0, 1), rate=0.5)
    assert e.kind == EventType.BIRTH
    assert e.clone_id == (0, 1)
    assert e.rate == 0.5


def test_event_is_dataclass():
    e1 = Event(kind=EventType.DEATH, clone_id=(1,), rate=0.2)
    e2 = Event(kind=EventType.DEATH, clone_id=(1,), rate=0.2)
    assert e1 == e2


def test_event_clone_id_is_tuple():
    e = Event(kind=EventType.MUTATION, clone_id=(0, 1, 2), rate=0.1)
    assert isinstance(e.clone_id, tuple)

