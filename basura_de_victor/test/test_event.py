from src.event import Event


# --- Event ---

def test_event_creation():
    e = Event(kind="birth", clone_id=(0, 1), rate=0.5)
    assert e.kind == "birth"
    assert e.clone_id == (0, 1)
    assert e.rate == 0.5


def test_event_is_dataclass():
    e1 = Event(kind="death", clone_id=(1,), rate=0.2)
    e2 = Event(kind="death", clone_id=(1,), rate=0.2)
    assert e1 == e2


def test_event_clone_id_is_tuple():
    e = Event(kind="mutation", clone_id=(0, 1, 2), rate=0.1)
    assert isinstance(e.clone_id, tuple)

